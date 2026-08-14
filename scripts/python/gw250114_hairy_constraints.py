"""Two-dimensional GW250114 constraints for the Zhen-Li hairy QNM atlas.

This calculation has two public-data branches:

1. the pyRing 220+221 deviation posterior in
   (Mf, af, domega_221, dtau_221), and
2. the direct RINGDOWN deviation posterior in
   (log f_220, log f_221, delta f_221), and
3. the RINGDOWN free-frequency time scans in
   (f_220, f_221, gamma_220, gamma_221).

The pyRing and free-spectrum products have constant native priors in their
fitted coordinates.  The RINGDOWN deviation run instead has a prior flat in
(Mf, chif, delta f_221); its induced density in the plotted log-frequency
coordinates is removed with the analytic transformation Jacobian.  No result
from this script is described as a newly sampled detector-strain likelihood.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import sys
from pathlib import Path


# Preload the installed copy before the optional local dependency path to
# avoid shadowing it with an incomplete local package copy.
import six  # noqa: F401


ROOT = Path(__file__).resolve().parents[2]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "results" / ".matplotlib"))
LOCAL_DEPS = ROOT / ".deps" / "python"
if LOCAL_DEPS.exists():
    sys.path.insert(0, str(LOCAL_DEPS))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import h5py
import joblib
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.special import logsumexp
from sklearn.mixture import GaussianMixture

from gw250114_bayesian_spectral_eft import kerr_tables, log_density
from gw250114_hairy_bayesian_benchmark import pseudo_observables_from_modes


SOLAR_MASS_TIME_SECONDS = 4.925490947e-6
DEFAULT_CONFIG = ROOT / "config" / "hairy_gw250114_publication.json"
DEFAULT_GRID = ROOT / "results" / "hairy_qnm_production_grid" / "hairy_qnm_production_grid.npz"
DEFAULT_PYRING_GMM = ROOT / "results" / "gw250114_bayesian_spectral_eft" / "spectral_likelihood_gmms.joblib"
DEFAULT_TIMESCAN = (
    ROOT
    / "data"
    / "raw"
    / "GW250114_data_release"
    / "data"
    / "220+221_amps_fs_gammas_merged_timescans.hdf5"
)
DEFAULT_RINGDOWN_DEVIATION = (
    ROOT
    / "data"
    / "raw"
    / "GW250114_data_release"
    / "data"
    / "220+221+df221+dg221_6M_f220meas_f221meas_df221meas_120Ksamps.hdf5"
)
DEFAULT_PRE_PRIOR = (
    ROOT / "results" / "article1_independent_prior" / "pre_minus40M_remnant_model_crosscheck_samples.csv"
)
DEFAULT_FULL_PRIOR = (
    ROOT / "results" / "gw250114_posterior_calibration" / "nrSur7dq4_selected_posterior_samples.csv"
)
DEFAULT_OUTPUT = ROOT / "results" / "gw250114_hairy_constraints"


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table: {path}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def quadrature_weights(axis: np.ndarray) -> np.ndarray:
    widths = np.diff(axis)
    output = np.empty(len(axis), dtype=float)
    output[0] = 0.5 * widths[0]
    output[-1] = 0.5 * widths[-1]
    output[1:-1] = 0.5 * (widths[:-1] + widths[1:])
    return output / np.sum(output)


def weighted_quantile(axis: np.ndarray, weights: np.ndarray, probabilities: list[float]) -> list[float]:
    # ``weights`` are nodal probability masses that already include the
    # trapezoidal quadrature weights.  Recover the nodal density before
    # constructing its continuous CDF; treating the masses as point masses
    # shifts every reported quantile by roughly half a grid cell (a uniform
    # [0, 0.5] prior would otherwise give q90=0.4375 on the 0.025 grid).
    axis = np.asarray(axis, dtype=float)
    masses = np.asarray(weights, dtype=float)
    density = masses / quadrature_weights(axis)
    intervals = 0.5 * (density[:-1] + density[1:]) * np.diff(axis)
    cumulative = np.concatenate(([0.0], np.cumsum(intervals)))
    if not np.isfinite(cumulative[-1]) or cumulative[-1] <= 0.0:
        raise ValueError("weighted quantile received zero or non-finite probability mass")
    cumulative /= cumulative[-1]
    return [float(np.interp(probability, cumulative, axis)) for probability in probabilities]


def load_prior_scenarios(
    pre_path: Path,
    full_path: Path,
    spin_min: float,
    spin_max: float,
    sample_count: int,
    seed: int,
) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    pre = np.genfromtxt(pre_path, names=True, delimiter=",")
    full = np.genfromtxt(full_path, names=True, delimiter=",")
    raw = {
        "pre_minus40M_BMR_HBR": (pre["final_mass_BMR2012"], pre["final_spin_HBR2016"]),
        "pre_minus40M_NRSur3dq8": (
            pre["final_mass_NRSur3dq8"],
            pre["final_spin_NRSur3dq8_projected"],
        ),
        "full_IMR_NRSur7dq4": (full["final_mass"], full["final_spin"]),
    }
    output: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for index, (name, (mass, spin)) in enumerate(raw.items()):
        finite = (
            np.isfinite(mass)
            & np.isfinite(spin)
            & (mass > 0.0)
            & (spin >= spin_min)
            & (spin <= spin_max)
        )
        mass = np.asarray(mass[finite], dtype=float)
        spin = np.asarray(spin[finite], dtype=float)
        rng = np.random.default_rng(seed + 1009 * index)
        if len(mass) > sample_count:
            selected = rng.choice(len(mass), sample_count, replace=False)
            mass = mass[selected]
            spin = spin[selected]
        output[name] = (mass, spin)
    return output


def read_timescan(path: Path) -> tuple[list[str], np.ndarray]:
    with h5py.File(path, "r") as handle:
        group = handle["samples"]
        columns = [value.decode() for value in group["axis0"][...]]
        values = np.asarray(group["block0_values"][...], dtype=float)
    return columns, values


def fit_timescan_gmm(
    samples: np.ndarray,
    component_counts: list[int],
    seed: int,
    maximum_samples: int | None = None,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    rng = np.random.default_rng(seed)
    if maximum_samples is not None and len(samples) > maximum_samples:
        samples = samples[rng.choice(len(samples), maximum_samples, replace=False)]
    order = rng.permutation(len(samples))
    split = max(1, int(0.8 * len(samples)))
    training = samples[order[:split]]
    test = samples[order[split:]]
    center = np.mean(training, axis=0)
    scale = np.std(training, axis=0, ddof=1)
    training_z = (training - center) / scale
    test_z = (test - center) / scale
    diagnostics: list[dict[str, object]] = []
    models: dict[int, GaussianMixture] = {}
    for count in component_counts:
        model = GaussianMixture(
            n_components=count,
            covariance_type="full",
            reg_covar=1.0e-6,
            max_iter=600,
            n_init=3,
            random_state=seed + count,
        )
        model.fit(training_z)
        models[count] = model
        diagnostics.append(
            {
                "components": count,
                "converged": bool(model.converged_),
                "iterations": int(model.n_iter_),
                "training_samples": len(training),
                "test_samples": len(test),
                "train_mean_log_density": float(np.mean(model.score_samples(training_z))),
                "test_mean_log_density": float(np.mean(model.score_samples(test_z))),
                "bic_training": float(model.bic(training_z)),
            }
        )
    preferred = int(max(diagnostics, key=lambda row: float(row["test_mean_log_density"]))["components"])
    return {
        "model": models[preferred],
        "center": center,
        "scale": scale,
        "preferred_components": preferred,
    }, diagnostics


def interpolate_modes_at_prior_spins(
    omega: np.ndarray,
    spins: np.ndarray,
    prior_spin: np.ndarray,
    overtone_indices: tuple[int, int],
) -> tuple[np.ndarray, np.ndarray]:
    shape = omega.shape[:2] + (len(prior_spin),)
    mode_220 = np.empty(shape, dtype=np.complex128)
    mode_221 = np.empty(shape, dtype=np.complex128)
    for ih in range(omega.shape[0]):
        for ia in range(omega.shape[1]):
            mode_220[ih, ia] = CubicSpline(spins, omega[ih, ia, :, overtone_indices[0]])(prior_spin)
            mode_221[ih, ia] = CubicSpline(spins, omega[ih, ia, :, overtone_indices[1]])(prior_spin)
    return mode_220, mode_221


def evaluate_pyring_surface(
    likelihood: dict[str, object],
    mode_220: np.ndarray,
    mode_221: np.ndarray,
    prior_mass: np.ndarray,
) -> tuple[np.ndarray, float, float]:
    model = likelihood["model"]
    center = np.asarray(likelihood["center"])
    scale = np.asarray(likelihood["scale"])
    kerr = kerr_tables()
    log_surface = np.full(mode_220.shape[:2], -np.inf)
    valid_min = 1.0
    for ih in range(mode_220.shape[0]):
        for ia in range(mode_220.shape[1]):
            mass, spin, domega, dtau, valid = pseudo_observables_from_modes(
                prior_mass,
                mode_220[ih, ia],
                mode_221[ih, ia],
                kerr,
            )
            valid_min = min(valid_min, float(np.mean(valid)))
            if np.any(valid):
                points = np.column_stack((mass[valid], spin[valid], domega[valid], dtau[valid]))
                density = log_density(model, points, center, scale)
                log_surface[ih, ia] = float(logsumexp(density) - math.log(len(prior_mass)))
    gr_log = float(np.mean(log_surface[:, 0]))
    return log_surface, gr_log, valid_min


def evaluate_ringdown_surface(
    likelihood: dict[str, object],
    mode_220: np.ndarray,
    mode_221: np.ndarray,
    prior_mass: np.ndarray,
) -> tuple[np.ndarray, float]:
    model = likelihood["model"]
    center = np.asarray(likelihood["center"])
    scale = np.asarray(likelihood["scale"])
    mass_seconds = prior_mass * SOLAR_MASS_TIME_SECONDS
    log_surface = np.full(mode_220.shape[:2], -np.inf)
    for ih in range(mode_220.shape[0]):
        for ia in range(mode_220.shape[1]):
            om0 = mode_220[ih, ia]
            om1 = mode_221[ih, ia]
            points = np.column_stack(
                (
                    om0.real / (2.0 * np.pi * mass_seconds),
                    om1.real / (2.0 * np.pi * mass_seconds),
                    -om0.imag / mass_seconds,
                    -om1.imag / mass_seconds,
                )
            )
            density = log_density(model, points, center, scale)
            log_surface[ih, ia] = float(logsumexp(density) - math.log(len(prior_mass)))
    gr_log = float(np.mean(log_surface[:, 0]))
    return log_surface, gr_log


def evaluate_ringdown_deviation_surface(
    likelihood: dict[str, object],
    mode_220: np.ndarray,
    mode_221: np.ndarray,
    prior_mass: np.ndarray,
) -> tuple[np.ndarray, float, float]:
    """Evaluate the direct public RINGDOWN spectroscopy product.

    The public delta-f variable is log(f_221 measured / f_221 inferred
    from the Kerr 220 mode).  The same pseudo-Kerr map used for pyRing
    therefore supplies the model prediction without linearization.
    """
    model = likelihood["model"]
    center = np.asarray(likelihood["center"])
    scale = np.asarray(likelihood["scale"])
    kerr = kerr_tables()
    spin_grid = np.asarray(kerr["spin"], dtype=float)
    dlog_220 = np.gradient(np.log(np.asarray(kerr["omega_220"]).real), spin_grid)
    dlog_221 = np.gradient(np.log(np.asarray(kerr["omega_221"]).real), spin_grid)
    mass_seconds = prior_mass * SOLAR_MASS_TIME_SECONDS
    log_surface = np.full(mode_220.shape[:2], -np.inf)
    valid_min = 1.0
    for ih in range(mode_220.shape[0]):
        for ia in range(mode_220.shape[1]):
            pseudo_mass, pseudo_spin, domega, _, valid = pseudo_observables_from_modes(
                prior_mass,
                mode_220[ih, ia],
                mode_221[ih, ia],
                kerr,
            )
            valid &= domega > -1.0
            valid_min = min(valid_min, float(np.mean(valid)))
            if np.any(valid):
                points = np.column_stack(
                    (
                        np.log(mode_220[ih, ia][valid].real / (2.0 * np.pi * mass_seconds[valid])),
                        np.log(mode_221[ih, ia][valid].real / (2.0 * np.pi * mass_seconds[valid])),
                        np.log1p(domega[valid]),
                    )
                )
                posterior_density = log_density(model, points, center, scale)
                # The RINGDOWN deviation analysis samples a prior uniform in
                # (Mf, chif, delta-f), while the released projection used here
                # is expressed in (log f220, log f221, delta-f).  For
                # f_n = F_n(chi)/M (and f_221 multiplied by exp(delta-f)),
                # |d y / d theta| = |d_chi log F220 - d_chi log F221| / M.
                # Dividing the posterior density by the induced prior density
                # therefore amounts to adding log of this Jacobian.  Overall
                # prior-volume constants cancel in all within-branch ratios.
                derivative_difference = np.abs(
                    np.interp(pseudo_spin[valid], spin_grid, dlog_220 - dlog_221)
                )
                log_jacobian = np.log(derivative_difference) - np.log(pseudo_mass[valid])
                density = posterior_density + log_jacobian
                log_surface[ih, ia] = float(logsumexp(density) - math.log(len(prior_mass)))
    gr_log = float(np.mean(log_surface[:, 0]))
    return log_surface, gr_log, valid_min


def summarize_surface(
    log_surface: np.ndarray,
    gr_log: float,
    h0s: np.ndarray,
    alphas: np.ndarray,
) -> tuple[dict[str, object], np.ndarray, np.ndarray, np.ndarray]:
    prior_mass = quadrature_weights(h0s)[:, None] * quadrature_weights(alphas)[None, :]
    log_evidence = float(logsumexp(log_surface + np.log(prior_mass)))
    posterior_mass = np.exp(log_surface + np.log(prior_mass) - log_evidence)
    alpha_marginal = np.sum(posterior_mass, axis=0)
    h0_marginal = np.sum(posterior_mass, axis=1)
    alpha_q05, alpha_q50, alpha_q90, alpha_q95 = weighted_quantile(
        alphas, alpha_marginal, [0.05, 0.50, 0.90, 0.95]
    )
    h0_q05, h0_q50, h0_q90, h0_q95 = weighted_quantile(
        h0s, h0_marginal, [0.05, 0.50, 0.90, 0.95]
    )
    ratio = posterior_mass / prior_mass
    kl_bits = float(np.sum(posterior_mass * np.log2(np.maximum(ratio, 1.0e-300))))
    order = np.argsort(ratio.ravel())[::-1]
    cumulative = np.cumsum(posterior_mass.ravel()[order])
    hpd = np.empty(posterior_mass.size, dtype=float)
    hpd[order] = cumulative
    hpd = hpd.reshape(posterior_mass.shape)
    best = np.unravel_index(int(np.argmax(log_surface)), log_surface.shape)
    summary = {
        "log_evidence_relative_constant": log_evidence,
        "log_gr_evidence_relative_constant": gr_log,
        "log_bayes_factor_hairy_vs_gr": log_evidence - gr_log,
        "bayes_factor_hairy_vs_gr": float(np.exp(np.clip(log_evidence - gr_log, -700.0, 700.0))),
        "maximum_log_likelihood_ratio_vs_gr": float(np.max(log_surface) - gr_log),
        "best_h0": float(h0s[best[0]]),
        "best_alpha_hair": float(alphas[best[1]]),
        "alpha_q05": alpha_q05,
        "alpha_median": alpha_q50,
        "alpha_q90": alpha_q90,
        "alpha_q95": alpha_q95,
        "h0_q05": h0_q05,
        "h0_median": h0_q50,
        "h0_q90": h0_q90,
        "h0_q95": h0_q95,
        "kl_divergence_bits": kl_bits,
    }
    return summary, posterior_mass, prior_mass, hpd


def conditional_limits(
    log_surface: np.ndarray,
    h0s: np.ndarray,
    alphas: np.ndarray,
    branch: str,
    prior_scenario: str,
    start_time: float | None,
) -> list[dict[str, object]]:
    alpha_weight = quadrature_weights(alphas)
    rows: list[dict[str, object]] = []
    for ih, h0 in enumerate(h0s):
        weights = np.exp(log_surface[ih] - np.max(log_surface[ih])) * alpha_weight
        weights /= np.sum(weights)
        q50, q68, q90, q95 = weighted_quantile(alphas, weights, [0.50, 0.68, 0.90, 0.95])
        uniform_q90 = float(alphas[0] + 0.9 * (alphas[-1] - alphas[0]))
        rows.append(
            {
                "likelihood_branch": branch,
                "prior_scenario": prior_scenario,
                "start_time_M": "" if start_time is None else start_time,
                "h0": float(h0),
                "alpha_median_conditional": q50,
                "alpha_q68_upper_conditional": q68,
                "alpha_q90_upper_conditional": q90,
                "alpha_q95_upper_conditional": q95,
                "uniform_prior_q90": uniform_q90,
                "q90_shift_from_uniform_prior": q90 - uniform_q90,
                "log_likelihood_ratio_alpha_max_vs_gr": float(log_surface[ih, -1] - log_surface[ih, 0]),
            }
        )
    return rows


def surface_rows(
    log_surface: np.ndarray,
    posterior_mass: np.ndarray,
    prior_mass: np.ndarray,
    hpd: np.ndarray,
    h0s: np.ndarray,
    alphas: np.ndarray,
    gr_log: float,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for ih, h0 in enumerate(h0s):
        for ia, alpha in enumerate(alphas):
            rows.append(
                {
                    "h0": float(h0),
                    "alpha_hair": float(alpha),
                    "log_likelihood_relative_gr": float(log_surface[ih, ia] - gr_log),
                    "posterior_probability_mass": float(posterior_mass[ih, ia]),
                    "prior_probability_mass": float(prior_mass[ih, ia]),
                    "posterior_to_prior_density_ratio": float(posterior_mass[ih, ia] / prior_mass[ih, ia]),
                    "hpd_enclosed_probability": float(hpd[ih, ia]),
                    "inside_90pct_hpd": bool(hpd[ih, ia] <= 0.90),
                    "inside_95pct_hpd": bool(hpd[ih, ia] <= 0.95),
                }
            )
    return rows


def save_qnm_shift_figure(
    path_base: Path,
    omega: np.ndarray,
    spins: np.ndarray,
    h0s: np.ndarray,
    alphas: np.ndarray,
    overtones: np.ndarray,
    reference_spin: float = 0.68,
) -> None:
    idx0 = int(np.where(overtones == 0)[0][0])
    idx1 = int(np.where(overtones == 1)[0][0])
    modes: list[np.ndarray] = []
    for index in (idx0, idx1):
        values = np.empty((len(h0s), len(alphas)), dtype=np.complex128)
        for ih in range(len(h0s)):
            for ia in range(len(alphas)):
                values[ih, ia] = CubicSpline(spins, omega[ih, ia, :, index])(reference_spin)
        modes.append(values)
    fig, axes = plt.subplots(2, 2, figsize=(9.0, 6.8), constrained_layout=True)
    labels = [
        ("220 frequency shift [%]", 100.0 * (modes[0].real / modes[0][:, :1].real - 1.0)),
        ("221 frequency shift [%]", 100.0 * (modes[1].real / modes[1][:, :1].real - 1.0)),
        ("220 damping-rate shift [%]", 100.0 * (modes[0].imag / modes[0][:, :1].imag - 1.0)),
        ("221 damping-rate shift [%]", 100.0 * (modes[1].imag / modes[1][:, :1].imag - 1.0)),
    ]
    for ax, (title, values) in zip(axes.ravel(), labels):
        mesh = ax.pcolormesh(alphas, h0s, values, shading="auto", cmap="viridis")
        fig.colorbar(mesh, ax=ax, label="percent")
        ax.set_title(title)
        ax.set_xlabel(r"$\alpha$")
        ax.set_ylabel(r"$h_0$")
    fig.suptitle(f"Direct QNM atlas at a={reference_spin:.2f}")
    fig.savefig(path_base.with_suffix(".png"), dpi=220)
    fig.savefig(path_base.with_suffix(".pdf"))
    plt.close(fig)


def save_primary_posterior_figure(
    path_base: Path,
    h0s: np.ndarray,
    alphas: np.ndarray,
    posterior_mass: np.ndarray,
    prior_mass: np.ndarray,
    title: str,
) -> None:
    density_ratio = posterior_mass / prior_mass
    fig, ax = plt.subplots(figsize=(6.8, 5.1), constrained_layout=True)
    mesh = ax.pcolormesh(alphas, h0s, density_ratio, shading="auto", cmap="magma")
    colorbar = fig.colorbar(mesh, ax=ax)
    colorbar.set_label("posterior / prior density")
    if float(np.min(density_ratio)) < 1.0 < float(np.max(density_ratio)):
        contour = ax.contour(
            alphas,
            h0s,
            density_ratio,
            levels=[1.0],
            colors=["white"],
            linewidths=1.5,
        )
        ax.clabel(contour, fmt={1.0: "posterior = prior"}, fontsize=8)
    ax.set_xlabel(r"hair parameter $\alpha$")
    ax.set_ylabel(r"hair parameter $h_0$")
    ax.set_title(title)
    fig.savefig(path_base.with_suffix(".png"), dpi=240)
    fig.savefig(path_base.with_suffix(".pdf"))
    plt.close(fig)


def save_conditional_limit_figure(path_base: Path, rows: list[dict[str, object]]) -> None:
    selected = [
        row
        for row in rows
        if row["likelihood_branch"] in {"pyRing", "RINGDOWN_deviation"}
    ]
    fig, ax = plt.subplots(figsize=(7.2, 4.8), constrained_layout=True)
    groups: dict[tuple[str, str], list[dict[str, object]]] = {}
    for row in selected:
        key = (str(row["likelihood_branch"]), str(row["start_time_M"]))
        groups.setdefault(key, []).append(row)
    for (branch, start), group in groups.items():
        group.sort(key=lambda row: float(row["h0"]))
        display = {"pyRing": "pyRing", "RINGDOWN_deviation": "RINGDOWN"}.get(branch, branch)
        label = display if start == "" else f"{display}, t0={start} M"
        ax.plot(
            [float(row["h0"]) for row in group],
            [float(row["alpha_q90_upper_conditional"]) for row in group],
            marker="o",
            ms=3,
            label=label,
        )
    ax.axhline(0.45, color="black", ls=":", lw=1.2, label="90% quantile of uniform prior")
    ax.set_xlabel(r"$h_0$")
    ax.set_ylabel(r"conditional 90% upper quantile of $\alpha$")
    plotted = [float(row["alpha_q90_upper_conditional"]) for row in selected]
    lower = min(plotted + [0.45])
    upper = max(plotted + [0.45])
    padding = max(0.005, 0.12 * (upper - lower))
    ax.set_ylim(max(0.0, lower - padding), min(0.5, upper + padding))
    ax.legend(fontsize=8, ncol=2)
    ax.grid(alpha=0.25)
    fig.savefig(path_base.with_suffix(".png"), dpi=220)
    fig.savefig(path_base.with_suffix(".pdf"))
    plt.close(fig)


def save_start_time_figure(path_base: Path, rows: list[dict[str, object]]) -> None:
    rows = sorted(rows, key=lambda row: float(row["start_time_M"]))
    times = [float(row["start_time_M"]) for row in rows]
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 4.0), constrained_layout=True)
    axes[0].plot(times, [float(row["log_bayes_factor_hairy_vs_gr"]) for row in rows], marker="o")
    axes[0].axhline(0.0, color="black", lw=1.0)
    axes[0].set_ylabel(r"$\ln B_{\rm hair/GR}$")
    axes[1].plot(times, [float(row["kl_divergence_bits"]) for row in rows], marker="o", color="tab:orange")
    axes[1].set_ylabel("information gain [bits]")
    for ax in axes:
        ax.set_xlabel(r"ringdown start time [$M_f$]")
        ax.grid(alpha=0.25)
    fig.savefig(path_base.with_suffix(".png"), dpi=220)
    fig.savefig(path_base.with_suffix(".pdf"))
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--grid", type=Path, default=DEFAULT_GRID)
    parser.add_argument("--pyring-gmm", type=Path, default=DEFAULT_PYRING_GMM)
    parser.add_argument("--ringdown-deviation", type=Path, default=DEFAULT_RINGDOWN_DEVIATION)
    parser.add_argument("--timescan", type=Path, default=DEFAULT_TIMESCAN)
    parser.add_argument("--pre-prior", type=Path, default=DEFAULT_PRE_PRIOR)
    parser.add_argument("--full-prior", type=Path, default=DEFAULT_FULL_PRIOR)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--reuse-fitted-density-models",
        action="store_true",
        help="Reuse deterministic GMM bundles already present in the output directory.",
    )
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    config = json.loads(args.config.read_text(encoding="utf-8"))
    inference = config["inference"]

    grid = np.load(args.grid)
    spins = np.asarray(grid["spins"], dtype=float)
    alphas = np.asarray(grid["alphas"], dtype=float)
    h0s = np.asarray(grid["h0s"], dtype=float)
    overtones = np.asarray(grid["overtones"], dtype=int)
    omega = np.asarray(grid["omega"], dtype=np.complex128)
    index_220 = int(np.where(overtones == 0)[0][0])
    index_221 = int(np.where(overtones == 1)[0][0])

    priors = load_prior_scenarios(
        args.pre_prior,
        args.full_prior,
        float(spins[0]),
        float(spins[-1]),
        int(inference["prior_samples"]),
        int(inference["seed"]),
    )
    prior_rows = [
        {
            "prior_scenario": name,
            "samples": len(mass),
            "mass_q05": float(np.quantile(mass, 0.05)),
            "mass_median": float(np.median(mass)),
            "mass_q95": float(np.quantile(mass, 0.95)),
            "spin_q05": float(np.quantile(spin, 0.05)),
            "spin_median": float(np.median(spin)),
            "spin_q95": float(np.quantile(spin, 0.95)),
            "independent_of_ringdown": not name.startswith("full_IMR"),
        }
        for name, (mass, spin) in priors.items()
    ]
    write_csv(args.output / "remnant_prior_scenarios.csv", prior_rows)

    py_bundle = joblib.load(args.pyring_gmm)
    py_likelihoods = {
        int(count): {"model": model, "center": py_bundle["center"], "scale": py_bundle["scale"]}
        for count, model in py_bundle["models"].items()
    }
    preferred_py_components = max(int(value) for value in inference["pyRing_gmm_components"])
    py_likelihood = py_likelihoods[preferred_py_components]

    with h5py.File(args.ringdown_deviation, "r") as handle:
        direct_ringdown_samples = np.column_stack(
            (
                np.log(np.asarray(handle["f_220"], dtype=float)),
                np.log(np.asarray(handle["f_221"], dtype=float)),
                np.asarray(handle["df_221"], dtype=float),
            )
        )
    direct_component_counts = [
        int(value) for value in inference["ringdown_deviation_gmm_components"]
    ]
    direct_model_path = args.output / "ringdown_deviation_likelihoods.joblib"
    direct_diagnostic_path = args.output / "ringdown_deviation_gmm_diagnostics.csv"
    if args.reuse_fitted_density_models and direct_model_path.exists() and direct_diagnostic_path.exists():
        direct_models = joblib.load(direct_model_path)
        direct_gmm_rows = list(csv.DictReader(direct_diagnostic_path.open(encoding="utf-8", newline="")))
    else:
        direct_models: dict[int, dict[str, object]] = {}
        direct_gmm_rows: list[dict[str, object]] = []
        # Fit every requested representation so that the nonlinear parameter
        # surface can be checked against density-estimator complexity.
        for components in direct_component_counts:
            bundle, diagnostics = fit_timescan_gmm(
                direct_ringdown_samples,
                [components],
                int(inference["seed"]) + 50000,
                int(inference["ringdown_deviation_fit_samples"]),
            )
            direct_models[components] = bundle
            direct_gmm_rows.extend(diagnostics)
        preferred_direct_components = int(
            max(direct_gmm_rows, key=lambda row: float(row["test_mean_log_density"]))["components"]
        )
        for row in direct_gmm_rows:
            row["selected"] = int(row["components"]) == preferred_direct_components
        write_csv(direct_diagnostic_path, direct_gmm_rows)
        joblib.dump(direct_models, direct_model_path)
    preferred_direct_components = int(
        max(direct_gmm_rows, key=lambda row: float(row["test_mean_log_density"]))["components"]
    )

    timescan_columns, timescan_values = read_timescan(args.timescan)
    column = {name: index for index, name in enumerate(timescan_columns)}
    timescan_model_path = args.output / "ringdown_timescan_likelihoods.joblib"
    timescan_diagnostic_path = args.output / "ringdown_timescan_gmm_diagnostics.csv"
    if args.reuse_fitted_density_models and timescan_model_path.exists() and timescan_diagnostic_path.exists():
        time_likelihoods = joblib.load(timescan_model_path)
    else:
        time_likelihoods: dict[float, dict[str, object]] = {}
        gmm_rows: list[dict[str, object]] = []
        for offset, start_time in enumerate(inference["ringdown_start_times_M"]):
            selected = timescan_values[np.isclose(timescan_values[:, column["start time [M]"]], float(start_time))]
            spectral = selected[:, [column["f_220"], column["f_221"], column["g_220"], column["g_221"]]]
            likelihood, diagnostics = fit_timescan_gmm(
                spectral,
                [int(value) for value in inference["ringdown_gmm_components"]],
                int(inference["seed"]) + 10000 + offset,
            )
            time_likelihoods[float(start_time)] = likelihood
            for row in diagnostics:
                gmm_rows.append({"start_time_M": start_time, **row, "selected": int(row["components"]) == likelihood["preferred_components"]})
        write_csv(timescan_diagnostic_path, gmm_rows)
        joblib.dump(time_likelihoods, timescan_model_path)

    # Use one coherent pre-peak remnant prescription for the headline result.
    # The analytic Barausse/Hofmann fit is retained as a separate robustness
    # calculation; averaging the two families would introduce an arbitrary
    # mixture weight.
    primary_prior_name = "pre_minus40M_NRSur3dq8"
    primary_mass, primary_spin = priors[primary_prior_name]
    primary_modes = interpolate_modes_at_prior_spins(omega, spins, primary_spin, (index_220, index_221))

    scenario_rows: list[dict[str, object]] = []
    conditional_rows: list[dict[str, object]] = []
    primary_products: tuple[np.ndarray, np.ndarray, np.ndarray] | None = None
    direct_products: tuple[np.ndarray, np.ndarray, np.ndarray] | None = None

    # GMM representation robustness for the preferred pyRing branch.
    for components in [int(value) for value in inference["pyRing_gmm_components"]]:
        surface, gr_log, valid_min = evaluate_pyring_surface(
            py_likelihoods[components], primary_modes[0], primary_modes[1], primary_mass
        )
        summary, posterior, prior_mass_grid, hpd = summarize_surface(surface, gr_log, h0s, alphas)
        scenario_rows.append(
            {
                "likelihood_branch": "pyRing",
                "start_time_M": 6.0,
                "prior_scenario": primary_prior_name,
                "gmm_components": components,
                "primary": components == preferred_py_components,
                "minimum_valid_fraction": valid_min,
                **summary,
            }
        )
        if components == preferred_py_components:
            primary_products = (posterior, prior_mass_grid, hpd)
            write_csv(
                args.output / "primary_posterior_surface.csv",
                surface_rows(surface, posterior, prior_mass_grid, hpd, h0s, alphas, gr_log),
            )
            np.savez_compressed(
                args.output / "primary_posterior_surface.npz",
                h0s=h0s,
                alphas=alphas,
                log_likelihood=surface,
                posterior_mass=posterior,
                prior_mass=prior_mass_grid,
                hpd_enclosed_probability=hpd,
                log_gr=gr_log,
            )
            conditional_rows.extend(
                conditional_limits(surface, h0s, alphas, "pyRing", primary_prior_name, None)
            )

    # Direct RINGDOWN deviation-posterior branch at the published 6M start.
    direct_primary_row: dict[str, object] | None = None
    for components in direct_component_counts:
        surface, gr_log, valid_min = evaluate_ringdown_deviation_surface(
            direct_models[components], primary_modes[0], primary_modes[1], primary_mass
        )
        summary, posterior, prior_mass_grid, hpd = summarize_surface(
            surface, gr_log, h0s, alphas
        )
        row = {
            "likelihood_branch": "RINGDOWN_deviation",
            "start_time_M": 6.0,
            "prior_scenario": primary_prior_name,
            "gmm_components": components,
            "primary": components == preferred_direct_components,
            "minimum_valid_fraction": valid_min,
            **summary,
        }
        scenario_rows.append(row)
        if components == preferred_direct_components:
            direct_primary_row = row
            direct_products = (posterior, prior_mass_grid, hpd)
            write_csv(
                args.output / "ringdown_deviation_posterior_surface.csv",
                surface_rows(
                    surface,
                    posterior,
                    prior_mass_grid,
                    hpd,
                    h0s,
                    alphas,
                    gr_log,
                ),
            )
            np.savez_compressed(
                args.output / "ringdown_deviation_posterior_surface.npz",
                h0s=h0s,
                alphas=alphas,
                log_likelihood=surface,
                posterior_mass=posterior,
                prior_mass=prior_mass_grid,
                hpd_enclosed_probability=hpd,
                log_gr=gr_log,
            )
            conditional_rows.extend(
                conditional_limits(
                    surface,
                    h0s,
                    alphas,
                    "RINGDOWN_deviation",
                    primary_prior_name,
                    6.0,
                )
            )

    # Remnant-prior robustness for both public likelihood branches at 6M.
    reference_time = float(inference["primary_start_time_M"])
    for prior_name, (prior_mass_values, prior_spin_values) in priors.items():
        modes = interpolate_modes_at_prior_spins(omega, spins, prior_spin_values, (index_220, index_221))
        for branch in ("pyRing", "RINGDOWN_deviation", "RINGDOWN_free_spectrum"):
            # The primary pyRing/prior combination, including all requested
            # GMM component counts, was already evaluated above.
            if branch in {"pyRing", "RINGDOWN_deviation", "RINGDOWN_free_spectrum"} and prior_name == primary_prior_name:
                continue
            if branch == "pyRing":
                surface, gr_log, valid_min = evaluate_pyring_surface(
                    py_likelihood, modes[0], modes[1], prior_mass_values
                )
                components = preferred_py_components
            elif branch == "RINGDOWN_deviation":
                surface, gr_log, valid_min = evaluate_ringdown_deviation_surface(
                    direct_models[preferred_direct_components],
                    modes[0],
                    modes[1],
                    prior_mass_values,
                )
                components = preferred_direct_components
            else:
                likelihood = time_likelihoods[reference_time]
                surface, gr_log = evaluate_ringdown_surface(
                    likelihood, modes[0], modes[1], prior_mass_values
                )
                valid_min = 1.0
                components = int(likelihood["preferred_components"])
            summary, _, _, _ = summarize_surface(surface, gr_log, h0s, alphas)
            scenario_rows.append(
                {
                    "likelihood_branch": branch,
                    "start_time_M": reference_time,
                    "prior_scenario": prior_name,
                    "gmm_components": components,
                    "primary": branch == "pyRing" and prior_name == primary_prior_name,
                    "minimum_valid_fraction": valid_min,
                    **summary,
                }
            )

    # Start-time robustness with the independent pre-merger remnant prior.
    start_rows: list[dict[str, object]] = []
    for start_time, likelihood in sorted(time_likelihoods.items()):
        surface, gr_log = evaluate_ringdown_surface(
            likelihood, primary_modes[0], primary_modes[1], primary_mass
        )
        summary, _, _, _ = summarize_surface(surface, gr_log, h0s, alphas)
        row = {
            "likelihood_branch": "RINGDOWN_free_spectrum",
            "start_time_M": start_time,
            "prior_scenario": primary_prior_name,
            "gmm_components": int(likelihood["preferred_components"]),
            "primary": start_time == reference_time,
            "minimum_valid_fraction": 1.0,
            **summary,
        }
        start_rows.append(row)
        scenario_rows.append(row)
        if start_time == reference_time:
            conditional_rows.extend(
                conditional_limits(
                    surface,
                    h0s,
                    alphas,
                    "RINGDOWN_free_spectrum",
                    primary_prior_name,
                    reference_time,
                )
            )

    # Hair-prior sensitivity reuses the primary pyRing likelihood surface.
    primary_surface_npz = np.load(args.output / "primary_posterior_surface.npz")
    primary_surface = np.asarray(primary_surface_npz["log_likelihood"])
    primary_gr_log = float(primary_surface_npz["log_gr"])
    hair_prior_rows: list[dict[str, object]] = []
    for alpha_max in config["hair_prior_robustness"]["alpha_max_values"]:
        for h0_max in config["hair_prior_robustness"]["h0_max_values"]:
            ia = np.where(alphas <= float(alpha_max) + 1.0e-12)[0]
            ih = np.where(h0s <= float(h0_max) + 1.0e-12)[0]
            sub_surface = primary_surface[np.ix_(ih, ia)]
            summary, _, _, _ = summarize_surface(
                sub_surface, primary_gr_log, h0s[ih], alphas[ia]
            )
            hair_prior_rows.append(
                {
                    "alpha_prior_min": float(alphas[ia[0]]),
                    "alpha_prior_max": float(alphas[ia[-1]]),
                    "h0_prior_min": float(h0s[ih[0]]),
                    "h0_prior_max": float(h0s[ih[-1]]),
                    **summary,
                }
            )

    write_csv(args.output / "constraint_scenario_summary.csv", scenario_rows)
    write_csv(args.output / "start_time_robustness.csv", start_rows)
    write_csv(args.output / "conditional_alpha_limits.csv", conditional_rows)
    write_csv(args.output / "hair_prior_robustness.csv", hair_prior_rows)

    save_qnm_shift_figure(args.output / "hairy_qnm_shift_atlas_a0p68", omega, spins, h0s, alphas, overtones)
    if primary_products is None:
        raise RuntimeError("primary posterior products were not produced")
    save_primary_posterior_figure(
        args.output / "gw250114_hairy_primary_posterior",
        h0s,
        alphas,
        primary_products[0],
        primary_products[1],
        "GW250114: pyRing hairy-BH projection",
    )
    if direct_products is None:
        raise RuntimeError("direct RINGDOWN posterior products were not produced")
    save_primary_posterior_figure(
        args.output / "gw250114_hairy_ringdown_deviation_posterior",
        h0s,
        alphas,
        direct_products[0],
        direct_products[1],
        "GW250114: RINGDOWN hairy-BH projection",
    )
    figure_limits = [
        row
        for row in conditional_rows
        if row["prior_scenario"] == primary_prior_name
        and (row["likelihood_branch"] == "pyRing" or float(row["start_time_M"]) == reference_time)
    ]
    save_conditional_limit_figure(args.output / "conditional_alpha_90_limits", figure_limits)
    save_start_time_figure(args.output / "ringdown_start_time_robustness", start_rows)

    primary_summary = next(
        row
        for row in scenario_rows
        if row["likelihood_branch"] == "pyRing"
        and row["prior_scenario"] == primary_prior_name
        and int(row["gmm_components"]) == preferred_py_components
    )
    reference_timescan = next(
        row
        for row in start_rows
        if float(row["start_time_M"]) == reference_time
    )
    if direct_primary_row is None:
        raise RuntimeError("direct RINGDOWN deviation result was not produced")
    output_files = sorted(path for path in args.output.iterdir() if path.is_file())
    summary = {
        "event": config["event"],
        "model": config["model"],
        "analysis_scope": "public_posterior_spectral_likelihood_not_new_detector_strain_sampling",
        "qnm_grid": str(args.grid.relative_to(ROOT)),
        "primary_likelihood": "pyRing_220_221_delta_posterior",
        "independent_crosscheck": "RINGDOWN_direct_deviation_posterior",
        "start_time_crosscheck": "RINGDOWN_free_frequency_and_damping_timescan",
        "primary_remnant_prior": primary_prior_name,
        "primary_result": primary_summary,
        "reference_RINGDOWN_result": direct_primary_row,
        "reference_RINGDOWN_timescan_result": reference_timescan,
        "start_times_tested_M": sorted(time_likelihoods),
        "prior_scenarios": list(priors),
        "publication_outputs_ready": True,
        "files": {path.name: sha256(path) for path in output_files},
    }
    (args.output / "constraint_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps({key: value for key, value in summary.items() if key != "files"}, indent=2))


if __name__ == "__main__":
    main()
