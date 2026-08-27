"""Controlled evolving-Kerr injection recovered with stationary hairy QNMs.

The injected spectrum follows the NR-calibrated mass and spin relaxation used
elsewhere in the project.  The recovery profiles independent sine/cosine mode
amplitudes and the remnant mass and spin at every point of a two-dimensional
hair grid.  This is a model-misspecification control, not an exact dynamical
spacetime waveform.
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


import six  # noqa: F401


ROOT = Path(__file__).resolve().parents[2]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "results" / ".matplotlib"))
LOCAL_DEPS = ROOT / ".deps" / "python"
if LOCAL_DEPS.exists():
    sys.path.insert(0, str(LOCAL_DEPS))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import minimize
from scipy.special import logsumexp

from evolving_kerr_bias import cumulative_trapezoid, time_axis, whiten_columns, whiten_vector
from gw250114_hairy_constraints import quadrature_weights, weighted_quantile


SOLAR_MASS_TIME_SECONDS = 4.925490947e-6
DEFAULT_CONFIG = ROOT / "config" / "hairy_gw250114_publication.json"
DEFAULT_GRID = ROOT / "results" / "hairy_qnm_production_grid" / "hairy_qnm_production_grid.npz"
DEFAULT_PRIOR_SAMPLES = (
    ROOT / "results" / "article1_independent_prior" / "pre_minus40M_remnant_model_crosscheck_samples.csv"
)
DEFAULT_OUTPUT = ROOT / "results" / "hairy_evolving_kerr_control"


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
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


def load_primary_nr_sur_covariance(path: Path, mass0: float, spin0: float) -> np.ndarray:
    """Gaussian control covariance for the primary NRSur3dq8 remnant map."""
    table = np.genfromtxt(path, names=True, delimiter=",")
    mass = np.asarray(table["final_mass_NRSur3dq8"], dtype=float)
    spin = np.asarray(table["final_spin_NRSur3dq8_projected"], dtype=float)
    finite = np.isfinite(mass) & np.isfinite(spin) & (mass > 0.0) & (spin > 0.0) & (spin < 1.0)
    offsets = np.column_stack((np.log(mass[finite] / mass0), spin[finite] - spin0))
    return np.cov(offsets, rowvar=False, ddof=1)


def build_interpolators(grid: np.lib.npyio.NpzFile) -> dict[str, RegularGridInterpolator]:
    spins = np.asarray(grid["spins"], dtype=float)
    alphas = np.asarray(grid["alphas"], dtype=float)
    h0s = np.asarray(grid["h0s"], dtype=float)
    overtones = np.asarray(grid["overtones"], dtype=int)
    omega = np.asarray(grid["omega"], dtype=np.complex128)
    output: dict[str, RegularGridInterpolator] = {}
    for overtone, mode in ((0, "220"), (1, "221")):
        index = int(np.where(overtones == overtone)[0][0])
        output[mode] = RegularGridInterpolator(
            (h0s, alphas, spins),
            omega[:, :, :, index],
            method="cubic",
            bounds_error=True,
        )
    return output


def mode_omega(
    interpolator: RegularGridInterpolator,
    h0: float,
    alpha: float,
    spin: np.ndarray | float,
) -> np.ndarray:
    spin_array = np.asarray(spin, dtype=float)
    points = np.column_stack(
        (
            np.full(spin_array.size, h0),
            np.full(spin_array.size, alpha),
            spin_array.ravel(),
        )
    )
    return np.asarray(interpolator(points), dtype=np.complex128).reshape(spin_array.shape)


def stationary_hairy_basis(
    modes: list[str],
    times: np.ndarray,
    mass: float,
    spin: float,
    alpha: float,
    h0: float,
    interpolators: dict[str, RegularGridInterpolator],
) -> np.ndarray:
    columns: list[np.ndarray] = []
    mass_seconds = mass * SOLAR_MASS_TIME_SECONDS
    for mode in modes:
        omega = complex(mode_omega(interpolators[mode], h0, alpha, np.asarray([spin]))[0])
        frequency = omega.real / (2.0 * np.pi * mass_seconds)
        tau = -mass_seconds / omega.imag
        decay = np.exp(-times / tau)
        phase = 2.0 * np.pi * frequency * times
        columns.extend((decay * np.cos(phase), decay * np.sin(phase)))
    return np.column_stack(columns)


def evolving_kerr_basis(
    modes: list[str],
    segment_times: np.ndarray,
    analysis_start_seconds: float,
    mass_final: float,
    spin_final: float,
    epsilon_mass: float,
    epsilon_spin: float,
    tau_mass_M: float,
    tau_spin_M: float,
    interpolators: dict[str, RegularGridInterpolator],
) -> np.ndarray:
    final_mass_time = mass_final * SOLAR_MASS_TIME_SECONDS
    base_step = float(segment_times[1] - segment_times[0])
    maximum_time = analysis_start_seconds + float(segment_times[-1])
    integration_times = np.arange(0.0, maximum_time + base_step, base_step)
    mass = mass_final * (1.0 + epsilon_mass * np.exp(-integration_times / (tau_mass_M * final_mass_time)))
    spin = spin_final * (1.0 + epsilon_spin * np.exp(-integration_times / (tau_spin_M * final_mass_time)))
    target = analysis_start_seconds + segment_times
    columns: list[np.ndarray] = []
    for mode in modes:
        omega = mode_omega(interpolators[mode], 1.0, 0.0, spin)
        frequency = omega.real / (2.0 * np.pi * mass * SOLAR_MASS_TIME_SECONDS)
        inverse_tau = -omega.imag / (mass * SOLAR_MASS_TIME_SECONDS)
        phase = cumulative_trapezoid(2.0 * np.pi * frequency, integration_times)
        log_decay = cumulative_trapezoid(inverse_tau, integration_times)
        phase_target = np.interp(target, integration_times, phase)
        decay_target = np.exp(-np.interp(target, integration_times, log_decay))
        columns.extend((decay_target * np.cos(phase_target), decay_target * np.sin(phase_target)))
    return np.column_stack(columns)


def profile_point(
    data_white: np.ndarray,
    modes: list[str],
    times: np.ndarray,
    sample_rate: float,
    mass0: float,
    spin0: float,
    alpha: float,
    h0: float,
    interpolators: dict[str, RegularGridInterpolator],
    inverse_prior_covariance: np.ndarray,
) -> tuple[float, float, float, bool]:
    def objective(parameters: np.ndarray) -> float:
        delta_ln_mass, delta_spin = parameters
        mass = mass0 * math.exp(float(delta_ln_mass))
        spin = spin0 + float(delta_spin)
        if not (0.1 <= spin <= 0.9):
            return 1.0e30
        basis = stationary_hairy_basis(modes, times, mass, spin, alpha, h0, interpolators)
        white_basis = whiten_columns(basis, sample_rate)
        coefficients, _, _, _ = np.linalg.lstsq(white_basis, data_white, rcond=1.0e-12)
        residual = data_white - white_basis @ coefficients
        offset = np.asarray([delta_ln_mass, delta_spin], dtype=float)
        return float(np.dot(residual, residual) + offset @ inverse_prior_covariance @ offset)

    fit = minimize(
        objective,
        np.zeros(2),
        method="L-BFGS-B",
        bounds=[(-0.12, 0.12), (-0.12, 0.12)],
        options={"ftol": 1.0e-12, "gtol": 1.0e-8, "maxiter": 250},
    )
    # L-BFGS-B occasionally reports an abnormal line-search termination when
    # the exact Kerr solution lies on a numerically flat ridge (alpha=0, or
    # h0 very close to 2).  Re-evaluate such points with a derivative-free
    # bounded method.  Keep the better finite solution and record convergence
    # only when one of the optimizers completed normally.
    if not fit.success:
        fallback = minimize(
            objective,
            np.asarray(fit.x, dtype=float),
            method="Powell",
            bounds=[(-0.12, 0.12), (-0.12, 0.12)],
            options={"xtol": 1.0e-9, "ftol": 1.0e-12, "maxiter": 500},
        )
        if np.isfinite(fallback.fun) and (
            not np.isfinite(fit.fun) or float(fallback.fun) <= float(fit.fun) + 1.0e-12
        ):
            fit = fallback
    finite = bool(np.isfinite(fit.fun) and np.all(np.isfinite(fit.x)))
    return (
        float(fit.fun),
        mass0 * math.exp(float(fit.x[0])),
        spin0 + float(fit.x[1]),
        bool(fit.success and finite),
    )


def summarize_grid(
    log_likelihood: np.ndarray,
    alphas: np.ndarray,
    h0s: np.ndarray,
) -> tuple[dict[str, object], np.ndarray]:
    prior = quadrature_weights(h0s)[:, None] * quadrature_weights(alphas)[None, :]
    log_evidence = float(logsumexp(log_likelihood + np.log(prior)))
    posterior = np.exp(log_likelihood + np.log(prior) - log_evidence)
    alpha_marginal = np.sum(posterior, axis=0)
    q50, q90, q95 = weighted_quantile(alphas, alpha_marginal, [0.50, 0.90, 0.95])
    gr_log = float(np.mean(log_likelihood[:, 0]))
    best = np.unravel_index(int(np.argmax(log_likelihood)), log_likelihood.shape)
    return {
        "log_bayes_factor_hairy_vs_gr": log_evidence - gr_log,
        "maximum_log_likelihood_ratio_vs_gr": float(np.max(log_likelihood) - gr_log),
        "best_alpha_hair": float(alphas[best[1]]),
        "best_h0": float(h0s[best[0]]),
        "alpha_median": q50,
        "alpha_q90": q90,
        "alpha_q95": q95,
        "posterior_probability_alpha_ge_0p1": float(np.sum(posterior[:, alphas >= 0.1])),
    }, posterior


def save_control_figure(path_base: Path, summaries: list[dict[str, object]]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 4.0), constrained_layout=True)
    for injection, color in (("stationary_Kerr", "black"), ("NR_calibrated_evolving_Kerr", "tab:red")):
        rows = sorted(
            [row for row in summaries if row["injection"] == injection],
            key=lambda row: float(row["start_time_M"]),
        )
        times = [float(row["start_time_M"]) for row in rows]
        label = {
            "stationary_Kerr": "stationary Kerr",
            "NR_calibrated_evolving_Kerr": "NR-calibrated evolving Kerr",
        }[injection]
        axes[0].plot(times, [float(row["maximum_log_likelihood_ratio_vs_gr"]) for row in rows], marker="o", color=color, label=label)
        axes[1].plot(times, [float(row["alpha_q90"]) for row in rows], marker="o", color=color, label=label)
    axes[0].set_ylabel(r"maximum $\Delta\ln L$ relative to Kerr")
    axes[1].set_ylabel(r"90% upper quantile of $\alpha$")
    axes[1].axhline(0.45, color="grey", ls=":", label="uniform-prior value")
    for ax in axes:
        ax.set_xlabel(r"start time [$M_f$]")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8)
    fig.savefig(path_base.with_suffix(".png"), dpi=220)
    fig.savefig(path_base.with_suffix(".pdf"))
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--grid", type=Path, default=DEFAULT_GRID)
    parser.add_argument("--prior-samples", type=Path, default=DEFAULT_PRIOR_SAMPLES)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    config = json.loads(args.config.read_text(encoding="utf-8"))["evolving_kerr_control"]
    grid = np.load(args.grid)
    interpolators = build_interpolators(grid)
    mass0 = float(config["mass_detector_msun"])
    spin0 = float(config["spin"])
    covariance = load_primary_nr_sur_covariance(args.prior_samples, mass0, spin0)
    inverse_covariance = np.linalg.inv(covariance)
    sample_rate = float(config["sample_rate_hz"])
    times = time_axis(sample_rate, float(config["duration_s"]))
    modes = list(config["modes"])
    amplitudes = np.asarray([value for mode in modes for value in config["amplitudes"][mode]], dtype=float)
    alpha_axis = np.asarray(config["recovery_alpha_values"], dtype=float)
    h0_axis = np.asarray(config["recovery_h0_values"], dtype=float)
    final_mass_time = mass0 * SOLAR_MASS_TIME_SECONDS

    stationary_peak = stationary_hairy_basis(modes, times, mass0, spin0, 0.0, 1.0, interpolators) @ amplitudes
    scale = float(config["network_snr"]) / math.sqrt(float(np.dot(whiten_vector(stationary_peak, sample_rate), whiten_vector(stationary_peak, sample_rate))))
    amplitudes *= scale

    rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    for injection_name, epsilon_mass, epsilon_spin in (
        ("stationary_Kerr", 0.0, 0.0),
        ("NR_calibrated_evolving_Kerr", float(config["epsilon_mass"]), float(config["epsilon_spin"])),
    ):
        for start_time_M in config["start_times_M"]:
            start_seconds = float(start_time_M) * final_mass_time
            injection_basis = evolving_kerr_basis(
                modes,
                times,
                start_seconds,
                mass0,
                spin0,
                epsilon_mass,
                epsilon_spin,
                float(config["tau_mass_M"]),
                float(config["tau_spin_M"]),
                interpolators,
            )
            data_white = whiten_vector(injection_basis @ amplitudes, sample_rate)
            segment_snr = math.sqrt(float(np.dot(data_white, data_white)))
            chi2 = np.empty((len(h0_axis), len(alpha_axis)), dtype=float)
            for ih, h0 in enumerate(h0_axis):
                for ia, alpha in enumerate(alpha_axis):
                    value, fitted_mass, fitted_spin, success = profile_point(
                        data_white,
                        modes,
                        times,
                        sample_rate,
                        mass0,
                        spin0,
                        float(alpha),
                        float(h0),
                        interpolators,
                        inverse_covariance,
                    )
                    chi2[ih, ia] = value
                    rows.append(
                        {
                            "injection": injection_name,
                            "start_time_M": float(start_time_M),
                            "segment_snr": segment_snr,
                            "h0": float(h0),
                            "alpha_hair": float(alpha),
                            "profile_chi2": value,
                            "fitted_mass": fitted_mass,
                            "fitted_spin": fitted_spin,
                            "optimizer_success": success,
                        }
                    )
            log_likelihood = -0.5 * chi2
            summary, posterior = summarize_grid(log_likelihood, alpha_axis, h0_axis)
            summary_rows.append(
                {
                    "injection": injection_name,
                    "start_time_M": float(start_time_M),
                    "segment_snr": segment_snr,
                    **summary,
                }
            )
            np.savez_compressed(
                args.output / f"surface_{injection_name}_t{float(start_time_M):g}M.npz",
                h0s=h0_axis,
                alphas=alpha_axis,
                profile_chi2=chi2,
                log_likelihood=log_likelihood,
                posterior_mass=posterior,
            )

    write_csv(args.output / "evolving_kerr_hairy_profile_grid.csv", rows)
    write_csv(args.output / "evolving_kerr_hairy_summary.csv", summary_rows)
    save_control_figure(args.output / "evolving_kerr_false_hair_control", summary_rows)

    stationary = [row for row in summary_rows if row["injection"] == "stationary_Kerr"]
    dynamic = [row for row in summary_rows if row["injection"] == "NR_calibrated_evolving_Kerr"]
    summary = {
        "scope": "controlled_evolving_Kerr_spectral_injection_not_exact_dynamical_metric",
        "network_snr_at_peak": float(config["network_snr"]),
        "drift": {
            "epsilon_mass": float(config["epsilon_mass"]),
            "epsilon_spin": float(config["epsilon_spin"]),
            "tau_mass_M": float(config["tau_mass_M"]),
            "tau_spin_M": float(config["tau_spin_M"]),
        },
        "stationary_maximum_false_log_likelihood_gain": max(float(row["maximum_log_likelihood_ratio_vs_gr"]) for row in stationary),
        "evolving_maximum_false_log_likelihood_gain": max(float(row["maximum_log_likelihood_ratio_vs_gr"]) for row in dynamic),
        "evolving_maximum_false_hair_best_alpha": max(float(row["best_alpha_hair"]) for row in dynamic),
        "all_optimizers_successful": all(bool(row["optimizer_success"]) for row in rows),
    }
    files = sorted(path for path in args.output.iterdir() if path.is_file())
    summary["files"] = {path.name: sha256(path) for path in files}
    (args.output / "evolving_kerr_control_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps({key: value for key, value in summary.items() if key != "files"}, indent=2))


if __name__ == "__main__":
    main()
