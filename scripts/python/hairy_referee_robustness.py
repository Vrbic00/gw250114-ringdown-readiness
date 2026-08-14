"""Robustness calculations requested during referee-level review.

The script performs four calculations that must remain distinct:

* paired remnant-subsampling uncertainty for the event Bayes factors;
* propagation of the alternative horizon-index convention through the event;
* repeated held-out GMM tail diagnostics, including event-level propagation;
* an audit of the small-deformation and effective-charge conditions quoted by Li.

The alternative recurrence is a convention sensitivity, not a draw from a
calibrated physical error distribution.  It is therefore reported as a second
likelihood surface rather than converted into a statistical uncertainty band.
"""

from __future__ import annotations

import argparse
import csv
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

import h5py
import joblib
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import logsumexp
from scipy.stats import ks_2samp
from sklearn.mixture import GaussianMixture

from gw250114_bayesian_spectral_eft import kerr_tables, log_density
from gw250114_hairy_bayesian_benchmark import pseudo_observables_from_modes
from gw250114_hairy_constraints import (
    SOLAR_MASS_TIME_SECONDS,
    conditional_limits,
    evaluate_pyring_surface,
    evaluate_ringdown_deviation_surface,
    interpolate_modes_at_prior_spins,
    quadrature_weights,
    summarize_surface,
)
from hairy_continued_fraction import approximate_horizons


DEFAULT_CONFIG = ROOT / "config" / "hairy_gw250114_publication.json"
DEFAULT_GRID = ROOT / "results" / "hairy_qnm_production_grid" / "hairy_qnm_production_grid.npz"
DEFAULT_SYSTEMATICS = ROOT / "results" / "hairy_qnm_internal_systematics" / "internal_systematics.npz"
DEFAULT_PRE_PRIOR = ROOT / "results" / "article1_independent_prior" / "pre_minus40M_remnant_model_crosscheck_samples.csv"
DEFAULT_PYRING_POSTERIOR = ROOT / "data" / "raw" / "GW250114_data_release" / "data" / "pyring_220_221_delta_posterior_without_frequencies.dat"
DEFAULT_PYRING_GMM = ROOT / "results" / "gw250114_bayesian_spectral_eft" / "spectral_likelihood_gmms.joblib"
DEFAULT_RINGDOWN_POSTERIOR = ROOT / "data" / "raw" / "GW250114_data_release" / "data" / "220+221+df221+dg221_6M_f220meas_f221meas_df221meas_120Ksamps.hdf5"
DEFAULT_RINGDOWN_GMM = ROOT / "results" / "gw250114_hairy_constraints" / "ringdown_deviation_likelihoods.joblib"
DEFAULT_RINGDOWN_DIAGNOSTICS = ROOT / "results" / "gw250114_hairy_constraints" / "ringdown_deviation_gmm_diagnostics.csv"
DEFAULT_OUTPUT = ROOT / "results" / "hairy_referee_robustness"


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table: {path}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def read_primary_remnant_samples(
    path: Path,
    spin_min: float,
    spin_max: float,
    sample_count: int,
    seed: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    table = np.genfromtxt(path, names=True, delimiter=",")
    mass = np.asarray(table["final_mass_NRSur3dq8"], dtype=float)
    spin = np.asarray(table["final_spin_NRSur3dq8_projected"], dtype=float)
    finite = np.isfinite(mass) & np.isfinite(spin) & (mass > 0.0) & (spin >= spin_min) & (spin <= spin_max)
    mass = mass[finite]
    spin = spin[finite]
    # NRSur3dq8 is the second scenario in load_prior_scenarios().  Reusing its
    # offset reproduces the exact headline subsample.
    rng = np.random.default_rng(seed + 1009)
    selected = rng.choice(len(mass), min(sample_count, len(mass)), replace=False)
    return mass, spin, selected


def direct_preferred_components(path: Path) -> int:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    return int(max(rows, key=lambda row: float(row["test_mean_log_density"]))["components"])


def py_log_terms(
    likelihood: dict[str, object],
    mode_220: np.ndarray,
    mode_221: np.ndarray,
    prior_mass: np.ndarray,
    h0s: np.ndarray,
    alphas: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    model = likelihood["model"]
    center = np.asarray(likelihood["center"])
    scale = np.asarray(likelihood["scale"])
    kerr = kerr_tables()
    weights = quadrature_weights(h0s)[:, None] * quadrature_weights(alphas)[None, :]
    log_numerator = np.full(len(prior_mass), -np.inf)
    log_denominator = np.full(len(prior_mass), -np.inf)
    for ih in range(len(h0s)):
        for ia in range(len(alphas)):
            mass, spin, domega, dtau, valid = pseudo_observables_from_modes(
                prior_mass, mode_220[ih, ia], mode_221[ih, ia], kerr
            )
            terms = np.full(len(prior_mass), -np.inf)
            if np.any(valid):
                points = np.column_stack((mass[valid], spin[valid], domega[valid], dtau[valid]))
                terms[valid] = log_density(model, points, center, scale)
            log_numerator = np.logaddexp(log_numerator, terms + math.log(float(weights[ih, ia])))
            if ih == 0 and ia == 0:
                log_denominator = terms.copy()
    return log_numerator, log_denominator


def direct_log_terms(
    likelihood: dict[str, object],
    mode_220: np.ndarray,
    mode_221: np.ndarray,
    prior_mass: np.ndarray,
    h0s: np.ndarray,
    alphas: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    model = likelihood["model"]
    center = np.asarray(likelihood["center"])
    scale = np.asarray(likelihood["scale"])
    kerr = kerr_tables()
    spin_grid = np.asarray(kerr["spin"], dtype=float)
    dlog_220 = np.gradient(np.log(np.asarray(kerr["omega_220"]).real), spin_grid)
    dlog_221 = np.gradient(np.log(np.asarray(kerr["omega_221"]).real), spin_grid)
    mass_seconds = prior_mass * SOLAR_MASS_TIME_SECONDS
    weights = quadrature_weights(h0s)[:, None] * quadrature_weights(alphas)[None, :]
    log_numerator = np.full(len(prior_mass), -np.inf)
    log_denominator = np.full(len(prior_mass), -np.inf)
    for ih in range(len(h0s)):
        for ia in range(len(alphas)):
            pseudo_mass, pseudo_spin, domega, _, valid = pseudo_observables_from_modes(
                prior_mass, mode_220[ih, ia], mode_221[ih, ia], kerr
            )
            valid &= domega > -1.0
            terms = np.full(len(prior_mass), -np.inf)
            if np.any(valid):
                points = np.column_stack(
                    (
                        np.log(mode_220[ih, ia][valid].real / (2.0 * np.pi * mass_seconds[valid])),
                        np.log(mode_221[ih, ia][valid].real / (2.0 * np.pi * mass_seconds[valid])),
                        np.log1p(domega[valid]),
                    )
                )
                posterior_density = log_density(model, points, center, scale)
                derivative_difference = np.abs(
                    np.interp(pseudo_spin[valid], spin_grid, dlog_220 - dlog_221)
                )
                terms[valid] = posterior_density + np.log(derivative_difference) - np.log(pseudo_mass[valid])
            log_numerator = np.logaddexp(log_numerator, terms + math.log(float(weights[ih, ia])))
            if ih == 0 and ia == 0:
                log_denominator = terms.copy()
    return log_numerator, log_denominator


def log_ratio_for_indices(log_num: np.ndarray, log_den: np.ndarray, indices: np.ndarray) -> float:
    return float(
        logsumexp(log_num[indices])
        - math.log(len(indices))
        - logsumexp(log_den[indices])
        + math.log(len(indices))
    )


def subsampling_rows(
    branch: str,
    log_num: np.ndarray,
    log_den: np.ndarray,
    selected: np.ndarray,
    repeats: int,
    seed: int,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    rng = np.random.default_rng(seed)
    rows: list[dict[str, object]] = []
    for repeat in range(repeats):
        indices = rng.choice(len(log_num), len(selected), replace=False)
        rows.append(
            {
                "branch": branch,
                "repeat": repeat,
                "samples": len(indices),
                "log_bayes_factor_hairy_vs_gr": log_ratio_for_indices(log_num, log_den, indices),
            }
        )
    values = np.asarray([float(row["log_bayes_factor_hairy_vs_gr"]) for row in rows])
    summary = {
        "branch": branch,
        "population_samples": len(log_num),
        "subsample_size": len(selected),
        "repeats": repeats,
        "headline_subsample_log_bayes_factor": log_ratio_for_indices(log_num, log_den, selected),
        "full_empirical_population_log_bayes_factor": log_ratio_for_indices(
            log_num, log_den, np.arange(len(log_num))
        ),
        "subsample_mean": float(np.mean(values)),
        "subsample_standard_deviation": float(np.std(values, ddof=1)),
        "subsample_q05": float(np.quantile(values, 0.05)),
        "subsample_q95": float(np.quantile(values, 0.95)),
    }
    return rows, summary


def convention_row(
    branch: str,
    baseline_surface: np.ndarray,
    baseline_gr: float,
    alternative_surface: np.ndarray,
    alternative_gr: float,
    fine_log_bayes: float,
    h0s: np.ndarray,
    alphas: np.ndarray,
) -> tuple[dict[str, object], dict[str, np.ndarray]]:
    base, base_post, _, _ = summarize_surface(baseline_surface, baseline_gr, h0s, alphas)
    alt, alt_post, _, _ = summarize_surface(alternative_surface, alternative_gr, h0s, alphas)
    base_cond = conditional_limits(baseline_surface, h0s, alphas, branch, "pre_minus40M_NRSur3dq8", 6.0)
    alt_cond = conditional_limits(alternative_surface, h0s, alphas, branch, "pre_minus40M_NRSur3dq8", 6.0)
    base_h05 = next(row for row in base_cond if np.isclose(float(row["h0"]), 0.5))
    alt_h05 = next(row for row in alt_cond if np.isclose(float(row["h0"]), 0.5))
    base_relative = baseline_surface - baseline_gr
    alt_relative = alternative_surface - alternative_gr
    return {
        "branch": branch,
        "coarse_baseline_log_bayes_factor": float(base["log_bayes_factor_hairy_vs_gr"]),
        "alternative_convention_log_bayes_factor": float(alt["log_bayes_factor_hairy_vs_gr"]),
        "convention_delta_log_bayes_factor": float(alt["log_bayes_factor_hairy_vs_gr"] - base["log_bayes_factor_hairy_vs_gr"]),
        "coarse_baseline_alpha_q90": float(base["alpha_q90"]),
        "alternative_convention_alpha_q90": float(alt["alpha_q90"]),
        "coarse_baseline_alpha_q90_at_h0_0p5": float(base_h05["alpha_q90_upper_conditional"]),
        "alternative_convention_alpha_q90_at_h0_0p5": float(alt_h05["alpha_q90_upper_conditional"]),
        "maximum_abs_change_in_relative_log_likelihood": float(np.max(np.abs(alt_relative - base_relative))),
        "posterior_total_variation_distance": float(0.5 * np.sum(np.abs(alt_post - base_post))),
        "fine_grid_log_bayes_factor": float(fine_log_bayes),
    }, {
        "baseline_log_likelihood": baseline_surface,
        "alternative_log_likelihood": alternative_surface,
        "baseline_posterior": base_post,
        "alternative_posterior": alt_post,
    }


def fit_tail_refits(
    branch: str,
    samples: np.ndarray,
    components: int,
    fit_samples: int,
    repeats: int,
    tail_probability: float,
    seed: int,
    modes: tuple[np.ndarray, np.ndarray],
    prior_mass: np.ndarray,
    h0s: np.ndarray,
    alphas: np.ndarray,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for repeat in range(repeats):
        rng = np.random.default_rng(seed + 10007 * repeat)
        pool = samples
        if len(pool) > fit_samples:
            pool = pool[rng.choice(len(pool), fit_samples, replace=False)]
        order = rng.permutation(len(pool))
        split = int(0.8 * len(pool))
        train = pool[order[:split]]
        test = pool[order[split:]]
        center = np.mean(train, axis=0)
        scale = np.std(train, axis=0, ddof=1)
        train_z = (train - center) / scale
        test_z = (test - center) / scale
        model = GaussianMixture(
            n_components=components,
            covariance_type="full",
            reg_covar=1.0e-6,
            max_iter=600,
            n_init=2,
            random_state=seed + 10007 * repeat + components,
        )
        model.fit(train_z)
        generated_z, _ = model.sample(len(test_z))
        radius_test = np.linalg.norm(test_z, axis=1)
        radius_generated = np.linalg.norm(generated_z, axis=1)
        threshold = float(np.quantile(radius_test, 1.0 - tail_probability))
        marginal_ks = [float(ks_2samp(test_z[:, index], generated_z[:, index]).statistic) for index in range(test_z.shape[1])]
        radial_ks = float(ks_2samp(radius_test, radius_generated).statistic)
        likelihood = {"model": model, "center": center, "scale": scale}
        if branch == "pyRing":
            surface, gr_log, valid_min = evaluate_pyring_surface(
                likelihood, modes[0], modes[1], prior_mass
            )
        else:
            surface, gr_log, valid_min = evaluate_ringdown_deviation_surface(
                likelihood, modes[0], modes[1], prior_mass
            )
        summary, _, _, _ = summarize_surface(surface, gr_log, h0s, alphas)
        score_test = model.score_samples(test_z)
        tail_mask = radius_test >= threshold
        rows.append(
            {
                "branch": branch,
                "repeat": repeat,
                "components": components,
                "training_samples": len(train),
                "heldout_samples": len(test),
                "converged": bool(model.converged_),
                "heldout_mean_log_density_standardized": float(np.mean(score_test)),
                "heldout_tail_mean_log_density_standardized": float(np.mean(score_test[tail_mask])),
                "heldout_radial_tail_probability": tail_probability,
                "generated_fraction_beyond_heldout_radial_tail": float(np.mean(radius_generated >= threshold)),
                "maximum_marginal_ks_distance": max(marginal_ks),
                "radial_ks_distance": radial_ks,
                "minimum_valid_fraction": valid_min,
                "log_bayes_factor_hairy_vs_gr": float(summary["log_bayes_factor_hairy_vs_gr"]),
                "alpha_q90": float(summary["alpha_q90"]),
            }
        )
    return rows


def validity_audit(spins: np.ndarray, alphas: np.ndarray, h0s: np.ndarray) -> tuple[list[dict[str, object]], dict[str, object]]:
    rows: list[dict[str, object]] = []
    for h0 in h0s:
        denominator = 1.0 - 0.5 * float(h0)
        for alpha in alphas:
            for spin in spins:
                r_plus, r_minus = approximate_horizons(float(spin), float(alpha), float(h0))
                expansion = 0.0 if denominator <= 0.0 else float(alpha) * math.exp(-r_plus / denominator)
                mass_effective = 0.5 * (r_plus + r_minus)
                charge_squared = max(0.0, r_plus * r_minus - float(spin) ** 2)
                charge_ratio = math.sqrt(charge_squared) / mass_effective
                rows.append(
                    {
                        "h0": float(h0),
                        "alpha_hair": float(alpha),
                        "spin": float(spin),
                        "approximate_r_plus": r_plus,
                        "approximate_r_minus": r_minus,
                        "outer_horizon_expansion_parameter": expansion,
                        "effective_charge_to_mass_ratio": charge_ratio,
                    }
                )
    reference = min(
        rows,
        key=lambda row: abs(float(row["h0"]) - 1.0)
        + abs(float(row["alpha_hair"]) - 0.5)
        + abs(float(row["spin"]) - 0.7),
    )
    event_rows = [row for row in rows if 0.6 <= float(row["spin"]) <= 0.75]
    return rows, {
        "condition_from_Li": "alpha*exp[-r_plus/(1-h0/2)] much less than 1 and Q_effect/M_effect small",
        "maximum_outer_horizon_expansion_parameter_full_grid": max(float(row["outer_horizon_expansion_parameter"]) for row in rows),
        "maximum_effective_charge_to_mass_ratio_full_grid": max(float(row["effective_charge_to_mass_ratio"]) for row in rows),
        "maximum_outer_horizon_expansion_parameter_spin_0p6_to_0p75": max(float(row["outer_horizon_expansion_parameter"]) for row in event_rows),
        "maximum_effective_charge_to_mass_ratio_spin_0p6_to_0p75": max(float(row["effective_charge_to_mass_ratio"]) for row in event_rows),
        "reference_h0_1_alpha_0p5_spin_0p7": reference,
    }


def save_convention_figure(path_base: Path, rows: list[dict[str, object]]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.0), sharey=True, constrained_layout=True)
    for ax, branch in zip(axes, ("pyRing", "RINGDOWN_deviation")):
        for convention, color, linestyle, label in (
            ("table_calibrated", "tab:blue", "-", "table-calibrated"),
            ("ode_consistent", "tab:orange", "--", "ODE-consistent"),
        ):
            selected = sorted(
                [row for row in rows if row["branch"] == branch and row["convention"] == convention],
                key=lambda row: float(row["h0"]),
            )
            ax.plot(
                [float(row["h0"]) for row in selected],
                [float(row["alpha_q90_upper_conditional"]) for row in selected],
                marker="o",
                ms=3,
                color=color,
                linestyle=linestyle,
                label=label,
            )
        ax.axhline(0.45, color="black", ls=":", lw=1.0, label="uniform-prior value")
        ax.set_xlabel(r"fixed $h_0/M$")
        ax.set_title("pyRing" if branch == "pyRing" else "RINGDOWN")
        ax.grid(alpha=0.25)
    axes[0].set_ylabel(r"conditional 90% upper quantile of $\alpha$")
    axes[0].legend(fontsize=8)
    fig.savefig(path_base.with_suffix(".png"), dpi=240)
    fig.savefig(path_base.with_suffix(".pdf"))
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--grid", type=Path, default=DEFAULT_GRID)
    parser.add_argument("--systematics", type=Path, default=DEFAULT_SYSTEMATICS)
    parser.add_argument("--pre-prior", type=Path, default=DEFAULT_PRE_PRIOR)
    parser.add_argument("--pyring-posterior", type=Path, default=DEFAULT_PYRING_POSTERIOR)
    parser.add_argument("--pyring-gmm", type=Path, default=DEFAULT_PYRING_GMM)
    parser.add_argument("--ringdown-posterior", type=Path, default=DEFAULT_RINGDOWN_POSTERIOR)
    parser.add_argument("--ringdown-gmm", type=Path, default=DEFAULT_RINGDOWN_GMM)
    parser.add_argument("--ringdown-diagnostics", type=Path, default=DEFAULT_RINGDOWN_DIAGNOSTICS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--reuse-gmm-refits",
        action="store_true",
        help="Reuse the existing split/refit table while regenerating other outputs.",
    )
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    config_all = json.loads(args.config.read_text(encoding="utf-8"))
    config = config_all["inference"]
    robustness = config_all["referee_robustness"]
    production = np.load(args.grid)
    spins = np.asarray(production["spins"], dtype=float)
    alphas = np.asarray(production["alphas"], dtype=float)
    h0s = np.asarray(production["h0s"], dtype=float)
    overtones = np.asarray(production["overtones"], dtype=int)
    omega = np.asarray(production["omega"], dtype=np.complex128)
    i220 = int(np.where(overtones == 0)[0][0])
    i221 = int(np.where(overtones == 1)[0][0])

    mass_all, spin_all, selected = read_primary_remnant_samples(
        args.pre_prior, float(spins[0]), float(spins[-1]), int(config["prior_samples"]), int(config["seed"])
    )
    production_modes_all = interpolate_modes_at_prior_spins(omega, spins, spin_all, (i220, i221))
    production_modes_selected = (
        production_modes_all[0][:, :, selected],
        production_modes_all[1][:, :, selected],
    )
    mass_selected = mass_all[selected]

    py_bundle = joblib.load(args.pyring_gmm)
    py_components = max(int(value) for value in config["pyRing_gmm_components"])
    py_likelihood = {"model": py_bundle["models"][py_components], "center": py_bundle["center"], "scale": py_bundle["scale"]}
    direct_components = direct_preferred_components(args.ringdown_diagnostics)
    direct_models = joblib.load(args.ringdown_gmm)
    direct_likelihood = direct_models[direct_components]

    mc_rows: list[dict[str, object]] = []
    mc_summary: list[dict[str, object]] = []
    for offset, (branch, likelihood, calculator) in enumerate(
        (
            ("pyRing", py_likelihood, py_log_terms),
            ("RINGDOWN_deviation", direct_likelihood, direct_log_terms),
        )
    ):
        log_num, log_den = calculator(
            likelihood, production_modes_all[0], production_modes_all[1], mass_all, h0s, alphas
        )
        rows, summary = subsampling_rows(
            branch,
            log_num,
            log_den,
            selected,
            int(robustness["remnant_subsample_repeats"]),
            int(config["seed"]) + 70000 + offset,
        )
        mc_rows.extend(rows)
        mc_summary.append(summary)
    write_csv(args.output / "remnant_subsampling_log_bayes.csv", mc_rows)

    internal = np.load(args.systematics)
    coarse_spins = np.asarray(internal["spins"], dtype=float)
    coarse_alphas = np.asarray(internal["alphas"], dtype=float)
    coarse_h0s = np.asarray(internal["h0s"], dtype=float)
    coarse_overtones = np.asarray(internal["overtones"], dtype=int)
    coarse_i220 = int(np.where(coarse_overtones == 0)[0][0])
    coarse_i221 = int(np.where(coarse_overtones == 1)[0][0])
    base_modes = interpolate_modes_at_prior_spins(
        np.asarray(internal["baseline"]), coarse_spins, spin_all[selected], (coarse_i220, coarse_i221)
    )
    alternative_modes = interpolate_modes_at_prior_spins(
        np.asarray(internal["ode_approximate"]), coarse_spins, spin_all[selected], (coarse_i220, coarse_i221)
    )

    convention_rows: list[dict[str, object]] = []
    convention_conditional_rows: list[dict[str, object]] = []
    convention_arrays: dict[str, np.ndarray] = {
        "h0s": coarse_h0s,
        "alphas": coarse_alphas,
    }
    for branch, likelihood, evaluator in (
        ("pyRing", py_likelihood, evaluate_pyring_surface),
        ("RINGDOWN_deviation", direct_likelihood, evaluate_ringdown_deviation_surface),
    ):
        base_surface, base_gr, _ = evaluator(likelihood, base_modes[0], base_modes[1], mass_selected)
        alt_surface, alt_gr, _ = evaluator(likelihood, alternative_modes[0], alternative_modes[1], mass_selected)
        fine_surface, fine_gr, _ = evaluator(
            likelihood, production_modes_selected[0], production_modes_selected[1], mass_selected
        )
        fine_summary, _, _, _ = summarize_surface(fine_surface, fine_gr, h0s, alphas)
        row, arrays = convention_row(
            branch,
            base_surface,
            base_gr,
            alt_surface,
            alt_gr,
            float(fine_summary["log_bayes_factor_hairy_vs_gr"]),
            coarse_h0s,
            coarse_alphas,
        )
        convention_rows.append(row)
        for convention, surface in (
            ("table_calibrated", base_surface),
            ("ode_consistent", alt_surface),
        ):
            for conditional in conditional_limits(
                surface, coarse_h0s, coarse_alphas, branch, "pre_minus40M_NRSur3dq8", 6.0
            ):
                convention_conditional_rows.append(
                    {
                        "branch": branch,
                        "convention": convention,
                        "h0": conditional["h0"],
                        "alpha_median_conditional": conditional["alpha_median_conditional"],
                        "alpha_q90_upper_conditional": conditional["alpha_q90_upper_conditional"],
                        "alpha_q95_upper_conditional": conditional["alpha_q95_upper_conditional"],
                    }
                )
        for name, value in arrays.items():
            convention_arrays[f"{branch}_{name}"] = value
    write_csv(args.output / "horizon_index_event_robustness.csv", convention_rows)
    write_csv(args.output / "horizon_index_conditional_limits.csv", convention_conditional_rows)
    np.savez_compressed(args.output / "horizon_index_event_surfaces.npz", **convention_arrays)
    save_convention_figure(args.output / "hairy_horizon_convention_sensitivity", convention_conditional_rows)

    py_posterior = np.genfromtxt(args.pyring_posterior, names=True)
    py_samples = np.column_stack(
        (py_posterior["Mf"], py_posterior["af"], py_posterior["domega_221"], py_posterior["dtau_221"])
    )
    with h5py.File(args.ringdown_posterior, "r") as handle:
        direct_samples = np.column_stack(
            (
                np.log(np.asarray(handle["f_220"], dtype=float)),
                np.log(np.asarray(handle["f_221"], dtype=float)),
                np.asarray(handle["df_221"], dtype=float),
            )
        )
    gmm_path = args.output / "gmm_tail_refit_diagnostics.csv"
    if args.reuse_gmm_refits and gmm_path.exists():
        with gmm_path.open(encoding="utf-8", newline="") as handle:
            gmm_rows = list(csv.DictReader(handle))
    else:
        gmm_rows: list[dict[str, object]] = []
        gmm_rows.extend(fit_tail_refits(
            "pyRing",
            py_samples,
            py_components,
            int(robustness["gmm_fit_samples"]),
            int(robustness["gmm_refit_repeats"]),
            float(robustness["tail_probability"]),
            int(config["seed"]) + 80000,
            base_modes,
            mass_selected,
            coarse_h0s,
            coarse_alphas,
        ))
        gmm_rows.extend(fit_tail_refits(
            "RINGDOWN_deviation",
            direct_samples,
            direct_components,
            int(robustness["gmm_fit_samples"]),
            int(robustness["gmm_refit_repeats"]),
            float(robustness["tail_probability"]),
            int(config["seed"]) + 90000,
            base_modes,
            mass_selected,
            coarse_h0s,
            coarse_alphas,
        ))
        write_csv(gmm_path, gmm_rows)
    gmm_summary: list[dict[str, object]] = []
    for branch in ("pyRing", "RINGDOWN_deviation"):
        rows = [row for row in gmm_rows if row["branch"] == branch]
        log_bayes = np.asarray([float(row["log_bayes_factor_hairy_vs_gr"]) for row in rows])
        gmm_summary.append(
            {
                "branch": branch,
                "refits": len(rows),
                "all_converged": all(str(row["converged"]).lower() == "true" for row in rows),
                "log_bayes_mean": float(np.mean(log_bayes)),
                "log_bayes_standard_deviation": float(np.std(log_bayes, ddof=1)),
                "log_bayes_min": float(np.min(log_bayes)),
                "log_bayes_max": float(np.max(log_bayes)),
                "generated_radial_tail_fraction_min": min(float(row["generated_fraction_beyond_heldout_radial_tail"]) for row in rows),
                "generated_radial_tail_fraction_max": max(float(row["generated_fraction_beyond_heldout_radial_tail"]) for row in rows),
                "maximum_marginal_ks_distance": max(float(row["maximum_marginal_ks_distance"]) for row in rows),
                "maximum_radial_ks_distance": max(float(row["radial_ks_distance"]) for row in rows),
            }
        )

    validity_rows, validity_summary = validity_audit(spins, alphas, h0s)
    write_csv(args.output / "li_approximation_validity_grid.csv", validity_rows)
    summary = {
        "primary_remnant_prior": "pre_minus40M_NRSur3dq8",
        "mc_integration": mc_summary,
        "horizon_index_convention": convention_rows,
        "gmm_tail_refits": gmm_summary,
        "li_approximation_validity": validity_summary,
        "interpretation": {
            "mc": "random-subset integration variability over the finite pre-peak NRSur3dq8 remnant sample",
            "convention": "deterministic alternative recurrence, not a statistical physical-error distribution",
            "gmm": "held-out radial and marginal diagnostics plus event-level split/refit sensitivity",
        },
    }
    (args.output / "referee_robustness_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
