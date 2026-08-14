"""Conditional GW250114 Bayesian benchmark for the rotating-hairy grid.

Four internally computable implementations are propagated through the same
public pyRing spectral likelihood.  The comparison deliberately remains a
conditional benchmark because none of these variants quantifies the missing
full gravitational perturbation sector or the Dudley--Finley structural error.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
LOCAL_DEPS = ROOT / ".deps" / "python"
if LOCAL_DEPS.exists():
    sys.path.insert(0, str(LOCAL_DEPS))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import joblib
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.special import logsumexp

from gw250114_bayesian_spectral_eft import interp_complex, kerr_tables, log_density, weighted_quantile


DEFAULT_GMM = ROOT / "results" / "gw250114_bayesian_spectral_eft" / "spectral_likelihood_gmms.joblib"
DEFAULT_SYSTEMATICS = ROOT / "results" / "hairy_qnm_internal_systematics" / "internal_systematics.npz"
DEFAULT_PRIOR = ROOT / "results" / "article1_independent_prior" / "pre_minus40M_remnant_model_crosscheck_samples.csv"
DEFAULT_OUTPUT = ROOT / "results" / "gw250114_hairy_bayesian_benchmark"


def quadrature_weights(axis: np.ndarray) -> np.ndarray:
    """Normalized trapezoid weights on a possibly irregular axis."""
    widths = np.diff(axis)
    weights = np.empty(len(axis), dtype=float)
    weights[0] = 0.5 * widths[0]
    weights[-1] = 0.5 * widths[-1]
    weights[1:-1] = 0.5 * (widths[:-1] + widths[1:])
    return weights / np.sum(weights)


def pseudo_observables_from_modes(
    mass_true: np.ndarray,
    deformed_220: np.ndarray,
    deformed_221: np.ndarray,
    kerr: dict[str, np.ndarray],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    quality_grid = kerr["omega_220"].real / (-kerr["omega_220"].imag)
    quality_deformed = deformed_220.real / (-deformed_220.imag)
    pseudo_spin = np.interp(quality_deformed, quality_grid, kerr["spin"], left=np.nan, right=np.nan)
    pseudo_220 = interp_complex(pseudo_spin, kerr["spin"], kerr["omega_220"])
    pseudo_221 = interp_complex(pseudo_spin, kerr["spin"], kerr["omega_221"])
    mass_ratio_fit_to_true = pseudo_220.real / deformed_220.real
    pseudo_mass = mass_true * mass_ratio_fit_to_true
    domega = deformed_221.real / pseudo_221.real * mass_ratio_fit_to_true - 1.0
    dtau = (1.0 / mass_ratio_fit_to_true) * ((-pseudo_221.imag) / (-deformed_221.imag)) - 1.0
    valid = (
        np.isfinite(pseudo_mass)
        & np.isfinite(pseudo_spin)
        & np.isfinite(domega)
        & np.isfinite(dtau)
        & (pseudo_spin > 0.0)
        & (pseudo_spin < 0.99)
        & (domega > -1.0)
        & (domega < 1.0)
        & (dtau > -0.9)
        & (dtau < 1.0)
    )
    return pseudo_mass, pseudo_spin, domega, dtau, valid


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gmm", type=Path, default=DEFAULT_GMM)
    parser.add_argument("--systematics", type=Path, default=DEFAULT_SYSTEMATICS)
    parser.add_argument("--prior", type=Path, default=DEFAULT_PRIOR)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--prior-samples", type=int, default=4096)
    parser.add_argument("--seed", type=int, default=25011405)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    likelihood = joblib.load(args.gmm)
    models = likelihood["models"]
    center = np.asarray(likelihood["center"])
    scale = np.asarray(likelihood["scale"])
    grid = np.load(args.systematics)
    spins = grid["spins"]
    alphas = grid["alphas"]
    h0s = grid["h0s"]
    overtones = grid["overtones"]
    index_220 = int(np.where(overtones == 0)[0][0])
    index_221 = int(np.where(overtones == 1)[0][0])
    variants = ["baseline", "ode_approximate", "table_exact", "ode_exact"]

    prior_rows = np.genfromtxt(args.prior, names=True, delimiter=",")
    mass_all = np.concatenate((prior_rows["final_mass_BMR2012"], prior_rows["final_mass_NRSur3dq8"]))
    spin_all = np.concatenate((prior_rows["final_spin_HBR2016"], prior_rows["final_spin_NRSur3dq8_projected"]))
    finite = np.isfinite(mass_all) & np.isfinite(spin_all) & (spin_all >= spins[0]) & (spin_all <= spins[-1])
    mass_all = mass_all[finite]
    spin_all = spin_all[finite]
    rng = np.random.default_rng(args.seed)
    selected = rng.choice(len(mass_all), size=min(args.prior_samples, len(mass_all)), replace=False)
    prior_mass = mass_all[selected]
    prior_spin = spin_all[selected]
    kerr = kerr_tables()

    alpha_weights = quadrature_weights(alphas)
    h0_weights = quadrature_weights(h0s)
    parameter_weights = h0_weights[:, None] * alpha_weights[None, :]
    rows: list[dict[str, object]] = []
    fixed_rows: list[dict[str, object]] = []

    for component_count, model in sorted(models.items()):
        gr_points = np.column_stack((prior_mass, prior_spin, np.zeros(len(prior_mass)), np.zeros(len(prior_mass))))
        gr_log = log_density(model, gr_points, center, scale)
        log_gr_evidence = float(logsumexp(gr_log) - math.log(len(gr_log)))
        for variant in variants:
            values = grid[variant]
            log_parameter = np.full((len(h0s), len(alphas)), -np.inf)
            fixed_lookup: dict[tuple[float, float], float] = {}
            valid_fractions: list[float] = []
            for ih, h0 in enumerate(h0s):
                for ia, alpha_hair in enumerate(alphas):
                    spline_220 = CubicSpline(spins, values[ih, ia, :, index_220])
                    spline_221 = CubicSpline(spins, values[ih, ia, :, index_221])
                    deformed_220 = np.asarray(spline_220(prior_spin), dtype=np.complex128)
                    deformed_221 = np.asarray(spline_221(prior_spin), dtype=np.complex128)
                    mass, spin, domega, dtau, valid = pseudo_observables_from_modes(
                        prior_mass, deformed_220, deformed_221, kerr
                    )
                    valid_fractions.append(float(np.mean(valid)))
                    if np.any(valid):
                        points = np.column_stack((mass[valid], spin[valid], domega[valid], dtau[valid]))
                        density = log_density(model, points, center, scale)
                        log_parameter[ih, ia] = logsumexp(density) - math.log(len(prior_mass))
                    fixed_lookup[(float(h0), float(alpha_hair))] = float(log_parameter[ih, ia])

            log_evidence = float(logsumexp(log_parameter + np.log(parameter_weights)))
            posterior_weights = np.exp(
                log_parameter + np.log(parameter_weights) - logsumexp(log_parameter + np.log(parameter_weights))
            )
            alpha_posterior = np.sum(posterior_weights, axis=0)
            h0_posterior = np.sum(posterior_weights, axis=1)
            alpha_q05, alpha_median, alpha_q95 = weighted_quantile(
                alphas, alpha_posterior, [0.05, 0.50, 0.95]
            )
            h0_q05, h0_median, h0_q95 = weighted_quantile(
                h0s, h0_posterior, [0.05, 0.50, 0.95]
            )
            rows.append(
                {
                    "gmm_components": int(component_count),
                    "variant": variant,
                    "alpha_prior_min": float(alphas[0]),
                    "alpha_prior_max": float(alphas[-1]),
                    "h0_prior_min": float(h0s[0]),
                    "h0_prior_max": float(h0s[-1]),
                    "log_bayes_factor_hairy_vs_gr": log_evidence - log_gr_evidence,
                    "bayes_factor_hairy_vs_gr": float(np.exp(log_evidence - log_gr_evidence)),
                    "alpha_q05": alpha_q05,
                    "alpha_median": alpha_median,
                    "alpha_q95": alpha_q95,
                    "h0_q05": h0_q05,
                    "h0_median": h0_median,
                    "h0_q95": h0_q95,
                    "minimum_valid_fraction": min(valid_fractions),
                    "maximum_log_likelihood_ratio_over_gr_evidence": float(
                        np.max(log_parameter) - log_gr_evidence
                    ),
                }
            )
            for fixed_alpha in (0.10, 0.25, 0.50):
                fixed_log = fixed_lookup[(1.0, fixed_alpha)]
                fixed_rows.append(
                    {
                        "gmm_components": int(component_count),
                        "variant": variant,
                        "alpha_hair": fixed_alpha,
                        "h0": 1.0,
                        "log_likelihood_ratio_fixed_hairy_vs_gr_marginalized_remnant": fixed_log
                        - log_gr_evidence,
                        "likelihood_ratio_fixed_hairy_vs_gr_marginalized_remnant": float(
                            np.exp(fixed_log - log_gr_evidence)
                        ),
                    }
                )

    def write_csv(path: Path, data: list[dict[str, object]]) -> None:
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(data[0]))
            writer.writeheader()
            writer.writerows(data)

    write_csv(args.output / "hairy_model_evidence.csv", rows)
    write_csv(args.output / "hairy_fixed_benchmarks.csv", fixed_rows)
    preferred = [row for row in rows if int(row["gmm_components"]) == max(models)]
    summary = {
        "event": "GW250114",
        "posterior_product": "pyRing_220_221_delta_weighted_posterior",
        "independent_remnant_prior_samples_used": len(prior_mass),
        "hair_prior": {"alpha_hair": [float(alphas[0]), float(alphas[-1])], "h0": [float(h0s[0]), float(h0s[-1])]},
        "variants": variants,
        "conditional_only": True,
        "exclusion_ready": False,
        "blocking_theory_error": "unquantified_Dudley_Finley_and_missing_full_perturbation_sector",
        "preferred_gmm_results": preferred,
    }
    (args.output / "benchmark_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps({key: value for key, value in summary.items() if key != "preferred_gmm_results"}, indent=2))
    print(f"evidence_rows={len(rows)}, fixed_rows={len(fixed_rows)}")


if __name__ == "__main__":
    main()
