"""Public-posterior Bayesian EFT calculation for GW250114.

The pyRing weighted posterior is converted into a smooth marginalized
spectral likelihood in (Mf, af, domega_221, dtau_221).  The conversion is
valid here because the four native priors are separable and uniform over the
interior sampled by the calculation; their constant density cancels in all
reported Bayes factors.  This remains a posterior-product likelihood, not a
new strain-level analysis.

The EFT mapping is nonlinear and accounts for the fact that pyRing defines
Mf and af through a Kerr 220 mode.  A deformed 220 is therefore first mapped
to its pseudo-Kerr remnant before the predicted 221 deviations are evaluated.
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

import joblib
import numpy as np
import qnm
from scipy.special import logsumexp
from sklearn.mixture import GaussianMixture


DEFAULT_POSTERIOR = (
    ROOT
    / "data"
    / "raw"
    / "GW250114_data_release"
    / "data"
    / "pyring_220_221_delta_posterior_without_frequencies.dat"
)
DEFAULT_PRIOR = (
    ROOT
    / "results"
    / "article1_independent_prior"
    / "pre_minus40M_remnant_model_crosscheck_samples.csv"
)
DEFAULT_COEFFICIENTS = ROOT / "data" / "beyond_kerr_qnm_selected_fits.csv"
DEFAULT_OUTPUT = ROOT / "results" / "gw250114_bayesian_spectral_eft"


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def read_coefficients(path: Path) -> dict[tuple[str, str, str], np.ndarray]:
    grouped: dict[tuple[str, str, str], list[tuple[int, complex]]] = {}
    for row in csv.DictReader(path.open(encoding="utf-8", newline="")):
        if row["mode"] not in {"220", "221"}:
            continue
        key = (row["operator"], row["branch"], row["mode"])
        grouped.setdefault(key, []).append(
            (int(row["k"]), complex(float(row["coefficient_re"]), float(row["coefficient_im"])))
        )
    result: dict[tuple[str, str, str], np.ndarray] = {}
    for key, values in grouped.items():
        maximum = max(index for index, _ in values)
        coefficients = np.zeros(maximum + 1, dtype=np.complex128)
        for index, value in values:
            coefficients[index] = value
        result[key] = coefficients
    return result


def polynomial(coefficients: np.ndarray, spin: np.ndarray) -> np.ndarray:
    return np.polynomial.polynomial.polyval(spin, coefficients)


def kerr_tables() -> dict[str, np.ndarray]:
    spins = np.linspace(0.20, 0.98, 1561)
    result: dict[str, np.ndarray] = {"spin": spins}
    for overtone in (0, 1):
        sequence = qnm.modes_cache(s=-2, l=2, m=2, n=overtone)
        result[f"omega_22{overtone}"] = np.asarray(
            [complex(sequence(a=float(spin))[0]) for spin in spins], dtype=np.complex128
        )
    return result


def interp_complex(spin: np.ndarray, grid: np.ndarray, values: np.ndarray) -> np.ndarray:
    return np.interp(spin, grid, values.real) + 1j * np.interp(spin, grid, values.imag)


def pseudo_kerr_observables(
    mass_true: np.ndarray,
    spin_true: np.ndarray,
    alpha: float,
    delta_220_coefficients: np.ndarray,
    delta_221_coefficients: np.ndarray,
    kerr: dict[str, np.ndarray],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    spin_grid = kerr["spin"]
    kerr_220 = interp_complex(spin_true, spin_grid, kerr["omega_220"])
    kerr_221 = interp_complex(spin_true, spin_grid, kerr["omega_221"])
    deformed_220 = kerr_220 + alpha * polynomial(delta_220_coefficients, spin_true)
    deformed_221 = kerr_221 + alpha * polynomial(delta_221_coefficients, spin_true)

    quality_grid = kerr["omega_220"].real / (-kerr["omega_220"].imag)
    quality_deformed = deformed_220.real / (-deformed_220.imag)
    pseudo_spin = np.interp(quality_deformed, quality_grid, spin_grid, left=np.nan, right=np.nan)
    pseudo_220 = interp_complex(pseudo_spin, spin_grid, kerr["omega_220"])
    pseudo_221 = interp_complex(pseudo_spin, spin_grid, kerr["omega_221"])
    mass_ratio_fit_to_true = pseudo_220.real / deformed_220.real
    pseudo_mass = mass_true * mass_ratio_fit_to_true
    domega_221 = deformed_221.real / pseudo_221.real * mass_ratio_fit_to_true - 1.0
    dtau_221 = (
        (1.0 / mass_ratio_fit_to_true)
        * ((-pseudo_221.imag) / (-deformed_221.imag))
        - 1.0
    )
    valid = (
        np.isfinite(pseudo_mass)
        & np.isfinite(pseudo_spin)
        & np.isfinite(domega_221)
        & np.isfinite(dtau_221)
        & (pseudo_spin > 0.0)
        & (pseudo_spin < 0.99)
        & (domega_221 > -1.0)
        & (domega_221 < 1.0)
        & (dtau_221 > -0.9)
        & (dtau_221 < 1.0)
    )
    return pseudo_mass, pseudo_spin, domega_221, dtau_221, valid


def weighted_quantile(grid: np.ndarray, weights: np.ndarray, probabilities: list[float]) -> list[float]:
    normalized = weights / np.sum(weights)
    cumulative = np.cumsum(normalized)
    return [float(np.interp(probability, cumulative, grid)) for probability in probabilities]


def fit_likelihood_models(samples: np.ndarray, components: list[int], seed: int) -> tuple[dict[int, GaussianMixture], list[dict[str, object]], np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    permutation = rng.permutation(len(samples))
    split = int(0.8 * len(samples))
    training = samples[permutation[:split]]
    test = samples[permutation[split:]]
    center = np.mean(training, axis=0)
    scale = np.std(training, axis=0, ddof=1)
    training_z = (training - center) / scale
    test_z = (test - center) / scale
    models: dict[int, GaussianMixture] = {}
    diagnostics: list[dict[str, object]] = []
    for count in components:
        model = GaussianMixture(
            n_components=count,
            covariance_type="full",
            reg_covar=1.0e-6,
            max_iter=500,
            n_init=2,
            random_state=seed + count,
        )
        model.fit(training_z)
        models[count] = model
        diagnostics.append(
            {
                "components": count,
                "converged": bool(model.converged_),
                "iterations": int(model.n_iter_),
                "train_mean_log_density": float(np.mean(model.score_samples(training_z))),
                "test_mean_log_density": float(np.mean(model.score_samples(test_z))),
                "bic_training": float(model.bic(training_z)),
                "aic_training": float(model.aic(training_z)),
            }
        )
    return models, diagnostics, center, scale


def log_density(model: GaussianMixture, values: np.ndarray, center: np.ndarray, scale: np.ndarray) -> np.ndarray:
    standardized = (values - center) / scale
    return model.score_samples(standardized) - float(np.sum(np.log(scale)))


def integrate_model(
    model: GaussianMixture,
    center: np.ndarray,
    scale: np.ndarray,
    prior_mass: np.ndarray,
    prior_spin: np.ndarray,
    alpha_grid: np.ndarray,
    coefficients_220: np.ndarray,
    coefficients_221: np.ndarray,
    kerr: dict[str, np.ndarray],
    gr_log_density: np.ndarray,
    bootstrap_seed: int,
) -> dict[str, object]:
    n_prior = len(prior_mass)
    n_alpha = len(alpha_grid)
    log_by_alpha = np.full(n_alpha, -np.inf)
    row_log_integrals = np.full((n_prior, n_alpha), -np.inf)
    valid_fractions = np.zeros(n_alpha)
    for index, alpha in enumerate(alpha_grid):
        mass, spin, domega, dtau, valid = pseudo_kerr_observables(
            prior_mass,
            prior_spin,
            float(alpha),
            coefficients_220,
            coefficients_221,
            kerr,
        )
        valid_fractions[index] = float(np.mean(valid))
        if np.any(valid):
            points = np.column_stack((mass[valid], spin[valid], domega[valid], dtau[valid]))
            density = log_density(model, points, center, scale)
            row_log_integrals[valid, index] = density
            log_by_alpha[index] = logsumexp(density) - math.log(n_prior)

    trapezoid = np.ones(n_alpha)
    trapezoid[[0, -1]] = 0.5
    trapezoid /= np.sum(trapezoid)
    log_evidence = float(logsumexp(log_by_alpha + np.log(trapezoid)))
    log_gr_evidence = float(logsumexp(gr_log_density) - math.log(n_prior))
    posterior_log_weights = log_by_alpha + np.log(trapezoid)
    posterior_weights = np.exp(posterior_log_weights - logsumexp(posterior_log_weights))
    q05, q16, median, q84, q95 = weighted_quantile(
        alpha_grid, posterior_weights, [0.05, 0.16, 0.50, 0.84, 0.95]
    )

    # Integrate alpha separately for every remnant-prior draw.  Paired
    # bootstrap resampling then preserves covariance with the GR evidence.
    row_log_evidence = logsumexp(row_log_integrals + np.log(trapezoid)[None, :], axis=1)
    common_offset = max(float(np.max(row_log_evidence)), float(np.max(gr_log_density)))
    numerator = np.exp(row_log_evidence - common_offset)
    denominator = np.exp(gr_log_density - common_offset)
    rng = np.random.default_rng(bootstrap_seed)
    boot_log_bayes: list[float] = []
    for _ in range(300):
        indices = rng.integers(0, n_prior, n_prior)
        boot_log_bayes.append(float(np.log(np.mean(numerator[indices]) / np.mean(denominator[indices]))))

    return {
        "alpha_prior_min": float(alpha_grid[0]),
        "alpha_prior_max": float(alpha_grid[-1]),
        "alpha_grid_points": n_alpha,
        "log_evidence_relative_constant": log_evidence,
        "log_gr_evidence_relative_constant": log_gr_evidence,
        "log_bayes_factor_eft_vs_gr": log_evidence - log_gr_evidence,
        "bayes_factor_eft_vs_gr": float(np.exp(np.clip(log_evidence - log_gr_evidence, -700, 700))),
        "bootstrap_log_bayes_std": float(np.std(boot_log_bayes, ddof=1)),
        "alpha_q05": q05,
        "alpha_q16": q16,
        "alpha_median": median,
        "alpha_q84": q84,
        "alpha_q95": q95,
        "posterior_probability_alpha_positive": float(np.sum(posterior_weights[alpha_grid > 0.0])),
        "minimum_valid_fraction": float(np.min(valid_fractions)),
        "maximum_log_likelihood_ratio_over_gr_evidence": float(np.max(log_by_alpha) - log_gr_evidence),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--posterior", type=Path, default=DEFAULT_POSTERIOR)
    parser.add_argument("--prior", type=Path, default=DEFAULT_PRIOR)
    parser.add_argument("--coefficients", type=Path, default=DEFAULT_COEFFICIENTS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--prior-samples", type=int, default=4096)
    parser.add_argument("--alpha-points", type=int, default=161)
    parser.add_argument("--seed", type=int, default=250114)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    posterior = np.genfromtxt(args.posterior, names=True)
    spectral_samples = np.column_stack(
        (posterior["Mf"], posterior["af"], posterior["domega_221"], posterior["dtau_221"])
    )
    component_counts = [8, 16, 32, 64]
    models, diagnostics, center, scale = fit_likelihood_models(
        spectral_samples, component_counts, args.seed
    )
    write_csv(args.output / "gmm_likelihood_diagnostics.csv", diagnostics)
    joblib.dump(
        {"models": models, "center": center, "scale": scale, "variables": ["Mf", "af", "domega_221", "dtau_221"]},
        args.output / "spectral_likelihood_gmms.joblib",
    )

    prior_rows = np.genfromtxt(args.prior, names=True, delimiter=",")
    mass_all = np.concatenate((prior_rows["final_mass_BMR2012"], prior_rows["final_mass_NRSur3dq8"]))
    spin_all = np.concatenate((prior_rows["final_spin_HBR2016"], prior_rows["final_spin_NRSur3dq8_projected"]))
    finite = np.isfinite(mass_all) & np.isfinite(spin_all) & (spin_all > 0.20) & (spin_all < 0.98)
    mass_all = mass_all[finite]
    spin_all = spin_all[finite]
    rng = np.random.default_rng(args.seed)
    selected = rng.choice(len(mass_all), size=min(args.prior_samples, len(mass_all)), replace=False)
    prior_mass = mass_all[selected]
    prior_spin = spin_all[selected]

    coefficients = read_coefficients(args.coefficients)
    kerr = kerr_tables()
    model_keys = sorted({(operator, branch) for operator, branch, mode in coefficients if mode == "220"})
    evidence_rows: list[dict[str, object]] = []

    for component_count in component_counts:
        model = models[component_count]
        gr_points = np.column_stack((prior_mass, prior_spin, np.zeros(len(prior_mass)), np.zeros(len(prior_mass))))
        gr_log_density = log_density(model, gr_points, center, scale)
        for model_index, (operator, branch) in enumerate(model_keys):
            coeff_220 = coefficients[(operator, branch, "220")]
            coeff_221 = coefficients[(operator, branch, "221")]
            spin_probe = np.linspace(float(np.quantile(prior_spin, 0.005)), float(np.quantile(prior_spin, 0.995)), 401)
            response_220 = polynomial(coeff_220, spin_probe) / interp_complex(
                spin_probe, kerr["spin"], kerr["omega_220"]
            )
            response_221 = polynomial(coeff_221, spin_probe) / interp_complex(
                spin_probe, kerr["spin"], kerr["omega_221"]
            )
            maximum_response = float(max(np.max(np.abs(response_220)), np.max(np.abs(response_221))))
            prior_widths = {
                "spectral_shift_5pct": min(0.5, 0.05 / maximum_response),
                "spectral_shift_10pct": min(0.5, 0.10 / maximum_response),
                "fixed_abs_alpha_0p1_stress": 0.10,
            }
            for prior_index, (prior_name, half_width) in enumerate(prior_widths.items()):
                alpha_grid = np.linspace(-half_width, half_width, args.alpha_points)
                result = integrate_model(
                    model,
                    center,
                    scale,
                    prior_mass,
                    prior_spin,
                    alpha_grid,
                    coeff_220,
                    coeff_221,
                    kerr,
                    gr_log_density,
                    args.seed + 10000 * component_count + 100 * model_index + prior_index,
                )
                evidence_rows.append(
                    {
                        "gmm_components": component_count,
                        "operator": operator,
                        "branch": branch,
                        "coupling_mass_power": 4 if operator.startswith("lambda") else 6,
                        "prior_name": prior_name,
                        "maximum_abs_complex_response_per_alpha": maximum_response,
                        **result,
                    }
                )

    write_csv(args.output / "eft_bayesian_evidence.csv", evidence_rows)
    preferred_components = int(
        max(diagnostics, key=lambda row: float(row["test_mean_log_density"]))["components"]
    )
    preferred = [row for row in evidence_rows if int(row["gmm_components"]) == preferred_components]
    summary = {
        "posterior_samples": int(len(spectral_samples)),
        "independent_remnant_prior_samples_total": int(len(mass_all)),
        "independent_remnant_prior_samples_used": int(len(prior_mass)),
        "posterior_product": "pyRing_220_221_delta_weighted_posterior",
        "likelihood_scope": "marginalized_public_posterior_product_not_strain_level",
        "uniform_native_prior_assumption": {
            "Mf": True,
            "af": True,
            "domega_221": True,
            "dtau_221": True,
        },
        "gmm_components_tested": component_counts,
        "preferred_components_for_summary": preferred_components,
        "all_gmms_converged": all(bool(row["converged"]) for row in diagnostics),
        "preferred_results": preferred,
    }
    (args.output / "bayesian_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps({key: value for key, value in summary.items() if key != "preferred_results"}, indent=2))
    print(f"evidence_rows={len(evidence_rows)}")


if __name__ == "__main__":
    main()
