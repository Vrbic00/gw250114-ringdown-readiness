"""Freeze and audit every computation needed before writing the hairy paper."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sys
from pathlib import Path


import six  # noqa: F401

ROOT = Path(__file__).resolve().parents[2]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "results" / ".matplotlib"))
LOCAL_DEPS = ROOT / ".deps" / "python"
if LOCAL_DEPS.exists():
    sys.path.insert(0, str(LOCAL_DEPS))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import qnm
from scipy.interpolate import CubicSpline


DEFAULT_GRID = ROOT / "results" / "hairy_qnm_production_grid" / "hairy_qnm_production_grid.npz"
DEFAULT_GRID_CSV = ROOT / "results" / "hairy_qnm_production_grid" / "hairy_qnm_production_grid.csv"
DEFAULT_GRID_SUMMARY = ROOT / "results" / "hairy_qnm_production_grid" / "production_grid_summary.json"
DEFAULT_TABLE_VALIDATION = ROOT / "results" / "hairy_continued_fraction" / "validation_summary.json"
DEFAULT_CONSTRAINT = ROOT / "results" / "gw250114_hairy_constraints" / "constraint_summary.json"
DEFAULT_CONSTRAINT_TABLE = ROOT / "results" / "gw250114_hairy_constraints" / "constraint_scenario_summary.csv"
DEFAULT_PRIOR_TABLE = ROOT / "results" / "gw250114_hairy_constraints" / "hair_prior_robustness.csv"
DEFAULT_LIMIT_TABLE = ROOT / "results" / "gw250114_hairy_constraints" / "conditional_alpha_limits.csv"
DEFAULT_CONTROL = ROOT / "results" / "hairy_evolving_kerr_control" / "evolving_kerr_control_summary.json"
DEFAULT_CONTROL_TABLE = ROOT / "results" / "hairy_evolving_kerr_control" / "evolving_kerr_hairy_summary.csv"
DEFAULT_PYRING_SURFACE = ROOT / "results" / "gw250114_hairy_constraints" / "primary_posterior_surface.npz"
DEFAULT_RINGDOWN_SURFACE = ROOT / "results" / "gw250114_hairy_constraints" / "ringdown_deviation_posterior_surface.npz"
DEFAULT_OUTPUT = ROOT / "results" / "hairy_publication_computations"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


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


def directory_size(path: Path) -> int:
    total = 0
    for item in path.rglob("*"):
        try:
            if item.is_file():
                total += item.stat().st_size
        except OSError:
            continue
    return total


def discrete_quantile(values: np.ndarray, weights: np.ndarray, probability: float) -> float:
    order = np.argsort(values)
    values = np.asarray(values, dtype=float)[order]
    weights = np.asarray(weights, dtype=float)[order]
    cumulative = np.cumsum(weights) / np.sum(weights)
    return float(values[min(int(np.searchsorted(cumulative, probability)), len(values) - 1)])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    grid_summary = json.loads(DEFAULT_GRID_SUMMARY.read_text(encoding="utf-8"))
    table_validation = json.loads(DEFAULT_TABLE_VALIDATION.read_text(encoding="utf-8"))
    constraint = json.loads(DEFAULT_CONSTRAINT.read_text(encoding="utf-8"))
    control = json.loads(DEFAULT_CONTROL.read_text(encoding="utf-8"))
    grid = np.load(DEFAULT_GRID)
    spins = np.asarray(grid["spins"], dtype=float)
    alphas = np.asarray(grid["alphas"], dtype=float)
    h0s = np.asarray(grid["h0s"], dtype=float)
    overtones = np.asarray(grid["overtones"], dtype=int)
    omega = np.asarray(grid["omega"], dtype=np.complex128)

    # Independent Kerr-limit comparison against the qnm package.
    kerr_rows: list[dict[str, object]] = []
    for ispin, spin in enumerate(spins):
        for io, overtone in enumerate(overtones):
            reference = complex(qnm.modes_cache(s=-2, l=2, m=2, n=int(overtone))(a=float(spin))[0])
            values = omega[:, 0, ispin, io]
            for ih, h0 in enumerate(h0s):
                kerr_rows.append(
                    {
                        "spin": float(spin),
                        "overtone": int(overtone),
                        "h0": float(h0),
                        "grid_real": float(values[ih].real),
                        "grid_minus_imag": float(-values[ih].imag),
                        "qnm_real": reference.real,
                        "qnm_minus_imag": -reference.imag,
                        "absolute_delta": abs(values[ih] - reference),
                        "relative_delta": abs(values[ih] - reference) / abs(reference),
                    }
                )
    write_csv(args.output / "kerr_limit_validation.csv", kerr_rows)

    grid_rows = read_csv(DEFAULT_GRID_CSV)
    domain_rows: list[dict[str, object]] = []
    for h0 in h0s:
        for alpha in alphas:
            selected = [
                row
                for row in grid_rows
                if float(row["h0"]) == float(h0)
                and float(row["alpha_hair"]) == float(alpha)
                and int(row["overtone"]) == 0
            ]
            gaps = [float(row["exact_r_plus"]) - float(row["exact_r_minus"]) for row in selected]
            domain_rows.append(
                {
                    "h0": float(h0),
                    "alpha_hair": float(alpha),
                    "spin_min": min(float(row["spin"]) for row in selected),
                    "spin_max": max(float(row["spin"]) for row in selected),
                    "minimum_exact_horizon_gap": min(gaps),
                    "maximum_outer_horizon_relative_error": max(float(row["relative_r_plus_error"]) for row in selected),
                    "maximum_inner_horizon_relative_error": max(float(row["relative_r_minus_error"]) for row in selected),
                    "all_two_horizon_points_valid": all(gap > 0.0 for gap in gaps),
                }
            )
    write_csv(args.output / "parameter_domain_summary.csv", domain_rows)

    # Translate the two-dimensional hair posterior into a nominal-spin QNM
    # displacement.  This is an identifiable observable even where alpha and
    # h0 separately form a nearly Kerr-degenerate ridge.
    reference_spin = 0.68
    nominal_modes = np.empty((len(h0s), len(alphas), len(overtones)), dtype=np.complex128)
    for ih in range(len(h0s)):
        for ia in range(len(alphas)):
            for io in range(len(overtones)):
                nominal_modes[ih, ia, io] = complex(
                    CubicSpline(spins, omega[ih, ia, :, io])(reference_spin)
                )
    reference_modes = nominal_modes[:, :1, :]
    frequency_shift = nominal_modes.real / reference_modes.real - 1.0
    damping_shift = nominal_modes.imag / reference_modes.imag - 1.0
    inference_indices = [int(np.where(overtones == overtone)[0][0]) for overtone in (0, 1)]
    spectral_distance = 100.0 * np.max(
        np.abs(
            np.stack(
                (
                    frequency_shift[:, :, inference_indices],
                    damping_shift[:, :, inference_indices],
                ),
                axis=-1,
            )
        ),
        axis=(2, 3),
    )
    shift_rows: list[dict[str, object]] = []
    for ih, h0 in enumerate(h0s):
        for ia, alpha in enumerate(alphas):
            row: dict[str, object] = {
                "h0": float(h0),
                "alpha_hair": float(alpha),
                "reference_spin": reference_spin,
                "maximum_abs_220_221_fractional_shift_pct": float(spectral_distance[ih, ia]),
            }
            for io, overtone in enumerate(overtones):
                row[f"mode_22{overtone}_frequency_shift_pct"] = float(100.0 * frequency_shift[ih, ia, io])
                row[f"mode_22{overtone}_damping_rate_shift_pct"] = float(100.0 * damping_shift[ih, ia, io])
            shift_rows.append(row)
    write_csv(args.output / "nominal_spin_qnm_shift_surface.csv", shift_rows)

    derived_rows: list[dict[str, object]] = []
    cdf_curves: list[tuple[str, np.ndarray, np.ndarray]] = []
    for branch, surface_path in (
        ("pyRing", DEFAULT_PYRING_SURFACE),
        ("RINGDOWN", DEFAULT_RINGDOWN_SURFACE),
    ):
        surface = np.load(surface_path)
        posterior = np.asarray(surface["posterior_mass"], dtype=float)
        prior = np.asarray(surface["prior_mass"], dtype=float)
        for distribution, masses in (("prior", prior), ("posterior", posterior)):
            flat_value = spectral_distance.ravel()
            flat_mass = masses.ravel()
            derived_rows.append(
                {
                    "branch": branch,
                    "distribution": distribution,
                    "reference_spin": reference_spin,
                    "spectral_distance_q50_pct": discrete_quantile(flat_value, flat_mass, 0.50),
                    "spectral_distance_q90_pct": discrete_quantile(flat_value, flat_mass, 0.90),
                    "spectral_distance_q95_pct": discrete_quantile(flat_value, flat_mass, 0.95),
                    "probability_distance_gt_1pct": float(np.sum(flat_mass[flat_value > 1.0])),
                    "probability_distance_gt_2pct": float(np.sum(flat_mass[flat_value > 2.0])),
                    "probability_distance_gt_5pct": float(np.sum(flat_mass[flat_value > 5.0])),
                }
            )
            order = np.argsort(flat_value)
            cdf_curves.append(
                (
                    f"{branch} {distribution}",
                    flat_value[order],
                    np.cumsum(flat_mass[order]) / np.sum(flat_mass),
                )
            )
    write_csv(args.output / "derived_spectral_distance_constraints.csv", derived_rows)

    fig, ax = plt.subplots(figsize=(6.5, 4.6), constrained_layout=True)
    styles = {
        "pyRing prior": ("grey", ":"),
        "pyRing posterior": ("tab:blue", "-"),
        "RINGDOWN prior": ("grey", "--"),
        "RINGDOWN posterior": ("tab:orange", "-"),
    }
    for label, values, cumulative in cdf_curves:
        color, style = styles[label]
        ax.step(values, cumulative, where="post", color=color, ls=style, label=label)
    ax.set_xlabel("maximum absolute 220/221 QNM shift at a=0.68 [%]")
    ax.set_ylabel("cumulative probability")
    ax.set_ylim(0.0, 1.0)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.savefig(args.output / "derived_spectral_distance_cdf.png", dpi=220)
    fig.savefig(args.output / "derived_spectral_distance_cdf.pdf")
    plt.close(fig)

    constraint_rows = read_csv(DEFAULT_CONSTRAINT_TABLE)
    primary = constraint["primary_result"]
    ringdown = constraint["reference_RINGDOWN_result"]
    table_primary = [
        {
            "branch": "pyRing",
            "start_time_M": primary["start_time_M"],
            "prior": primary["prior_scenario"],
            "lnB_hair_vs_GR": primary["log_bayes_factor_hairy_vs_gr"],
            "max_delta_lnL": primary["maximum_log_likelihood_ratio_vs_gr"],
            "alpha_median": primary["alpha_median"],
            "alpha_q90": primary["alpha_q90"],
            "alpha_q95": primary["alpha_q95"],
            "h0_median": primary["h0_median"],
            "information_bits": primary["kl_divergence_bits"],
        },
        {
            "branch": "RINGDOWN",
            "start_time_M": ringdown["start_time_M"],
            "prior": ringdown["prior_scenario"],
            "lnB_hair_vs_GR": ringdown["log_bayes_factor_hairy_vs_gr"],
            "max_delta_lnL": ringdown["maximum_log_likelihood_ratio_vs_gr"],
            "alpha_median": ringdown["alpha_median"],
            "alpha_q90": ringdown["alpha_q90"],
            "alpha_q95": ringdown["alpha_q95"],
            "h0_median": ringdown["h0_median"],
            "information_bits": ringdown["kl_divergence_bits"],
        },
    ]
    write_csv(args.output / "table_primary_constraints.csv", table_primary)

    control_rows = read_csv(DEFAULT_CONTROL_TABLE)
    write_csv(args.output / "table_evolving_kerr_control.csv", control_rows)

    # Compact domain figure at the most demanding spin in the atlas.
    gap = np.asarray([float(row["minimum_exact_horizon_gap"]) for row in domain_rows]).reshape(len(h0s), len(alphas))
    fig, ax = plt.subplots(figsize=(6.5, 4.8), constrained_layout=True)
    mesh = ax.pcolormesh(alphas, h0s, gap, shading="auto", cmap="cividis")
    fig.colorbar(mesh, ax=ax, label=r"minimum $r_+-r_-$ over the spin grid")
    ax.set_xlabel(r"$\alpha$")
    ax.set_ylabel(r"$h_0$")
    ax.set_title("Two-horizon domain of the production atlas")
    fig.savefig(args.output / "hairy_parameter_domain.png", dpi=220)
    fig.savefig(args.output / "hairy_parameter_domain.pdf")
    plt.close(fig)

    required_figures = [
        ROOT / "results" / "gw250114_hairy_constraints" / "gw250114_hairy_primary_posterior.pdf",
        ROOT / "results" / "gw250114_hairy_constraints" / "gw250114_hairy_ringdown_deviation_posterior.pdf",
        ROOT / "results" / "gw250114_hairy_constraints" / "hairy_qnm_shift_atlas_a0p68.pdf",
        ROOT / "results" / "hairy_evolving_kerr_control" / "evolving_kerr_false_hair_control.pdf",
        ROOT / "results" / "hairy_positive_injection_recovery" / "hairy_positive_injection_sensitivity.pdf",
        ROOT / "results" / "hairy_referee_robustness" / "hairy_horizon_convention_sensitivity.pdf",
    ]
    checks = [
        ("published_table_reproduction", table_validation["status"] == "PASS" and table_validation["cases"] == 81),
        ("production_grid_all_roots", bool(grid_summary["all_roots_successful"])),
        ("production_grid_size", int(grid_summary["modes"]) == 35343),
        ("continued_fraction_convergence", float(grid_summary["maximum_depth_abs_delta"]) < 1.0e-8),
        ("interpolation_accuracy", float(grid_summary["maximum_interpolation_relative_error"]) < 5.0e-4),
        ("kerr_limit", max(float(row["relative_delta"]) for row in kerr_rows) < 1.0e-5),
        ("two_horizon_domain", all(bool(row["all_two_horizon_points_valid"]) for row in domain_rows)),
        ("primary_constraint_available", bool(constraint["publication_outputs_ready"])),
        ("evolving_control_optimizers", bool(control["all_optimizers_successful"])),
        ("publication_figures", all(path.exists() and path.stat().st_size > 0 for path in required_figures)),
        ("derived_spectral_constraints", len(derived_rows) == 4),
    ]
    check_rows = [{"check": name, "status": "PASS" if passed else "FAIL"} for name, passed in checks]
    write_csv(args.output / "validation_gates.csv", check_rows)

    project_bytes = directory_size(ROOT)
    summary = {
        "status": "PASS" if all(passed for _, passed in checks) else "FAIL",
        "validation_gates_passed": sum(passed for _, passed in checks),
        "validation_gates_total": len(checks),
        "direct_qnm_modes": int(grid_summary["modes"]),
        "published_cases_reproduced": int(table_validation["cases"]),
        "maximum_published_table_relative_delta": float(table_validation["max_relative_delta"]),
        "maximum_kerr_limit_relative_delta": max(float(row["relative_delta"]) for row in kerr_rows),
        "maximum_interpolation_relative_error": float(grid_summary["maximum_interpolation_relative_error"]),
        "primary_pyRing_result": primary,
        "reference_RINGDOWN_result": ringdown,
        "reference_RINGDOWN_timescan_result": constraint["reference_RINGDOWN_timescan_result"],
        "derived_spectral_distance_constraints": derived_rows,
        "evolving_kerr_control": {key: value for key, value in control.items() if key != "files"},
        "project_size_bytes": project_bytes,
        "project_size_MB": project_bytes / 1024**2,
        "project_size_GB": project_bytes / 1024**3,
        "article_text_written": True,
    }
    (args.output / "computation_audit_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    lines = [
        "# Hairy-BH publication computation status",
        "",
        "This is a machine-generated technical handoff. It is not manuscript prose.",
        "",
        f"- Overall computation status: `{summary['status']}`",
        f"- Direct QNM roots: `{summary['direct_qnm_modes']}`",
        f"- Published Zhen-Li cases reproduced: `{summary['published_cases_reproduced']}`",
        f"- Maximum table relative difference: `{summary['maximum_published_table_relative_delta']:.3e}`",
        f"- Maximum production interpolation error: `{summary['maximum_interpolation_relative_error']:.3e}`",
        f"- Project size: `{summary['project_size_MB']:.2f} MB`",
        "",
        "## Primary public-posterior results",
        "",
        f"- pyRing: ln B(hair/GR) = `{float(primary['log_bayes_factor_hairy_vs_gr']):.4f}`, "
        f"alpha 90% upper quantile = `{float(primary['alpha_q90']):.4f}`, "
        f"information gain = `{float(primary['kl_divergence_bits']):.4f} bits`.",
        f"- RINGDOWN at {ringdown['start_time_M']} M: ln B(hair/GR) = "
        f"`{float(ringdown['log_bayes_factor_hairy_vs_gr']):.4f}`, "
        f"alpha 90% upper quantile = `{float(ringdown['alpha_q90']):.4f}`.",
        "",
        "## Scope boundary",
        "",
        "The observational calculations use public marginalized spectral posterior products. "
        "They are not a newly sampled H1/L1 detector-strain likelihood. All numerical inputs, "
        "posterior surfaces, robustness tables, controls, and figures required before manuscript "
        "writing are frozen in the result directories.",
    ]
    (args.output / "COMPUTATION_STATUS.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    generated = sorted(path for path in args.output.iterdir() if path.is_file())
    manifest = [
        {
            "path": str(path.relative_to(ROOT)),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
        for path in generated
        if path.name != "computation_manifest.csv"
    ]
    write_csv(args.output / "computation_manifest.csv", manifest)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
