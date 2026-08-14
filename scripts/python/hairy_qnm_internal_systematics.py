"""Quantify computable internal systematics of the rotating-hairy QNM grid.

This does not attempt to estimate the uncomputable error of replacing the
full coupled perturbation problem by the Dudley--Finley-like equation.  It
isolates two narrower effects that can be calculated from the published
setup: second-order versus numerical horizons, and the ambiguous horizon
indices in the continued-fraction reduction.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
LOCAL_DEPS = ROOT / ".deps" / "python"
if LOCAL_DEPS.exists():
    sys.path.insert(0, str(LOCAL_DEPS))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np

from build_hairy_qnm_grid import exact_horizons
from hairy_continued_fraction import kerr_seed, solve_mode


DEFAULT_GRID = ROOT / "results" / "hairy_qnm_grid" / "hairy_qnm_grid.npz"
DEFAULT_OUTPUT = ROOT / "results" / "hairy_qnm_internal_systematics"


def statistics(values: list[float]) -> dict[str, float]:
    array = np.asarray(values, dtype=float)
    return {
        "median": float(np.median(array)),
        "p90": float(np.quantile(array, 0.90)),
        "p95": float(np.quantile(array, 0.95)),
        "max": float(np.max(array)),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--grid", type=Path, default=DEFAULT_GRID)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--depth", type=int, default=320)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    data = np.load(args.grid)
    spins = data["spins"]
    alphas = data["alphas"]
    h0s = data["h0s"]
    overtones = data["overtones"]
    baseline = data["omega"]
    variants = {
        "ode_approximate": np.full_like(baseline, np.nan + 1j * np.nan),
        "table_exact": np.full_like(baseline, np.nan + 1j * np.nan),
        "ode_exact": np.full_like(baseline, np.nan + 1j * np.nan),
    }
    rows: list[dict[str, object]] = []

    for ih, h0 in enumerate(h0s):
        for ia, alpha_hair in enumerate(alphas):
            for ispin, spin in enumerate(spins):
                exact = exact_horizons(float(spin), float(alpha_hair), float(h0))
                for io, overtone in enumerate(overtones):
                    base = complex(baseline[ih, ia, ispin, io])
                    if alpha_hair == 0.0:
                        ode_approximate = base
                        table_exact = base
                        ode_exact = base
                    else:
                        ode_approximate = solve_mode(
                            float(spin),
                            float(alpha_hair),
                            float(h0),
                            int(overtone),
                            initial=base,
                            truncation=args.depth,
                            horizon_index_convention="ode_consistent",
                        ).omega
                        table_exact = solve_mode(
                            float(spin),
                            float(alpha_hair),
                            float(h0),
                            int(overtone),
                            initial=base,
                            truncation=args.depth,
                            horizons=exact,
                            horizon_source="numerical_original_delta_roots",
                        ).omega
                        ode_exact = solve_mode(
                            float(spin),
                            float(alpha_hair),
                            float(h0),
                            int(overtone),
                            initial=ode_approximate,
                            truncation=args.depth,
                            horizon_index_convention="ode_consistent",
                            horizons=exact,
                            horizon_source="numerical_original_delta_roots",
                        ).omega

                    variants["ode_approximate"][ih, ia, ispin, io] = ode_approximate
                    variants["table_exact"][ih, ia, ispin, io] = table_exact
                    variants["ode_exact"][ih, ia, ispin, io] = ode_exact
                    kerr = kerr_seed(float(spin), int(overtone))
                    hair_shift = abs(base - kerr)
                    index_delta = abs(ode_approximate - base)
                    horizon_delta = abs(table_exact - base)
                    combined_delta = abs(ode_exact - base)
                    rows.append(
                        {
                            "h0": float(h0),
                            "alpha_hair": float(alpha_hair),
                            "spin": float(spin),
                            "overtone": int(overtone),
                            "baseline_real": base.real,
                            "baseline_minus_imag": -base.imag,
                            "ode_approximate_real": ode_approximate.real,
                            "ode_approximate_minus_imag": -ode_approximate.imag,
                            "table_exact_real": table_exact.real,
                            "table_exact_minus_imag": -table_exact.imag,
                            "ode_exact_real": ode_exact.real,
                            "ode_exact_minus_imag": -ode_exact.imag,
                            "hair_shift_abs": hair_shift,
                            "index_convention_abs_delta": index_delta,
                            "horizon_source_abs_delta": horizon_delta,
                            "combined_abs_delta": combined_delta,
                            "index_delta_over_hair_shift": np.nan if hair_shift == 0.0 else index_delta / hair_shift,
                            "horizon_delta_over_hair_shift": np.nan if hair_shift == 0.0 else horizon_delta / hair_shift,
                            "combined_delta_over_hair_shift": np.nan if hair_shift == 0.0 else combined_delta / hair_shift,
                        }
                    )

    with (args.output / "internal_systematics.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    np.savez_compressed(
        args.output / "internal_systematics.npz",
        spins=spins,
        alphas=alphas,
        h0s=h0s,
        overtones=overtones,
        baseline=baseline,
        **variants,
    )

    deformed = [row for row in rows if float(row["alpha_hair"]) > 0.0]
    summary: dict[str, object] = {
        "cases": len(rows),
        "deformed_cases": len(deformed),
        "depth": args.depth,
        "scope_warning": (
            "Internal implementation effects only; excludes the structural error of the "
            "Dudley-Finley approximation and any missing perturbation sector."
        ),
    }
    for key in (
        "index_convention_abs_delta",
        "horizon_source_abs_delta",
        "combined_abs_delta",
        "index_delta_over_hair_shift",
        "horizon_delta_over_hair_shift",
        "combined_delta_over_hair_shift",
    ):
        summary[key] = statistics([float(row[key]) for row in deformed if np.isfinite(float(row[key]))])
    reference = [
        row
        for row in rows
        if float(row["spin"]) == 0.7
        and float(row["alpha_hair"]) == 0.5
        and float(row["h0"]) == 1.0
    ]
    summary["reference_a0p7_alpha0p5_h01"] = reference
    (args.output / "systematics_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
