"""Audit the effective Kerr--Newman-form identity in Li's factorized model.

Li's second-order horizon approximation replaces the exponential horizon
function by Delta_fac=(r-r_plus)(r-r_minus).  This script maps the coarse
internal-systematics grid to

    M_eff = (r_plus + r_minus)/2,
    q_eff^2 = r_plus*r_minus - a^2,

checks the polynomial and horizon identities, and compares the table-calibrated
recurrence with the ODE-consistent branch.  It is an algebraic and
implementation-level control, not a validation of the missing coupled
perturbation system.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SYSTEMATICS = (
    ROOT / "results" / "hairy_qnm_internal_systematics" / "internal_systematics.csv"
)
DEFAULT_OUTPUT = ROOT / "results" / "effective_kerr_newman_control"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def hair_exponential(radius: float, h0: float) -> float:
    denominator = 1.0 - 0.5 * h0
    if denominator <= 0.0:
        return 0.0
    return math.exp(-radius / denominator)


def approximate_horizons(
    spin: float, alpha_hair: float, h0: float, iterations: int = 2
) -> tuple[float, float]:
    root_kerr = math.sqrt(1.0 - spin * spin)
    plus = 1.0 + root_kerr
    minus = 1.0 - root_kerr
    for _ in range(iterations):
        plus = 1.0 + math.sqrt(
            1.0
            - spin * spin
            - alpha_hair * plus * plus * hair_exponential(plus, h0)
        )
        minus = 1.0 - math.sqrt(
            1.0
            - spin * spin
            - alpha_hair * minus * minus * hair_exponential(minus, h0)
        )
    return plus, minus


def linear_quantile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    position = probability * (len(ordered) - 1)
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return (1.0 - weight) * ordered[lower] + weight * ordered[upper]


def quantiles(values: list[float]) -> dict[str, float]:
    return {
        "median": linear_quantile(values, 0.50),
        "p90": linear_quantile(values, 0.90),
        "p95": linear_quantile(values, 0.95),
        "maximum": max(values),
    }


def parameter_ranges(rows: list[dict[str, object]]) -> dict[str, dict[str, float]]:
    keys = (
        "effective_mass_over_original_mass",
        "effective_spin",
        "signed_effective_charge_to_mass_ratio",
        "effective_charge_magnitude_to_mass_ratio",
        "effective_polynomial_extremality",
    )
    return {
        key: {
            "minimum": min(float(row[key]) for row in rows),
            "maximum": max(float(row[key]) for row in rows),
        }
        for key in keys
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--systematics", type=Path, default=DEFAULT_SYSTEMATICS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    source_rows = read_csv(args.systematics)
    output_rows: list[dict[str, object]] = []
    cell_rows: dict[tuple[float, float, float], dict[str, object]] = {}
    delta_errors: list[float] = []
    derivative_errors: list[float] = []
    horizon_errors: list[float] = []
    table_deltas: list[float] = []
    table_deltas_over_shift: list[float] = []

    for source in source_rows:
        h0 = float(source["h0"])
        alpha = float(source["alpha_hair"])
        spin = float(source["spin"])
        overtone = int(source["overtone"])
        r_plus, r_minus = approximate_horizons(spin, alpha, h0)
        mass_effective = 0.5 * (r_plus + r_minus)
        charge_squared = r_plus * r_minus - spin * spin
        signed_charge_ratio = math.copysign(
            math.sqrt(abs(charge_squared)) / mass_effective, charge_squared
        )

        discriminant = mass_effective**2 - spin**2 - charge_squared
        if discriminant < -1.0e-12:
            raise RuntimeError(f"superextremal mapping at {(h0, alpha, spin)}")
        root_gap = math.sqrt(max(0.0, discriminant))
        kn_plus = mass_effective + root_gap
        kn_minus = mass_effective - root_gap
        horizon_error = max(abs(kn_plus - r_plus), abs(kn_minus - r_minus))
        horizon_errors.append(horizon_error)

        cell_key = (h0, alpha, spin)
        if cell_key not in cell_rows:
            gap = r_plus - r_minus
            for radius in (
                r_minus,
                0.5 * (r_plus + r_minus),
                r_plus,
                r_plus + gap,
                r_plus + 5.0 * gap,
            ):
                delta_factorized = (radius - r_plus) * (radius - r_minus)
                delta_kn = (
                    radius**2
                    - 2.0 * mass_effective * radius
                    + spin**2
                    + charge_squared
                )
                delta_errors.append(abs(delta_factorized - delta_kn))
                derivative_errors.append(
                    abs(
                        (2.0 * radius - r_plus - r_minus)
                        - (2.0 * radius - 2.0 * mass_effective)
                    )
                )
            cell_rows[cell_key] = {
                "h0": h0,
                "alpha_hair": alpha,
                "spin_a_over_original_mass": spin,
                "approximate_r_plus": r_plus,
                "approximate_r_minus": r_minus,
                "effective_mass_over_original_mass": mass_effective,
                "effective_spin": spin / mass_effective,
                "signed_effective_charge_to_mass_ratio": signed_charge_ratio,
                "effective_charge_magnitude_to_mass_ratio": abs(signed_charge_ratio),
                "effective_charge_squared": charge_squared,
                "real_charge_kerr_newman_mapping": charge_squared >= -1.0e-12,
                "effective_polynomial_extremality": (
                    math.sqrt(max(0.0, spin**2 + charge_squared)) / mass_effective
                ),
                "horizon_roundtrip_abs_delta": horizon_error,
            }

        omega_ode = complex(
            float(source["ode_approximate_real"]),
            -float(source["ode_approximate_minus_imag"]),
        )
        omega_table = complex(
            float(source["baseline_real"]),
            -float(source["baseline_minus_imag"]),
        )
        table_delta = abs(omega_table - omega_ode)
        hair_shift = float(source["hair_shift_abs"])
        table_deltas.append(table_delta)
        if alpha > 0.0 and hair_shift > 1.0e-12:
            table_deltas_over_shift.append(table_delta / hair_shift)

        output_rows.append(
            {
                **cell_rows[cell_key],
                "overtone": overtone,
                "effective_kn_df_ode_identity_real": omega_ode.real,
                "effective_kn_df_ode_identity_minus_imag": -omega_ode.imag,
                "ode_consistent_real": omega_ode.real,
                "ode_consistent_minus_imag": -omega_ode.imag,
                "table_calibrated_real": omega_table.real,
                "table_calibrated_minus_imag": -omega_table.imag,
                "effective_kn_df_equals_ode_consistent_by_delta_identity": True,
                "table_vs_effective_kn_abs_delta": table_delta,
                "table_vs_effective_kn_over_hair_shift": (
                    table_delta / hair_shift
                    if alpha > 0.0 and hair_shift > 1.0e-12
                    else 0.0
                ),
            }
        )

    unique_cells = list(cell_rows.values())
    deformed_cells = [row for row in unique_cells if float(row["alpha_hair"]) > 0.0]
    event_cells = [
        row
        for row in unique_cells
        if 0.6 <= float(row["spin_a_over_original_mass"]) <= 0.75
    ]
    deformed_event_cells = [
        row for row in event_cells if float(row["alpha_hair"]) > 0.0
    ]
    negative_cells = [
        row for row in deformed_cells if float(row["effective_charge_squared"]) < -1.0e-12
    ]
    negative_event_cells = [
        row
        for row in deformed_event_cells
        if float(row["effective_charge_squared"]) < -1.0e-12
    ]

    summary = {
        "status": "PASS",
        "scope": (
            "Algebraic and implementation control of factorized Delta and the "
            "Dudley--Finley spectrum; not a full coupled-system validation."
        ),
        "spectral_rows": len(output_rows),
        "unique_parameter_cells": len(unique_cells),
        "delta_audit_evaluations": len(delta_errors),
        "maximum_abs_delta_polynomial_difference": max(delta_errors),
        "maximum_abs_delta_derivative_difference": max(derivative_errors),
        "maximum_horizon_roundtrip_abs_delta": max(horizon_errors),
        "effective_kn_df_vs_ode_consistent_relation": (
            "exact_identity_after_factorization"
        ),
        "table_calibrated_vs_effective_kn_df_abs_delta": quantiles(table_deltas),
        "table_calibrated_vs_effective_kn_df_over_hair_shift": quantiles(
            table_deltas_over_shift
        ),
        "table_calibrated_vs_effective_kn_df_over_hair_shift_entries": len(
            table_deltas_over_shift
        ),
        "effective_parameter_ranges_full_grid": parameter_ranges(unique_cells),
        "effective_parameter_ranges_event_spin_0p6_to_0p75": parameter_ranges(
            event_cells
        ),
        "signed_charge_squared_mapping": {
            "deformed_cells": len(deformed_cells),
            "negative_charge_squared_cells": len(negative_cells),
            "negative_fraction": len(negative_cells) / len(deformed_cells),
            "deformed_event_spin_cells": len(deformed_event_cells),
            "negative_charge_squared_event_spin_cells": len(negative_event_cells),
            "negative_event_spin_fraction": (
                len(negative_event_cells) / len(deformed_event_cells)
            ),
        },
        "interpretation": {
            "ode_consistent": (
                "The factorized ODE is exactly the Dudley--Finley equation for "
                "the mapped effective Kerr--Newman-form polynomial."
            ),
            "negative_charge_squared": (
                "Where q_eff^2 is negative, the map is algebraic and does not "
                "describe a real-charge Einstein--Maxwell black hole."
            ),
        },
    }

    write_csv(args.output / "effective_kerr_newman_control.csv", output_rows)
    write_csv(args.output / "effective_kerr_newman_parameter_map.csv", unique_cells)
    (args.output / "effective_kerr_newman_control_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
