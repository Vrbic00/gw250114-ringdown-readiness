"""Build and validate a regular rotating-hairy QNM grid.

The grid is restricted to the parameter region for which Zhen Li quotes the
second-order horizon approximation as controlled: alpha <= 0.5, h0 >= 0.5,
and a <= 0.9.  Every stored mode is solved at two continued-fraction depths.
An independent midpoint sample measures trilinear interpolation error.
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
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import brentq

from hairy_continued_fraction import approximate_horizons, kerr_seed, solve_mode


DEFAULT_OUTPUT = ROOT / "results" / "hairy_qnm_grid"


def exact_horizons(spin: float, alpha_hair: float, h0: float) -> tuple[float, float]:
    """Numerically locate the two positive roots of the original Delta."""
    denominator = 1.0 - 0.5 * h0

    def delta(radius: float) -> float:
        exponential = 0.0 if denominator <= 0.0 else np.exp(-radius / denominator)
        return radius * radius + spin * spin - 2.0 * radius + alpha_hair * radius * radius * exponential

    # Throughout the controlled production domain the inner and outer roots
    # bracket r=1.  Direct bracketing is more than two orders of magnitude
    # faster than scanning thousands of radii for every grid point.  Retain a
    # dense fallback so that a future extension of the domain fails safely
    # instead of silently returning the wrong pair of roots.
    roots: list[float] = []
    left_value = delta(1.0e-10)
    middle_value = delta(1.0)
    right_value = delta(4.0)
    if left_value * middle_value < 0.0 and middle_value * right_value < 0.0:
        roots = [
            float(brentq(delta, 1.0e-10, 1.0, xtol=1.0e-14)),
            float(brentq(delta, 1.0, 4.0, xtol=1.0e-14)),
        ]
    else:
        radii = np.linspace(1.0e-10, 4.0, 8001)
        values = np.asarray([delta(float(radius)) for radius in radii])
        for left, right, f_left, f_right in zip(radii[:-1], radii[1:], values[:-1], values[1:]):
            if f_left == 0.0:
                roots.append(float(left))
            elif f_left * f_right < 0.0:
                roots.append(float(brentq(delta, float(left), float(right), xtol=1.0e-14)))
    unique = sorted({round(value, 13) for value in roots})
    if len(unique) != 2:
        raise ValueError(f"expected two exact horizons, found {unique}")
    return float(unique[1]), float(unique[0])


def solve_grid(
    spins: np.ndarray,
    alphas: np.ndarray,
    h0s: np.ndarray,
    overtones: np.ndarray,
    low_depth: int,
    high_depth: int,
) -> tuple[np.ndarray, list[dict[str, object]]]:
    shape = (len(h0s), len(alphas), len(spins), len(overtones))
    omega = np.full(shape, np.nan + 1j * np.nan, dtype=np.complex128)
    rows: list[dict[str, object]] = []

    for ih, h0 in enumerate(h0s):
        for ia, alpha_hair in enumerate(alphas):
            for ispin, spin in enumerate(spins):
                approximate_plus, approximate_minus = approximate_horizons(float(spin), float(alpha_hair), float(h0))
                exact_plus, exact_minus = exact_horizons(float(spin), float(alpha_hair), float(h0))
                for io, overtone in enumerate(overtones):
                    seed = kerr_seed(float(spin), int(overtone))
                    if ia > 0 and np.isfinite(omega[ih, ia - 1, ispin, io].real):
                        seed = complex(omega[ih, ia - 1, ispin, io])
                    low = solve_mode(
                        float(spin),
                        float(alpha_hair),
                        float(h0),
                        int(overtone),
                        initial=seed,
                        truncation=low_depth,
                    )
                    high = solve_mode(
                        float(spin),
                        float(alpha_hair),
                        float(h0),
                        int(overtone),
                        initial=low.omega,
                        truncation=high_depth,
                    )
                    omega[ih, ia, ispin, io] = high.omega
                    rows.append(
                        {
                            "h0": float(h0),
                            "alpha_hair": float(alpha_hair),
                            "spin": float(spin),
                            "overtone": int(overtone),
                            "omega_real": high.omega.real,
                            "minus_omega_imag": -high.omega.imag,
                            "low_depth": low_depth,
                            "high_depth": high_depth,
                            "depth_abs_delta": abs(high.omega - low.omega),
                            "residual_norm": high.residual_norm,
                            "root_success": high.root_success,
                            "approximate_r_plus": approximate_plus,
                            "approximate_r_minus": approximate_minus,
                            "exact_r_plus": exact_plus,
                            "exact_r_minus": exact_minus,
                            "relative_r_plus_error": abs(approximate_plus - exact_plus) / exact_plus,
                            "relative_r_minus_error": abs(approximate_minus - exact_minus) / exact_minus,
                            "horizon_index_convention": high.horizon_index_convention,
                        }
                    )
    return omega, rows


def validate_interpolator(
    spins: np.ndarray,
    alphas: np.ndarray,
    h0s: np.ndarray,
    overtones: np.ndarray,
    omega: np.ndarray,
    high_depth: int,
    samples: int,
) -> list[dict[str, object]]:
    rng = np.random.default_rng(221208112)
    interpolators = [
        RegularGridInterpolator(
            (h0s, alphas, spins),
            omega[:, :, :, io],
            method="cubic",
            bounds_error=True,
        )
        for io in range(len(overtones))
    ]
    rows: list[dict[str, object]] = []
    for sample in range(samples):
        ih = int(rng.integers(0, len(h0s) - 1))
        ia = int(rng.integers(0, len(alphas) - 1))
        ispin = int(rng.integers(0, len(spins) - 1))
        point = np.asarray(
            [
                0.5 * (h0s[ih] + h0s[ih + 1]),
                0.5 * (alphas[ia] + alphas[ia + 1]),
                0.5 * (spins[ispin] + spins[ispin + 1]),
            ],
            dtype=float,
        )
        for io, overtone in enumerate(overtones):
            predicted = complex(interpolators[io](point).item())
            direct = solve_mode(
                spin=float(point[2]),
                alpha_hair=float(point[1]),
                h0=float(point[0]),
                overtone=int(overtone),
                initial=predicted,
                truncation=high_depth,
            )
            rows.append(
                {
                    "sample": sample,
                    "h0": point[0],
                    "alpha_hair": point[1],
                    "spin": point[2],
                    "overtone": int(overtone),
                    "interpolated_real": predicted.real,
                    "interpolated_minus_imag": -predicted.imag,
                    "direct_real": direct.omega.real,
                    "direct_minus_imag": -direct.omega.imag,
                    "absolute_error": abs(predicted - direct.omega),
                    "relative_error": abs(predicted - direct.omega) / abs(direct.omega),
                    "direct_root_success": direct.root_success,
                    "interpolation_method": "tensor_product_cubic",
                }
            )
    return rows


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--low-depth", type=int, default=160)
    parser.add_argument("--high-depth", type=int, default=320)
    parser.add_argument("--validation-samples", type=int, default=128)
    args = parser.parse_args()

    spins = np.round(np.arange(0.10, 0.9001, 0.05), 10)
    alphas = np.round(np.arange(0.0, 0.5001, 0.05), 10)
    h0s = np.asarray([0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 1.99])
    overtones = np.asarray([0, 1, 2], dtype=int)
    args.output.mkdir(parents=True, exist_ok=True)

    omega, rows = solve_grid(spins, alphas, h0s, overtones, args.low_depth, args.high_depth)
    write_csv(args.output / "hairy_qnm_grid.csv", rows)
    np.savez_compressed(
        args.output / "hairy_qnm_grid.npz",
        spins=spins,
        alphas=alphas,
        h0s=h0s,
        overtones=overtones,
        omega=omega,
    )

    interpolation_rows = validate_interpolator(
        spins,
        alphas,
        h0s,
        overtones,
        omega,
        args.high_depth,
        args.validation_samples,
    )
    write_csv(args.output / "interpolation_validation.csv", interpolation_rows)

    summary = {
        "grid_shape": list(omega.shape),
        "modes": int(omega.size),
        "all_roots_successful": all(bool(row["root_success"]) for row in rows),
        "max_depth_abs_delta": max(float(row["depth_abs_delta"]) for row in rows),
        "median_depth_abs_delta": float(np.median([float(row["depth_abs_delta"]) for row in rows])),
        "max_outer_horizon_relative_error": max(float(row["relative_r_plus_error"]) for row in rows),
        "max_inner_horizon_relative_error": max(float(row["relative_r_minus_error"]) for row in rows),
        "interpolation_cases": len(interpolation_rows),
        "max_interpolation_relative_error": max(float(row["relative_error"]) for row in interpolation_rows),
        "median_interpolation_relative_error": float(
            np.median([float(row["relative_error"]) for row in interpolation_rows])
        ),
        "domain": {
            "spin": [float(spins[0]), float(spins[-1])],
            "alpha_hair": [float(alphas[0]), float(alphas[-1])],
            "h0": [float(h0s[0]), float(h0s[-1])],
            "overtones": overtones.tolist(),
        },
        "continued_fraction_depths": [args.low_depth, args.high_depth],
        "horizon_index_convention": "table_calibrated",
        "interpolation_method": "tensor_product_cubic",
    }
    (args.output / "grid_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
