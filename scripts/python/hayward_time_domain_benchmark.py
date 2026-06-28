"""Time-domain axial gravitational benchmark for the Hayward spacetime.

Reference:
    Bolokhov & Skvortsova, "Gravitational quasinormal modes of the Hayward
    spacetime", arXiv:2508.19989.

Metric:
    f(r) = 1 - 2 r^2 / (r^3 + gamma)

Axial gravitational potential:
    V = f [2 f/r^2 - f'/r + ((ell+2)(ell-1))/r^2]

The comparison uses the paper's Table V Prony fits for ell=2, n=0.
"""

from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "python"))

from static_master_potential_time_domain import (  # noqa: E402
    characteristic_evolve_with_potential,
    fit_windows,
    pct_delta,
)


REFERENCE_PRONY_L2_N0 = {
    0.10: (0.375316, 0.088223),
    0.35: (0.379646, 0.086115),
    0.60: (0.384318, 0.083517),
    0.85: (0.389335, 0.080214),
    1.10: (0.394567, 0.075869),
    1.18: (0.396224, 0.074198),
}


def f_lapse(r: float | np.ndarray, gamma: float) -> float | np.ndarray:
    return 1.0 - 2.0 * r * r / (r**3 + gamma)


def f_prime(r: float | np.ndarray, gamma: float) -> float | np.ndarray:
    return 2.0 * (r**4 - 2.0 * gamma * r) / (r**3 + gamma) ** 2


def axial_potential_r(r: float | np.ndarray, gamma: float, ell: int = 2) -> float | np.ndarray:
    f = f_lapse(r, gamma)
    return f * (2.0 * f / (r * r) - f_prime(r, gamma) / r + ((ell + 2.0) * (ell - 1.0)) / (r * r))


def outer_horizon(gamma: float) -> float:
    if gamma == 0.0:
        return 2.0

    xs = np.linspace(1e-6, 10.0, 30000)
    fs = f_lapse(xs, gamma)
    brackets: list[tuple[float, float]] = []
    for left, right, f_left, f_right in zip(xs[:-1], xs[1:], fs[:-1], fs[1:]):
        if f_left == 0.0:
            brackets.append((float(left), float(left)))
        elif f_left * f_right < 0.0:
            brackets.append((float(left), float(right)))
    if not brackets:
        raise ValueError(f"no Hayward horizon found for gamma={gamma}")

    low, high = brackets[-1]
    for _ in range(100):
        mid = 0.5 * (low + high)
        if f_lapse(low, gamma) * f_lapse(mid, gamma) <= 0.0:
            high = mid
        else:
            low = mid
    return 0.5 * (low + high)


def peak_radius(gamma: float, ell: int = 2) -> float:
    rh = outer_horizon(gamma)
    left = rh + 1e-7
    right = 30.0
    inv_phi = (math.sqrt(5.0) - 1.0) / 2.0
    inv_phi2 = (3.0 - math.sqrt(5.0)) / 2.0
    h = right - left
    c = left + inv_phi2 * h
    d = left + inv_phi * h
    yc = axial_potential_r(c, gamma, ell)
    yd = axial_potential_r(d, gamma, ell)
    for _ in range(120):
        if yc < yd:
            left = c
            c = d
            yc = yd
            h = right - left
            d = left + inv_phi * h
            yd = axial_potential_r(d, gamma, ell)
        else:
            right = d
            d = c
            yd = yc
            h = right - left
            c = left + inv_phi2 * h
            yc = axial_potential_r(c, gamma, ell)
    return 0.5 * (left + right)


def shifted_potential_factory(gamma: float):
    rh = outer_horizon(gamma)
    r_peak = peak_radius(gamma)

    r_max = 700.0
    rho = np.geomspace(1e-8, r_max - rh, 240000)
    r_grid = rh + rho
    f_grid = f_lapse(r_grid, gamma)
    inv_f = 1.0 / f_grid
    dx = 0.5 * (inv_f[1:] + inv_f[:-1]) * np.diff(r_grid)
    x_raw = np.concatenate([[0.0], np.cumsum(dx)])
    x_peak = float(np.interp(r_peak, r_grid, x_raw))
    x_grid = x_raw - x_peak
    v_grid = axial_potential_r(r_grid, gamma)

    if x_grid[-1] < 300.0:
        raise ValueError(
            f"insufficient tortoise range for gamma={gamma}: {x_grid[0]} to {x_grid[-1]}"
        )

    def potential(x_values: np.ndarray) -> np.ndarray:
        return np.interp(x_values, x_grid, v_grid, left=0.0, right=v_grid[-1])

    return potential, rh, r_peak, x_grid[0], x_grid[-1]


def best_fit_for_waveform(times: np.ndarray, waveform: np.ndarray) -> dict[str, float]:
    windows = fit_windows(
        times,
        waveform,
        t_start_min=40.0,
        t_start_max=220.0,
        widths=(60.0, 70.0, 80.0, 100.0, 120.0, 140.0),
        step=10.0,
    )
    finite = [
        row
        for row in windows
        if math.isfinite(row["td_Momega_real"])
        and math.isfinite(row["td_Momega_imag_abs"])
        and math.isfinite(row["relative_residual"])
    ]
    if not finite:
        return {
            "fit_t_start": math.nan,
            "fit_t_stop": math.nan,
            "fit_width": math.nan,
            "td_Momega_real": math.nan,
            "td_Momega_imag_abs": math.nan,
            "relative_residual": math.inf,
        }
    return min(finite, key=lambda row: row["relative_residual"])


def main() -> None:
    out_dir = ROOT / "results" / "hayward_time_domain"
    out_dir.mkdir(parents=True, exist_ok=True)

    summary_rows: list[dict[str, float | str]] = []
    for gamma, (ref_real, ref_imag) in REFERENCE_PRONY_L2_N0.items():
        potential, rh, r_peak, x_min, x_max = shifted_potential_factory(gamma)
        times, waveform = characteristic_evolve_with_potential(
            potential,
            h=0.2,
            n_grid=2700,
            observer_rstar=20.0,
            gaussian_center=0.0,
            gaussian_width=3.0,
        )
        np.savetxt(
            out_dir / f"hayward_gamma_{gamma:g}_axial_l2_waveform.csv",
            np.column_stack([times, waveform]),
            delimiter=",",
            header="t_over_M,psi",
            comments="",
        )
        best = best_fit_for_waveform(times, waveform)
        summary_rows.append(
            {
                "case": f"hayward_gamma_{gamma:g}_axial_l2_n0",
                "gamma": gamma,
                "ell": 2,
                "n": 0,
                "horizon_outer": rh,
                "r_peak": r_peak,
                "x_min": x_min,
                "x_max": x_max,
                "reference_source": "Bolokhov_Skvortsova_2025_Table_V_Prony",
                "reference_Momega_real": ref_real,
                "reference_Momega_imag_abs": ref_imag,
                "fit_t_start": best["fit_t_start"],
                "fit_t_stop": best["fit_t_stop"],
                "fit_width": best["fit_width"],
                "td_Momega_real": best["td_Momega_real"],
                "td_Momega_imag_abs": best["td_Momega_imag_abs"],
                "delta_real_pct": pct_delta(float(best["td_Momega_real"]), ref_real),
                "delta_imag_pct": pct_delta(float(best["td_Momega_imag_abs"]), ref_imag),
                "relative_residual": best["relative_residual"],
            }
        )

    summary_csv = out_dir / "hayward_time_domain_summary.csv"
    with summary_csv.open("w", newline="", encoding="utf-8") as handle:
        fields = list(summary_rows[0].keys())
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(summary_rows)

    lines = [
        "# Hayward Axial Gravitational Time-Domain Benchmark",
        "",
        "Reference: Bolokhov & Skvortsova 2025, arXiv:2508.19989, Table V.",
        "",
        "Metric and potential:",
        "",
        "```text",
        "f(r) = 1 - 2 r^2 / (r^3 + gamma)",
        "V = f [2 f/r^2 - f'/r + ((ell+2)(ell-1))/r^2]",
        "```",
        "",
        "| gamma | r_h | r_peak | reference Re | TD Re | delta Re [%] | reference Abs(Im) | TD Abs(Im) | delta Abs(Im) [%] | window | residual |",
        "| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in summary_rows:
        lines.append(
            "| {gamma:.2f} | {horizon_outer:.6f} | {r_peak:.6f} | {reference_Momega_real:.6f} | {td_Momega_real:.6f} | {delta_real_pct:.3f} | {reference_Momega_imag_abs:.6f} | {td_Momega_imag_abs:.6f} | {delta_imag_pct:.3f} | {fit_t_start:.0f}-{fit_t_stop:.0f} | {relative_residual:.4f} |".format(
                **row
            )
        )
    lines.extend(
        [
            "",
            "Interpretation:",
            "",
            "- This benchmark reproduces the published fundamental-mode time-domain/Prony values for the Hayward axial gravitational potential.",
            "- The first overtone is intentionally left for a separate multi-exponential extraction layer.",
            "- The result is a static readiness/stress test, not a rotating GW250114 remnant constraint.",
            "",
        ]
    )
    report_path = out_dir / "hayward_time_domain_report.md"
    report_path.write_text("\n".join(lines), encoding="utf-8")

    print(f"Summary: {summary_csv}")
    print(f"Report: {report_path}")
    for row in summary_rows:
        print(
            "gamma={gamma:.2f}\tTD={td_Momega_real:.6f}-{td_Momega_imag_abs:.6f}i\t"
            "ref={reference_Momega_real:.6f}-{reference_Momega_imag_abs:.6f}i\t"
            "delta=({delta_real_pct:.2f}%, {delta_imag_pct:.2f}%)".format(**row)
        )


if __name__ == "__main__":
    main()
