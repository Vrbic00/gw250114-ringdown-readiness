"""Time-domain axial gravitational benchmark for the Bardeen black hole.

Reference:
    Ulhoa, "On Quasinormal Modes for Gravitational Perturbations of Bardeen
    Black Hole", arXiv:1303.3143.

The metric is

    f(r) = 1 - 2 r^2 / (r^2 + alpha^2)^(3/2)

for m=1. The axial gravitational potential is the source-dependent
nonlinear-electrodynamics potential used in the paper,

    V = f [(l(l+1)+2(f-1))/r^2 + f'/r + f'' + 2 k L]

with k=8 pi and L = 3 alpha^2 / (r^2 + alpha^2)^(5/2). This is intentionally
not a metric-only Regge-Wheeler guess.
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


ALPHA_EXTREMAL = 4.0 / math.sqrt(27.0)

REFERENCE_L2_N0 = {
    0.0: (0.3731620888, 0.08921749033),
    0.3: (0.4066400107, 0.08797961013),
    ALPHA_EXTREMAL: (0.7828924561, 0.07622286290),
}


def f_lapse(r: float | np.ndarray, alpha: float) -> float | np.ndarray:
    return 1.0 - 2.0 * r**2 / (r * r + alpha * alpha) ** 1.5


def f_prime(r: float | np.ndarray, alpha: float) -> float | np.ndarray:
    s = r * r + alpha * alpha
    return 2.0 * r * (r * r - 2.0 * alpha * alpha) / s**2.5


def f_second(r: float | np.ndarray, alpha: float) -> float | np.ndarray:
    s = r * r + alpha * alpha
    return 2.0 * (-2.0 * r**4 + 11.0 * alpha * alpha * r * r - 2.0 * alpha**4) / s**3.5


def ned_lagrangian(r: float | np.ndarray, alpha: float) -> float | np.ndarray:
    if alpha == 0.0:
        return np.zeros_like(r, dtype=float) if isinstance(r, np.ndarray) else 0.0
    return 3.0 * alpha * alpha / (r * r + alpha * alpha) ** 2.5


def axial_potential_r(r: float | np.ndarray, alpha: float, ell: int = 2) -> float | np.ndarray:
    f = f_lapse(r, alpha)
    source_term = 16.0 * math.pi * ned_lagrangian(r, alpha)
    return f * (
        (ell * (ell + 1.0) + 2.0 * (f - 1.0)) / (r * r)
        + f_prime(r, alpha) / r
        + f_second(r, alpha)
        + source_term
    )


def outer_horizon(alpha: float) -> float:
    if alpha == 0.0:
        return 2.0
    if abs(alpha - ALPHA_EXTREMAL) < 1e-10:
        return math.sqrt(2.0) * alpha

    xs = np.linspace(1e-6, 10.0, 20000)
    fs = f_lapse(xs, alpha)
    brackets: list[tuple[float, float]] = []
    for left, right, f_left, f_right in zip(xs[:-1], xs[1:], fs[:-1], fs[1:]):
        if f_left == 0.0:
            brackets.append((float(left), float(left)))
        elif f_left * f_right < 0.0:
            brackets.append((float(left), float(right)))
    if not brackets:
        raise ValueError(f"no Bardeen horizon found for alpha={alpha}")

    low, high = brackets[-1]
    for _ in range(100):
        mid = 0.5 * (low + high)
        if f_lapse(low, alpha) * f_lapse(mid, alpha) <= 0.0:
            high = mid
        else:
            low = mid
    return 0.5 * (low + high)


def peak_radius(alpha: float, ell: int = 2) -> float:
    rh = outer_horizon(alpha)
    left = rh + (1e-5 if abs(alpha - ALPHA_EXTREMAL) < 1e-10 else 1e-7)
    right = 30.0
    inv_phi = (math.sqrt(5.0) - 1.0) / 2.0
    inv_phi2 = (3.0 - math.sqrt(5.0)) / 2.0
    h = right - left
    c = left + inv_phi2 * h
    d = left + inv_phi * h
    yc = axial_potential_r(c, alpha, ell)
    yd = axial_potential_r(d, alpha, ell)
    for _ in range(120):
        if yc < yd:
            left = c
            c = d
            yc = yd
            h = right - left
            d = left + inv_phi * h
            yd = axial_potential_r(d, alpha, ell)
        else:
            right = d
            d = c
            yd = yc
            h = right - left
            c = left + inv_phi2 * h
            yc = axial_potential_r(c, alpha, ell)
    return 0.5 * (left + right)


def shifted_potential_factory(alpha: float):
    rh = outer_horizon(alpha)
    r_peak = peak_radius(alpha)

    rho_min = 1e-6 if abs(alpha - ALPHA_EXTREMAL) < 1e-10 else 1e-8
    r_max = 700.0
    rho = np.geomspace(rho_min, r_max - rh, 240000)
    r_grid = rh + rho
    f_grid = f_lapse(r_grid, alpha)
    inv_f = 1.0 / f_grid
    dx = 0.5 * (inv_f[1:] + inv_f[:-1]) * np.diff(r_grid)
    x_raw = np.concatenate([[0.0], np.cumsum(dx)])
    x_peak = float(np.interp(r_peak, r_grid, x_raw))
    x_grid = x_raw - x_peak
    v_grid = axial_potential_r(r_grid, alpha)

    if x_grid[-1] < 300.0:
        raise ValueError(
            f"insufficient tortoise range for alpha={alpha}: {x_grid[0]} to {x_grid[-1]}"
        )

    def potential(x_values: np.ndarray) -> np.ndarray:
        return np.interp(x_values, x_grid, v_grid, left=0.0, right=v_grid[-1])

    return potential, rh, r_peak, x_grid[0], x_grid[-1]


def best_fit_for_waveform(times: np.ndarray, waveform: np.ndarray) -> dict[str, float]:
    windows = fit_windows(
        times,
        waveform,
        t_start_min=40.0,
        t_start_max=200.0,
        widths=(60.0, 70.0, 80.0, 100.0, 120.0),
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
    out_dir = ROOT / "results" / "bardeen_time_domain"
    out_dir.mkdir(parents=True, exist_ok=True)

    summary_rows: list[dict[str, float | str]] = []
    for alpha, (ref_real, ref_imag) in REFERENCE_L2_N0.items():
        potential, rh, r_peak, x_min, x_max = shifted_potential_factory(alpha)
        times, waveform = characteristic_evolve_with_potential(
            potential,
            h=0.2,
            n_grid=2600,
            observer_rstar=20.0,
            gaussian_center=0.0,
            gaussian_width=3.0,
        )
        tag = "extremal" if abs(alpha - ALPHA_EXTREMAL) < 1e-10 else f"{alpha:g}"
        np.savetxt(
            out_dir / f"bardeen_alpha_{tag}_axial_l2_waveform.csv",
            np.column_stack([times, waveform]),
            delimiter=",",
            header="t_over_M,psi",
            comments="",
        )
        best = best_fit_for_waveform(times, waveform)
        summary_rows.append(
            {
                "case": f"bardeen_alpha_{tag}_axial_l2_n0",
                "alpha": alpha,
                "ell": 2,
                "n": 0,
                "horizon_outer": rh,
                "r_peak": r_peak,
                "x_min": x_min,
                "x_max": x_max,
                "reference_source": "Ulhoa_2013_ThirdOrderWKB",
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

    summary_csv = out_dir / "bardeen_time_domain_summary.csv"
    with summary_csv.open("w", newline="", encoding="utf-8") as handle:
        fields = list(summary_rows[0].keys())
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(summary_rows)

    lines = [
        "# Bardeen Axial Gravitational Time-Domain Benchmark",
        "",
        "Reference: Ulhoa 2013, arXiv:1303.3143, Tables I-III.",
        "",
        "Metric and potential:",
        "",
        "```text",
        "f(r) = 1 - 2 r^2 / (r^2 + alpha^2)^(3/2)",
        "V = f [(l(l+1)+2(f-1))/r^2 + f'/r + f'' + 2 k L]",
        "k = 8 pi,  L = 3 alpha^2 / (r^2 + alpha^2)^(5/2)",
        "```",
        "",
        "| alpha | r_h | r_peak | reference Re | TD Re | delta Re [%] | reference Abs(Im) | TD Abs(Im) | delta Abs(Im) [%] | window | residual |",
        "| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in summary_rows:
        lines.append(
            "| {alpha:.9f} | {horizon_outer:.6f} | {r_peak:.6f} | {reference_Momega_real:.6f} | {td_Momega_real:.6f} | {delta_real_pct:.3f} | {reference_Momega_imag_abs:.6f} | {td_Momega_imag_abs:.6f} | {delta_imag_pct:.3f} | {fit_t_start:.0f}-{fit_t_stop:.0f} | {relative_residual:.4f} |".format(
                **row
            )
        )
    lines.extend(
        [
            "",
            "Interpretation:",
            "",
            "- This is a source-dependent nonlinear-electrodynamics Bardeen axial gravitational potential, not a metric-only proxy.",
            "- The comparison is against the paper's third-order WKB tables, so sub-percent agreement is not expected in every case.",
            "- The extremal case is numerically delicate because the horizon is degenerate and the tortoise coordinate is more singular.",
            "",
        ]
    )
    report_path = out_dir / "bardeen_time_domain_report.md"
    report_path.write_text("\n".join(lines), encoding="utf-8")

    print(f"Summary: {summary_csv}")
    print(f"Report: {report_path}")
    for row in summary_rows:
        print(
            "alpha={alpha:.9f}\tTD={td_Momega_real:.6f}-{td_Momega_imag_abs:.6f}i\t"
            "ref={reference_Momega_real:.6f}-{reference_Momega_imag_abs:.6f}i\t"
            "delta=({delta_real_pct:.2f}%, {delta_imag_pct:.2f}%)".format(**row)
        )


if __name__ == "__main__":
    main()
