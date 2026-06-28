"""Time-domain axial gravitational benchmark for braneworld tidal charge.

Reference:
    Toshmatov et al., "Quasinormal frequencies of black hole in the
    braneworld", arXiv:1605.02058, Table II.

The metric is

    f(r) = 1 - 2/r - q/r^2

where q = Q*/M^2 > 0 is the positive tidal-charge parameter used in that
paper. The axial gravitational supplied potential is

    V = f [l(l+1)/r^2 - 2(3r + 2q)/r^4]

for M=1. The tortoise origin is shifted so that the potential peak sits near
r_* = 0 for each q; this keeps the generic time-domain extraction stable.
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


REFERENCE_AXIAL_L2_N0 = {
    0.1: (0.3673, 0.0883),
    0.4: (0.3508, 0.0867),
    0.7: (0.3370, 0.0850),
    1.0: (0.3250, 0.0835),
    2.0: (0.2944, 0.0788),
}


def roots(q: float) -> tuple[float, float]:
    s = math.sqrt(1.0 + q)
    return 1.0 + s, 1.0 - s


def f_lapse(r: float, q: float) -> float:
    return 1.0 - 2.0 / r - q / (r * r)


def axial_potential_r(r: float, q: float, ell: int = 2) -> float:
    f = f_lapse(r, q)
    return f * (ell * (ell + 1.0) / (r * r) - 2.0 * (3.0 * r + 2.0 * q) / (r**4))


def tortoise_raw(r: float, q: float) -> float:
    rp, rm = roots(q)
    denom = rp - rm
    a_plus = rp * rp / denom
    a_minus = -rm * rm / denom
    return r + a_plus * math.log(r - rp) + a_minus * math.log(r - rm)


def peak_radius(q: float, ell: int = 2) -> float:
    rp, _rm = roots(q)
    left = rp + 1e-7
    right = 30.0
    inv_phi = (math.sqrt(5.0) - 1.0) / 2.0
    inv_phi2 = (3.0 - math.sqrt(5.0)) / 2.0

    h = right - left
    c = left + inv_phi2 * h
    d = left + inv_phi * h
    yc = axial_potential_r(c, q, ell)
    yd = axial_potential_r(d, q, ell)
    for _ in range(120):
        if yc < yd:
            left = c
            c = d
            yc = yd
            h = right - left
            d = left + inv_phi * h
            yd = axial_potential_r(d, q, ell)
        else:
            right = d
            d = c
            yd = yc
            h = right - left
            c = left + inv_phi2 * h
            yc = axial_potential_r(c, q, ell)
    return 0.5 * (left + right)


def invert_tortoise(raw_x: float, q: float) -> float:
    rp, _rm = roots(q)
    low = rp + 1e-13
    high = max(rp + 2.0, raw_x + 10.0)

    while tortoise_raw(high, q) < raw_x:
        high *= 2.0

    for _ in range(90):
        mid = 0.5 * (low + high)
        if tortoise_raw(mid, q) < raw_x:
            low = mid
        else:
            high = mid
    return 0.5 * (low + high)


def shifted_potential_factory(q: float):
    r_peak = peak_radius(q)
    x_peak = tortoise_raw(r_peak, q)

    def potential(x_values: np.ndarray) -> np.ndarray:
        values = []
        for x in x_values:
            r = invert_tortoise(float(x) + x_peak, q)
            values.append(axial_potential_r(r, q))
        return np.array(values, dtype=float)

    return potential, r_peak, x_peak


def best_fit_for_waveform(times: np.ndarray, waveform: np.ndarray) -> dict[str, float]:
    windows = fit_windows(
        times,
        waveform,
        t_start_min=40.0,
        t_start_max=180.0,
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
    out_dir = ROOT / "results" / "tidal_charge_time_domain"
    out_dir.mkdir(parents=True, exist_ok=True)

    summary_rows: list[dict[str, float | str]] = []
    for q, (ref_real, ref_imag) in REFERENCE_AXIAL_L2_N0.items():
        potential, r_peak, x_peak = shifted_potential_factory(q)
        times, waveform = characteristic_evolve_with_potential(
            potential,
            h=0.2,
            n_grid=2400,
            observer_rstar=20.0,
            gaussian_center=0.0,
            gaussian_width=3.0,
        )
        np.savetxt(
            out_dir / f"tidal_charge_q{q:g}_axial_l2_waveform.csv",
            np.column_stack([times, waveform]),
            delimiter=",",
            header="t_over_M,psi",
            comments="",
        )
        best = best_fit_for_waveform(times, waveform)
        summary_rows.append(
            {
                "case": f"tidal_charge_q{q:g}_axial_l2_n0",
                "q_tidal": q,
                "ell": 2,
                "n": 0,
                "r_peak": r_peak,
                "x_peak_raw": x_peak,
                "reference_source": "Toshmatov_2016_Table_II_6th_order_WKB",
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

    summary_csv = out_dir / "tidal_charge_time_domain_summary.csv"
    with summary_csv.open("w", newline="", encoding="utf-8") as handle:
        fields = list(summary_rows[0].keys())
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(summary_rows)

    lines = [
        "# Braneworld Tidal-Charge Time-Domain Benchmark",
        "",
        "Reference: Toshmatov et al. 2016, arXiv:1605.02058, Table II.",
        "",
        "Metric and potential:",
        "",
        "```text",
        "f(r) = 1 - 2/r - q/r^2",
        "V_axial = f [l(l+1)/r^2 - 2(3r + 2q)/r^4]",
        "```",
        "",
        "The tortoise coordinate is shifted so that the axial potential peak is near `r* = 0` for each `q`.",
        "",
        "| q=Q*/M^2 | r_peak | reference Re | TD Re | delta Re [%] | reference Abs(Im) | TD Abs(Im) | delta Abs(Im) [%] | window | residual |",
        "| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in summary_rows:
        lines.append(
            "| {q_tidal:.1f} | {r_peak:.6f} | {reference_Momega_real:.6f} | {td_Momega_real:.6f} | {delta_real_pct:.3f} | {reference_Momega_imag_abs:.6f} | {td_Momega_imag_abs:.6f} | {delta_imag_pct:.3f} | {fit_t_start:.0f}-{fit_t_stop:.0f} | {relative_residual:.4f} |".format(
                **row
            )
        )
    lines.extend(
        [
            "",
            "Interpretation:",
            "",
            "- This is the first non-Schwarzschild supplied gravitational-potential benchmark in the static branch.",
            "- The comparison is against the paper's sixth-order WKB table, not an exact Leaver spectrum.",
            "- Agreement at the percent level is sufficient for adopting the case as a validation/stress-test target before moving to Bardeen and Hayward.",
            "",
        ]
    )
    report_path = out_dir / "tidal_charge_time_domain_report.md"
    report_path.write_text("\n".join(lines), encoding="utf-8")

    print(f"Summary: {summary_csv}")
    print(f"Report: {report_path}")
    for row in summary_rows:
        print(
            "q={q_tidal:.1f}\tTD={td_Momega_real:.6f}-{td_Momega_imag_abs:.6f}i\t"
            "ref={reference_Momega_real:.6f}-{reference_Momega_imag_abs:.6f}i\t"
            "delta=({delta_real_pct:.2f}%, {delta_imag_pct:.2f}%)".format(**row)
        )


if __name__ == "__main__":
    main()
