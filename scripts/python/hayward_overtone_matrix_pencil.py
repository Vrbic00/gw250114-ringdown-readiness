"""Matrix-pencil extraction of the Hayward first overtone.

This script extends the Hayward static benchmark to the first overtone
(`ell=2, n=1`). The source paper reports that the overtone is more sensitive
to the Hayward quantum parameter gamma. A single damped-sinusoid Prony fit is
not enough for this mode, so we use a rank-8 matrix-pencil extraction on a
fixed early ringdown window.

Reference values are the WKB-8 Padé entries in Table III of:
    Bolokhov & Skvortsova, arXiv:2508.19989.
"""

from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "python"))

from hayward_time_domain_benchmark import (  # noqa: E402
    REFERENCE_PRONY_L2_N0,
    shifted_potential_factory,
)
from static_master_potential_time_domain import (  # noqa: E402
    characteristic_evolve_with_potential,
    pct_delta,
)


REFERENCE_WKB8_L2_N1 = {
    0.00: (0.346003, 0.273555),
    0.02: (0.346515, 0.272975),
    0.04: (0.347024, 0.272300),
    0.06: (0.347609, 0.271405),
    0.08: (0.348888, 0.270209),
    0.10: (0.350735, 0.270779),
    1.18: (0.370157, 0.225937),
}


def matrix_pencil_modes(
    times: np.ndarray,
    waveform: np.ndarray,
    *,
    rank: int = 8,
    fit_t_start: float = 60.0,
    fit_t_stop: float = 120.0,
) -> list[dict[str, float]]:
    mask = (times >= fit_t_start) & (times <= fit_t_stop)
    t = times[mask]
    y = waveform[mask].astype(float).copy()
    y -= float(np.mean(y))
    scale = float(np.max(np.abs(y)))
    if scale <= 0.0 or not math.isfinite(scale):
        return []
    y /= scale

    n_samples = len(y)
    pencil_rows = n_samples // 2
    pencil_cols = n_samples - pencil_rows
    if min(pencil_rows, pencil_cols) <= rank + 2:
        return []

    y0 = np.column_stack([y[i : i + pencil_rows] for i in range(pencil_cols)])
    y1 = np.column_stack([y[i + 1 : i + pencil_rows + 1] for i in range(pencil_cols)])
    u, s, vh = np.linalg.svd(y0, full_matrices=False)
    ur = u[:, :rank]
    sr = s[:rank]
    vr = vh[:rank, :].conj().T
    reduced = (ur.conj().T @ y1 @ vr) @ np.diag(1.0 / sr)
    roots = np.linalg.eigvals(reduced)

    dt = float(np.median(np.diff(t)))
    rows: list[dict[str, float]] = []
    for root in roots:
        if abs(root) <= 0.0 or abs(root) >= 1.2:
            continue
        exponent = np.log(root) / dt
        omega = abs(float(np.imag(exponent)))
        gamma = -float(np.real(exponent))
        if omega > 0.02 and gamma > 0.0 and gamma < 1.0:
            rows.append(
                {
                    "td_Momega_real": omega,
                    "td_Momega_imag_abs": gamma,
                    "root_abs": abs(root),
                }
            )

    unique: list[dict[str, float]] = []
    for row in sorted(rows, key=lambda item: item["td_Momega_real"]):
        if not any(
            abs(row["td_Momega_real"] - old["td_Momega_real"]) < 1e-5
            and abs(row["td_Momega_imag_abs"] - old["td_Momega_imag_abs"]) < 1e-5
            for old in unique
        ):
            unique.append(row)
    return unique


def select_overtone(
    modes: list[dict[str, float]],
    fundamental_real: float,
    fundamental_imag_abs: float,
) -> dict[str, float] | None:
    candidates = [
        mode
        for mode in modes
        if 0.6 * fundamental_real < mode["td_Momega_real"] < 0.98 * fundamental_real
        and 1.8 * fundamental_imag_abs < mode["td_Momega_imag_abs"] < 5.0 * fundamental_imag_abs
    ]
    if not candidates:
        return None

    # Prefer the mode closest to the expected overtone damping hierarchy,
    # without using the reference overtone frequency itself.
    return min(
        candidates,
        key=lambda mode: abs(mode["td_Momega_imag_abs"] / fundamental_imag_abs - 3.0)
        + 0.25 * abs(mode["td_Momega_real"] / fundamental_real - 0.93),
    )


def main() -> None:
    out_dir = ROOT / "results" / "hayward_overtone_matrix_pencil"
    out_dir.mkdir(parents=True, exist_ok=True)

    all_mode_rows: list[dict[str, float | str]] = []
    summary_rows: list[dict[str, float | str]] = []
    for gamma, (ref_real, ref_imag) in REFERENCE_WKB8_L2_N1.items():
        potential, rh, r_peak, _x_min, _x_max = shifted_potential_factory(gamma)
        times, waveform = characteristic_evolve_with_potential(
            potential,
            h=0.2,
            n_grid=2700,
            observer_rstar=20.0,
            gaussian_center=0.0,
            gaussian_width=3.0,
        )
        modes = matrix_pencil_modes(times, waveform)

        # Use the published/literature-consistent fundamental scale only for
        # mode classification, not for fitting the overtone frequency.
        if gamma in REFERENCE_PRONY_L2_N0:
            fundamental_real, fundamental_imag = REFERENCE_PRONY_L2_N0[gamma]
        else:
            # For small gamma values not in Table V, interpolate from the
            # smooth WKB-8 fundamental trend in Table I.
            fundamental_real = 0.373669 + (0.375306 - 0.373669) * (gamma / 0.10)
            fundamental_imag = 0.088972 + (0.088227 - 0.088972) * (gamma / 0.10)

        for mode in modes:
            all_mode_rows.append(
                {
                    "gamma": gamma,
                    "fit_t_start": 60.0,
                    "fit_t_stop": 120.0,
                    "rank": 8,
                    **mode,
                }
            )

        overtone = select_overtone(modes, fundamental_real, fundamental_imag)
        if overtone is None:
            td_real = math.nan
            td_imag = math.nan
            delta_real = math.nan
            delta_imag = math.nan
            status = "FAIL:no_overtone_candidate"
        else:
            td_real = overtone["td_Momega_real"]
            td_imag = overtone["td_Momega_imag_abs"]
            delta_real = pct_delta(td_real, ref_real)
            delta_imag = pct_delta(td_imag, ref_imag)
            max_delta = max(abs(delta_real), abs(delta_imag))
            status = "PASS:sub_percent_overtone" if max_delta <= 1.0 else "WARN:few_percent_overtone"

        summary_rows.append(
            {
                "case": f"hayward_gamma_{gamma:g}_axial_l2_n1",
                "gamma": gamma,
                "ell": 2,
                "n": 1,
                "horizon_outer": rh,
                "r_peak": r_peak,
                "fit_t_start": 60.0,
                "fit_t_stop": 120.0,
                "rank": 8,
                "reference_source": "Bolokhov_Skvortsova_2025_Table_III_WKB8",
                "reference_Momega_real": ref_real,
                "reference_Momega_imag_abs": ref_imag,
                "td_Momega_real": td_real,
                "td_Momega_imag_abs": td_imag,
                "delta_real_pct": delta_real,
                "delta_imag_pct": delta_imag,
                "status": status,
            }
        )

    modes_csv = out_dir / "hayward_overtone_matrix_pencil_modes.csv"
    with modes_csv.open("w", newline="", encoding="utf-8") as handle:
        fields = list(all_mode_rows[0].keys())
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(all_mode_rows)

    summary_csv = out_dir / "hayward_overtone_matrix_pencil_summary.csv"
    with summary_csv.open("w", newline="", encoding="utf-8") as handle:
        fields = list(summary_rows[0].keys())
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(summary_rows)

    lines = [
        "# Hayward First-Overtone Matrix-Pencil Check",
        "",
        "Reference: Bolokhov & Skvortsova 2025, arXiv:2508.19989, Table III WKB-8 entries.",
        "",
        "Extraction settings: rank-8 matrix pencil, fixed window `60 <= t/M <= 120`, `h=0.2 M`.",
        "",
        "| gamma | reference Re | TD Re | delta Re [%] | reference Abs(Im) | TD Abs(Im) | delta Abs(Im) [%] | status |",
        "| ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for row in summary_rows:
        lines.append(
            "| {gamma:.2f} | {reference_Momega_real:.6f} | {td_Momega_real:.6f} | {delta_real_pct:.3f} | {reference_Momega_imag_abs:.6f} | {td_Momega_imag_abs:.6f} | {delta_imag_pct:.3f} | {status} |".format(
                **row
            )
        )
    lines.extend(
        [
            "",
            "Interpretation:",
            "",
            "- The first overtone is extracted from the same time-domain waveform as the fundamental mode, but requires a multi-mode matrix-pencil fit.",
            "- The selection rule uses the fundamental-mode damping hierarchy, not the reference overtone frequency itself.",
            "- This is a useful validation of overtone sensitivity, but a final precision overtone solver should still use Leaver, high-order WKB/Pade, or a more controlled multi-mode fit.",
            "",
        ]
    )
    report_path = out_dir / "hayward_overtone_matrix_pencil_report.md"
    report_path.write_text("\n".join(lines), encoding="utf-8")

    print(f"Overtone summary: {summary_csv}")
    print(f"Overtone modes: {modes_csv}")
    print(f"Report: {report_path}")
    for row in summary_rows:
        print(
            "gamma={gamma:.2f}\tTD={td_Momega_real:.6f}-{td_Momega_imag_abs:.6f}i\t"
            "ref={reference_Momega_real:.6f}-{reference_Momega_imag_abs:.6f}i\t"
            "delta=({delta_real_pct:.2f}%, {delta_imag_pct:.2f}%)\t{status}".format(
                **row
            )
        )


if __name__ == "__main__":
    main()
