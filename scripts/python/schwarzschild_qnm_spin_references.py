r"""Schwarzschild spin-s QNM references from the `qnm` package.

This compares first- and third-order WKB values for scalar, electromagnetic, and
gravitational Schwarzschild potentials against cached `qnm` references at
`a=0`. It is a Level 5 benchmark for known perturbation equations, not a
generic alternative-gravity solver.

Run:
    C:\Users\vrb0015\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe scripts/python/schwarzschild_qnm_spin_references.py
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / ".deps" / "python"))

import qnm  # noqa: E402


CASES = [
    {
        "name": "Schwarzschild_scalar_l2",
        "spin_weight": 0,
        "l": 2,
        "m": 0,
        "n": 0,
        "wkb_row": "Schwarzschild_scalar_l2",
    },
    {
        "name": "Schwarzschild_EM_l2",
        "spin_weight": -1,
        "l": 2,
        "m": 1,
        "n": 0,
        "wkb_row": "Schwarzschild_EM_l2",
    },
    {
        "name": "Schwarzschild_gravitational_l2",
        "spin_weight": -2,
        "l": 2,
        "m": 2,
        "n": 0,
        "wkb_row": "Schwarzschild_ReggeWheeler_l2",
    },
]


def pct_difference(value: float, reference: float) -> float:
    return 100.0 * (value - reference) / reference


def read_wkb_rows(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return {row["name"]: row for row in csv.DictReader(handle)}


def main() -> None:
    result_dir = ROOT / "results" / "master_potential_wkb"
    wkb_rows = read_wkb_rows(result_dir / "master_potential_wkb.csv")
    out_csv = result_dir / "schwarzschild_qnm_spin_references.csv"
    out_report = result_dir / "schwarzschild_qnm_spin_references.md"

    rows = []
    for case in CASES:
        seq = qnm.modes_cache(
            s=case["spin_weight"],
            l=case["l"],
            m=case["m"],
            n=case["n"],
        )
        omega, _angular_constant, _mixing = seq(a=0.0)
        reference_real = float(omega.real)
        reference_imag_abs = abs(float(omega.imag))
        wkb = wkb_rows[case["wkb_row"]]
        wkb1_real = float(wkb.get("WKB1_Momega_real", wkb["WKB_Momega_real"]))
        wkb1_imag_abs = float(wkb.get("WKB1_Momega_imag_abs", wkb["WKB_Momega_imag_abs"]))
        wkb3_real = float(wkb.get("WKB3_Momega_real", wkb["WKB_Momega_real"]))
        wkb3_imag_abs = float(wkb.get("WKB3_Momega_imag_abs", wkb["WKB_Momega_imag_abs"]))
        rows.append(
            {
                "name": case["name"],
                "spin_weight": case["spin_weight"],
                "l": case["l"],
                "m": case["m"],
                "n": case["n"],
                "qnm_real": reference_real,
                "qnm_imag_abs": reference_imag_abs,
                "wkb_row": case["wkb_row"],
                "wkb1_real": wkb1_real,
                "wkb1_imag_abs": wkb1_imag_abs,
                "wkb1_delta_real_pct": pct_difference(wkb1_real, reference_real),
                "wkb1_delta_imag_pct": pct_difference(wkb1_imag_abs, reference_imag_abs),
                "wkb3_real": wkb3_real,
                "wkb3_imag_abs": wkb3_imag_abs,
                "wkb3_delta_real_pct": pct_difference(wkb3_real, reference_real),
                "wkb3_delta_imag_pct": pct_difference(wkb3_imag_abs, reference_imag_abs),
            }
        )

    fields = [
        "name",
        "spin_weight",
        "l",
        "m",
        "n",
        "qnm_real",
        "qnm_imag_abs",
        "wkb_row",
        "wkb1_real",
        "wkb1_imag_abs",
        "wkb1_delta_real_pct",
        "wkb1_delta_imag_pct",
        "wkb3_real",
        "wkb3_imag_abs",
        "wkb3_delta_real_pct",
        "wkb3_delta_imag_pct",
    ]
    with out_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    report = [
        "# Schwarzschild Spin-s QNM References",
        "",
        "References computed with the Python `qnm` package at `a=0`.",
        "These are benchmark values for known Schwarzschild perturbation problems.",
        "",
        "| case | s | qnm Re(M omega) | WKB Re(M omega) | delta Re [%] | qnm Abs(Im) | WKB Abs(Im) | delta Abs(Im) [%] |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in rows:
        report.append(
            "| {name} | {spin_weight} | {qnm_real:.12f} | {wkb1_real:.12f} | {wkb1_delta_real_pct:.3f} | {qnm_imag_abs:.12f} | {wkb1_imag_abs:.12f} | {wkb1_delta_imag_pct:.3f} |".format(
                **row
            )
        )
    report.extend(
        [
            "",
            "## Third-Order WKB",
            "",
            "| case | s | qnm Re(M omega) | WKB3 Re(M omega) | delta Re [%] | qnm Abs(Im) | WKB3 Abs(Im) | delta Abs(Im) [%] |",
            "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    for row in rows:
        report.append(
            "| {name} | {spin_weight} | {qnm_real:.12f} | {wkb3_real:.12f} | {wkb3_delta_real_pct:.3f} | {qnm_imag_abs:.12f} | {wkb3_imag_abs:.12f} | {wkb3_delta_imag_pct:.3f} |".format(
                **row
            )
        )
    report.extend(
        [
            "",
            "Interpretation:",
            "",
            "- This is the current Level 5 benchmark set for Schwarzschild.",
            "- First-order WKB is useful for screening but noticeably overestimates the real frequency for low `l=2` modes.",
            "- Third-order WKB removes most of that real-frequency bias in these Schwarzschild benchmark cases.",
            "- The next numerical upgrade should be a direct solver or a higher-order WKB implementation for supplied potentials.",
            "",
        ]
    )
    out_report.write_text("\n".join(report), encoding="utf-8")

    print(f"Spin-s reference CSV: {out_csv}")
    print(f"Spin-s reference report: {out_report}")
    for row in rows:
        print(
            "{name}: qnm={qnm_real:.12f}-{qnm_imag_abs:.12f}i, WKB1 delta Re={wkb1_delta_real_pct:.3f}%, WKB3 delta Re={wkb3_delta_real_pct:.3f}%".format(
                **row
            )
        )


if __name__ == "__main__":
    main()
