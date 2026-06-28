r"""Schwarzschild gravitational QNM reference from the `qnm` package.

This gives a reference value for the l=2, n=0 gravitational mode at a=0,
to compare against the first-order WKB Regge-Wheeler/Zerilli estimates.

Run:
    C:\Users\vrb0015\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe scripts/python/schwarzschild_qnm_reference.py
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / ".deps" / "python"))

import qnm  # noqa: E402


def pct_difference(value: float, reference: float) -> float:
    return 100.0 * (value - reference) / reference


def main() -> None:
    result_dir = ROOT / "results" / "master_potential_wkb"
    wkb_csv = result_dir / "master_potential_wkb.csv"
    out_report = result_dir / "schwarzschild_qnm_reference.md"

    seq = qnm.modes_cache(s=-2, l=2, m=2, n=0)
    omega, _angular_constant, _mixing = seq(a=0.0)
    exact_real = float(omega.real)
    exact_imag_abs = abs(float(omega.imag))

    with wkb_csv.open(newline="", encoding="utf-8") as handle:
        rows = {row["name"]: row for row in csv.DictReader(handle)}

    comparison_rows = []
    for name in ("Schwarzschild_ReggeWheeler_l2", "Schwarzschild_Zerilli_l2"):
        row = rows[name]
        wkb_real = float(row["WKB_Momega_real"])
        wkb_imag_abs = float(row["WKB_Momega_imag_abs"])
        comparison_rows.append(
            {
                "name": name,
                "wkb_real": wkb_real,
                "wkb_imag_abs": wkb_imag_abs,
                "delta_real_pct": pct_difference(wkb_real, exact_real),
                "delta_imag_pct": pct_difference(wkb_imag_abs, exact_imag_abs),
            }
        )

    lines = [
        "# Schwarzschild Gravitational QNM Reference",
        "",
        "Reference computed with the Python `qnm` package at `a=0` for `s=-2, l=2, m=2, n=0`.",
        "",
        f"- `Re(M omega) = {exact_real:.12f}`",
        f"- `Abs(Im(M omega)) = {exact_imag_abs:.12f}`",
        "",
        "| WKB potential | Re(M omega) WKB | delta Re [%] | Abs(Im) WKB | delta Abs(Im) [%] |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for row in comparison_rows:
        lines.append(
            "| {name} | {wkb_real:.12f} | {delta_real_pct:.3f} | {wkb_imag_abs:.12f} | {delta_imag_pct:.3f} |".format(
                **row
            )
        )
    lines.extend(
        [
            "",
            "Interpretation:",
            "",
            "- RW and Zerilli WKB estimates agree with each other, as expected for the Schwarzschild gravitational sectors.",
            "- First-order WKB at `l=2` is only a screening approximation; the real part is several percent high relative to the `qnm` reference.",
            "- A higher-order WKB or direct QNM solver is needed before making precision claims.",
            "",
        ]
    )
    out_report.write_text("\n".join(lines), encoding="utf-8")

    print(f"Schwarzschild qnm reference report: {out_report}")
    print(f"qnm Re(M omega): {exact_real:.12f}")
    print(f"qnm Abs(Im(M omega)): {exact_imag_abs:.12f}")


if __name__ == "__main__":
    main()
