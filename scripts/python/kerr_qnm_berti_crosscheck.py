r"""Pure-Python cross-check for the Wolfram Kerr QNM baseline.

This does not replace the numerical `qnm` package check. It independently
recomputes the same Berti-Cardoso-Will fitting formulas with Python stdlib
only and compares them against the Wolfram-generated CSV.

Run:
    C:\Users\vrb0015\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe scripts/python/kerr_qnm_berti_crosscheck.py
"""

from __future__ import annotations

import csv
import math
from pathlib import Path


SOLAR_MASS_TIME_SECONDS = 4.925490947e-6

EVENT = {
    "mass_detector_msun": 68.1,
    "mass_detector_minus": 0.9,
    "mass_detector_plus": 0.8,
    "spin": 0.68,
    "spin_minus": 0.01,
    "spin_plus": 0.01,
}

COEFFICIENTS = {
    "220": {"l": 2, "m": 2, "n": 0, "f1": 1.5251, "f2": -1.1568, "f3": 0.1292, "q1": 0.7000, "q2": 1.4187, "q3": -0.4990},
    "221": {"l": 2, "m": 2, "n": 1, "f1": 1.3673, "f2": -1.0260, "f3": 0.1628, "q1": 0.1000, "q2": 0.5436, "q3": -0.4731},
    "222": {"l": 2, "m": 2, "n": 2, "f1": 1.3223, "f2": -1.0257, "f3": 0.1860, "q1": -0.1000, "q2": 0.4206, "q3": -0.4256},
    "440": {"l": 4, "m": 4, "n": 0, "f1": 2.3000, "f2": -1.5056, "f3": 0.2244, "q1": 1.1929, "q2": 3.1191, "q3": -0.4825},
}


def qnm_values(mode: str, mass_solar: float, spin: float) -> dict[str, float | str | int]:
    coeff = COEFFICIENTS[mode]
    m_seconds = mass_solar * SOLAR_MASS_TIME_SECONDS
    momega = coeff["f1"] + coeff["f2"] * (1.0 - spin) ** coeff["f3"]
    omega_si = momega / m_seconds
    quality = coeff["q1"] + coeff["q2"] * (1.0 - spin) ** coeff["q3"]
    tau_s = 2.0 * quality / omega_si
    return {
        "mode": mode,
        "l": coeff["l"],
        "m": coeff["m"],
        "n": coeff["n"],
        "Momega": momega,
        "f_Hz": omega_si / (2.0 * math.pi),
        "Q": quality,
        "tau_ms": 1000.0 * tau_s,
    }


def qnm_with_uncertainty(mode: str) -> dict[str, float | str | int]:
    mass = EVENT["mass_detector_msun"]
    spin = EVENT["spin"]
    masses = [
        mass - EVENT["mass_detector_minus"],
        mass,
        mass + EVENT["mass_detector_plus"],
    ]
    spins = [
        max(0.0, spin - EVENT["spin_minus"]),
        spin,
        min(0.999, spin + EVENT["spin_plus"]),
    ]
    samples = [qnm_values(mode, m, s) for m in masses for s in spins]
    central = qnm_values(mode, mass, spin)
    f_values = [float(row["f_Hz"]) for row in samples]
    tau_values = [float(row["tau_ms"]) for row in samples]
    central.update(
        {
            "f_Hz_min": min(f_values),
            "f_Hz_max": max(f_values),
            "tau_ms_min": min(tau_values),
            "tau_ms_max": max(tau_values),
        }
    )
    return central


def read_wolfram_csv(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return {row["mode"]: row for row in csv.DictReader(handle)}


def max_abs_difference(py_row: dict[str, float | str | int], wl_row: dict[str, str]) -> float:
    fields = [
        "Momega",
        "f_Hz",
        "f_Hz_min",
        "f_Hz_max",
        "Q",
        "tau_ms",
        "tau_ms_min",
        "tau_ms_max",
    ]
    return max(abs(float(py_row[field]) - float(wl_row[field])) for field in fields)


def main() -> None:
    root = Path(__file__).resolve().parents[2]
    result_dir = root / "results" / "gw250114_kerr_qnm"
    wolfram_csv = result_dir / "gw250114_kerr_qnm_berti.csv"
    crosscheck_csv = result_dir / "python_berti_crosscheck.csv"
    crosscheck_report = result_dir / "python_crosscheck_report.md"

    wolfram_rows = read_wolfram_csv(wolfram_csv)
    python_rows = [qnm_with_uncertainty(mode) for mode in ("220", "221", "222", "440")]
    comparisons = [
        {
            "mode": row["mode"],
            "max_abs_difference": max_abs_difference(row, wolfram_rows[str(row["mode"])]),
        }
        for row in python_rows
    ]

    with crosscheck_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["mode", "max_abs_difference"])
        writer.writeheader()
        writer.writerows(comparisons)

    max_diff = max(row["max_abs_difference"] for row in comparisons)
    status = "PASS" if max_diff < 1e-10 else "CHECK"

    report = "\n".join(
        [
            "# Python Berti-Formula Cross-Check",
            "",
            "This is a pure-Python stdlib check of the same Berti-Cardoso-Will formulas used by the Wolfram script.",
            "It is not the independent `qnm` package calculation; that is handled by `scripts/python/qnm_solver_crosscheck.py`.",
            "",
            f"Status: `{status}`",
            f"Maximum absolute difference vs Wolfram CSV: `{max_diff:.3e}`",
            "",
            "| mode | max abs difference |",
            "| --- | ---: |",
            *[
                f"| {row['mode']} | {row['max_abs_difference']:.3e} |"
                for row in comparisons
            ],
            "",
            "The true Python `qnm` solver check is handled separately by `scripts/python/qnm_solver_crosscheck.py`.",
            "",
        ]
    )
    crosscheck_report.write_text(report, encoding="utf-8")

    print(f"Python cross-check status: {status}")
    print(f"Maximum absolute difference: {max_diff:.3e}")
    print(f"CSV: {crosscheck_csv}")
    print(f"Report: {crosscheck_report}")


if __name__ == "__main__":
    main()
