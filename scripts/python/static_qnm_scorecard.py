"""Build a paper-facing scorecard for static supplied-potential QNM tests."""

from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def as_float(value: str) -> float:
    try:
        return float(value)
    except Exception:
        return math.nan


def verdict(delta_re: float, delta_im: float) -> str:
    max_abs = max(abs(delta_re), abs(delta_im))
    if max_abs <= 1.0:
        return "PASS:sub_percent_validation"
    if max_abs <= 3.0:
        return "WARN:few_percent_validation"
    return "FAIL:check_solver_or_reference"


def main() -> None:
    out_dir = ROOT / "results" / "static_qnm_scorecard"
    out_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, object]] = []

    static_path = ROOT / "results" / "static_master_potential_time_domain" / "time_domain_qnm_summary.csv"
    for row in read_csv(static_path):
        delta_re = as_float(row["delta_real_pct"])
        delta_im = as_float(row["delta_imag_pct"])
        rows.append(
            {
                "family": "Schwarzschild",
                "case": row["name"],
                "readiness_level": "validation_anchor",
                "sector": row["sector"],
                "ell": row["ell"],
                "n": "0",
                "parameter_name": "none",
                "parameter_value": "0",
                "reference_source": row["reference_source"],
                "reference_Momega_real": row["reference_Momega_real"],
                "td_Momega_real": row["td_Momega_real"],
                "delta_real_pct": delta_re,
                "reference_Momega_imag_abs": row["reference_Momega_imag_abs"],
                "td_Momega_imag_abs": row["td_Momega_imag_abs"],
                "delta_imag_pct": delta_im,
                "verdict": verdict(delta_re, delta_im),
            }
        )

    tidal_path = ROOT / "results" / "tidal_charge_time_domain" / "tidal_charge_time_domain_summary.csv"
    for row in read_csv(tidal_path):
        delta_re = as_float(row["delta_real_pct"])
        delta_im = as_float(row["delta_imag_pct"])
        rows.append(
            {
                "family": "Braneworld_tidal_charge",
                "case": row["case"],
                "readiness_level": "validated_supplied_gravitational_potential",
                "sector": "axial_gravitational",
                "ell": row["ell"],
                "n": row["n"],
                "parameter_name": "Qstar_over_M2",
                "parameter_value": row["q_tidal"],
                "reference_source": row["reference_source"],
                "reference_Momega_real": row["reference_Momega_real"],
                "td_Momega_real": row["td_Momega_real"],
                "delta_real_pct": delta_re,
                "reference_Momega_imag_abs": row["reference_Momega_imag_abs"],
                "td_Momega_imag_abs": row["td_Momega_imag_abs"],
                "delta_imag_pct": delta_im,
                "verdict": verdict(delta_re, delta_im),
            }
        )

    bardeen_path = ROOT / "results" / "bardeen_time_domain" / "bardeen_time_domain_summary.csv"
    for row in read_csv(bardeen_path):
        delta_re = as_float(row["delta_real_pct"])
        delta_im = as_float(row["delta_imag_pct"])
        rows.append(
            {
                "family": "Bardeen_NED",
                "case": row["case"],
                "readiness_level": "validated_source_dependent_gravitational_potential",
                "sector": "axial_gravitational",
                "ell": row["ell"],
                "n": row["n"],
                "parameter_name": "alpha",
                "parameter_value": row["alpha"],
                "reference_source": row["reference_source"],
                "reference_Momega_real": row["reference_Momega_real"],
                "td_Momega_real": row["td_Momega_real"],
                "delta_real_pct": delta_re,
                "reference_Momega_imag_abs": row["reference_Momega_imag_abs"],
                "td_Momega_imag_abs": row["td_Momega_imag_abs"],
                "delta_imag_pct": delta_im,
                "verdict": verdict(delta_re, delta_im),
            }
        )

    hayward_path = ROOT / "results" / "hayward_time_domain" / "hayward_time_domain_summary.csv"
    for row in read_csv(hayward_path):
        delta_re = as_float(row["delta_real_pct"])
        delta_im = as_float(row["delta_imag_pct"])
        rows.append(
            {
                "family": "Hayward",
                "case": row["case"],
                "readiness_level": "validated_axial_gravitational_potential_fundamental",
                "sector": "axial_gravitational",
                "ell": row["ell"],
                "n": row["n"],
                "parameter_name": "gamma",
                "parameter_value": row["gamma"],
                "reference_source": row["reference_source"],
                "reference_Momega_real": row["reference_Momega_real"],
                "td_Momega_real": row["td_Momega_real"],
                "delta_real_pct": delta_re,
                "reference_Momega_imag_abs": row["reference_Momega_imag_abs"],
                "td_Momega_imag_abs": row["td_Momega_imag_abs"],
                "delta_imag_pct": delta_im,
                "verdict": verdict(delta_re, delta_im),
            }
        )

    hayward_overtone_path = (
        ROOT
        / "results"
        / "hayward_overtone_matrix_pencil"
        / "hayward_overtone_matrix_pencil_summary.csv"
    )
    if hayward_overtone_path.exists():
        for row in read_csv(hayward_overtone_path):
            delta_re = as_float(row["delta_real_pct"])
            delta_im = as_float(row["delta_imag_pct"])
            rows.append(
                {
                    "family": "Hayward_overtone",
                    "case": row["case"],
                    "readiness_level": "validated_axial_gravitational_potential_first_overtone",
                    "sector": "axial_gravitational",
                    "ell": row["ell"],
                    "n": row["n"],
                    "parameter_name": "gamma",
                    "parameter_value": row["gamma"],
                    "reference_source": row["reference_source"],
                    "reference_Momega_real": row["reference_Momega_real"],
                    "td_Momega_real": row["td_Momega_real"],
                    "delta_real_pct": delta_re,
                    "reference_Momega_imag_abs": row["reference_Momega_imag_abs"],
                    "td_Momega_imag_abs": row["td_Momega_imag_abs"],
                    "delta_imag_pct": delta_im,
                    "verdict": verdict(delta_re, delta_im),
                }
            )

    fields = [
        "family",
        "case",
        "readiness_level",
        "sector",
        "ell",
        "n",
        "parameter_name",
        "parameter_value",
        "reference_source",
        "reference_Momega_real",
        "td_Momega_real",
        "delta_real_pct",
        "reference_Momega_imag_abs",
        "td_Momega_imag_abs",
        "delta_imag_pct",
        "verdict",
    ]
    csv_path = out_dir / "static_qnm_validation_summary.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    by_family: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        by_family[str(row["family"])].append(row)

    lines = [
        "# Static QNM Readiness Scorecard",
        "",
        "This table joins the local time-domain validations for supplied static master potentials.",
        "",
        "| family | rows | max abs delta Re [%] | max abs delta Abs(Im) [%] | verdicts |",
        "| --- | ---: | ---: | ---: | --- |",
    ]
    for family in sorted(by_family):
        fam_rows = by_family[family]
        max_re = max(abs(float(row["delta_real_pct"])) for row in fam_rows)
        max_im = max(abs(float(row["delta_imag_pct"])) for row in fam_rows)
        verdicts = ", ".join(sorted({str(row["verdict"]) for row in fam_rows}))
        lines.append(f"| {family} | {len(fam_rows)} | {max_re:.3f} | {max_im:.3f} | {verdicts} |")

    lines.extend(
        [
            "",
            "Interpretation:",
            "",
            "- `PASS:sub_percent_validation` means the local time-domain extraction reproduces the reference table to within 1 percent in both real frequency and damping rate.",
            "- These are static supplied-potential validations. They do not replace rotating-remnant GW250114 constraints.",
            "- The Hayward overtone rows use a multi-mode matrix-pencil extraction and should be read as an overtone-sensitivity check, not as a final precision spectroscopy solver.",
            "- The scorecard is intended as the paper-facing readiness layer for metrics often used in geodesic, QPO, and shadow phenomenology.",
            "",
        ]
    )
    report_path = out_dir / "static_qnm_scorecard.md"
    report_path.write_text("\n".join(lines), encoding="utf-8")

    print(f"Scorecard CSV: {csv_path}")
    print(f"Scorecard report: {report_path}")
    for family in sorted(by_family):
        fam_rows = by_family[family]
        max_re = max(abs(float(row["delta_real_pct"])) for row in fam_rows)
        max_im = max(abs(float(row["delta_imag_pct"])) for row in fam_rows)
        print(f"{family}: rows={len(fam_rows)} max_delta=({max_re:.3f}%, {max_im:.3f}%)")


if __name__ == "__main__":
    main()
