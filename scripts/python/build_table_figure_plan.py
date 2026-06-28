"""Build main/supplement table and figure plan for the manuscript package."""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = ROOT / "results" / "manuscript_package"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def as_float(value: object, default: float = 0.0) -> float:
    try:
        return float(value)
    except Exception:
        return default


def as_int(value: object, default: int = 0) -> int:
    try:
        return int(float(str(value)))
    except Exception:
        return default


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT)).replace("\\", "/")


def build_public_projection_main_table() -> list[dict[str, object]]:
    summary_rows = read_csv(
        ROOT / "results" / "gw250114_constraints_comparison" / "projection_constraints_summary.csv"
    )
    consistency = read_csv(
        ROOT / "results" / "gw250114_constraints_comparison" / "projection_consistency_by_operator.csv"
    )
    max_difference_row = max(consistency, key=lambda row: as_float(row["normalized_projection_difference"]))
    max_difference = as_float(max_difference_row["normalized_projection_difference"])

    rows: list[dict[str, object]] = []
    for row in summary_rows:
        rows.append(
            {
                "projection": row["projection"],
                "observable_set": (
                    "{log f_220, log f_221, df_221}"
                    if row["projection"] == "RINGDOWN"
                    else "{log(1+domega_221), log(1+dtau_221)}"
                ),
                "tested_rows": as_int(row["rows"]),
                "zero_outside_90pct": as_int(row["zero_outside_90pct_count"]),
                "max_sigma_from_zero": round(as_float(row["max_sigma_from_zero"]), 3),
                "max_sigma_row": f"{row['max_sigma_operator']} {row['max_sigma_polarization']}",
                "tightest_alpha_sigma": as_float(row["tightest_alpha_sigma"]),
                "note": "Public-product projection; not a strain-level EFT likelihood.",
            }
        )

    rows.append(
        {
            "projection": "RINGDOWN_vs_PYRING",
            "observable_set": "operator-by-operator comparison",
            "tested_rows": len(consistency),
            "zero_outside_90pct": 0,
            "max_sigma_from_zero": round(max_difference, 3),
            "max_sigma_row": f"{max_difference_row['operator']} {max_difference_row['polarization']}",
            "tightest_alpha_sigma": "",
            "note": "Maximum normalized projection difference; branches compared, not combined.",
        }
    )
    return rows


def scorecard_by_family() -> dict[str, dict[str, object]]:
    rows = read_csv(ROOT / "results" / "static_qnm_scorecard" / "static_qnm_validation_summary.csv")
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row["family"]].append(row)
    out: dict[str, dict[str, object]] = {}
    for family, fam_rows in grouped.items():
        out[family] = {
            "validation_rows": len(fam_rows),
            "max_validation_delta_pct": max(
                max(abs(as_float(row["delta_real_pct"])), abs(as_float(row["delta_imag_pct"])))
                for row in fam_rows
            ),
            "modes": "; ".join(sorted({f"l={row['ell']},n={row['n']}" for row in fam_rows})),
        }
    return out


def physical_by_family() -> dict[str, dict[str, object]]:
    rows = read_csv(
        ROOT / "results" / "static_qnm_physical_deviations" / "static_qnm_physical_deviations.csv"
    )
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    family_map = {
        "Hayward_fundamental": "Hayward",
        "Hayward_overtone": "Hayward_overtone",
        "Bardeen_NED": "Bardeen_NED",
        "Braneworld_tidal_charge": "Braneworld_tidal_charge",
    }
    for row in rows:
        family = family_map.get(row["family"])
        if family:
            grouped[family].append(row)
    out: dict[str, dict[str, object]] = {}
    for family, fam_rows in grouped.items():
        max_row = max(fam_rows, key=lambda row: as_float(row["max_abs_physical_delta_pct"]))
        out[family] = {
            "max_physical_delta_pct": as_float(max_row["max_abs_physical_delta_pct"]),
            "max_physical_parameter": f"{max_row['parameter_name']}={as_float(max_row['parameter_value']):g}",
            "physical_baseline": max_row["baseline_label"],
        }
    return out


def readiness_by_candidate() -> dict[str, dict[str, str]]:
    rows = read_csv(ROOT / "results" / "static_qnm_readiness_audit" / "static_metric_readiness_audit.csv")
    return {row["candidate"]: row for row in rows}


def build_static_main_table() -> list[dict[str, object]]:
    scorecard = scorecard_by_family()
    physical = physical_by_family()
    family_to_candidate = {
        "Schwarzschild": "schwarzschild_rw_zerilli",
        "Braneworld_tidal_charge": "braneworld_tidal_charge_rn_like",
        "Bardeen_NED": "bardeen_ned_axial",
        "Hayward": "hayward_axial_gravitational",
        "Hayward_overtone": "hayward_axial_gravitational",
    }
    readiness = readiness_by_candidate()
    rows: list[dict[str, object]] = []
    for family in sorted(scorecard):
        candidate = family_to_candidate.get(family, "")
        ready_row = readiness.get(candidate, {})
        phys = physical.get(family, {})
        rows.append(
            {
                "family": family,
                "readiness_status": ready_row.get("audit_status", "validation_anchor"),
                "validated_modes": scorecard[family]["modes"],
                "validation_rows": scorecard[family]["validation_rows"],
                "max_validation_delta_pct": round(as_float(scorecard[family]["max_validation_delta_pct"]), 3),
                "max_sampled_physical_delta_pct": (
                    round(as_float(phys["max_physical_delta_pct"]), 3) if phys else ""
                ),
                "largest_sampled_parameter": phys.get("max_physical_parameter", "not_applicable"),
                "interpretation_boundary": "Static supplied-potential benchmark; not a GW250114 rotating-remnant exclusion.",
            }
        )
    return rows


def build_decisions() -> list[dict[str, object]]:
    return [
        {
            "item": "Main Table 1",
            "placement": "main_text",
            "source_path": rel(OUT_DIR / "main_table1_public_projection_summary.csv"),
            "caption": (
                "Public GW250114 projected higher-derivative constraints. The table reports separate "
                "RINGDOWN and pyRing projections and their maximum normalized comparison difference; "
                "the two public-product branches are compared but not statistically combined."
            ),
            "referee_guardrail": "Use projection/consistency language, not combined bound or strain-level likelihood.",
        },
        {
            "item": "Main Table 2",
            "placement": "main_text",
            "source_path": rel(OUT_DIR / "main_table2_static_readiness_summary.csv"),
            "caption": (
                "Static supplied-potential QNM readiness summary. Validation errors quantify "
                "reproduction of published reference spectra, while physical deviations quantify "
                "sampled shifts from Schwarzschild or zero-parameter baselines."
            ),
            "referee_guardrail": "State that these are static benchmarks, not observational exclusions of rotating remnants.",
        },
        {
            "item": "Supplement Table S1",
            "placement": "supplement",
            "source_path": "results/gw250114_paper_tables/table1_main_projected_constraints.csv",
            "caption": "Full row-level public-product projected EFT intervals.",
            "referee_guardrail": "Machine-readable detail supporting Main Table 1.",
        },
        {
            "item": "Supplement Table S2",
            "placement": "supplement",
            "source_path": "results/gw250114_paper_tables/table2_pipeline_consistency.csv",
            "caption": "Operator-by-operator RINGDOWN versus pyRing comparison.",
            "referee_guardrail": "Comparison only; no statistical combination.",
        },
        {
            "item": "Supplement Table S3",
            "placement": "supplement",
            "source_path": "results/gw250114_paper_tables/table3_pyring_filter_robustness.csv",
            "caption": "pyRing lower-tail filter robustness scenarios.",
            "referee_guardrail": "Robustness diagnostic, not independent evidence.",
        },
        {
            "item": "Supplement Table S4",
            "placement": "supplement",
            "source_path": "results/gw250114_paper_tables/table4_linearized_posterior_check.csv",
            "caption": "Empirical linearized posterior-sample projection check.",
            "referee_guardrail": "Non-Gaussian sanity check, not a separate likelihood.",
        },
        {
            "item": "Supplement Table S5",
            "placement": "supplement",
            "source_path": "results/gw250114_paper_tables/table5_public_ringdown_observables.csv",
            "caption": "Public ringdown observables used in the projection.",
            "referee_guardrail": "Defines inputs; keep variable naming tied to public products.",
        },
        {
            "item": "Supplement Table S6",
            "placement": "supplement",
            "source_path": "results/static_qnm_scorecard/static_qnm_validation_summary.csv",
            "caption": "Full static supplied-potential validation scorecard.",
            "referee_guardrail": "Validation against reference tables at their approximation level.",
        },
        {
            "item": "Supplement Table S7",
            "placement": "supplement",
            "source_path": "results/static_qnm_readiness_audit/static_metric_readiness_audit.csv",
            "caption": "Candidate metric readiness audit and negative-control classifications.",
            "referee_guardrail": "Classification standard, not a complete literature census.",
        },
        {
            "item": "Supplement Table S8",
            "placement": "supplement",
            "source_path": "results/static_qnm_physical_deviations/static_qnm_physical_thresholds.csv",
            "caption": "Sparse-grid static QNM physical-deviation threshold crossings.",
            "referee_guardrail": "Diagnostic thresholds only; no observational exclusion intervals.",
        },
        {
            "item": "Figure 1",
            "placement": "main_text",
            "source_path": "results/gw250114_constraints_comparison/projection_sigma_from_zero_comparison.png",
            "caption": (
                "Nominal displacement of one-at-a-time public-product projected couplings from alpha=0. "
                "Both RINGDOWN and pyRing branches retain alpha=0 within the configured 90 percent intervals."
            ),
            "referee_guardrail": "Do not imply detection or independent combination.",
        },
        {
            "item": "Figure 2",
            "placement": "main_text",
            "source_path": "results/static_qnm_readiness_audit/static_metric_readiness_ladder.svg",
            "caption": (
                "Ringdown-readiness ladder for metric models. A line element or geodesic observable is "
                "not sufficient for gravitational ringdown use without perturbation physics and a validated spectrum."
            ),
            "referee_guardrail": "This figure carries the conceptual unification of the static branch.",
        },
        {
            "item": "Figure 3",
            "placement": "main_text_or_short_appendix",
            "source_path": "results/static_qnm_readiness_audit/static_qnm_validation_scorecard.svg",
            "caption": (
                "Static supplied-potential validation errors for the implemented benchmark families. "
                "Sub-percent reproduction validates the local implementation, not the physical realism of the metric."
            ),
            "referee_guardrail": "Keep wording on numerical validation, not observational agreement.",
        },
        {
            "item": "Figure S1",
            "placement": "supplement",
            "source_path": "results/static_qnm_physical_deviations/static_qnm_physical_deviations.svg",
            "caption": (
                "Sampled physical QNM deviations from Schwarzschild or zero-parameter baselines. "
                "The large shifts show that the static tests can be discriminating once a physical tolerance is declared."
            ),
            "referee_guardrail": "Diagnostic only; physical parameter priors are required before constraints.",
        },
    ]


def write_markdown(
    path: Path,
    public_rows: list[dict[str, object]],
    static_rows: list[dict[str, object]],
    decisions: list[dict[str, object]],
) -> None:
    lines = [
        "# Table And Figure Plan",
        "",
        "This plan is deliberately conservative. The main text should carry only the compact results needed for the central claim; row-level and diagnostic material goes to the supplement.",
        "",
        "## Main Table 1: Public Projection Summary",
        "",
        "| projection | observable set | rows | zero outside 90% | max sigma/statistic | row | note |",
        "| --- | --- | ---: | ---: | ---: | --- | --- |",
    ]
    for row in public_rows:
        lines.append(
            f"| {row['projection']} | {row['observable_set']} | {row['tested_rows']} | "
            f"{row['zero_outside_90pct']} | {row['max_sigma_from_zero']} | "
            f"{row['max_sigma_row']} | {row['note']} |"
        )

    lines.extend(
        [
            "",
            "## Main Table 2: Static Readiness Summary",
            "",
            "| family | readiness | modes | validation rows | max validation delta [%] | max physical delta [%] | largest sampled parameter |",
            "| --- | --- | --- | ---: | ---: | ---: | --- |",
        ]
    )
    for row in static_rows:
        lines.append(
            f"| {row['family']} | {row['readiness_status']} | {row['validated_modes']} | "
            f"{row['validation_rows']} | {row['max_validation_delta_pct']} | "
            f"{row['max_sampled_physical_delta_pct']} | {row['largest_sampled_parameter']} |"
        )

    lines.extend(
        [
            "",
            "## Placement And Captions",
            "",
            "| item | placement | source | caption draft | guardrail |",
            "| --- | --- | --- | --- | --- |",
        ]
    )
    for row in decisions:
        lines.append(
            f"| {row['item']} | {row['placement']} | {row['source_path']} | "
            f"{row['caption']} | {row['referee_guardrail']} |"
        )

    lines.extend(
        [
            "",
            "## Referee-Safe Structural Choice",
            "",
            "- Main text should contain no more than two compact tables and two or three figures.",
            "- The public GW250114 projection remains the observational anchor.",
            "- The static branch is justified as a ringdown-readiness standard and negative-control audit.",
            "- Detailed robustness and sparse-threshold rows should be supplementary to avoid over-selling diagnostics as constraints.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    public_rows = build_public_projection_main_table()
    static_rows = build_static_main_table()
    decisions = build_decisions()

    table1_path = OUT_DIR / "main_table1_public_projection_summary.csv"
    write_csv(
        table1_path,
        public_rows,
        [
            "projection",
            "observable_set",
            "tested_rows",
            "zero_outside_90pct",
            "max_sigma_from_zero",
            "max_sigma_row",
            "tightest_alpha_sigma",
            "note",
        ],
    )

    table2_path = OUT_DIR / "main_table2_static_readiness_summary.csv"
    write_csv(
        table2_path,
        static_rows,
        [
            "family",
            "readiness_status",
            "validated_modes",
            "validation_rows",
            "max_validation_delta_pct",
            "max_sampled_physical_delta_pct",
            "largest_sampled_parameter",
            "interpretation_boundary",
        ],
    )

    decisions_path = OUT_DIR / "table_figure_decisions.csv"
    write_csv(
        decisions_path,
        decisions,
        ["item", "placement", "source_path", "caption", "referee_guardrail"],
    )

    plan_path = OUT_DIR / "table_figure_plan.md"
    write_markdown(plan_path, public_rows, static_rows, decisions)

    print(f"Main table 1: {table1_path}")
    print(f"Main table 2: {table2_path}")
    print(f"Decisions: {decisions_path}")
    print(f"Plan: {plan_path}")
    print(f"Main/public rows: {len(public_rows)}; static rows: {len(static_rows)}; decisions: {len(decisions)}")


if __name__ == "__main__":
    main()
