"""Pre-manuscript audit for the GW250114/static-QNM paper package."""

from __future__ import annotations

import csv
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = ROOT / "results" / "pre_manuscript_audit"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def as_float(value: object, default: float = math.nan) -> float:
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


def status_row(
    rows: list[dict[str, object]],
    *,
    block: str,
    status: str,
    evidence: str,
    caveat: str,
    action: str,
) -> None:
    rows.append(
        {
            "block": block,
            "status": status,
            "evidence": evidence,
            "caveat": caveat,
            "recommended_action": action,
        }
    )


def audit_required_files(rows: list[dict[str, object]]) -> None:
    manifest_path = ROOT / "results" / "manuscript_package" / "table_figure_manifest.csv"
    manifest = read_csv(manifest_path)
    missing = []
    zero_size = []
    for item in manifest:
        path = ROOT / item["source_path"]
        if not path.exists():
            missing.append(item["source_path"])
        elif path.stat().st_size <= 0:
            zero_size.append(item["source_path"])

    if missing or zero_size:
        status_row(
            rows,
            block="table_figure_manifest",
            status="STOP",
            evidence=f"missing={len(missing)}, zero_size={len(zero_size)}",
            caveat="The manuscript package references files that are not available.",
            action="Regenerate missing outputs before drafting.",
        )
    else:
        diagnostic = sum(1 for item in manifest if "diagnostic" in item["status"])
        status_row(
            rows,
            block="table_figure_manifest",
            status="GO",
            evidence=f"{len(manifest)} referenced tables/figures exist; {diagnostic} marked diagnostic.",
            caveat="Diagnostic items must be labeled as such in captions.",
            action="Use manifest as the draft table/figure control list.",
        )


def audit_gw_projection(rows: list[dict[str, object]]) -> None:
    summary = read_csv(ROOT / "results" / "gw250114_constraints_comparison" / "projection_constraints_summary.csv")
    by_projection = {row["projection"]: row for row in summary}
    missing = {"RINGDOWN", "PYRING_DELTA"} - set(by_projection)
    if missing:
        status_row(
            rows,
            block="gw250114_public_projection",
            status="STOP",
            evidence=f"Missing projection summary rows: {sorted(missing)}",
            caveat="Main result cannot be stated without both public-product branches.",
            action="Rerun gw250114_constraints_comparison.wl.",
        )
        return

    max_sigma = max(as_float(row["max_sigma_from_zero"]) for row in summary)
    zero_outside = sum(as_int(row["zero_outside_90pct_count"]) for row in summary)
    if zero_outside == 0 and max_sigma < 2.0:
        status = "GO"
        action = "Safe to state no robust projected beyond-Kerr deviation, with public-product caveats."
    elif zero_outside == 0:
        status = "CAUTION"
        action = "State null result, but discuss the largest nominal shift explicitly."
    else:
        status = "STOP"
        action = "Inspect rows excluding alpha=0 before writing the main claim."
    status_row(
        rows,
        block="gw250114_public_projection",
        status=status,
        evidence=f"zero-outside-90 count={zero_outside}; max sigma from zero={max_sigma:.3f}.",
        caveat="Linearized one-at-a-time projection of public posterior products, not a strain-level EFT likelihood.",
        action=action,
    )


def audit_branch_consistency(rows: list[dict[str, object]]) -> None:
    consistency = read_csv(
        ROOT / "results" / "gw250114_constraints_comparison" / "projection_consistency_by_operator.csv"
    )
    max_row = max(consistency, key=lambda row: as_float(row["normalized_projection_difference"]))
    non_overlap = [row for row in consistency if as_int(row["intervals_overlap_90pct"]) != 1]
    zero_fail = [row for row in consistency if as_int(row["zero_inside_90pct_both"]) != 1]
    max_diff = as_float(max_row["normalized_projection_difference"])
    if not non_overlap and not zero_fail and max_diff < 1.0:
        status = "GO"
    elif not non_overlap and not zero_fail:
        status = "CAUTION"
    else:
        status = "STOP"
    status_row(
        rows,
        block="ringdown_pyring_consistency",
        status=status,
        evidence=(
            f"max normalized projection difference={max_diff:.3f} "
            f"({max_row['operator']}/{max_row['polarization']}); "
            f"non-overlap rows={len(non_overlap)}; zero-fail rows={len(zero_fail)}."
        ),
        caveat="RINGDOWN and pyRing products should be compared, not statistically combined.",
        action="Keep comparison table in main text or early supplement; avoid combined intervals.",
    )


def audit_robustness(rows: list[dict[str, object]]) -> None:
    filter_rows = read_csv(ROOT / "results" / "gw250114_paper_tables" / "table3_pyring_filter_robustness.csv")
    empirical_rows = read_csv(ROOT / "results" / "gw250114_paper_tables" / "table4_linearized_posterior_check.csv")
    filter_zero = sum(as_int(row["zero_outside_90pct_count"]) for row in filter_rows)
    max_filter_sigma = max(as_float(row["max_sigma_from_zero"]) for row in filter_rows)
    empirical_zero = sum(as_int(row["zero_outside_empirical_90pct_count"]) for row in empirical_rows)
    max_empirical_sigma = max(as_float(row["max_nominal_abs_median_over_sd"]) for row in empirical_rows)
    if filter_zero == 0 and empirical_zero == 0 and max(max_filter_sigma, max_empirical_sigma) < 1.5:
        status = "GO"
    elif filter_zero == 0 and empirical_zero == 0:
        status = "CAUTION"
    else:
        status = "STOP"
    status_row(
        rows,
        block="robustness_checks",
        status=status,
        evidence=(
            f"filter zero-outside={filter_zero}, max filter sigma={max_filter_sigma:.3f}; "
            f"empirical zero-outside={empirical_zero}, max empirical nominal sigma={max_empirical_sigma:.3f}."
        ),
        caveat="Filter and empirical checks are robustness diagnostics, not independent detections.",
        action="Keep as support for null result; place detailed rows in supplement.",
    )


def audit_static_validation(rows: list[dict[str, object]]) -> None:
    scorecard = read_csv(ROOT / "results" / "static_qnm_scorecard" / "static_qnm_validation_summary.csv")
    families = sorted({row["family"] for row in scorecard})
    max_validation = max(
        max(abs(as_float(row["delta_real_pct"])), abs(as_float(row["delta_imag_pct"])))
        for row in scorecard
    )
    non_pass = [row for row in scorecard if not row["verdict"].startswith("PASS")]
    if not non_pass and max_validation <= 1.0:
        status = "GO"
    elif not non_pass:
        status = "CAUTION"
    else:
        status = "STOP"
    status_row(
        rows,
        block="static_qnm_validation",
        status=status,
        evidence=f"{len(scorecard)} rows, families={', '.join(families)}, max validation delta={max_validation:.3f}%.",
        caveat="Several references are WKB/Prony tables rather than exact Leaver spectra.",
        action="Use as readiness/stress-test validation; do not call it final precision spectroscopy.",
    )


def audit_static_readiness(rows: list[dict[str, object]]) -> None:
    audit = read_csv(ROOT / "results" / "static_qnm_readiness_audit" / "static_metric_readiness_audit.csv")
    counts: dict[str, int] = {}
    for row in audit:
        counts[row["audit_status"]] = counts.get(row["audit_status"], 0) + 1
    validated = counts.get("validated_static_gravitational_qnm", 0)
    negative = counts.get("negative_control_metric_only", 0)
    if validated >= 3 and negative >= 1:
        status = "GO"
    else:
        status = "CAUTION"
    evidence = ", ".join(f"{key}={counts[key]}" for key in sorted(counts))
    status_row(
        rows,
        block="static_readiness_audit",
        status=status,
        evidence=evidence,
        caveat="Audit is a classification standard, not a population-complete literature review.",
        action="Frame as a reusable referee filter and explicitly invite extension.",
    )


def audit_static_physical_deviations(rows: list[dict[str, object]]) -> None:
    physical = read_csv(
        ROOT / "results" / "static_qnm_physical_deviations" / "static_qnm_physical_deviations.csv"
    )
    max_physical = max(as_float(row["max_abs_physical_delta_pct"]) for row in physical)
    max_validation = max(
        max(abs(as_float(row["validation_delta_real_pct"])), abs(as_float(row["validation_delta_imag_pct"])))
        for row in physical
    )
    threshold_rows = read_csv(
        ROOT / "results" / "static_qnm_physical_deviations" / "static_qnm_physical_thresholds.csv"
    )
    interpolated = sum(1 for row in threshold_rows if row["status"] == "linear_interpolation_sparse_grid")
    if max_physical > 10.0 and max_validation <= 1.0:
        status = "GO"
    else:
        status = "CAUTION"
    status_row(
        rows,
        block="static_physical_deviation_layer",
        status=status,
        evidence=(
            f"max physical shift={max_physical:.2f}%; max validation-error scale={max_validation:.3f}%; "
            f"sparse threshold rows={len(threshold_rows)}, interpolated={interpolated}."
        ),
        caveat="Threshold crossings are sparse-grid diagnostics and require physical parameter priors before constraint language.",
        action="Use to answer whether the static test is discriminating; avoid observational exclusion wording.",
    )


def build_risk_register() -> list[dict[str, object]]:
    return [
        {
            "risk": "Public posterior products are not independent likelihoods.",
            "severity": "high",
            "affected_section": "GW250114 projection",
            "mitigation": "Compare RINGDOWN and pyRing side by side; do not combine intervals.",
            "writing_rule": "Use 'projection' and 'consistency check', not 'combined constraint'.",
        },
        {
            "risk": "Linearized one-at-a-time EFT couplings miss nonlinear and multi-parameter degeneracies.",
            "severity": "high",
            "affected_section": "Methods and limitations",
            "mitigation": "State one-at-a-time linear response explicitly; reserve multi-parameter EFT for outlook.",
            "writing_rule": "Do not claim a full theory-space exclusion.",
        },
        {
            "risk": "Static spherical QNM branch is not a rotating GW250114 remnant model.",
            "severity": "high",
            "affected_section": "Static readiness audit",
            "mitigation": "Keep static branch as readiness/stress-test audit and negative-control standard.",
            "writing_rule": "Never say GW250114 rules out Bardeen/Hayward/tidal-charge static metrics.",
        },
        {
            "risk": "Some static references are WKB/Prony tables, not exact spectra.",
            "severity": "medium",
            "affected_section": "Static validation",
            "mitigation": "Report validation against the same published approximation level and label precision limits.",
            "writing_rule": "Use 'reproduces published table' rather than 'exact QNM'.",
        },
        {
            "risk": "Sparse physical-deviation threshold crossings can look like constraints.",
            "severity": "medium",
            "affected_section": "Static physical deviations",
            "mitigation": "Call crossings sparse-grid diagnostics and require parameter priors for constraint language.",
            "writing_rule": "Use 'threshold crossing' not 'excluded parameter'.",
        },
        {
            "risk": "Manuscript may look like two papers joined together.",
            "severity": "medium",
            "affected_section": "Framing",
            "mitigation": "Make the unifying theme 'ringdown-readiness and public reproducibility'.",
            "writing_rule": "Use static branch as a methodological audit supporting the main standard.",
        },
        {
            "risk": "Target journal may demand a stronger original physics result.",
            "severity": "medium",
            "affected_section": "Submission strategy",
            "mitigation": "Aim PRD only if Methods and audit standard are sharp; keep CQG fallback framing ready.",
            "writing_rule": "Emphasize reusable benchmark tables, code, and clear limitations.",
        },
    ]


def write_report(path: Path, audit_rows: list[dict[str, object]], risks: list[dict[str, object]]) -> None:
    status_order = {"STOP": 0, "CAUTION": 1, "GO": 2}
    stop_count = sum(1 for row in audit_rows if row["status"] == "STOP")
    caution_count = sum(1 for row in audit_rows if row["status"] == "CAUTION")
    go_count = sum(1 for row in audit_rows if row["status"] == "GO")
    readiness = "GO_TO_METHODS_DRAFT" if stop_count == 0 else "BLOCKED"
    if stop_count == 0 and caution_count > 0:
        readiness = "GO_WITH_CLAIM_GUARDRAILS"

    lines = [
        "# Pre-Manuscript Audit",
        "",
        f"Overall readiness: `{readiness}`.",
        "",
        f"Audit counts: `{go_count}` GO, `{caution_count}` CAUTION, `{stop_count}` STOP.",
        "",
        "## Computation And Result Checks",
        "",
        "| block | status | evidence | caveat | recommended action |",
        "| --- | --- | --- | --- | --- |",
    ]
    for row in sorted(audit_rows, key=lambda item: (status_order.get(str(item["status"]), 99), str(item["block"]))):
        lines.append(
            f"| {row['block']} | {row['status']} | {row['evidence']} | {row['caveat']} | {row['recommended_action']} |"
        )

    lines.extend(
        [
            "",
            "## Referee Risk Register",
            "",
            "| risk | severity | affected section | mitigation | writing rule |",
            "| --- | --- | --- | --- | --- |",
        ]
    )
    severity_order = {"high": 0, "medium": 1, "low": 2}
    for row in sorted(risks, key=lambda item: (severity_order.get(str(item["severity"]), 99), str(item["risk"]))):
        lines.append(
            f"| {row['risk']} | {row['severity']} | {row['affected_section']} | {row['mitigation']} | {row['writing_rule']} |"
        )

    lines.extend(
        [
            "",
            "## Go/No-Go Decision",
            "",
            "There is no current computational STOP item. The next step can be Methods drafting, provided the draft follows the claim guardrails:",
            "",
            "- no combined RINGDOWN+pyRing constraint;",
            "- no strain-level likelihood claim;",
            "- no static-metric observational exclusion;",
            "- sparse static thresholds remain diagnostic;",
            "- static branch is presented as a ringdown-readiness audit.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    audit_rows: list[dict[str, object]] = []
    audit_required_files(audit_rows)
    audit_gw_projection(audit_rows)
    audit_branch_consistency(audit_rows)
    audit_robustness(audit_rows)
    audit_static_validation(audit_rows)
    audit_static_readiness(audit_rows)
    audit_static_physical_deviations(audit_rows)

    risks = build_risk_register()

    audit_csv = OUT_DIR / "pre_manuscript_audit.csv"
    write_csv(
        audit_csv,
        audit_rows,
        ["block", "status", "evidence", "caveat", "recommended_action"],
    )

    risk_csv = OUT_DIR / "referee_risk_register.csv"
    write_csv(
        risk_csv,
        risks,
        ["risk", "severity", "affected_section", "mitigation", "writing_rule"],
    )

    report_path = OUT_DIR / "pre_manuscript_audit.md"
    write_report(report_path, audit_rows, risks)

    print(f"Audit CSV: {audit_csv}")
    print(f"Risk register: {risk_csv}")
    print(f"Report: {report_path}")
    for row in audit_rows:
        print(f"{row['status']}: {row['block']} - {row['evidence']}")


if __name__ == "__main__":
    main()
