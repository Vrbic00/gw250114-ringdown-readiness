"""Assemble a draft manuscript package from current project outputs."""

from __future__ import annotations

import csv
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


def as_float(value: object) -> float:
    try:
        return float(value)
    except Exception:
        return float("nan")


def path_text(path: Path) -> str:
    return str(path.relative_to(ROOT)).replace("\\", "/")


def constraints_summary() -> dict[str, dict[str, str]]:
    rows = read_csv(ROOT / "results" / "gw250114_constraints_comparison" / "projection_constraints_summary.csv")
    return {row["projection"]: row for row in rows}


def max_consistency_difference() -> tuple[float, str, str]:
    rows = read_csv(ROOT / "results" / "gw250114_constraints_comparison" / "projection_consistency_by_operator.csv")
    row = max(rows, key=lambda item: as_float(item["normalized_projection_difference"]))
    return (
        as_float(row["normalized_projection_difference"]),
        row["operator"],
        row["polarization"],
    )


def static_validation_summary() -> tuple[int, str]:
    rows = read_csv(ROOT / "results" / "static_qnm_scorecard" / "static_qnm_validation_summary.csv")
    families = sorted({row["family"] for row in rows})
    return len(rows), ", ".join(families)


def readiness_counts() -> dict[str, int]:
    rows = read_csv(ROOT / "results" / "static_qnm_readiness_audit" / "static_metric_readiness_audit.csv")
    counts: dict[str, int] = {}
    for row in rows:
        counts[row["audit_status"]] = counts.get(row["audit_status"], 0) + 1
    return counts


def largest_physical_deviations() -> list[dict[str, str]]:
    rows = read_csv(ROOT / "results" / "static_qnm_physical_deviations" / "static_qnm_physical_deviations.csv")
    best: dict[str, dict[str, str]] = {}
    for row in rows:
        key = f"{row['family']} n={row['n']}"
        if key not in best or as_float(row["max_abs_physical_delta_pct"]) > as_float(best[key]["max_abs_physical_delta_pct"]):
            best[key] = row
    return [best[key] for key in sorted(best)]


def build_manifest() -> list[dict[str, object]]:
    rows = [
        {
            "item": "Table 1",
            "role": "main projected EFT constraints",
            "source_path": path_text(ROOT / "results" / "gw250114_paper_tables" / "table1_main_projected_constraints.csv"),
            "manuscript_section": "Results: public GW250114 projection",
            "status": "ready",
        },
        {
            "item": "Table 2",
            "role": "RINGDOWN versus pyRing consistency",
            "source_path": path_text(ROOT / "results" / "gw250114_paper_tables" / "table2_pipeline_consistency.csv"),
            "manuscript_section": "Results: consistency and robustness",
            "status": "ready",
        },
        {
            "item": "Table 3",
            "role": "pyRing lower-tail filter robustness",
            "source_path": path_text(ROOT / "results" / "gw250114_paper_tables" / "table3_pyring_filter_robustness.csv"),
            "manuscript_section": "Supplement: robustness",
            "status": "ready",
        },
        {
            "item": "Table 4",
            "role": "empirical linearized posterior check",
            "source_path": path_text(ROOT / "results" / "gw250114_paper_tables" / "table4_linearized_posterior_check.csv"),
            "manuscript_section": "Supplement: non-Gaussian check",
            "status": "ready",
        },
        {
            "item": "Table 5",
            "role": "public ringdown observables",
            "source_path": path_text(ROOT / "results" / "gw250114_paper_tables" / "table5_public_ringdown_observables.csv"),
            "manuscript_section": "Data inputs",
            "status": "ready",
        },
        {
            "item": "Table 6",
            "role": "static QNM validation scorecard",
            "source_path": path_text(ROOT / "results" / "static_qnm_scorecard" / "static_qnm_validation_summary.csv"),
            "manuscript_section": "Static-QNM readiness audit",
            "status": "ready",
        },
        {
            "item": "Table 7",
            "role": "static metric readiness audit",
            "source_path": path_text(ROOT / "results" / "static_qnm_readiness_audit" / "static_metric_readiness_audit.csv"),
            "manuscript_section": "Static-QNM readiness audit",
            "status": "ready",
        },
        {
            "item": "Table 8",
            "role": "static physical QNM deviations and sparse thresholds",
            "source_path": path_text(ROOT / "results" / "static_qnm_physical_deviations" / "static_qnm_physical_thresholds.csv"),
            "manuscript_section": "Static physical-deviation stress tests",
            "status": "ready_diagnostic_not_constraint",
        },
        {
            "item": "Figure 1",
            "role": "projected EFT sigma from zero comparison",
            "source_path": path_text(ROOT / "results" / "gw250114_constraints_comparison" / "projection_sigma_from_zero_comparison.png"),
            "manuscript_section": "Results: public GW250114 projection",
            "status": "ready",
        },
        {
            "item": "Figure 2",
            "role": "static QNM validation scorecard",
            "source_path": path_text(ROOT / "results" / "static_qnm_readiness_audit" / "static_qnm_validation_scorecard.svg"),
            "manuscript_section": "Static-QNM readiness audit",
            "status": "ready",
        },
        {
            "item": "Figure 3",
            "role": "static metric ringdown-readiness ladder",
            "source_path": path_text(ROOT / "results" / "static_qnm_readiness_audit" / "static_metric_readiness_ladder.svg"),
            "manuscript_section": "Methods or limitations",
            "status": "ready",
        },
        {
            "item": "Figure 4",
            "role": "physical static-QNM deviations from baseline",
            "source_path": path_text(ROOT / "results" / "static_qnm_physical_deviations" / "static_qnm_physical_deviations.svg"),
            "manuscript_section": "Static physical-deviation stress tests",
            "status": "ready_diagnostic_not_constraint",
        },
    ]
    return rows


def build_skeleton() -> str:
    projection = constraints_summary()
    consistency_value, consistency_operator, consistency_polarization = max_consistency_difference()
    static_rows, static_families = static_validation_summary()
    counts = readiness_counts()
    physical = largest_physical_deviations()

    ringdown = projection["RINGDOWN"]
    pyring = projection["PYRING_DELTA"]
    count_text = ", ".join(f"{key}: {counts[key]}" for key in sorted(counts))

    physical_lines = []
    for row in physical:
        physical_lines.append(
            f"- {row['family']} `n={row['n']}` reaches max sampled shift "
            f"`{as_float(row['max_abs_physical_delta_pct']):.2f}%` at "
            f"`{row['parameter_name']}={as_float(row['parameter_value']):g}`."
        )

    return "\n".join(
        [
            "# Manuscript Skeleton",
            "",
            "Working title:",
            "",
            "```text",
            "Public-data projected higher-derivative ringdown constraints from GW250114 and a static-QNM readiness audit",
            "```",
            "",
            "## Draft Abstract",
            "",
            "We present a reproducible public-data projection of higher-derivative Kerr-ringdown fingerprints using the public GW250114 spectroscopy products. The RINGDOWN and pyRing projection branches both retain the GR value `alpha=0` inside every one-at-a-time 90 percent interval. The largest nominal displacement from zero is `"
            f"{as_float(ringdown['max_sigma_from_zero']):.3f} sigma` in the RINGDOWN projection and `"
            f"{as_float(pyring['max_sigma_from_zero']):.3f} sigma` in the pyRing deviation projection. The two public-product projections are mutually consistent, with maximum normalized projection difference `"
            f"{consistency_value:.3f} sigma` for `{consistency_operator}/{consistency_polarization}`. We also provide a static spherical QNM readiness audit for metric models often used in QPO and shadow phenomenology, validating `"
            f"{static_rows}` static supplied-potential rows across `{static_families}`. The static branch is not a rotating-remnant constraint, but it demonstrates that validated master-potential spectra can produce percent-to-tens-of-percent physical QNM shifts while metric-only models fail the ringdown-readiness gate.",
            "",
            "## Central Claim",
            "",
            "Using public GW250114 ringdown/spectroscopy posterior products and published higher-derivative Kerr QNM fingerprints, the linearized one-at-a-time projection finds no robust beyond-Kerr deviation. In parallel, static spherical metrics can be audited for ringdown readiness only when they supply gravitational perturbation physics and reproducible QNM spectra.",
            "",
            "## Must Not Claim",
            "",
            "- This is not a full LVK-style strain-level EFT likelihood.",
            "- The RINGDOWN and pyRing products are not independent likelihoods and should not be statistically combined.",
            "- The static spherical branch does not rule out rotating remnant alternatives to GW250114.",
            "- Metric-only, QPO-only, or shadow-only models should not be called ringdown-constrained.",
            "- Sparse static threshold crossings are diagnostic interpolation results, not observational exclusion intervals.",
            "",
            "## Suggested Paper Structure",
            "",
            "1. Introduction: GW250114 as public ringdown benchmark and why reproducible projections matter.",
            "2. Public data and theory inputs: RINGDOWN, pyRing, NRSur7dq4 remnant calibration, and Cano/BeyondKerrQNM fingerprints.",
            "3. Projection method: linearized EFT fingerprints, mass/spin profiling, separate RINGDOWN and pyRing branches.",
            "4. Public-data results: projected constraints, branch consistency, robustness, empirical posterior check.",
            "5. Static-QNM readiness audit: candidate registry, validation scorecard, readiness ladder.",
            "6. Static physical-deviation stress tests: distinguish validation error from physical spectral shifts.",
            "7. Limitations and outlook: strain likelihood, multi-parameter EFT, future events, theory-backed rotating alternatives.",
            "",
            "## Key Numerical Points",
            "",
            f"- RINGDOWN projection rows: `{ringdown['rows']}`, zero outside 90 percent: `{ringdown['zero_outside_90pct_count']}`, max sigma from zero: `{as_float(ringdown['max_sigma_from_zero']):.3f}`.",
            f"- pyRing projection rows: `{pyring['rows']}`, zero outside 90 percent: `{pyring['zero_outside_90pct_count']}`, max sigma from zero: `{as_float(pyring['max_sigma_from_zero']):.3f}`.",
            f"- Largest RINGDOWN/pyRing normalized projection difference: `{consistency_value:.3f} sigma`.",
            f"- Static readiness counts: {count_text}.",
            *physical_lines,
            "",
            "## Table And Figure Plan",
            "",
            "See `table_figure_manifest.csv` in this package.",
            "",
            "## Immediate Writing Tasks",
            "",
            "1. Convert this skeleton into a LaTeX or Markdown manuscript draft.",
            "2. Write the Methods section with exact public-data variable definitions.",
            "3. Add short derivations for the Gaussian projection and nuisance profiling.",
            "4. Polish figure captions so every limitation is explicit.",
            "5. Decide target journal and tune the emphasis: methods/reproducibility versus phenomenological audit.",
            "",
        ]
    )


def build_guardrails() -> str:
    return "\n".join(
        [
            "# Claim Guardrails",
            "",
            "## Safe Main-Text Claims",
            "",
            "- Public GW250114 RINGDOWN and pyRing products give mutually consistent linearized one-at-a-time EFT projections.",
            "- No projected coupling excludes `alpha=0` at the configured 90 percent interval level.",
            "- The static-QNM branch is a readiness and stress-test audit for static spherical metrics, not a direct GW250114 alternative-metric constraint.",
            "- Validated static supplied-potential examples show physical QNM shifts far larger than numerical validation errors.",
            "",
            "## Unsafe Or Overstated Claims",
            "",
            "- Do not say that GW250114 rules out Bardeen, Hayward, tidal-charge, or other static metrics.",
            "- Do not say the public-product projection is equivalent to a full strain-level likelihood.",
            "- Do not combine RINGDOWN and pyRing as independent constraints.",
            "- Do not treat scalar, electromagnetic, shadow, or QPO calculations as gravitational ringdown evidence unless the gravitational perturbation system is supplied and validated.",
            "",
            "## Referee-Resistant Framing",
            "",
            "The paper is strongest as a reproducible benchmark plus a model-readiness filter. Its critical edge is not an overclaim of exclusion, but a clear standard: a metric used for black-hole phenomenology is not ringdown-ready until its perturbation physics and QNM spectrum are explicit and reproducible.",
            "",
        ]
    )


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    manifest = build_manifest()
    manifest_csv = OUT_DIR / "table_figure_manifest.csv"
    write_csv(manifest_csv, manifest, ["item", "role", "source_path", "manuscript_section", "status"])

    skeleton_path = OUT_DIR / "manuscript_skeleton.md"
    skeleton_path.write_text(build_skeleton(), encoding="utf-8")

    guardrails_path = OUT_DIR / "claim_guardrails.md"
    guardrails_path.write_text(build_guardrails(), encoding="utf-8")

    print(f"Manifest: {manifest_csv}")
    print(f"Skeleton: {skeleton_path}")
    print(f"Guardrails: {guardrails_path}")
    print(f"Items: {len(manifest)}")


if __name__ == "__main__":
    main()
