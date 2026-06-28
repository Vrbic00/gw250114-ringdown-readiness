"""Build a paper-facing readiness audit for static spherical QNM candidates."""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = ROOT / "results" / "static_qnm_readiness_audit"


VALIDATION_FAMILIES = {
    "schwarzschild_rw_zerilli": ["Schwarzschild"],
    "braneworld_tidal_charge_rn_like": ["Braneworld_tidal_charge"],
    "bardeen_ned_axial": ["Bardeen_NED"],
    "hayward_axial_gravitational": ["Hayward", "Hayward_overtone"],
}


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


def registry_status(row: dict[str, str], has_validation: bool) -> tuple[str, str]:
    readiness = row["readiness"]
    status = row["master_potential_status"].lower()

    if has_validation:
        return (
            "validated_static_gravitational_qnm",
            "May be used as a static supplied-potential benchmark, with the non-rotating limitation explicit.",
        )

    if readiness.startswith("A"):
        return (
            "ready_for_next_implementation",
            "Good candidate for the next reproducible benchmark before any observational language.",
        )

    if "source" in status or "ambiguity" in status:
        return (
            "source_model_caution",
            "Useful as a conceptual caution: the same metric can imply different perturbation physics.",
        )

    if readiness.startswith("C"):
        return (
            "future_theory_specific_project",
            "Do not fold into the current paper unless the full perturbation setup is reproduced.",
        )

    return (
        "negative_control_metric_only",
        "Do not promote to gravitational ringdown constraints; cite as a readiness failure or outlook item.",
    )


def scorecard_by_family(rows: list[dict[str, str]]) -> dict[str, dict[str, object]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row["family"]].append(row)

    result: dict[str, dict[str, object]] = {}
    for family, fam_rows in grouped.items():
        max_re = max(abs(as_float(row["delta_real_pct"])) for row in fam_rows)
        max_im = max(abs(as_float(row["delta_imag_pct"])) for row in fam_rows)
        result[family] = {
            "rows": len(fam_rows),
            "max_abs_delta_real_pct": max_re,
            "max_abs_delta_imag_pct": max_im,
            "verdicts": "; ".join(sorted({row["verdict"] for row in fam_rows})),
            "modes": "; ".join(sorted({f"l={row['ell']},n={row['n']}" for row in fam_rows})),
        }
    return result


def svg_escape(text: object) -> str:
    return (
        str(text)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def write_scorecard_svg(path: Path, families: dict[str, dict[str, object]]) -> None:
    ordered = sorted(families.items())
    width = 980
    row_h = 54
    top = 70
    label_w = 210
    plot_w = 620
    height = top + row_h * len(ordered) + 80
    scale_max = 1.05

    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<style>text{font-family:Arial,Helvetica,sans-serif;font-size:15px;fill:#1f2933}.small{font-size:12px;fill:#52606d}.title{font-size:22px;font-weight:700}.axis{stroke:#9aa5b1;stroke-width:1}.grid{stroke:#d9e2ec;stroke-width:1}.threshold{stroke:#d64545;stroke-width:2;stroke-dasharray:5 5}</style>',
        '<text x="30" y="34" class="title">Static QNM validation scorecard</text>',
        '<text x="30" y="56" class="small">Maximum absolute deviation from reference tables. Red dashed line: 1 percent validation target.</text>',
    ]

    x0 = label_w + 40
    y0 = top
    for tick in [0.0, 0.25, 0.5, 0.75, 1.0]:
        x = x0 + tick / scale_max * plot_w
        lines.append(f'<line x1="{x:.1f}" y1="{y0 - 10}" x2="{x:.1f}" y2="{height - 56}" class="grid"/>')
        lines.append(f'<text x="{x - 8:.1f}" y="{height - 34}" class="small">{tick:g}</text>')
    threshold_x = x0 + 1.0 / scale_max * plot_w
    lines.append(f'<line x1="{threshold_x:.1f}" y1="{y0 - 18}" x2="{threshold_x:.1f}" y2="{height - 54}" class="threshold"/>')

    for idx, (family, stats) in enumerate(ordered):
        y = y0 + idx * row_h
        max_re = as_float(stats["max_abs_delta_real_pct"])
        max_im = as_float(stats["max_abs_delta_imag_pct"])
        re_w = min(max_re / scale_max * plot_w, plot_w)
        im_w = min(max_im / scale_max * plot_w, plot_w)
        lines.append(f'<text x="30" y="{y + 20}" font-weight="700">{svg_escape(family)}</text>')
        lines.append(f'<text x="30" y="{y + 39}" class="small">{svg_escape(stats["rows"])} rows, {svg_escape(stats["modes"])}</text>')
        lines.append(f'<rect x="{x0:.1f}" y="{y + 5}" width="{re_w:.1f}" height="17" fill="#2f80ed"/>')
        lines.append(f'<rect x="{x0:.1f}" y="{y + 28}" width="{im_w:.1f}" height="17" fill="#f2994a"/>')
        lines.append(f'<text x="{x0 + re_w + 6:.1f}" y="{y + 18}" class="small">Re {max_re:.3f}%</text>')
        lines.append(f'<text x="{x0 + im_w + 6:.1f}" y="{y + 41}" class="small">|Im| {max_im:.3f}%</text>')

    lines.append(f'<rect x="{x0:.1f}" y="{height - 24}" width="16" height="10" fill="#2f80ed"/><text x="{x0 + 22:.1f}" y="{height - 15}" class="small">Re(omega)</text>')
    lines.append(f'<rect x="{x0 + 120:.1f}" y="{height - 24}" width="16" height="10" fill="#f2994a"/><text x="{x0 + 142:.1f}" y="{height - 15}" class="small">Abs(Im omega)</text>')
    lines.append("</svg>")
    path.write_text("\n".join(lines), encoding="utf-8")


def write_ladder_svg(path: Path) -> None:
    steps = [
        ("Metric/geodesic only", "horizon, photon sphere, ISCO, QPO, shadow"),
        ("Test-field QNM", "scalar or electromagnetic fields on fixed background"),
        ("Supplied master potential", "axial/polar gravitational potential and boundary conditions"),
        ("Validated static QNM", "local solver reproduces published spectra below target error"),
        ("Theory-backed rotating QNM", "action, perturbations, Kerr limit, mode fits for event projection"),
    ]
    width = 1100
    height = 250
    box_w = 190
    gap = 22
    x0 = 28
    y0 = 72
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<style>text{font-family:Arial,Helvetica,sans-serif;font-size:14px;fill:#1f2933}.small{font-size:12px;fill:#52606d}.title{font-size:22px;font-weight:700}.box{fill:#f5f7fa;stroke:#bcccdc;stroke-width:1.5}.ready{fill:#e3fcef;stroke:#2f9e44}.warn{fill:#fff8e1;stroke:#f2a900}.future{fill:#e7f0ff;stroke:#2f80ed}</style>',
        '<text x="28" y="34" class="title">Static metric ringdown-readiness ladder</text>',
        '<text x="28" y="56" class="small">A line element alone is not a gravitational ringdown model.</text>',
    ]
    for idx, (title, body) in enumerate(steps):
        x = x0 + idx * (box_w + gap)
        cls = "box"
        if idx == 2:
            cls = "box warn"
        if idx == 3:
            cls = "box ready"
        if idx == 4:
            cls = "box future"
        lines.append(f'<rect x="{x}" y="{y0}" width="{box_w}" height="115" rx="6" class="{cls}"/>')
        lines.append(f'<text x="{x + 12}" y="{y0 + 26}" font-weight="700">{svg_escape(title)}</text>')
        words = body.split()
        line = ""
        yy = y0 + 52
        for word in words:
            candidate = f"{line} {word}".strip()
            if len(candidate) > 24:
                lines.append(f'<text x="{x + 12}" y="{yy}" class="small">{svg_escape(line)}</text>')
                yy += 18
                line = word
            else:
                line = candidate
        if line:
            lines.append(f'<text x="{x + 12}" y="{yy}" class="small">{svg_escape(line)}</text>')
        if idx < len(steps) - 1:
            ax = x + box_w + 5
            ay = y0 + 58
            lines.append(f'<line x1="{ax}" y1="{ay}" x2="{ax + gap - 10}" y2="{ay}" stroke="#52606d" stroke-width="1.6"/>')
            lines.append(f'<polygon points="{ax + gap - 10},{ay - 5} {ax + gap - 10},{ay + 5} {ax + gap - 2},{ay}" fill="#52606d"/>')
    lines.append("</svg>")
    path.write_text("\n".join(lines), encoding="utf-8")


def write_markdown(path: Path, rows: list[dict[str, object]], family_stats: dict[str, dict[str, object]]) -> None:
    counts: dict[str, int] = defaultdict(int)
    for row in rows:
        counts[str(row["audit_status"])] += 1

    lines = [
        "# Static Metric QNM Readiness Audit",
        "",
        "This report converts the candidate registry and validation scorecard into a paper-facing referee filter.",
        "",
        "Core rule: a static line element, geodesic observable, QPO fit, or shadow calculation is not enough for a gravitational ringdown claim. A candidate needs at least a supplied gravitational master potential plus a reproduced QNM spectrum.",
        "",
        "## Status Counts",
        "",
        "| audit status | count |",
        "| --- | ---: |",
    ]
    for status in sorted(counts):
        lines.append(f"| {status} | {counts[status]} |")

    lines.extend(
        [
            "",
            "## Validated Families",
            "",
            "| family | rows | max abs delta Re [%] | max abs delta Abs(Im) [%] | modes |",
            "| --- | ---: | ---: | ---: | --- |",
        ]
    )
    for family in sorted(family_stats):
        stats = family_stats[family]
        lines.append(
            f"| {family} | {stats['rows']} | "
            f"{as_float(stats['max_abs_delta_real_pct']):.3f} | "
            f"{as_float(stats['max_abs_delta_imag_pct']):.3f} | "
            f"{stats['modes']} |"
        )

    lines.extend(
        [
            "",
            "## Candidate Audit",
            "",
            "| candidate | readiness | audit status | validation families | referee action |",
            "| --- | --- | --- | --- | --- |",
        ]
    )
    for row in rows:
        lines.append(
            f"| {row['candidate']} | {row['readiness']} | {row['audit_status']} | "
            f"{row['validation_families'] or 'none'} | {row['referee_action']} |"
        )

    lines.extend(
        [
            "",
            "## Manuscript Use",
            "",
            "- Main text: use the validated families as positive examples of a static supplied-potential readiness audit.",
            "- Appendix or outlook: list ready-but-not-yet-implemented candidates such as RN/gravito-electromagnetic perturbations.",
            "- Negative-control paragraph: metric-only regular black holes, shadow/QPO-only models, and source-ambiguous examples cannot be promoted to gravitational ringdown constraints without perturbation physics.",
            "- GW250114 connection: keep this branch separate from rotating-remnant constraints; its role is community-facing model triage.",
            "",
            "Figures generated with this audit:",
            "",
            "- `static_qnm_validation_scorecard.svg`",
            "- `static_metric_readiness_ladder.svg`",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    registry = read_csv(ROOT / "data" / "static_master_potential_candidate_registry.csv")
    scorecard = read_csv(ROOT / "results" / "static_qnm_scorecard" / "static_qnm_validation_summary.csv")
    family_stats = scorecard_by_family(scorecard)

    audit_rows: list[dict[str, object]] = []
    for row in registry:
        families = VALIDATION_FAMILIES.get(row["candidate"], [])
        matched = [family for family in families if family in family_stats]
        has_validation = bool(matched)
        audit_status, referee_action = registry_status(row, has_validation)
        matched_stats = [family_stats[family] for family in matched]
        max_re = max((as_float(stats["max_abs_delta_real_pct"]) for stats in matched_stats), default="")
        max_im = max((as_float(stats["max_abs_delta_imag_pct"]) for stats in matched_stats), default="")
        validation_rows = sum(int(stats["rows"]) for stats in matched_stats)
        audit_rows.append(
            {
                "candidate": row["candidate"],
                "readiness": row["readiness"],
                "audit_status": audit_status,
                "physics_level": row["physics_level"],
                "validation_families": "; ".join(matched),
                "validation_rows": validation_rows,
                "max_abs_delta_real_pct": max_re,
                "max_abs_delta_imag_pct": max_im,
                "recommended_role": row["recommended_role"],
                "main_risk": row["main_risk"],
                "referee_action": referee_action,
                "primary_sources": row["primary_sources"],
            }
        )

    fields = [
        "candidate",
        "readiness",
        "audit_status",
        "physics_level",
        "validation_families",
        "validation_rows",
        "max_abs_delta_real_pct",
        "max_abs_delta_imag_pct",
        "recommended_role",
        "main_risk",
        "referee_action",
        "primary_sources",
    ]
    csv_path = OUT_DIR / "static_metric_readiness_audit.csv"
    write_csv(csv_path, audit_rows, fields)

    md_path = OUT_DIR / "static_metric_readiness_audit.md"
    write_markdown(md_path, audit_rows, family_stats)

    scorecard_svg = OUT_DIR / "static_qnm_validation_scorecard.svg"
    write_scorecard_svg(scorecard_svg, family_stats)

    ladder_svg = OUT_DIR / "static_metric_readiness_ladder.svg"
    write_ladder_svg(ladder_svg)

    print(f"Audit CSV: {csv_path}")
    print(f"Audit report: {md_path}")
    print(f"Validation figure: {scorecard_svg}")
    print(f"Readiness ladder: {ladder_svg}")
    for row in audit_rows:
        print(f"{row['candidate']}: {row['audit_status']} [{row['validation_families'] or 'no validation'}]")


if __name__ == "__main__":
    main()
