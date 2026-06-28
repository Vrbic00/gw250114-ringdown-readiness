"""Scan manuscript wording for known referee-risk phrases."""

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = ROOT / "results" / "manuscript_package"
DEFAULT_MANUSCRIPT = OUT_DIR / "manuscript_v0.md"


RISK_RULES = [
    {
        "rule": "direct_static_exclusion",
        "pattern": r"rule[s]? out .*?(Bardeen|Hayward|tidal|static)",
        "severity": "check_context",
        "safe_context": "Allowed only when explicitly negated, e.g. 'do not rule out'.",
    },
    {
        "rule": "combined_constraint_language",
        "pattern": r"combined .*?(constraint|bound|likelihood)",
        "severity": "check_context",
        "safe_context": "Allowed only when saying no combined likelihood/constraint is constructed.",
    },
    {
        "rule": "strain_likelihood_language",
        "pattern": r"full strain-level|strain-level beyond-Kerr likelihood|strain-level EFT",
        "severity": "check_context",
        "safe_context": "Allowed only in limitation/future-work language.",
    },
    {
        "rule": "constraint_without_projection",
        "pattern": r"(?<!projected )(?<!diagnostic )(?<!observational )constraint[s]?",
        "severity": "medium",
        "safe_context": "Prefer 'projected constraint', 'diagnostic threshold', or explicit limitation.",
    },
    {
        "rule": "excluded_language",
        "pattern": r"\b(excluded|excludes|exclusion)\b",
        "severity": "high",
        "safe_context": "Avoid unless explicitly saying no exclusion is claimed.",
    },
    {
        "rule": "independent_likelihood_language",
        "pattern": r"independent likelihood[s]?",
        "severity": "check_context",
        "safe_context": "Allowed only when saying RINGDOWN and pyRing are not independent likelihoods.",
    },
]


def context_for(lines: list[str], line_no: int, radius: int = 1) -> str:
    start = max(0, line_no - 1 - radius)
    stop = min(len(lines), line_no + radius)
    return " ".join(line.strip() for line in lines[start:stop])


def is_obviously_negated(context: str) -> bool:
    lowered = context.lower()
    negators = [
        "not ",
        "no ",
        "do not",
        "does not",
        "should not",
        "rather than",
        "without",
    ]
    return any(token in lowered for token in negators)


def main() -> None:
    manuscript = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_MANUSCRIPT
    if not manuscript.is_absolute():
        manuscript = ROOT / manuscript
    text = manuscript.read_text(encoding="utf-8")
    lines = text.splitlines()
    rows: list[dict[str, object]] = []

    for rule in RISK_RULES:
        pattern = re.compile(str(rule["pattern"]), re.IGNORECASE)
        for idx, line in enumerate(lines, start=1):
            for match in pattern.finditer(line):
                context = context_for(lines, idx)
                negated = is_obviously_negated(context)
                status = "OK_NEGATED_OR_LIMITED" if negated else "REVIEW"
                if rule["severity"] == "medium" and "projected" in context.lower():
                    status = "OK_PROJECTED_CONTEXT"
                rows.append(
                    {
                        "rule": rule["rule"],
                        "severity": rule["severity"],
                        "line": idx,
                        "match": match.group(0),
                        "status": status,
                        "safe_context": rule["safe_context"],
                        "context": context,
                    }
                )

    stem = manuscript.stem
    csv_path = OUT_DIR / f"{stem}_text_audit.csv"
    fields = ["rule", "severity", "line", "match", "status", "safe_context", "context"]
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    review_rows = [row for row in rows if row["status"] == "REVIEW"]
    report_path = OUT_DIR / f"{stem}_text_audit.md"
    report_lines = [
        "# Manuscript Text Audit",
        "",
        f"Scanned file: `{manuscript.relative_to(ROOT)}`",
        "",
        f"Risk hits: `{len(rows)}` total, `{len(review_rows)}` requiring review.",
        "",
        "| status | rule | line | match | context |",
        "| --- | --- | ---: | --- | --- |",
    ]
    for row in rows:
        report_lines.append(
            f"| {row['status']} | {row['rule']} | {row['line']} | {row['match']} | {row['context']} |"
        )
    report_lines.extend(
        [
            "",
            "Interpretation:",
            "",
            "- `OK_NEGATED_OR_LIMITED` means the risky phrase appears in explicit limitation language.",
            "- `OK_PROJECTED_CONTEXT` means constraint language is tied to projected constraints.",
            "- `REVIEW` rows should be edited before journal submission.",
            "",
        ]
    )
    report_path.write_text("\n".join(report_lines), encoding="utf-8")

    print(f"Text audit CSV: {csv_path}")
    print(f"Text audit report: {report_path}")
    print(f"Total hits: {len(rows)}")
    print(f"Review hits: {len(review_rows)}")
    for row in review_rows:
        print(f"REVIEW line {row['line']}: {row['rule']} -> {row['match']}")


if __name__ == "__main__":
    main()
