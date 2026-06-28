from __future__ import annotations

import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

REQUIRED_PATHS = [
    "README.md",
    "DATA_SOURCES.md",
    "REPRODUCIBILITY.md",
    "CITATION.cff",
    "LICENSE",
    "requirements.txt",
    "paper/main.tex",
    "paper/references.bib",
    "results/gw250114_constraints_comparison/projection_constraints_long.csv",
    "results/static_qnm_scorecard/static_qnm_validation_summary.csv",
    "results/static_qnm_readiness_audit/static_metric_readiness_audit.csv",
]

FORBIDDEN_PREFIXES = [
    "data/raw/",
    ".deps/",
    ".deps_local/",
    ".venv/",
    "venv/",
    "__pycache__/",
]

FORBIDDEN_SUFFIXES = [
    ".h5",
    ".hdf5",
    ".tar.gz",
    ".zip",
    ".aux",
    ".bbl",
    ".blg",
    ".log",
    ".tmp",
]

MAX_TRACKED_FILE_MB = 10.0


def tracked_files() -> list[Path]:
    try:
        output = subprocess.check_output(
            ["git", "ls-files"],
            cwd=ROOT,
            text=True,
            stderr=subprocess.DEVNULL,
        )
    except Exception:
        return [p.relative_to(ROOT) for p in ROOT.rglob("*") if p.is_file() and ".git" not in p.parts]

    return [Path(line.strip()) for line in output.splitlines() if line.strip()]


def fail(message: str) -> None:
    print(f"ERROR: {message}", file=sys.stderr)
    raise SystemExit(1)


def main() -> None:
    files = tracked_files()
    file_strings = [path.as_posix() for path in files]

    missing = [path for path in REQUIRED_PATHS if not (ROOT / path).exists()]
    if missing:
        fail("Missing required release files: " + ", ".join(missing))

    forbidden = []
    oversized = []
    for rel, rel_string in zip(files, file_strings):
        if any(rel_string.startswith(prefix) for prefix in FORBIDDEN_PREFIXES):
            forbidden.append(rel_string)
        if any(rel_string.endswith(suffix) for suffix in FORBIDDEN_SUFFIXES):
            forbidden.append(rel_string)

        size_mb = (ROOT / rel).stat().st_size / (1024 * 1024)
        if size_mb > MAX_TRACKED_FILE_MB:
            oversized.append(f"{rel_string} ({size_mb:.2f} MB)")

    if forbidden:
        fail("Forbidden tracked files found: " + ", ".join(sorted(set(forbidden))))
    if oversized:
        fail("Tracked files exceed size limit: " + ", ".join(oversized))

    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    required_phrases = [
        "observational constraints",
        "not committed",
        "DATA_SOURCES.md",
    ]
    missing_phrases = [phrase for phrase in required_phrases if phrase not in readme]
    if missing_phrases:
        fail("README is missing required guardrail phrases: " + ", ".join(missing_phrases))

    total_mb = sum((ROOT / rel).stat().st_size for rel in files) / (1024 * 1024)
    print(f"Release package check passed: {len(files)} tracked files, {total_mb:.2f} MB.")


if __name__ == "__main__":
    main()
