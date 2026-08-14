from __future__ import annotations

import ast
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
REQUIRED = {
    "README.md",
    "DATA_SOURCES.md",
    "REPRODUCIBILITY.md",
    "LICENSE",
    "CITATION.cff",
    "requirements.txt",
    "config/hairy_gw250114_publication.json",
    "data/ovalle_rotating_hairy_qnm_zhen_li_2022.csv",
    "scripts/python/hairy_continued_fraction.py",
    "scripts/python/build_hairy_qnm_production_grid.py",
    "scripts/python/gw250114_hairy_constraints.py",
}
FORBIDDEN_EXTENSIONS = {".tex", ".bib", ".png", ".jpg", ".jpeg", ".pdf", ".svg"}
SENSITIVE_PATTERNS = {
    "absolute Windows user path": re.compile(r"[A-Za-z]:\\Users\\", re.IGNORECASE),
    "private key": re.compile(r"BEGIN (?:RSA|OPENSSH|EC) PRIVATE KEY", re.IGNORECASE),
    "common secret assignment": re.compile(
        r"(?:password|passwd|api[_-]?key|access[_-]?token|client[_-]?secret)\s*[:=]\s*['\"][^'\"]+",
        re.IGNORECASE,
    ),
}


def main() -> int:
    failures: list[str] = []
    files = [path for path in ROOT.rglob("*") if path.is_file()]
    relative = {path.relative_to(ROOT).as_posix() for path in files}

    for required in sorted(REQUIRED - relative):
        failures.append(f"missing required file: {required}")

    for path in files:
        rel = path.relative_to(ROOT).as_posix()
        if path.suffix.lower() in FORBIDDEN_EXTENSIONS or path.name.lower().startswith("manuscript"):
            failures.append(f"publication artifact present: {rel}")
        if path.suffix.lower() in {".py", ".md", ".json", ".wl", ".cff", ".txt"}:
            text = path.read_text(encoding="utf-8", errors="replace")
            for label, pattern in SENSITIVE_PATTERNS.items():
                if pattern.search(text):
                    failures.append(f"{label}: {rel}")
        if path.suffix.lower() == ".py":
            try:
                ast.parse(path.read_text(encoding="utf-8"), filename=rel)
            except SyntaxError as error:
                failures.append(f"Python syntax error in {rel}: {error}")

    if failures:
        for failure in failures:
            print(f"FAIL: {failure}")
        return 1

    print(f"PASS: {len(files)} files checked; code-only snapshot is structurally ready.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
