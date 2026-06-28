"""Consolidate selected BeyondKerrQNM polynomial fits into one CSV.

The source directory is the public GPL-3.0 repository associated with
arXiv:2409.04517:
https://github.com/pacmn91/BeyondKerrQNM

Only modes used by the current GW250114 project are selected.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


THEORIES = {
    "cubic_even": ("lambda_ev", 6, "preserving"),
    "cubic_odd": ("lambda_odd", 6, "breaking"),
    "quartic_1": ("epsilon1", 8, "preserving"),
    "quartic_2": ("epsilon2", 8, "preserving"),
    "quartic_3": ("epsilon3", 8, "breaking"),
}

MODES = {
    "220": (2, 2, 0),
    "221": (2, 2, 1),
    "222": (2, 2, 2),
    "330": (3, 3, 0),
    "440": (4, 4, 0),
}

REPOSITORY_URL = "https://github.com/pacmn91/BeyondKerrQNM"
PAPER_URL = "https://arxiv.org/abs/2409.04517"
LICENSE = "GPL-3.0"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("fits_dir", type=Path)
    parser.add_argument("output_csv", type=Path)
    parser.add_argument("--commit", required=True)
    return parser.parse_args()


def read_fit(path: Path) -> list[tuple[int, float, float]]:
    coefficients: list[tuple[int, float, float]] = []
    with path.open(encoding="utf-8-sig") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            order, real, imaginary = stripped.split()
            coefficients.append((int(order), float(real), float(imaginary)))
    return coefficients


def main() -> None:
    args = parse_args()
    rows: list[dict[str, object]] = []

    for mode, (ell, emm, overtone) in MODES.items():
        mode_token = f"l{ell}_m{emm}_n{overtone}"
        for theory, (operator, derivative_order, parity) in THEORIES.items():
            for branch in ("plus", "minus"):
                filename = f"{theory}_{branch}_{mode_token}.txt"
                path = args.fits_dir / filename
                if not path.exists():
                    raise FileNotFoundError(path)
                coefficients = read_fit(path)
                maximum_order = max(order for order, _, _ in coefficients)
                for order, real, imaginary in coefficients:
                    rows.append(
                        {
                            "mode": mode,
                            "l": ell,
                            "m": emm,
                            "n": overtone,
                            "theory": theory,
                            "operator": operator,
                            "derivative_order": derivative_order,
                            "parity": parity,
                            "branch": branch,
                            "k": order,
                            "coefficient_re": real,
                            "coefficient_im": imaginary,
                            "maximum_fit_order": maximum_order,
                            "source_file": filename,
                            "source_commit": args.commit,
                            "source_repository": REPOSITORY_URL,
                            "source_paper": PAPER_URL,
                            "license": LICENSE,
                        }
                    )

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0])
    with args.output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} coefficient rows to {args.output_csv}")
    print(f"Modes: {', '.join(MODES)}")
    print(f"Source commit: {args.commit}")


if __name__ == "__main__":
    main()
