"""Build the production QNM atlas for the GW250114 hairy-BH analysis.

The stored grid consists entirely of direct continued-fraction roots.  Dense
posterior surfaces are evaluated later by interpolation, but the publication
atlas itself is not an interpolated visualization product.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
LOCAL_DEPS = ROOT / ".deps" / "python"
if LOCAL_DEPS.exists():
    sys.path.insert(0, str(LOCAL_DEPS))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np

from build_hairy_qnm_grid import solve_grid, validate_interpolator


DEFAULT_CONFIG = ROOT / "config" / "hairy_gw250114_publication.json"
DEFAULT_OUTPUT = ROOT / "results" / "hairy_qnm_production_grid"


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    config = json.loads(args.config.read_text(encoding="utf-8"))["qnm_grid"]
    args.output.mkdir(parents=True, exist_ok=True)

    spins = np.round(
        np.arange(config["spin_min"], config["spin_max"] + 0.5 * config["spin_step"], config["spin_step"]),
        10,
    )
    alphas = np.round(
        np.arange(config["alpha_min"], config["alpha_max"] + 0.5 * config["alpha_step"], config["alpha_step"]),
        10,
    )
    h0s = np.asarray(config["h0_values"], dtype=float)
    overtones = np.asarray(config["overtones"], dtype=int)
    low_depth = int(config["low_depth"])
    high_depth = int(config["high_depth"])

    started = time.perf_counter()
    omega, rows = solve_grid(spins, alphas, h0s, overtones, low_depth, high_depth)
    elapsed_grid = time.perf_counter() - started
    csv_path = args.output / "hairy_qnm_production_grid.csv"
    npz_path = args.output / "hairy_qnm_production_grid.npz"
    write_csv(csv_path, rows)
    np.savez_compressed(
        npz_path,
        spins=spins,
        alphas=alphas,
        h0s=h0s,
        overtones=overtones,
        omega=omega,
    )

    validation_rows = validate_interpolator(
        spins,
        alphas,
        h0s,
        overtones,
        omega,
        high_depth,
        int(config["interpolation_validation_samples"]),
    )
    validation_path = args.output / "production_interpolation_validation.csv"
    write_csv(validation_path, validation_rows)

    deformed = [row for row in rows if float(row["alpha_hair"]) > 0.0]
    summary = {
        "model": "Zhen_Li_rotating_hairy_black_hole",
        "grid_kind": "direct_continued_fraction_roots",
        "shape": list(omega.shape),
        "modes": int(omega.size),
        "axis_counts": {
            "h0": len(h0s),
            "alpha_hair": len(alphas),
            "spin": len(spins),
            "overtone": len(overtones),
        },
        "domain": {
            "h0": [float(h0s[0]), float(h0s[-1])],
            "alpha_hair": [float(alphas[0]), float(alphas[-1])],
            "spin": [float(spins[0]), float(spins[-1])],
            "overtones": overtones.tolist(),
        },
        "all_roots_successful": all(bool(row["root_success"]) for row in rows),
        "maximum_residual_norm": max(float(row["residual_norm"]) for row in rows),
        "maximum_depth_abs_delta": max(float(row["depth_abs_delta"]) for row in rows),
        "median_depth_abs_delta": float(np.median([float(row["depth_abs_delta"]) for row in rows])),
        "maximum_outer_horizon_relative_error": max(float(row["relative_r_plus_error"]) for row in deformed),
        "maximum_inner_horizon_relative_error": max(float(row["relative_r_minus_error"]) for row in deformed),
        "interpolation_validation_cases": len(validation_rows),
        "maximum_interpolation_relative_error": max(float(row["relative_error"]) for row in validation_rows),
        "p95_interpolation_relative_error": float(
            np.quantile([float(row["relative_error"]) for row in validation_rows], 0.95)
        ),
        "median_interpolation_relative_error": float(
            np.median([float(row["relative_error"]) for row in validation_rows])
        ),
        "continued_fraction_depths": [low_depth, high_depth],
        "elapsed_grid_seconds": elapsed_grid,
        "files": {
            csv_path.name: sha256(csv_path),
            npz_path.name: sha256(npz_path),
            validation_path.name: sha256(validation_path),
        },
    }
    (args.output / "production_grid_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
