"""Construct a ringdown-data-independent remnant prior for GW250114.

Only the public pre_-40M posterior is used.  Component parameters are mapped to
the remnant with two published fit families:

* BMR2012 final mass + HBR2016 generic-precessing final spin (precession),
* NRSur3dq8 aligned remnant fit augmented by the standard leading in-plane
  spin projection used by PESummary.

An equal-weight mixture of the paired predictions propagates both progenitor
posterior uncertainty and the small between-fit spread.  The result is
independent of the analysed ringdown samples, but remains conditional on GR
remnant fitting formulae; the report states this explicitly.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import types
import warnings
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
LOCAL_DEPS = ROOT / ".deps" / "python"
if LOCAL_DEPS.exists():
    sys.path.insert(0, str(LOCAL_DEPS))

import numpy as np
import precession


def load_surfinbh():
    """Load the aligned fit without requiring unused LAL spin-evolution code."""
    if "lal" not in sys.modules:
        lal = types.ModuleType("lal")
        lal.MSUN_SI = 1.9884099021470416e30
        lal.MTSUN_SI = 4.9254909476412675e-6
        lal.PC_SI = 3.085677581491367e16
        lal.C_SI = 299792458.0
        sys.modules["lal"] = lal
    if "lalsimulation" not in sys.modules:
        sys.modules["lalsimulation"] = types.ModuleType("lalsimulation")
    import surfinBH

    return surfinBH.LoadFits("NRSur3dq8Remnant")


def write_rows(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def covariance_rows(covariance: np.ndarray) -> list[dict]:
    names = ["delta_lnM", "delta_chi"]
    sigma = np.sqrt(np.diag(covariance))
    rows = []
    for i, row_name in enumerate(names):
        for j, column_name in enumerate(names):
            rows.append(
                {
                    "row": row_name,
                    "column": column_name,
                    "covariance": float(covariance[i, j]),
                    "correlation": float(covariance[i, j] / (sigma[i] * sigma[j])),
                }
            )
    return rows


def summary_row(name: str, values: np.ndarray) -> dict:
    q05, q16, q50, q84, q95 = np.quantile(values, [0.05, 0.16, 0.5, 0.84, 0.95])
    return {
        "parameter": name,
        "mean": float(np.mean(values)),
        "standard_deviation": float(np.std(values, ddof=1)),
        "q05": float(q05),
        "q16": float(q16),
        "median": float(q50),
        "q84": float(q84),
        "q95": float(q95),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config", type=Path, default=ROOT / "config" / "article1_nr_validation.json"
    )
    parser.add_argument(
        "--output", type=Path, default=ROOT / "results" / "article1_independent_prior"
    )
    args = parser.parse_args()
    config = json.loads(args.config.read_text(encoding="utf-8"))
    event = config["event"]
    settings = config["independent_prior"]
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)

    samples = np.genfromtxt(ROOT / settings["samples"], names=True)
    mass_total = np.asarray(samples["total_mass"], dtype=float)
    mass_ratio = np.asarray(samples["mass_ratio"], dtype=float)  # m2/m1 <= 1
    spin1_mag = np.asarray(samples["spin1_magnitude"], dtype=float)
    spin2_mag = np.asarray(samples["spin2_magnitude"], dtype=float)
    tilt1 = np.asarray(samples["tilt1"], dtype=float)
    tilt2 = np.asarray(samples["tilt2"], dtype=float)
    phi12 = np.asarray(samples["phi12"], dtype=float)

    # Generic-spin fit family, vectorized over all posterior samples.
    bmr_mass_fraction = np.asarray(
        precession.remnantmass(tilt1, tilt2, mass_ratio, spin1_mag, spin2_mag),
        dtype=float,
    )
    hbr_spin = np.asarray(
        precession.remnantspin(
            tilt1, tilt2, phi12, mass_ratio, spin1_mag, spin2_mag,
            which="HBR16_34corr",
        ),
        dtype=float,
    )
    bmr_mass = mass_total * bmr_mass_fraction

    # Independent implementation/family. NRSur3dq8 itself is aligned-spin;
    # the leading in-plane angular-momentum contribution follows the public
    # PESummary projected-precessing prescription.
    fit = load_surfinbh()
    surfin_mass = np.empty(len(samples), dtype=float)
    surfin_spin_aligned = np.empty(len(samples), dtype=float)
    surfin_spin_projected = np.empty(len(samples), dtype=float)
    surfin_mass_error = np.empty(len(samples), dtype=float)
    surfin_spin_error = np.empty(len(samples), dtype=float)
    warnings.filterwarnings("ignore", module="surfinBH")
    for index, sample in enumerate(samples):
        q_surfin = 1.0 / float(sample["mass_ratio"])
        result = fit.all(
            q_surfin,
            [0.0, 0.0, float(sample["spin1_z"])],
            [0.0, 0.0, float(sample["spin2_z"])],
            allow_extrap=True,
        )
        mass_fraction = float(result[0])
        aligned_spin = abs(float(result[1][2]))
        w1 = (q_surfin / (1.0 + q_surfin)) ** 2
        w2 = (1.0 / (1.0 + q_surfin)) ** 2
        spin_perp_x = w1 * float(sample["spin1_x"]) + w2 * float(sample["spin2_x"])
        spin_perp_y = w1 * float(sample["spin1_y"]) + w2 * float(sample["spin2_y"])
        spin_perp = math.hypot(spin_perp_x, spin_perp_y)
        surfin_mass[index] = float(sample["total_mass"]) * mass_fraction
        surfin_spin_aligned[index] = aligned_spin
        surfin_spin_projected[index] = math.sqrt(aligned_spin**2 + spin_perp**2)
        surfin_mass_error[index] = float(sample["total_mass"]) * float(result[3])
        surfin_spin_error[index] = abs(float(result[4][2]))

    finite = (
        np.isfinite(bmr_mass)
        & np.isfinite(hbr_spin)
        & np.isfinite(surfin_mass)
        & np.isfinite(surfin_spin_projected)
        & (hbr_spin < 1.0)
        & (surfin_spin_projected < 1.0)
    )
    if not np.all(finite):
        warnings.warn(f"Discarding {np.count_nonzero(~finite)} nonphysical remnant samples")
    bmr_mass = bmr_mass[finite]
    hbr_spin = hbr_spin[finite]
    surfin_mass = surfin_mass[finite]
    surfin_spin_aligned = surfin_spin_aligned[finite]
    surfin_spin_projected = surfin_spin_projected[finite]
    surfin_mass_error = surfin_mass_error[finite]
    surfin_spin_error = surfin_spin_error[finite]

    mass0 = float(event["mass_detector_msun"])
    spin0 = float(event["spin"])
    # Equal model weights make the finite between-fit spread part of the prior.
    mixture_mass = np.concatenate((bmr_mass, surfin_mass))
    mixture_spin = np.concatenate((hbr_spin, surfin_spin_projected))
    delta_ln_mass = np.log(mixture_mass / mass0)
    delta_spin = mixture_spin - spin0
    nuisance = np.column_stack((delta_ln_mass, delta_spin))
    covariance = np.cov(nuisance, rowvar=False, ddof=1)
    mean_offset = np.mean(nuisance, axis=0)

    write_rows(output / "pre_minus40M_remnant_prior_covariance.csv", covariance_rows(covariance))
    write_rows(
        output / "pre_minus40M_remnant_prior_mean.csv",
        [
            {"parameter": "delta_lnM", "mean_offset": float(mean_offset[0])},
            {"parameter": "delta_chi", "mean_offset": float(mean_offset[1])},
        ],
    )
    summary = [
        summary_row("final_mass_detector_msun_model_mixture", mixture_mass),
        summary_row("final_spin_model_mixture", mixture_spin),
        summary_row("delta_lnM", delta_ln_mass),
        summary_row("delta_chi", delta_spin),
        summary_row("final_mass_detector_msun_BMR2012", bmr_mass),
        summary_row("final_spin_HBR2016", hbr_spin),
        summary_row("final_mass_detector_msun_NRSur3dq8", surfin_mass),
        summary_row("final_spin_NRSur3dq8_projected", surfin_spin_projected),
        summary_row("final_spin_NRSur3dq8_aligned_only", surfin_spin_aligned),
    ]
    write_rows(output / "pre_minus40M_remnant_prior_summary.csv", summary)

    sample_rows = []
    for index in range(len(bmr_mass)):
        sample_rows.append(
            {
                "sample_index": index,
                "final_mass_BMR2012": float(bmr_mass[index]),
                "final_spin_HBR2016": float(hbr_spin[index]),
                "final_mass_NRSur3dq8": float(surfin_mass[index]),
                "final_spin_NRSur3dq8_projected": float(surfin_spin_projected[index]),
                "final_spin_NRSur3dq8_aligned": float(surfin_spin_aligned[index]),
                "NRSur3dq8_mass_fit_error_msun": float(surfin_mass_error[index]),
                "NRSur3dq8_spin_fit_error": float(surfin_spin_error[index]),
            }
        )
    write_rows(output / "pre_minus40M_remnant_model_crosscheck_samples.csv", sample_rows)

    sigma = np.sqrt(np.diag(covariance))
    corr = covariance[0, 1] / (sigma[0] * sigma[1])
    mass_model_median_difference = float(np.median(surfin_mass - bmr_mass))
    spin_model_median_difference = float(np.median(surfin_spin_projected - hbr_spin))
    metadata = {
        "input_file": str(ROOT / settings["samples"]),
        "input_sample_count": int(len(samples)),
        "retained_sample_count_per_model": int(len(bmr_mass)),
        "mixture_sample_count": int(len(mixture_mass)),
        "data_cut_time_M": float(settings["cut_time_M"]),
        "ringdown_data_used": False,
        "conditional_on_GR_remnant_fits": True,
        "models": [
            "Barausse_Morozova_Rezzolla_2012_mass_plus_Hofmann_Barausse_Rezzolla_2016_spin",
            "NRSur3dq8Remnant_plus_PESummary_leading_in_plane_projection",
        ],
        "baseline_mass_detector_msun": mass0,
        "baseline_spin": spin0,
        "mean_offset_delta_lnM_delta_chi": mean_offset.tolist(),
        "covariance": covariance.tolist(),
        "correlation": float(corr),
        "median_between_model_mass_difference_msun": mass_model_median_difference,
        "median_between_model_spin_difference": spin_model_median_difference,
    }
    (output / "independent_prior_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )

    mass_summary = summary[0]
    spin_summary = summary[1]
    report = [
        "# Pre−40M independent remnant prior",
        "",
        "This prior uses only the public GW250114 posterior obtained after removing the signal later than −40M. No ringdown posterior or post-peak strain enters the construction.",
        "",
        f"- retained progenitor samples per fit family: `{len(bmr_mass)}`",
        f"- detector-frame final mass: `{mass_summary['median']:.3f}` Msun (90% `{mass_summary['q05']:.3f}–{mass_summary['q95']:.3f}` Msun)",
        f"- final spin: `{spin_summary['median']:.4f}` (90% `{spin_summary['q05']:.4f}–{spin_summary['q95']:.4f}`)",
        f"- sigma(delta ln M): `{sigma[0]:.5f}`",
        f"- sigma(delta chi): `{sigma[1]:.5f}`",
        f"- correlation: `{corr:.4f}`",
        f"- measured prior centre relative to (68.1 Msun, 0.68): `({mean_offset[0]:+.5f}, {mean_offset[1]:+.5f})`",
        f"- median NR-fit-family difference: `{mass_model_median_difference:+.3f}` Msun in mass and `{spin_model_median_difference:+.4f}` in spin",
        "",
        "The prior is statistically independent of the ringdown data but not theory-agnostic: both mappings assume GR binary evolution and GR-calibrated remnant fits. It is therefore appropriate for preventing data reuse while testing ringdown spectral deformations conditional on a GR progenitor map; it must not be advertised as a fully model-independent no-hair test.",
    ]
    (output / "INDEPENDENT_PRIOR_REPORT.md").write_text(
        "\n".join(report) + "\n", encoding="utf-8"
    )
    print("\n".join(report[4:15]))
    print(f"Output: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
