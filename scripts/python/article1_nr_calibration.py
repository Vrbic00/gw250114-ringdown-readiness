"""Calibrate the evolving-remnant envelope on SXS:BBH:3617.

The waveform peak fixes t=0.  Common-horizon Christodoulou mass and spin
magnitude are divided by their late-time medians and fitted with a signed
exponential.  Multiple fit starts are retained so that the quoted calibration
does not hide early-time coordinate or relaxation sensitivity.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
LOCAL_DEPS = ROOT / ".deps" / "python"
if LOCAL_DEPS.exists():
    sys.path.insert(0, str(LOCAL_DEPS))
sys.path.insert(0, str(ROOT / "scripts" / "python"))

import h5py
import numpy as np
from scipy.optimize import least_squares

import evolving_kerr_bias as ekb


def fit_exponential(
    time_mf: np.ndarray,
    fractional_offset: np.ndarray,
    start_mf: float,
    end_mf: float,
    *,
    fixed_amplitude: float | None = None,
) -> dict[str, float]:
    keep = (time_mf >= start_mf) & (time_mf <= end_mf)
    x = np.asarray(time_mf[keep], dtype=float)
    y = np.asarray(fractional_offset[keep], dtype=float)
    if len(x) < 10:
        raise ValueError("Too few horizon samples in the requested fit window")

    if fixed_amplitude is None:
        initial_amplitude = float(y[0] * math.exp(start_mf / 6.0))

        def residual(parameters: np.ndarray) -> np.ndarray:
            amplitude, tau, offset = parameters
            return amplitude * np.exp(-x / tau) + offset - y

        fit = least_squares(
            residual,
            np.asarray([initial_amplitude, 6.0, 0.0]),
            bounds=([-0.1, 0.1, -0.01], [0.1, 100.0, 0.01]),
            ftol=1e-14,
            xtol=1e-14,
            gtol=1e-14,
            max_nfev=5000,
        )
        amplitude, tau, offset = map(float, fit.x)
    else:

        def residual(parameters: np.ndarray) -> np.ndarray:
            tau = parameters[0]
            return fixed_amplitude * np.exp(-x / tau) - y

        fit = least_squares(
            residual,
            np.asarray([6.0]),
            bounds=([0.1], [100.0]),
            ftol=1e-14,
            xtol=1e-14,
            gtol=1e-14,
            max_nfev=5000,
        )
        amplitude = float(fixed_amplitude)
        tau = float(fit.x[0])
        offset = 0.0

    model = amplitude * np.exp(-x / tau) + offset
    rms = float(np.sqrt(np.mean((model - y) ** 2)))
    maximum = float(np.max(np.abs(model - y)))
    return {
        "fit_start_Mf": float(start_mf),
        "fit_end_Mf": float(end_mf),
        "amplitude_fractional": amplitude,
        "tau_Mf": tau,
        "constant_offset_fractional": offset,
        "rms_fractional": rms,
        "max_abs_residual_fractional": maximum,
        "n_samples": int(len(x)),
        "success": int(bool(fit.success)),
    }


def write_rows(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config", type=Path, default=ROOT / "config" / "article1_nr_validation.json"
    )
    parser.add_argument(
        "--output", type=Path, default=ROOT / "results" / "article1_nr_calibration"
    )
    args = parser.parse_args()
    config = json.loads(args.config.read_text(encoding="utf-8"))
    sxs_config = config["sxs"]
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)

    metadata = json.loads((ROOT / sxs_config["metadata"]).read_text(encoding="utf-8"))
    with h5py.File(ROOT / sxs_config["horizons"], "r") as horizons:
        mass_data = np.asarray(horizons["AhC.dir/ChristodoulouMass.dat"])
        spin_data = np.asarray(horizons["AhC.dir/chiMagInertial.dat"])

    # The packed strain file is decoded by sxs.  Importing lazily keeps the
    # horizon-only part of this script transparent and easy to audit.
    import sxs

    waveform = sxs.load(str(ROOT / sxs_config["strain"]))
    h22 = np.asarray(waveform[:, waveform.index(2, 2)])
    waveform_time = np.asarray(waveform.t, dtype=float)
    peak_index = int(np.argmax(np.abs(h22)))
    peak_time = float(waveform_time[peak_index])

    mass_final = float(np.median(mass_data[-500:, 1]))
    spin_final = float(np.median(spin_data[-500:, 1]))
    time_mass_mf = (mass_data[:, 0] - peak_time) / mass_final
    time_spin_mf = (spin_data[:, 0] - peak_time) / mass_final
    mass_offset = mass_data[:, 1] / mass_final - 1.0
    spin_offset = spin_data[:, 1] / spin_final - 1.0
    epsilon_mass_peak = float(np.interp(0.0, time_mass_mf, mass_offset))
    epsilon_spin_peak = float(np.interp(0.0, time_spin_mf, spin_offset))

    fit_start, fit_end = map(float, sxs_config["fit_window_Mf"])
    primary_mass_fit = fit_exponential(
        time_mass_mf,
        mass_offset,
        fit_start,
        fit_end,
        fixed_amplitude=epsilon_mass_peak,
    )
    primary_spin_fit = fit_exponential(
        time_spin_mf,
        spin_offset,
        fit_start,
        fit_end,
        fixed_amplitude=epsilon_spin_peak,
    )
    fits: list[dict] = []
    for quantity, time_values, offsets in (
        ("mass", time_mass_mf, mass_offset),
        ("spin", time_spin_mf, spin_offset),
    ):
        for start in sxs_config["fit_start_checks_Mf"]:
            row = fit_exponential(time_values, offsets, float(start), fit_end)
            row["quantity"] = quantity
            row["fit_kind"] = "free_amplitude_plus_offset"
            fits.append(row)
    for quantity, row in (("mass", primary_mass_fit), ("spin", primary_spin_fit)):
        fixed = dict(row)
        fixed["quantity"] = quantity
        fixed["fit_kind"] = "peak_amplitude_fixed_zero_offset"
        fits.append(fixed)
    fits.sort(key=lambda row: (row["quantity"], row["fit_kind"], row["fit_start_Mf"]))
    write_rows(output / "sxs3617_exponential_fits.csv", fits)

    end_track = float(sxs_config["track_end_Mf"])
    common_time = np.linspace(0.0, end_track, 601)
    track_rows = []
    for time_value in common_time:
        track_rows.append(
            {
                "time_from_h22_peak_Mf": float(time_value),
                "mass_fraction_of_initial_total": float(
                    np.interp(time_value, time_mass_mf, mass_data[:, 1])
                ),
                "spin_magnitude": float(np.interp(time_value, time_spin_mf, spin_data[:, 1])),
                "epsilon_mass_fractional": float(
                    np.interp(time_value, time_mass_mf, mass_offset)
                ),
                "epsilon_spin_fractional": float(
                    np.interp(time_value, time_spin_mf, spin_offset)
                ),
                "mass_exponential_fit": float(
                    epsilon_mass_peak * np.exp(-time_value / primary_mass_fit["tau_Mf"])
                ),
                "spin_exponential_fit": float(
                    epsilon_spin_peak * np.exp(-time_value / primary_spin_fit["tau_Mf"])
                ),
            }
        )
    write_rows(output / "sxs3617_horizon_track.csv", track_rows)

    calibrated = {
        "simulation": sxs_config["simulation"],
        "level": sxs_config["level"],
        "waveform_extrapolation_order": 2,
        "time_origin": "peak_abs_h22_at_future_null_infinity",
        "waveform_peak_coordinate_time_Mtotal": peak_time,
        "common_horizon_time_Mtotal": float(metadata["common_horizon_time"]),
        "peak_minus_common_horizon_time_Mtotal": peak_time
        - float(metadata["common_horizon_time"]),
        "mass_final_fraction_of_initial_total_late_median": mass_final,
        "mass_final_fraction_metadata": float(metadata["remnant_mass"]),
        "spin_final_late_median": spin_final,
        "spin_final_metadata": float(metadata["remnant_dimensionless_spin"][2]),
        "epsilon_mass_at_h22_peak": epsilon_mass_peak,
        "epsilon_spin_at_h22_peak": epsilon_spin_peak,
        "tau_mass_Mf": primary_mass_fit["tau_Mf"],
        "tau_spin_Mf": primary_spin_fit["tau_Mf"],
        "fit_window_Mf": [fit_start, fit_end],
        "calibration_scope": "equal_mass_nonspinning_SXS_system_representative_of_GW250114",
    }
    (output / "nr_calibrated_drift.json").write_text(
        json.dumps(calibrated, indent=2) + "\n", encoding="utf-8"
    )

    plotted_rows = track_rows[::50]
    series = {
        "NR mass": [(row["time_from_h22_peak_Mf"], 100.0 * row["epsilon_mass_fractional"]) for row in plotted_rows],
        "mass exponential": [(row["time_from_h22_peak_Mf"], 100.0 * row["mass_exponential_fit"]) for row in plotted_rows],
        "NR spin": [(row["time_from_h22_peak_Mf"], 100.0 * row["epsilon_spin_fractional"]) for row in plotted_rows],
        "spin exponential": [(row["time_from_h22_peak_Mf"], 100.0 * row["spin_exponential_fit"]) for row in plotted_rows],
    }
    ekb._save_line_plot(
        output / "sxs3617_remnant_relaxation.png",
        series,
        title="SXS:BBH:3617 common-horizon relaxation",
        xlabel="time from |h22| peak / Mf",
        ylabel="fractional offset from late remnant [%]",
        thresholds=(0.0,),
    )

    report = [
        "# NR calibration of the evolving-remnant envelope",
        "",
        "The calibration uses the common horizon of SXS:BBH:3617 (Lev3) and sets time zero at the peak of the extrapolated N=2 |h22| waveform.",
        "",
        f"- late-time mass fraction: `{mass_final:.10f}` (metadata `{metadata['remnant_mass']:.10f}`)",
        f"- late-time spin: `{spin_final:.10f}` (metadata `{metadata['remnant_dimensionless_spin'][2]:.10f}`)",
        f"- mass offset at the waveform peak: `{100.0 * epsilon_mass_peak:.4f}%`",
        f"- spin offset at the waveform peak: `{100.0 * epsilon_spin_peak:.4f}%`",
        f"- fixed-amplitude mass relaxation time: `{primary_mass_fit['tau_Mf']:.3f} Mf`",
        f"- fixed-amplitude spin relaxation time: `{primary_spin_fit['tau_Mf']:.3f} Mf`",
        f"- horizon formation precedes the waveform peak by `{calibrated['peak_minus_common_horizon_time_Mtotal']:.3f} Mtotal`",
        "",
        "The result calibrates one representative binary, not a population-wide envelope. The full NR-waveform injection is therefore the primary model-misspecification test; the fitted exponential is retained as an interpretable controlled surrogate.",
    ]
    (output / "NR_CALIBRATION_REPORT.md").write_text(
        "\n".join(report) + "\n", encoding="utf-8"
    )
    print(json.dumps(calibrated, indent=2))
    print(f"Output: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
