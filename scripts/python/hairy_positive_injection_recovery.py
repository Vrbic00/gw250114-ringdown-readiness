"""Noise-free positive injection--recovery test for the hairy QNM grid.

The calculation asks an intentionally narrow question: at the spectral SNR of
the public 6M ringdown segment, how large must an internally consistent hairy
two-mode signal be before the stationary recovery distinguishes it from Kerr?
It is an Asimov sensitivity benchmark, not a detection-efficiency study over
detector noise or an exact waveform of the underlying field theory.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
from pathlib import Path


import six  # noqa: F401


ROOT = Path(__file__).resolve().parents[2]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "results" / ".matplotlib"))
LOCAL_DEPS = ROOT / ".deps" / "python"
if LOCAL_DEPS.exists():
    sys.path.insert(0, str(LOCAL_DEPS))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from evolving_kerr_bias import time_axis, whiten_vector
from hairy_evolving_kerr_control import (
    build_interpolators,
    load_primary_nr_sur_covariance,
    mode_omega,
    profile_point,
    stationary_hairy_basis,
    summarize_grid,
)


DEFAULT_CONFIG = ROOT / "config" / "hairy_gw250114_publication.json"
DEFAULT_GRID = ROOT / "results" / "hairy_qnm_production_grid" / "hairy_qnm_production_grid.npz"
DEFAULT_PRIOR_SAMPLES = ROOT / "results" / "article1_independent_prior" / "pre_minus40M_remnant_model_crosscheck_samples.csv"
DEFAULT_OUTPUT = ROOT / "results" / "hairy_positive_injection_recovery"


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table: {path}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def spectral_distance_percent(
    interpolators: dict[str, object],
    spin: float,
    alpha: float,
    h0: float,
) -> float:
    shifts: list[float] = []
    for mode in ("220", "221"):
        deformed = complex(mode_omega(interpolators[mode], h0, alpha, np.asarray([spin]))[0])
        kerr = complex(mode_omega(interpolators[mode], h0, 0.0, np.asarray([spin]))[0])
        shifts.extend(
            (
                abs(deformed.real / kerr.real - 1.0),
                abs(deformed.imag / kerr.imag - 1.0),
            )
        )
    return 100.0 * max(shifts)


def save_figure(path_base: Path, rows: list[dict[str, object]]) -> None:
    selected = [row for row in rows if row["case_kind"] == "GW250114_SNR_scan"]
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.1), constrained_layout=True)
    for h0, color, marker in ((0.5, "tab:blue", "o"), (1.0, "tab:green", "X")):
        group = sorted(
            [row for row in selected if np.isclose(float(row["h0_injected"]), h0)],
            key=lambda row: float(row["alpha_injected"]),
        )
        x = np.asarray([float(row["spectral_distance_percent"]) for row in group])
        expected = np.asarray([float(row["maximum_log_likelihood_ratio_vs_gr"]) for row in group])
        log_bayes = np.asarray([float(row["log_bayes_factor_hairy_vs_gr"]) for row in group])
        linestyle = "-" if len(group) > 1 else "none"
        label = rf"$h_0/M={h0:g}$"
        axes[0].plot(x, expected, marker=marker, linestyle=linestyle, color=color, label=label)
        axes[1].plot(x, log_bayes, marker=marker, linestyle=linestyle, color=color, label=label)
    axes[0].set_ylabel(r"expected max. $\Delta\ln\mathcal{L}$ vs. Kerr")
    axes[1].axhline(0.0, color="black", lw=1.0)
    axes[1].set_ylabel(r"Asimov $\ln B_{\rm Li}^{\rm K}$")
    event_log_bayes = np.asarray(
        [float(row["log_bayes_factor_hairy_vs_gr"]) for row in selected], dtype=float
    )
    axes[1].set_ylim(float(np.min(event_log_bayes)) - 0.01, float(np.max(event_log_bayes)) + 0.02)
    axes[1].text(
        0.98,
        0.95,
        r"$\ln 10=2.30$ (off scale)",
        transform=axes[1].transAxes,
        ha="right",
        va="top",
        color="0.35",
        fontsize=8,
    )
    axes[0].legend(fontsize=8)
    axes[1].legend(fontsize=8, loc="lower right")
    for ax in axes:
        ax.set_xlabel(r"injected two-mode spectral distance $d_{\rm QNM}$ [%]")
        ax.grid(alpha=0.25)
    fig.savefig(path_base.with_suffix(".png"), dpi=240)
    fig.savefig(path_base.with_suffix(".pdf"))
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--grid", type=Path, default=DEFAULT_GRID)
    parser.add_argument("--prior-samples", type=Path, default=DEFAULT_PRIOR_SAMPLES)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--figure-only",
        action="store_true",
        help="Regenerate the figure from the existing summary table.",
    )
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    if args.figure_only:
        summary_path = args.output / "positive_injection_summary.csv"
        with summary_path.open(encoding="utf-8", newline="") as handle:
            existing_rows = list(csv.DictReader(handle))
        save_figure(args.output / "hairy_positive_injection_sensitivity", existing_rows)
        print(f"regenerated {args.output / 'hairy_positive_injection_sensitivity.pdf'}")
        return

    config_all = json.loads(args.config.read_text(encoding="utf-8"))
    config = config_all["positive_hair_control"]
    grid = np.load(args.grid)
    interpolators = build_interpolators(grid)
    alpha_axis = np.asarray(grid["alphas"], dtype=float)
    h0_axis = np.asarray(grid["h0s"], dtype=float)
    mass0 = float(config["mass_detector_msun"])
    spin0 = float(config["spin"])
    sample_rate = float(config["sample_rate_hz"])
    times = time_axis(sample_rate, float(config["duration_s"]))
    modes = list(config["modes"])
    base_amplitudes = np.asarray(
        [value for mode in modes for value in config["amplitudes"][mode]], dtype=float
    )
    covariance = load_primary_nr_sur_covariance(args.prior_samples, mass0, spin0)
    inverse_covariance = np.linalg.inv(covariance)

    grid_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    cases = [
        {**injection, "segment_snr": float(config["segment_snr"]), "case_kind": "GW250114_SNR_scan"}
        for injection in config["injections"]
    ]
    ladder = config["snr_ladder"]
    cases.extend(
        {
            "h0": float(ladder["h0"]),
            "alpha": float(ladder["alpha"]),
            "segment_snr": float(snr),
            "case_kind": "SNR_ladder",
        }
        for snr in ladder["values"]
    )
    for injection_index, injection in enumerate(cases):
        injected_h0 = float(injection["h0"])
        injected_alpha = float(injection["alpha"])
        target_snr = float(injection["segment_snr"])
        injection_basis = stationary_hairy_basis(
            modes,
            times,
            mass0,
            spin0,
            injected_alpha,
            injected_h0,
            interpolators,
        )
        raw = injection_basis @ base_amplitudes
        raw_white = whiten_vector(raw, sample_rate)
        scale = target_snr / math.sqrt(float(np.dot(raw_white, raw_white)))
        data_white = whiten_vector(raw * scale, sample_rate)
        segment_snr = math.sqrt(float(np.dot(data_white, data_white)))
        chi2 = np.empty((len(h0_axis), len(alpha_axis)), dtype=float)
        for ih, h0 in enumerate(h0_axis):
            for ia, alpha in enumerate(alpha_axis):
                value, fitted_mass, fitted_spin, success = profile_point(
                    data_white,
                    modes,
                    times,
                    sample_rate,
                    mass0,
                    spin0,
                    float(alpha),
                    float(h0),
                    interpolators,
                    inverse_covariance,
                )
                chi2[ih, ia] = value
                grid_rows.append(
                    {
                        "injection_index": injection_index,
                        "case_kind": str(injection["case_kind"]),
                        "h0_injected": injected_h0,
                        "alpha_injected": injected_alpha,
                        "h0_recovery": float(h0),
                        "alpha_recovery": float(alpha),
                        "profile_chi2": value,
                        "fitted_mass": fitted_mass,
                        "fitted_spin": fitted_spin,
                        "optimizer_success": success,
                    }
                )
        log_likelihood = -0.5 * chi2
        result, posterior = summarize_grid(log_likelihood, alpha_axis, h0_axis)
        best = np.unravel_index(int(np.argmax(log_likelihood)), log_likelihood.shape)
        summary_rows.append(
            {
                "injection_index": injection_index,
                "case_kind": str(injection["case_kind"]),
                "h0_injected": injected_h0,
                "alpha_injected": injected_alpha,
                "spectral_distance_percent": spectral_distance_percent(
                    interpolators, spin0, injected_alpha, injected_h0
                ),
                "segment_snr": segment_snr,
                "recovered_best_h0": float(h0_axis[best[0]]),
                "recovered_best_alpha": float(alpha_axis[best[1]]),
                **result,
            }
        )
        np.savez_compressed(
            args.output / f"positive_injection_{injection_index:02d}.npz",
            h0s=h0_axis,
            alphas=alpha_axis,
            profile_chi2=chi2,
            log_likelihood=log_likelihood,
            posterior_mass=posterior,
        )

    write_csv(args.output / "positive_injection_profile_grid.csv", grid_rows)
    write_csv(args.output / "positive_injection_summary.csv", summary_rows)
    save_figure(args.output / "hairy_positive_injection_sensitivity", summary_rows)
    event_rows = [row for row in summary_rows if row["case_kind"] == "GW250114_SNR_scan"]
    ordered = sorted(event_rows, key=lambda row: float(row["spectral_distance_percent"]))
    positive_evidence = [row for row in ordered if float(row["log_bayes_factor_hairy_vs_gr"]) > 0.0]
    ten_to_one = [row for row in ordered if float(row["log_bayes_factor_hairy_vs_gr"]) >= math.log(10.0)]
    output = {
        "scope": "noise_free_Asimov_two_mode_spectral_injection_recovery",
        "segment_snr": float(config["segment_snr"]),
        "primary_remnant_covariance": "pre_minus40M_NRSur3dq8",
        "injections": summary_rows,
        "smallest_tested_spectral_distance_with_positive_log_bayes": None if not positive_evidence else float(positive_evidence[0]["spectral_distance_percent"]),
        "smallest_tested_spectral_distance_with_bayes_factor_at_least_10": None if not ten_to_one else float(ten_to_one[0]["spectral_distance_percent"]),
        "snr_ladder": [row for row in summary_rows if row["case_kind"] == "SNR_ladder"],
        "all_optimizers_successful": all(bool(row["optimizer_success"]) for row in grid_rows),
        "limitation": "No detector-noise realizations and no theory-specific excitation or nonlinear waveform systematics.",
    }
    (args.output / "positive_injection_summary.json").write_text(
        json.dumps(output, indent=2), encoding="utf-8"
    )
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
