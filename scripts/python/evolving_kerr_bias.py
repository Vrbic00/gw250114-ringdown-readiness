"""Evolving-Kerr-spectrum misspecification study for article 1.

This script generates controlled GR injections whose instantaneous Kerr QNM
spectrum follows relaxing effective mass and spin functions.  It then recovers
those injections with a stationary Kerr+EFT model and measures the false EFT
bias after profiling linear mode amplitudes, remnant mass, and remnant spin.

The injected model is an adiabatic waveform surrogate.  It is not presented as
an exact time-dependent Kerr solution or as a numerical-relativity waveform.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[2]
LOCAL_DEPS = ROOT / ".deps" / "python"
if LOCAL_DEPS.exists():
    sys.path.insert(0, str(LOCAL_DEPS))

import numpy as np
from PIL import Image, ImageDraw, ImageFont
from scipy.optimize import least_squares, minimize


SOLAR_MASS_TIME_SECONDS = 4.925490947e-6

QNM_COEFFICIENTS = {
    "220": (1.5251, -1.1568, 0.1292, 0.7000, 1.4187, -0.4990),
    "221": (1.3673, -1.0260, 0.1628, 0.1000, 0.5436, -0.4731),
    "222": (1.3223, -1.0257, 0.1860, -0.1000, 0.4206, -0.4256),
    "330": (1.8956, -1.3043, 0.1818, 0.9000, 2.3430, -0.4810),
    "440": (2.3000, -1.5056, 0.2244, 1.1929, 3.1191, -0.4825),
}


@dataclass(frozen=True)
class Drift:
    name: str
    epsilon_mass: float
    epsilon_spin: float
    tau_mass_m: float
    tau_spin_m: float


@dataclass
class LinearResult:
    identifiable: bool
    alpha_hat: float
    sigma_alpha: float
    bias_sigma: float
    alpha_information: float
    retained_information_fraction: float
    chi2_nuisance: float
    chi2_best: float
    segment_snr: float
    nuisance_rank: int
    nuisance_condition: float


def load_config(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def load_anchors(path: Path, mass0: float, spin0: float) -> dict[str, complex]:
    anchors: dict[str, complex] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            if abs(float(row["mass_detector_msun"]) - mass0) > 1e-9:
                continue
            if abs(float(row["spin"]) - spin0) > 1e-9:
                continue
            mode = row["mode"]
            numerical = complex(
                float(row["qnm_Momega_real"]), float(row["qnm_Momega_imag"])
            )
            anchors[mode] = numerical - berti_momega(mode, spin0)
    return anchors


def load_eft_coefficients(path: Path) -> dict[tuple[str, str, str], list[tuple[int, complex]]]:
    coefficients: dict[tuple[str, str, str], list[tuple[int, complex]]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            key = (row["mode"], row["operator"], row["branch"])
            coefficients.setdefault(key, []).append(
                (
                    int(row["k"]),
                    complex(float(row["coefficient_re"]), float(row["coefficient_im"])),
                )
            )
    for values in coefficients.values():
        values.sort(key=lambda item: item[0])
    return coefficients


def load_remnant_covariance(path: Path) -> np.ndarray:
    labels = {"delta_lnM": 0, "delta_chi": 1}
    covariance = np.zeros((2, 2), dtype=float)
    with path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            covariance[labels[row["row"]], labels[row["column"]]] = float(row["covariance"])
    if not np.allclose(covariance, covariance.T, rtol=0.0, atol=1e-15):
        raise ValueError("Remnant covariance is not symmetric")
    if np.min(np.linalg.eigvalsh(covariance)) <= 0.0:
        raise ValueError("Remnant covariance is not positive definite")
    return covariance


def remnant_prior_covariances(config: dict) -> dict[str, np.ndarray | None]:
    relative_path = Path(config["remnant_prior_covariance_file"])
    base = load_remnant_covariance(ROOT / relative_path)
    output: dict[str, np.ndarray | None] = {}
    for name, width_scale in config["remnant_prior_scenarios"].items():
        output[name] = None if width_scale is None else base * float(width_scale) ** 2
    return output


def berti_momega(mode: str, spin: float) -> complex:
    f1, f2, f3, q1, q2, q3 = QNM_COEFFICIENTS[mode]
    omega_r = f1 + f2 * (1.0 - spin) ** f3
    quality = q1 + q2 * (1.0 - spin) ** q3
    omega_i = -omega_r / (2.0 * quality)
    return complex(omega_r, omega_i)


def kerr_momega(mode: str, spin: float, anchors: dict[str, complex]) -> complex:
    return berti_momega(mode, spin) + anchors.get(mode, 0.0j)


def eft_shift(
    mode: str,
    operator: str,
    branch: str,
    spin: float,
    coefficients: dict[tuple[str, str, str], list[tuple[int, complex]]],
) -> complex:
    values = coefficients.get((mode, operator, branch))
    if not values:
        raise KeyError(f"Missing EFT coefficients for {(mode, operator, branch)}")
    return sum(coefficient * spin**order for order, coefficient in values)


def mode_parameters(
    mode: str,
    mass_solar: np.ndarray | float,
    spin: np.ndarray | float,
    anchors: dict[str, complex],
    coefficients: dict[tuple[str, str, str], list[tuple[int, complex]]],
    operator: str = "kerr",
    branch: str = "none",
    alpha: float = 0.0,
) -> tuple[np.ndarray, np.ndarray]:
    spin_array = np.asarray(spin, dtype=float)
    mass_array = np.asarray(mass_solar, dtype=float)
    f1, f2, f3, q1, q2, q3 = QNM_COEFFICIENTS[mode]
    omega_r = f1 + f2 * (1.0 - spin_array) ** f3 + anchors.get(mode, 0.0j).real
    quality = q1 + q2 * (1.0 - spin_array) ** q3
    omega_i = -(f1 + f2 * (1.0 - spin_array) ** f3) / (2.0 * quality)
    omega_i = omega_i + anchors.get(mode, 0.0j).imag
    if operator != "kerr" and alpha != 0.0:
        shift_values = np.zeros_like(spin_array, dtype=complex)
        for order, coefficient in coefficients[(mode, operator, branch)]:
            shift_values = shift_values + coefficient * spin_array**order
        omega_r = omega_r + alpha * shift_values.real
        omega_i = omega_i + alpha * shift_values.imag
    if np.any(mass_array <= 0.0) or np.any(spin_array <= 0.0) or np.any(spin_array >= 0.99):
        raise ValueError("Mass or spin outside the configured physical domain")
    if np.any(omega_r <= 0.0) or np.any(omega_i >= 0.0):
        raise ValueError("Invalid QNM frequency or damping sign")
    mass_seconds = mass_array * SOLAR_MASS_TIME_SECONDS
    frequency_hz = omega_r / (2.0 * np.pi * mass_seconds)
    tau_s = -mass_seconds / omega_i
    return np.asarray(frequency_hz, dtype=float), np.asarray(tau_s, dtype=float)


def time_axis(sample_rate: float, duration: float) -> np.ndarray:
    count = int(round(duration * sample_rate)) + 1
    return np.arange(count, dtype=float) / sample_rate


def stationary_basis(
    modes: list[str],
    times: np.ndarray,
    mass: float,
    spin: float,
    anchors: dict[str, complex],
    coefficients: dict[tuple[str, str, str], list[tuple[int, complex]]],
    operator: str = "kerr",
    branch: str = "none",
    alpha: float = 0.0,
) -> np.ndarray:
    columns: list[np.ndarray] = []
    for mode in modes:
        frequency, tau = mode_parameters(
            mode,
            mass,
            spin,
            anchors,
            coefficients,
            operator,
            branch,
            alpha,
        )
        decay = np.exp(-times / float(tau))
        phase = 2.0 * np.pi * float(frequency) * times
        columns.extend((decay * np.cos(phase), decay * np.sin(phase)))
    return np.column_stack(columns)


def coefficient_vector(modes: list[str], amplitudes: dict[str, list[float]]) -> np.ndarray:
    return np.asarray(
        [value for mode in modes for value in amplitudes.get(mode, [0.0, 0.0])],
        dtype=float,
    )


def cumulative_trapezoid(values: np.ndarray, times: np.ndarray) -> np.ndarray:
    output = np.zeros_like(values, dtype=float)
    if len(values) > 1:
        output[1:] = np.cumsum(0.5 * (values[:-1] + values[1:]) * np.diff(times))
    return output


def dynamic_basis(
    modes: list[str],
    segment_times: np.ndarray,
    analysis_start_s: float,
    mass_final: float,
    spin_final: float,
    drift: Drift,
    anchors: dict[str, complex],
    coefficients: dict[tuple[str, str, str], list[tuple[int, complex]]],
) -> np.ndarray:
    final_mass_time = mass_final * SOLAR_MASS_TIME_SECONDS
    maximum_time = analysis_start_s + float(segment_times[-1])
    base_step = float(segment_times[1] - segment_times[0])
    integration_times = np.arange(0.0, maximum_time + base_step, base_step)
    tau_mass = drift.tau_mass_m * final_mass_time
    tau_spin = drift.tau_spin_m * final_mass_time
    mass_values = mass_final * (
        1.0 + drift.epsilon_mass * np.exp(-integration_times / tau_mass)
    )
    spin_values = spin_final * (
        1.0 + drift.epsilon_spin * np.exp(-integration_times / tau_spin)
    )
    target_times = analysis_start_s + segment_times
    columns: list[np.ndarray] = []
    for mode in modes:
        frequency, tau = mode_parameters(
            mode,
            mass_values,
            spin_values,
            anchors,
            coefficients,
        )
        phase = cumulative_trapezoid(2.0 * np.pi * frequency, integration_times)
        log_decay = cumulative_trapezoid(1.0 / tau, integration_times)
        phase_target = np.interp(target_times, integration_times, phase)
        decay_target = np.exp(-np.interp(target_times, integration_times, log_decay))
        columns.extend(
            (decay_target * np.cos(phase_target), decay_target * np.sin(phase_target))
        )
    return np.column_stack(columns)


def noise_psd(frequency: np.ndarray) -> np.ndarray:
    safe = np.maximum(frequency, 10.0)
    return 1.0 + (70.0 / safe) ** 4 + (frequency / 500.0) ** 2


def whiten_vector(series: np.ndarray, sample_rate: float) -> np.ndarray:
    transformed = np.fft.fft(np.asarray(series, dtype=float))
    count = len(series)
    positive = np.arange(1, count // 2 + 1)
    frequencies = sample_rate * positive / count
    factor = np.sqrt(4.0 * (sample_rate / count) / noise_psd(frequencies))
    weighted = factor * transformed[positive]
    return np.concatenate((weighted.real, weighted.imag))


def whiten_columns(columns: np.ndarray, sample_rate: float) -> np.ndarray:
    return np.column_stack(
        [whiten_vector(columns[:, index], sample_rate) for index in range(columns.shape[1])]
    )


def project_out(vector: np.ndarray, columns: np.ndarray) -> tuple[np.ndarray, int, float]:
    if columns.size == 0:
        return vector.copy(), 0, 1.0
    coefficients, _, rank, singular = np.linalg.lstsq(columns, vector, rcond=1e-12)
    residual = vector - columns @ coefficients
    condition = math.inf
    if len(singular) and singular[-1] > 0.0:
        condition = float(singular[0] / singular[-1])
    return residual, int(rank), condition


def scale_for_peak_snr(
    modes: list[str],
    times: np.ndarray,
    amplitudes: dict[str, list[float]],
    mass: float,
    spin: float,
    target_snr: float,
    sample_rate: float,
    anchors: dict[str, complex],
    coefficients: dict[tuple[str, str, str], list[tuple[int, complex]]],
) -> float:
    basis = stationary_basis(modes, times, mass, spin, anchors, coefficients)
    waveform = basis @ coefficient_vector(modes, amplitudes)
    norm = float(np.dot(whiten_vector(waveform, sample_rate), whiten_vector(waveform, sample_rate)))
    return target_snr / math.sqrt(norm)


def linear_bias_result(
    modes: list[str],
    amplitudes: dict[str, list[float]],
    times: np.ndarray,
    sample_rate: float,
    scale: float,
    mass: float,
    spin: float,
    drift: Drift,
    analysis_start_s: float,
    operator: str,
    branch: str,
    derivative_steps: dict[str, float],
    relative_tolerance: float,
    anchors: dict[str, complex],
    coefficients: dict[tuple[str, str, str], list[tuple[int, complex]]],
    added_whitened_noise: np.ndarray | None = None,
    prior_covariance: np.ndarray | None = None,
    prior_mean_offset: np.ndarray | None = None,
) -> LinearResult:
    injection_coefficients = scale * coefficient_vector(modes, amplitudes)
    base_basis = stationary_basis(modes, times, mass, spin, anchors, coefficients)
    dynamic_waveform = dynamic_basis(
        modes,
        times,
        analysis_start_s,
        mass,
        spin,
        drift,
        anchors,
        coefficients,
    ) @ injection_coefficients
    dynamic_white = whiten_vector(dynamic_waveform, sample_rate)
    amplitude_columns = whiten_columns(base_basis, sample_rate)
    coefficient_values, _, _, _ = np.linalg.lstsq(
        amplitude_columns, dynamic_white, rcond=1e-12
    )
    base_white = amplitude_columns @ coefficient_values
    residual = dynamic_white - base_white
    if added_whitened_noise is not None:
        residual = residual + added_whitened_noise

    h_mass = derivative_steps["ln_mass"]
    h_spin = derivative_steps["spin"]
    h_alpha = derivative_steps["alpha"]

    mass_plus = stationary_basis(
        modes, times, mass * math.exp(h_mass), spin, anchors, coefficients
    ) @ coefficient_values
    mass_minus = stationary_basis(
        modes, times, mass * math.exp(-h_mass), spin, anchors, coefficients
    ) @ coefficient_values
    derivative_mass = (
        whiten_vector(mass_plus, sample_rate) - whiten_vector(mass_minus, sample_rate)
    ) / (2.0 * h_mass)

    spin_plus = stationary_basis(
        modes, times, mass, spin + h_spin, anchors, coefficients
    ) @ coefficient_values
    spin_minus = stationary_basis(
        modes, times, mass, spin - h_spin, anchors, coefficients
    ) @ coefficient_values
    derivative_spin = (
        whiten_vector(spin_plus, sample_rate) - whiten_vector(spin_minus, sample_rate)
    ) / (2.0 * h_spin)

    alpha_plus = stationary_basis(
        modes,
        times,
        mass,
        spin,
        anchors,
        coefficients,
        operator,
        branch,
        h_alpha,
    ) @ coefficient_values
    alpha_minus = stationary_basis(
        modes,
        times,
        mass,
        spin,
        anchors,
        coefficients,
        operator,
        branch,
        -h_alpha,
    ) @ coefficient_values
    derivative_alpha = (
        whiten_vector(alpha_plus, sample_rate) - whiten_vector(alpha_minus, sample_rate)
    ) / (2.0 * h_alpha)

    if prior_covariance is not None:
        precision = np.linalg.inv(np.asarray(prior_covariance, dtype=float))
        prior_whitener = np.linalg.cholesky(precision).T
        prior_mean = (
            np.zeros(2, dtype=float)
            if prior_mean_offset is None
            else np.asarray(prior_mean_offset, dtype=float)
        )
        if prior_mean.shape != (2,):
            raise ValueError("prior_mean_offset must have shape (2,)")
        # The two pseudo-observations encode a Gaussian prior in
        # (delta ln M, delta chi).  A non-zero vector permits a genuinely
        # independent pre-merger posterior to retain its measured centre,
        # rather than silently recentering it on the injected remnant.
        residual = np.concatenate((residual, prior_whitener @ prior_mean))
        amplitude_columns = np.vstack(
            (amplitude_columns, np.zeros((2, amplitude_columns.shape[1])))
        )
        derivative_mass = np.concatenate((derivative_mass, prior_whitener[:, 0]))
        derivative_spin = np.concatenate((derivative_spin, prior_whitener[:, 1]))
        derivative_alpha = np.concatenate((derivative_alpha, np.zeros(2)))
    nuisance = np.column_stack((amplitude_columns, derivative_mass, derivative_spin))
    residual_perp, rank, condition = project_out(residual, nuisance)
    alpha_perp, _, _ = project_out(derivative_alpha, nuisance)
    information = float(np.dot(alpha_perp, alpha_perp))
    raw_information = float(np.dot(derivative_alpha, derivative_alpha))
    retained = information / raw_information if raw_information > 0.0 else 0.0
    identifiable = information > max(1e-12, relative_tolerance * raw_information)
    segment_snr = float(np.linalg.norm(dynamic_white))
    chi2_nuisance = float(np.dot(residual_perp, residual_perp))
    if not identifiable:
        return LinearResult(
            False,
            math.nan,
            math.inf,
            math.nan,
            information,
            retained,
            chi2_nuisance,
            chi2_nuisance,
            segment_snr,
            rank,
            condition,
        )
    alpha_hat = float(np.dot(alpha_perp, residual_perp) / information)
    sigma_alpha = 1.0 / math.sqrt(information)
    bias_sigma = alpha_hat / sigma_alpha
    best_residual = residual_perp - alpha_hat * alpha_perp
    chi2_best = float(np.dot(best_residual, best_residual))
    return LinearResult(
        True,
        alpha_hat,
        sigma_alpha,
        bias_sigma,
        information,
        retained,
        chi2_nuisance,
        chi2_best,
        segment_snr,
        rank,
        condition,
    )


def drift_scenarios(config: dict) -> list[tuple[str, float, float, float]]:
    scenarios = [("stationary", 0.0, 0.0, 0.0)]
    for family, multipliers in config["drift_families"].items():
        for strength in config["drift_strengths_fractional"]:
            scenarios.append(
                (
                    family,
                    float(strength),
                    float(strength) * float(multipliers[0]),
                    float(strength) * float(multipliers[1]),
                )
            )
    return scenarios


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"No rows available for {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def result_row(
    result: LinearResult,
    *,
    family: str,
    strength: float,
    drift: Drift,
    relaxation_name: str,
    start_m: float,
    mode_set: str,
    amplitude_scenario: str,
    operator: str,
    branch: str,
    prior_scenario: str = "free_remnant",
) -> dict:
    return {
        "drift_family": family,
        "drift_strength_fractional": strength,
        "epsilon_mass": drift.epsilon_mass,
        "epsilon_spin_fractional": drift.epsilon_spin,
        "relaxation_pair": relaxation_name,
        "tau_mass_M": drift.tau_mass_m,
        "tau_spin_M": drift.tau_spin_m,
        "start_time_M": start_m,
        "mode_set": mode_set,
        "amplitude_scenario": amplitude_scenario,
        "operator": operator,
        "branch": branch,
        "prior_scenario": prior_scenario,
        "identifiable": int(result.identifiable),
        "alpha_hat": result.alpha_hat,
        "sigma_alpha": result.sigma_alpha,
        "bias_sigma": result.bias_sigma,
        "alpha_information": result.alpha_information,
        "retained_information_fraction": result.retained_information_fraction,
        "chi2_nuisance": result.chi2_nuisance,
        "chi2_best": result.chi2_best,
        "delta_chi2_eft": result.chi2_nuisance - result.chi2_best,
        "segment_snr": result.segment_snr,
        "nuisance_rank": result.nuisance_rank,
        "nuisance_condition": result.nuisance_condition,
    }


def run_main_scan(
    config: dict,
    anchors: dict[str, complex],
    coefficients: dict[tuple[str, str, str], list[tuple[int, complex]]],
    prior_scenario: str = "free_remnant",
    prior_covariance: np.ndarray | None = None,
) -> list[dict]:
    event = config["event"]
    mass = float(event["mass_detector_msun"])
    spin = float(event["spin"])
    mass_time = mass * SOLAR_MASS_TIME_SECONDS
    sample_rate = float(config["sample_rate_hz"])
    times = time_axis(sample_rate, float(config["duration_s"]))
    operator = config["representative_operator"]
    branch = config["representative_branch"]
    rows: list[dict] = []
    scales: dict[tuple[str, str], float] = {}
    for mode_set, modes in config["mode_sets"].items():
        for amplitude_name, amplitudes in config["amplitude_scenarios"].items():
            scales[(mode_set, amplitude_name)] = scale_for_peak_snr(
                modes,
                times,
                amplitudes,
                mass,
                spin,
                float(config["peak_snr"]),
                sample_rate,
                anchors,
                coefficients,
            )
    for family, strength, epsilon_mass, epsilon_spin in drift_scenarios(config):
        relaxation_items: Iterable[tuple[str, list[float]]]
        if family == "stationary":
            relaxation_items = [("reference", config["relaxation_pairs_M"]["reference"])]
        else:
            relaxation_items = config["relaxation_pairs_M"].items()
        for relaxation_name, relaxation_pair in relaxation_items:
            drift = Drift(
                family,
                epsilon_mass,
                epsilon_spin,
                float(relaxation_pair[0]),
                float(relaxation_pair[1]),
            )
            for start_m in config["start_times_M"]:
                start_s = float(start_m) * mass_time
                for mode_set, modes in config["mode_sets"].items():
                    for amplitude_name, amplitudes in config["amplitude_scenarios"].items():
                        result = linear_bias_result(
                            modes,
                            amplitudes,
                            times,
                            sample_rate,
                            scales[(mode_set, amplitude_name)],
                            mass,
                            spin,
                            drift,
                            start_s,
                            operator,
                            branch,
                            config["derivative_steps"],
                            float(config["identifiability_relative_tolerance"]),
                            anchors,
                            coefficients,
                            prior_covariance=prior_covariance,
                        )
                        rows.append(
                            result_row(
                                result,
                                family=family,
                                strength=strength,
                                drift=drift,
                                relaxation_name=relaxation_name,
                                start_m=float(start_m),
                                mode_set=mode_set,
                                amplitude_scenario=amplitude_name,
                                operator=operator,
                                branch=branch,
                                prior_scenario=prior_scenario,
                            )
                        )
    return rows


def run_snr_scaling(
    config: dict,
    anchors: dict[str, complex],
    coefficients: dict[tuple[str, str, str], list[tuple[int, complex]]],
    priors: dict[str, np.ndarray | None],
) -> list[dict]:
    event = config["event"]
    mass = float(event["mass_detector_msun"])
    spin = float(event["spin"])
    mass_time = mass * SOLAR_MASS_TIME_SECONDS
    sample_rate = float(config["sample_rate_hz"])
    times = time_axis(sample_rate, float(config["duration_s"]))
    amplitude_name = "moderate_221"
    amplitudes = config["amplitude_scenarios"][amplitude_name]
    operator = config["representative_operator"]
    branch = config["representative_branch"]
    scenarios: list[tuple[str, float, Drift, str]] = [
        ("stationary", 0.0, Drift("stationary", 0.0, 0.0, 10.0, 10.0), "reference")
    ]
    for family, relaxation_name in (
        ("mass_up_spin_down", "reference"),
        ("mass_down_spin_up", "slow_both"),
    ):
        multipliers = config["drift_families"][family]
        relaxation = config["relaxation_pairs_M"][relaxation_name]
        for strength in config["drift_strengths_fractional"]:
            scenarios.append(
                (
                    family,
                    float(strength),
                    Drift(
                        family,
                        float(strength) * float(multipliers[0]),
                        float(strength) * float(multipliers[1]),
                        float(relaxation[0]),
                        float(relaxation[1]),
                    ),
                    relaxation_name,
                )
            )
    rows: list[dict] = []
    for peak_snr in config["snr_scaling_values"]:
        for mode_set, modes in config["mode_sets"].items():
            scale = scale_for_peak_snr(
                modes,
                times,
                amplitudes,
                mass,
                spin,
                float(peak_snr),
                sample_rate,
                anchors,
                coefficients,
            )
            for prior_name, prior_covariance in priors.items():
                for family, strength, drift, relaxation_name in scenarios:
                    for start_m in (0.0, 6.0, 12.0):
                        result = linear_bias_result(
                            modes,
                            amplitudes,
                            times,
                            sample_rate,
                            scale,
                            mass,
                            spin,
                            drift,
                            start_m * mass_time,
                            operator,
                            branch,
                            config["derivative_steps"],
                            float(config["identifiability_relative_tolerance"]),
                            anchors,
                            coefficients,
                            prior_covariance=prior_covariance,
                        )
                        row = result_row(
                            result,
                            family=family,
                            strength=strength,
                            drift=drift,
                            relaxation_name=relaxation_name,
                            start_m=start_m,
                            mode_set=mode_set,
                            amplitude_scenario=amplitude_name,
                            operator=operator,
                            branch=branch,
                            prior_scenario=prior_name,
                        )
                        row["peak_snr"] = float(peak_snr)
                        rows.append(row)
    return rows


def make_prior_report(
    output: Path, config: dict, prior_rows: list[dict], snr_rows: list[dict]
) -> None:
    lines = [
        "# Remnant-prior and SNR sensitivity of the evolving-Kerr bias",
        "",
        "The IMR-prior branches import a GR waveform estimate of the final remnant. They are sensitivity studies, not independent ringdown-only tests.",
        "",
    ]
    for prior_name in sorted({row["prior_scenario"] for row in prior_rows}):
        selected = [
            row
            for row in prior_rows
            if row["prior_scenario"] == prior_name
            and row["identifiable"] == 1
            and float(row["drift_strength_fractional"]) > 0.0
        ]
        maximum = max(selected, key=lambda row: abs(float(row["bias_sigma"])))
        lines.extend(
            (
                f"## {prior_name}",
                "",
                f"- maximum absolute bias at reference peak SNR {config['peak_snr']}: `{abs(float(maximum['bias_sigma'])):.3f} sigma`",
                f"- family / strength: `{maximum['drift_family']}` / `{100.0 * float(maximum['drift_strength_fractional']):.2f}%`",
                f"- mode set / start time: `{maximum['mode_set']}` / `{maximum['start_time_M']} M_f`",
                "",
            )
        )
    scaling_dynamic = [
        row
        for row in snr_rows
        if row["identifiable"] == 1 and float(row["drift_strength_fractional"]) > 0.0
    ]
    maximum_scaling = max(scaling_dynamic, key=lambda row: abs(float(row["bias_sigma"])))
    lines.extend(
        (
            "## SNR scaling envelope",
            "",
            f"- largest configured absolute bias: `{abs(float(maximum_scaling['bias_sigma'])):.3f} sigma`",
            f"- peak SNR: `{maximum_scaling['peak_snr']}`",
            f"- prior / mode set: `{maximum_scaling['prior_scenario']}` / `{maximum_scaling['mode_set']}`",
            f"- family / strength / start: `{maximum_scaling['drift_family']}` / `{100.0 * float(maximum_scaling['drift_strength_fractional']):.2f}%` / `{maximum_scaling['start_time_M']} M_f`",
            "",
            "All values are controlled-surrogate diagnostics. In particular, the tight-prior results must not be described as ringdown-only evidence.",
        )
    )
    (output / "evolving_kerr_prior_snr_report.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def run_operator_robustness(
    config: dict,
    anchors: dict[str, complex],
    coefficients: dict[tuple[str, str, str], list[tuple[int, complex]]],
) -> list[dict]:
    event = config["event"]
    mass = float(event["mass_detector_msun"])
    spin = float(event["spin"])
    mass_time = mass * SOLAR_MASS_TIME_SECONDS
    sample_rate = float(config["sample_rate_hz"])
    times = time_axis(sample_rate, float(config["duration_s"]))
    family = "mass_up_spin_down"
    multipliers = config["drift_families"][family]
    strength = 0.025
    relaxation_name = "reference"
    relaxation_pair = config["relaxation_pairs_M"][relaxation_name]
    drift = Drift(
        family,
        strength * float(multipliers[0]),
        strength * float(multipliers[1]),
        float(relaxation_pair[0]),
        float(relaxation_pair[1]),
    )
    rows: list[dict] = []
    for mode_set, modes in config["mode_sets"].items():
        for amplitude_name, amplitudes in config["amplitude_scenarios"].items():
            scale = scale_for_peak_snr(
                modes,
                times,
                amplitudes,
                mass,
                spin,
                float(config["peak_snr"]),
                sample_rate,
                anchors,
                coefficients,
            )
            for operator in config["all_operators"]:
                for branch in config["all_branches"]:
                    for start_m in config["start_times_M"]:
                        result = linear_bias_result(
                            modes,
                            amplitudes,
                            times,
                            sample_rate,
                            scale,
                            mass,
                            spin,
                            drift,
                            float(start_m) * mass_time,
                            operator,
                            branch,
                            config["derivative_steps"],
                            float(config["identifiability_relative_tolerance"]),
                            anchors,
                            coefficients,
                        )
                        rows.append(
                            result_row(
                                result,
                                family=family,
                                strength=strength,
                                drift=drift,
                                relaxation_name=relaxation_name,
                                start_m=float(start_m),
                                mode_set=mode_set,
                                amplitude_scenario=amplitude_name,
                                operator=operator,
                                branch=branch,
                            )
                        )
    return rows


def profile_stationary_chi2(
    theta: np.ndarray,
    target_white: np.ndarray,
    modes: list[str],
    times: np.ndarray,
    sample_rate: float,
    mass0: float,
    spin0: float,
    operator: str,
    branch: str,
    anchors: dict[str, complex],
    coefficients: dict[tuple[str, str, str], list[tuple[int, complex]]],
) -> float:
    mass = mass0 * math.exp(float(theta[0]))
    spin = spin0 + float(theta[1])
    alpha = float(theta[2])
    try:
        basis = stationary_basis(
            modes,
            times,
            mass,
            spin,
            anchors,
            coefficients,
            operator,
            branch,
            alpha,
        )
    except (ValueError, KeyError, FloatingPointError):
        return 1e100
    white_basis = whiten_columns(basis, sample_rate)
    fitted, _, _, _ = np.linalg.lstsq(white_basis, target_white, rcond=1e-12)
    residual = target_white - white_basis @ fitted
    return float(np.dot(residual, residual))


def run_profile_validation(
    config: dict,
    anchors: dict[str, complex],
    coefficients: dict[tuple[str, str, str], list[tuple[int, complex]]],
) -> list[dict]:
    settings = config["profile_validation"]
    event = config["event"]
    mass = float(event["mass_detector_msun"])
    spin = float(event["spin"])
    mass_time = mass * SOLAR_MASS_TIME_SECONDS
    sample_rate = float(config["sample_rate_hz"])
    times = time_axis(sample_rate, float(config["duration_s"]))
    modes = config["mode_sets"][settings["mode_set"]]
    amplitudes = config["amplitude_scenarios"][settings["amplitude_scenario"]]
    scale = scale_for_peak_snr(
        modes,
        times,
        amplitudes,
        mass,
        spin,
        float(config["peak_snr"]),
        sample_rate,
        anchors,
        coefficients,
    )
    coefficient_values = scale * coefficient_vector(modes, amplitudes)
    operator = config["representative_operator"]
    branch = config["representative_branch"]
    family = settings["scenario_family"]
    multipliers = config["drift_families"][family]
    relaxation_name = settings["relaxation_pair"]
    relaxation_pair = config["relaxation_pairs_M"][relaxation_name]
    bounds = [
        tuple(settings["ln_mass_bounds"]),
        tuple(settings["spin_offset_bounds"]),
        tuple(settings["alpha_bounds"]),
    ]
    rows: list[dict] = []
    for strength in settings["strengths_fractional"]:
        drift = Drift(
            family if strength else "stationary",
            float(strength) * float(multipliers[0]),
            float(strength) * float(multipliers[1]),
            float(relaxation_pair[0]),
            float(relaxation_pair[1]),
        )
        for start_m in settings["start_times_M"]:
            start_s = float(start_m) * mass_time
            target = dynamic_basis(
                modes,
                times,
                start_s,
                mass,
                spin,
                drift,
                anchors,
                coefficients,
            ) @ coefficient_values
            target_white = whiten_vector(target, sample_rate)
            linear = linear_bias_result(
                modes,
                amplitudes,
                times,
                sample_rate,
                scale,
                mass,
                spin,
                drift,
                start_s,
                operator,
                branch,
                config["derivative_steps"],
                float(config["identifiability_relative_tolerance"]),
                anchors,
                coefficients,
            )
            initial_alpha = 0.0 if not np.isfinite(linear.alpha_hat) else linear.alpha_hat
            initial_spectral = np.asarray(
                [0.0, 0.0, np.clip(initial_alpha, *bounds[2])], dtype=float
            )
            initial_basis = stationary_basis(
                modes,
                times,
                mass,
                spin,
                anchors,
                coefficients,
                operator,
                branch,
                float(initial_spectral[2]),
            )
            initial_white_basis = whiten_columns(initial_basis, sample_rate)
            initial_amplitudes, _, _, _ = np.linalg.lstsq(
                initial_white_basis, target_white, rcond=1e-12
            )
            initial = np.concatenate((initial_spectral, initial_amplitudes))
            lower = np.asarray([bounds[0][0], bounds[1][0], bounds[2][0]] + [-np.inf] * len(initial_amplitudes))
            upper = np.asarray([bounds[0][1], bounds[1][1], bounds[2][1]] + [np.inf] * len(initial_amplitudes))

            def full_residual(theta: np.ndarray) -> np.ndarray:
                candidate_mass = mass * math.exp(float(theta[0]))
                candidate_spin = spin + float(theta[1])
                candidate_alpha = float(theta[2])
                try:
                    candidate_basis = stationary_basis(
                        modes,
                        times,
                        candidate_mass,
                        candidate_spin,
                        anchors,
                        coefficients,
                        operator,
                        branch,
                        candidate_alpha,
                    )
                except (ValueError, KeyError, FloatingPointError):
                    return np.full_like(target_white, 1e50)
                candidate = whiten_columns(candidate_basis, sample_rate) @ theta[3:]
                return candidate - target_white

            fit = least_squares(
                full_residual,
                initial,
                bounds=(lower, upper),
                method="trf",
                x_scale="jac",
                ftol=1e-12,
                xtol=1e-12,
                gtol=1e-12,
                max_nfev=3000,
            )
            best = fit.x
            alpha_derivative = fit.jac[:, 2]
            nuisance_indices = [0, 1] + list(range(3, fit.jac.shape[1]))
            alpha_perp, local_rank, local_condition = project_out(
                alpha_derivative, fit.jac[:, nuisance_indices]
            )
            local_information = float(np.dot(alpha_perp, alpha_perp))
            nonlinear_sigma = (
                1.0 / math.sqrt(local_information) if local_information > 1e-12 else math.inf
            )
            rows.append(
                {
                    "drift_family": drift.name,
                    "drift_strength_fractional": float(strength),
                    "epsilon_mass": drift.epsilon_mass,
                    "epsilon_spin_fractional": drift.epsilon_spin,
                    "relaxation_pair": relaxation_name,
                    "start_time_M": float(start_m),
                    "linear_alpha_hat": linear.alpha_hat,
                    "linear_sigma_alpha": linear.sigma_alpha,
                    "linear_bias_sigma": linear.bias_sigma,
                    "nonlinear_success": int(bool(fit.success)),
                    "nonlinear_alpha_hat": float(best[2]),
                    "nonlinear_sigma_alpha_local": nonlinear_sigma,
                    "nonlinear_bias_sigma_local": (
                        float(best[2]) / nonlinear_sigma
                        if np.isfinite(nonlinear_sigma) and nonlinear_sigma > 0.0
                        else math.nan
                    ),
                    "best_delta_ln_mass": float(best[0]),
                    "best_delta_spin": float(best[1]),
                    "chi2_best": float(2.0 * fit.cost),
                    "optimizer_message": str(fit.message),
                    "local_profile_information": local_information,
                    "local_nuisance_rank": local_rank,
                    "local_nuisance_condition": local_condition,
                }
            )
    return rows


def run_noise_ensemble(
    config: dict,
    anchors: dict[str, complex],
    coefficients: dict[tuple[str, str, str], list[tuple[int, complex]]],
) -> tuple[list[dict], list[dict]]:
    settings = config["noise"]
    event = config["event"]
    mass = float(event["mass_detector_msun"])
    spin = float(event["spin"])
    mass_time = mass * SOLAR_MASS_TIME_SECONDS
    sample_rate = float(config["sample_rate_hz"])
    times = time_axis(sample_rate, float(config["duration_s"]))
    modes = config["mode_sets"][settings["mode_set"]]
    amplitudes = config["amplitude_scenarios"][settings["amplitude_scenario"]]
    scale = scale_for_peak_snr(
        modes,
        times,
        amplitudes,
        mass,
        spin,
        float(config["peak_snr"]),
        sample_rate,
        anchors,
        coefficients,
    )
    family = settings["scenario_family"]
    multipliers = config["drift_families"][family]
    relaxation_name = settings["relaxation_pair"]
    relaxation_pair = config["relaxation_pairs_M"][relaxation_name]
    operator = config["representative_operator"]
    branch = config["representative_branch"]
    rng = np.random.default_rng(int(settings["random_seed"]))
    rows: list[dict] = []
    summary: list[dict] = []
    white_dimension = 2 * (len(times) // 2)
    for strength in settings["strengths_fractional"]:
        drift = Drift(
            family if strength else "stationary",
            float(strength) * float(multipliers[0]),
            float(strength) * float(multipliers[1]),
            float(relaxation_pair[0]),
            float(relaxation_pair[1]),
        )
        for start_m in settings["start_times_M"]:
            scenario_rows: list[dict] = []
            for draw in range(int(settings["draws"])):
                noise = rng.normal(0.0, 1.0, white_dimension)
                result = linear_bias_result(
                    modes,
                    amplitudes,
                    times,
                    sample_rate,
                    scale,
                    mass,
                    spin,
                    drift,
                    float(start_m) * mass_time,
                    operator,
                    branch,
                    config["derivative_steps"],
                    float(config["identifiability_relative_tolerance"]),
                    anchors,
                    coefficients,
                    added_whitened_noise=noise,
                )
                row = {
                    "draw": draw,
                    "drift_family": drift.name,
                    "drift_strength_fractional": float(strength),
                    "start_time_M": float(start_m),
                    "alpha_hat": result.alpha_hat,
                    "sigma_alpha": result.sigma_alpha,
                    "z_alpha": result.bias_sigma,
                    "segment_snr": result.segment_snr,
                    "random_seed": int(settings["random_seed"]),
                }
                rows.append(row)
                scenario_rows.append(row)
            z_values = np.asarray([row["z_alpha"] for row in scenario_rows], dtype=float)
            count_90 = int(np.count_nonzero(np.abs(z_values) > 1.645))
            count_95 = int(np.count_nonzero(np.abs(z_values) > 1.96))
            count_3 = int(np.count_nonzero(np.abs(z_values) > 3.0))
            low_90, high_90 = wilson_interval(count_90, len(z_values))
            low_95, high_95 = wilson_interval(count_95, len(z_values))
            low_3, high_3 = wilson_interval(count_3, len(z_values))
            summary.append(
                {
                    "drift_family": drift.name,
                    "drift_strength_fractional": float(strength),
                    "start_time_M": float(start_m),
                    "draws": len(z_values),
                    "mean_z": float(np.mean(z_values)),
                    "sd_z": float(np.std(z_values, ddof=1)),
                    "median_z": float(np.median(z_values)),
                    "q05_z": float(np.quantile(z_values, 0.05)),
                    "q95_z": float(np.quantile(z_values, 0.95)),
                    "false_positive_abs_z_gt_1p645": count_90 / len(z_values),
                    "false_positive_abs_z_gt_1p645_wilson95_low": low_90,
                    "false_positive_abs_z_gt_1p645_wilson95_high": high_90,
                    "false_positive_abs_z_gt_1p96": count_95 / len(z_values),
                    "false_positive_abs_z_gt_1p96_wilson95_low": low_95,
                    "false_positive_abs_z_gt_1p96_wilson95_high": high_95,
                    "false_positive_abs_z_gt_3": count_3 / len(z_values),
                    "false_positive_abs_z_gt_3_wilson95_low": low_3,
                    "false_positive_abs_z_gt_3_wilson95_high": high_3,
                }
            )
    return rows, summary


def wilson_interval(successes: int, trials: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if trials <= 0:
        return math.nan, math.nan
    fraction = successes / trials
    denominator = 1.0 + z**2 / trials
    centre = (fraction + z**2 / (2.0 * trials)) / denominator
    half_width = (
        z
        * math.sqrt(fraction * (1.0 - fraction) / trials + z**2 / (4.0 * trials**2))
        / denominator
    )
    return max(0.0, centre - half_width), min(1.0, centre + half_width)


def run_convergence(
    config: dict,
    anchors: dict[str, complex],
    coefficients: dict[tuple[str, str, str], list[tuple[int, complex]]],
) -> list[dict]:
    settings = config["convergence"]
    event = config["event"]
    mass = float(event["mass_detector_msun"])
    spin = float(event["spin"])
    mass_time = mass * SOLAR_MASS_TIME_SECONDS
    modes = config["mode_sets"][settings["mode_set"]]
    amplitudes = config["amplitude_scenarios"][settings["amplitude_scenario"]]
    family = settings["scenario_family"]
    multipliers = config["drift_families"][family]
    strength = float(settings["strength_fractional"])
    relaxation_name = settings["relaxation_pair"]
    relaxation_pair = config["relaxation_pairs_M"][relaxation_name]
    drift = Drift(
        family,
        strength * float(multipliers[0]),
        strength * float(multipliers[1]),
        float(relaxation_pair[0]),
        float(relaxation_pair[1]),
    )
    rows: list[dict] = []
    for sample_rate in settings["sample_rates_hz"]:
        for duration in settings["durations_s"]:
            times = time_axis(float(sample_rate), float(duration))
            scale = scale_for_peak_snr(
                modes,
                times,
                amplitudes,
                mass,
                spin,
                float(config["peak_snr"]),
                float(sample_rate),
                anchors,
                coefficients,
            )
            result = linear_bias_result(
                modes,
                amplitudes,
                times,
                float(sample_rate),
                scale,
                mass,
                spin,
                drift,
                float(settings["start_time_M"]) * mass_time,
                config["representative_operator"],
                config["representative_branch"],
                config["derivative_steps"],
                float(config["identifiability_relative_tolerance"]),
                anchors,
                coefficients,
            )
            rows.append(
                {
                    "sample_rate_hz": float(sample_rate),
                    "duration_s": float(duration),
                    "samples": len(times),
                    "alpha_hat": result.alpha_hat,
                    "sigma_alpha": result.sigma_alpha,
                    "bias_sigma": result.bias_sigma,
                    "segment_snr": result.segment_snr,
                    "retained_information_fraction": result.retained_information_fraction,
                }
            )
    return rows


def validation_rows(
    config: dict,
    anchors: dict[str, complex],
    coefficients: dict[tuple[str, str, str], list[tuple[int, complex]]],
) -> list[dict]:
    event = config["event"]
    mass = float(event["mass_detector_msun"])
    spin = float(event["spin"])
    sample_rate = float(config["sample_rate_hz"])
    times = time_axis(sample_rate, float(config["duration_s"]))
    rows: list[dict] = []
    crosscheck_path = ROOT / "results" / "gw250114_kerr_qnm" / "qnm_solver_crosscheck.csv"
    references: dict[str, complex] = {}
    with crosscheck_path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            references[row["mode"]] = complex(
                float(row["qnm_Momega_real"]), float(row["qnm_Momega_imag"])
            )
    for mode in ("220", "221"):
        calculated = kerr_momega(mode, spin, anchors)
        difference = calculated - references[mode]
        rows.append(
            {
                "test": f"anchored_kerr_{mode}",
                "metric": "complex_Momega_abs_difference",
                "value": abs(difference),
                "tolerance": 1e-12,
                "pass": int(abs(difference) < 1e-12),
            }
        )
    modes = config["mode_sets"]["220_221"]
    amplitudes = config["amplitude_scenarios"]["moderate_221"]
    coefficients_vector = coefficient_vector(modes, amplitudes)
    stationary = stationary_basis(modes, times, mass, spin, anchors, coefficients)
    zero_drift = Drift("stationary", 0.0, 0.0, 10.0, 10.0)
    dynamic = dynamic_basis(
        modes, times, 0.0, mass, spin, zero_drift, anchors, coefficients
    )
    difference = np.max(np.abs(stationary @ coefficients_vector - dynamic @ coefficients_vector))
    rows.append(
        {
            "test": "stationary_waveform_limit",
            "metric": "max_abs_time_series_difference",
            "value": float(difference),
            "tolerance": 2e-12,
            "pass": int(difference < 2e-12),
        }
    )
    scale = scale_for_peak_snr(
        modes,
        times,
        amplitudes,
        mass,
        spin,
        float(config["peak_snr"]),
        sample_rate,
        anchors,
        coefficients,
    )
    result = linear_bias_result(
        modes,
        amplitudes,
        times,
        sample_rate,
        scale,
        mass,
        spin,
        zero_drift,
        0.0,
        config["representative_operator"],
        config["representative_branch"],
        config["derivative_steps"],
        float(config["identifiability_relative_tolerance"]),
        anchors,
        coefficients,
    )
    rows.append(
        {
            "test": "stationary_false_eft_bias",
            "metric": "abs_bias_sigma",
            "value": abs(result.bias_sigma),
            "tolerance": 1e-8,
            "pass": int(abs(result.bias_sigma) < 1e-8),
        }
    )
    return rows


def _figure_font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    candidates = ["arialbd.ttf" if bold else "arial.ttf", "DejaVuSans-Bold.ttf" if bold else "DejaVuSans.ttf"]
    for candidate in candidates:
        try:
            return ImageFont.truetype(candidate, size)
        except OSError:
            continue
    return ImageFont.load_default()


def _rgb_interpolate(left: tuple[int, int, int], right: tuple[int, int, int], fraction: float) -> tuple[int, int, int]:
    fraction = min(1.0, max(0.0, fraction))
    return tuple(round(a + fraction * (b - a)) for a, b in zip(left, right))


def _diverging_colour(value: float, limit: float) -> tuple[int, int, int]:
    if not np.isfinite(value):
        return (220, 220, 220)
    scaled = min(1.0, abs(value) / max(limit, 1e-15))
    endpoint = (180, 28, 42) if value >= 0.0 else (37, 87, 167)
    return _rgb_interpolate((250, 250, 248), endpoint, scaled)


def _sequential_colour(value: float) -> tuple[int, int, int]:
    if not np.isfinite(value):
        return (220, 220, 220)
    value = min(1.0, max(0.0, value))
    if value <= 0.5:
        return _rgb_interpolate((38, 68, 120), (49, 173, 176), 2.0 * value)
    return _rgb_interpolate((49, 173, 176), (241, 214, 83), 2.0 * value - 1.0)


def _save_heatmap(
    path: Path,
    matrix: np.ndarray,
    xlabels: list[str],
    ylabels: list[str],
    title: str,
    xlabel: str,
    ylabel: str,
    colourbar_label: str,
    diverging: bool,
) -> None:
    width, height = 1600, 900
    left, top, right, bottom = 230, 135, 1310, 735
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    title_font = _figure_font(38, bold=True)
    label_font = _figure_font(30)
    tick_font = _figure_font(25)
    small_font = _figure_font(23)
    draw.text((width // 2, 42), title, fill="black", font=title_font, anchor="ma")
    rows, columns = matrix.shape
    cell_width = (right - left) / max(columns, 1)
    cell_height = (bottom - top) / max(rows, 1)
    finite = np.abs(matrix[np.isfinite(matrix)])
    limit = max(float(np.max(finite)) if finite.size else 1.0, 1e-12)
    for row in range(rows):
        for column in range(columns):
            x0 = round(left + column * cell_width)
            x1 = round(left + (column + 1) * cell_width)
            y0 = round(bottom - (row + 1) * cell_height)
            y1 = round(bottom - row * cell_height)
            value = float(matrix[row, column])
            colour = _diverging_colour(value, limit) if diverging else _sequential_colour(value)
            draw.rectangle((x0, y0, x1, y1), fill=colour, outline=(245, 245, 245), width=2)
    draw.rectangle((left, top, right, bottom), outline="black", width=3)
    for column, label in enumerate(xlabels):
        x = left + (column + 0.5) * cell_width
        draw.text((x, bottom + 18), label, fill="black", font=tick_font, anchor="ma")
    for row, label in enumerate(ylabels):
        y = bottom - (row + 0.5) * cell_height
        draw.text((left - 22, y), label, fill="black", font=tick_font, anchor="rm")
    draw.text(((left + right) / 2, height - 68), xlabel, fill="black", font=label_font, anchor="mm")
    y_label_layer = Image.new("RGBA", (520, 60), (255, 255, 255, 0))
    y_draw = ImageDraw.Draw(y_label_layer)
    y_draw.text((260, 30), ylabel, fill="black", font=label_font, anchor="mm")
    y_label_layer = y_label_layer.rotate(90, expand=True)
    image.paste(y_label_layer, (42, round((height - y_label_layer.height) / 2)), y_label_layer)
    bar_left, bar_right = 1370, 1420
    for y in range(top, bottom):
        fraction = 1.0 - (y - top) / max(bottom - top - 1, 1)
        if diverging:
            value = (2.0 * fraction - 1.0) * limit
            colour = _diverging_colour(value, limit)
        else:
            colour = _sequential_colour(fraction)
        draw.line((bar_left, y, bar_right, y), fill=colour)
    draw.rectangle((bar_left, top, bar_right, bottom), outline="black", width=2)
    if diverging:
        labels = [(top, f"{limit:.2g}"), ((top + bottom) / 2, "0"), (bottom, f"{-limit:.2g}")]
    else:
        labels = [(top, "1.0"), ((top + bottom) / 2, "0.5"), (bottom, "0.0")]
    for y, label in labels:
        draw.text((bar_right + 14, y), label, fill="black", font=small_font, anchor="lm")
    bar_layer = Image.new("RGBA", (600, 55), (255, 255, 255, 0))
    bar_draw = ImageDraw.Draw(bar_layer)
    bar_draw.text((300, 27), colourbar_label, fill="black", font=small_font, anchor="mm")
    bar_layer = bar_layer.rotate(90, expand=True)
    image.paste(bar_layer, (1510, round((height - bar_layer.height) / 2)), bar_layer)
    image.save(path, dpi=(180, 180))


def _save_line_plot(
    path: Path,
    series: dict[str, list[tuple[float, float]]],
    title: str = "Overtone-amplitude dependence at 5% initial drift",
    xlabel: str = "analysis start time t0 / Mf",
    ylabel: str = "normalized false EFT bias Bα",
    thresholds: tuple[float, ...] = (-1.645, 0.0, 1.645),
) -> None:
    width, height = 1500, 900
    left, top, right, bottom = 180, 130, 1400, 740
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    title_font = _figure_font(38, bold=True)
    label_font = _figure_font(29)
    tick_font = _figure_font(24)
    draw.text((width // 2, 42), title, fill="black", font=title_font, anchor="ma")
    all_points = [point for points in series.values() for point in points]
    x_values = [point[0] for point in all_points]
    y_values = [point[1] for point in all_points] + list(thresholds)
    xmin, xmax = min(x_values), max(x_values)
    ymin, ymax = min(y_values), max(y_values)
    padding = max(0.2, 0.08 * (ymax - ymin))
    ymin -= padding
    ymax += padding
    to_x = lambda value: left + (value - xmin) * (right - left) / max(xmax - xmin, 1e-15)
    to_y = lambda value: bottom - (value - ymin) * (bottom - top) / max(ymax - ymin, 1e-15)
    for value in np.linspace(ymin, ymax, 7):
        y = to_y(float(value))
        draw.line((left, y, right, y), fill=(225, 225, 225), width=1)
        draw.text((left - 14, y), f"{value:.1f}", fill="black", font=tick_font, anchor="rm")
    for value in sorted(set(x_values)):
        x = to_x(value)
        draw.line((x, top, x, bottom), fill=(238, 238, 238), width=1)
        draw.text((x, bottom + 17), f"{value:g}", fill="black", font=tick_font, anchor="ma")
    for threshold in thresholds:
        colour = (50, 50, 50) if threshold == 0.0 else (115, 115, 115)
        draw.line((left, to_y(threshold), right, to_y(threshold)), fill=colour, width=3)
    colours = [(35, 95, 160), (196, 60, 45), (44, 145, 91), (135, 77, 157)]
    for index, (label, points) in enumerate(sorted(series.items())):
        colour = colours[index % len(colours)]
        points = sorted(points)
        pixels = [(to_x(x), to_y(y)) for x, y in points]
        if len(pixels) > 1:
            draw.line(pixels, fill=colour, width=5, joint="curve")
        for x, y in pixels:
            draw.ellipse((x - 7, y - 7, x + 7, y + 7), fill=colour, outline="white", width=2)
        legend_y = top + 35 + 42 * index
        draw.line((right - 355, legend_y, right - 285, legend_y), fill=colour, width=5)
        draw.text((right - 270, legend_y), label.replace("_", " "), fill="black", font=tick_font, anchor="lm")
    draw.rectangle((left, top, right, bottom), outline="black", width=3)
    draw.text(((left + right) / 2, height - 70), xlabel, fill="black", font=label_font, anchor="mm")
    y_layer = Image.new("RGBA", (650, 60), (255, 255, 255, 0))
    y_draw = ImageDraw.Draw(y_layer)
    y_draw.text((325, 30), ylabel, fill="black", font=label_font, anchor="mm")
    y_layer = y_layer.rotate(90, expand=True)
    image.paste(y_layer, (32, round((height - y_layer.height) / 2)), y_layer)
    image.save(path, dpi=(180, 180))


def make_figures(output: Path, main_rows: list[dict], noise_summary: list[dict]) -> None:
    filtered = [
        row
        for row in main_rows
        if row["drift_family"] in {"stationary", "mass_up_spin_down"}
        and row["relaxation_pair"] == "reference"
        and row["mode_set"] == "220_221"
        and row["amplitude_scenario"] == "moderate_221"
        and row["identifiable"] == 1
    ]
    strengths = sorted({float(row["drift_strength_fractional"]) for row in filtered})
    starts = sorted({float(row["start_time_M"]) for row in filtered})
    matrix = np.full((len(strengths), len(starts)), np.nan)
    for row in filtered:
        i = strengths.index(float(row["drift_strength_fractional"]))
        j = starts.index(float(row["start_time_M"]))
        matrix[i, j] = float(row["bias_sigma"])
    _save_heatmap(
        output / "evolving_kerr_bias_heatmap.png",
        matrix,
        [f"{value:g}" for value in starts],
        [f"{100 * value:g}" for value in strengths],
        "False EFT bias: mass up, spin down; 220+221",
        "analysis start time t0 / Mf",
        "initial fractional drift [%]",
        "normalized bias Bα",
        True,
    )

    compare = [
        row
        for row in main_rows
        if row["drift_family"] == "mass_up_spin_down"
        and abs(float(row["drift_strength_fractional"]) - 0.05) < 1e-12
        and row["relaxation_pair"] == "reference"
        and row["mode_set"] == "220_221"
        and row["identifiable"] == 1
    ]
    series: dict[str, list[tuple[float, float]]] = {}
    for amplitude in sorted({row["amplitude_scenario"] for row in compare}):
        selected = sorted(
            [row for row in compare if row["amplitude_scenario"] == amplitude],
            key=lambda row: float(row["start_time_M"]),
        )
        series[amplitude] = [
            (float(row["start_time_M"]), float(row["bias_sigma"])) for row in selected
        ]
    _save_line_plot(output / "evolving_kerr_start_time_amplitude_comparison.png", series)

    strengths_noise = sorted(
        {float(row["drift_strength_fractional"]) for row in noise_summary}
    )
    starts_noise = sorted({float(row["start_time_M"]) for row in noise_summary})
    noise_matrix = np.full((len(strengths_noise), len(starts_noise)), np.nan)
    for row in noise_summary:
        i = strengths_noise.index(float(row["drift_strength_fractional"]))
        j = starts_noise.index(float(row["start_time_M"]))
        noise_matrix[i, j] = float(row["false_positive_abs_z_gt_1p645"])
    _save_heatmap(
        output / "evolving_kerr_noise_false_positive_heatmap.png",
        noise_matrix,
        [f"{value:g}" for value in starts_noise],
        [f"{100 * value:g}" for value in strengths_noise],
        "Empirical false-positive rate: |z_alpha| > 1.645",
        "analysis start time t0 / Mf",
        "initial fractional drift [%]",
        "fraction of Gaussian-noise realizations",
        False,
    )


def make_prior_figures(
    output: Path, main_rows: list[dict], prior_rows: list[dict], snr_rows: list[dict]
) -> None:
    combined = main_rows + prior_rows
    styles = (
        ("free_remnant", "220_221", "free remnant; 220+221"),
        ("nrSur7dq4_imr_prior", "220_only", "IMR prior; 220 only"),
        ("nrSur7dq4_imr_prior", "220_221", "IMR prior; 220+221"),
        ("loose_3x_imr_prior", "220_only", "3x IMR prior; 220 only"),
    )
    start_series: dict[str, list[tuple[float, float]]] = {}
    for prior_name, mode_set, label in styles:
        selected = sorted(
            (
                row
                for row in combined
                if row["prior_scenario"] == prior_name
                and row["mode_set"] == mode_set
                and row["drift_family"] == "mass_down_spin_up"
                and abs(float(row["drift_strength_fractional"]) - 0.05) < 1e-12
                and row["relaxation_pair"] == "slow_both"
                and row["amplitude_scenario"] == "moderate_221"
                and row["identifiable"] == 1
            ),
            key=lambda row: float(row["start_time_M"]),
        )
        start_series[label] = [
            (float(row["start_time_M"]), float(row["bias_sigma"])) for row in selected
        ]
    _save_line_plot(
        output / "evolving_kerr_prior_mode_start_time_comparison.png",
        start_series,
        title="Remnant-prior and mode-set dependence at 5% slow drift",
        xlabel="analysis start time t0 / Mf",
        ylabel="normalized false EFT bias Bα",
    )

    snr_styles = (
        ("free_remnant", "220_221", "free remnant; 220+221"),
        ("nrSur7dq4_imr_prior", "220_only", "IMR prior; 220 only"),
        ("nrSur7dq4_imr_prior", "220_221", "IMR prior; 220+221"),
        ("loose_3x_imr_prior", "220_only", "3x IMR prior; 220 only"),
    )
    snr_series: dict[str, list[tuple[float, float]]] = {}
    for prior_name, mode_set, label in snr_styles:
        selected = sorted(
            (
                row
                for row in snr_rows
                if row["prior_scenario"] == prior_name
                and row["mode_set"] == mode_set
                and row["drift_family"] == "mass_down_spin_up"
                and abs(float(row["drift_strength_fractional"]) - 0.05) < 1e-12
                and row["relaxation_pair"] == "slow_both"
                and float(row["start_time_M"]) == 0.0
                and row["identifiable"] == 1
            ),
            key=lambda row: float(row["peak_snr"]),
        )
        snr_series[label] = [
            (float(row["peak_snr"]), float(row["bias_sigma"])) for row in selected
        ]
    _save_line_plot(
        output / "evolving_kerr_prior_mode_snr_comparison.png",
        snr_series,
        title="SNR dependence of the false EFT bias at early start",
        xlabel="peak-time synthetic SNR",
        ylabel="normalized false EFT bias Bα",
        thresholds=(-3.0, -1.645, 0.0, 1.645, 3.0),
    )


def make_report(
    output: Path,
    config: dict,
    validation: list[dict],
    main_rows: list[dict],
    operator_rows: list[dict],
    profile_rows: list[dict],
    noise_summary: list[dict],
    convergence: list[dict],
) -> None:
    identifiable_main = [row for row in main_rows if row["identifiable"] == 1]
    dynamic_main = [
        row
        for row in identifiable_main
        if float(row["drift_strength_fractional"]) > 0.0
    ]
    maximum = max(dynamic_main, key=lambda row: abs(float(row["bias_sigma"])))
    stationary = [
        row
        for row in identifiable_main
        if float(row["drift_strength_fractional"]) == 0.0
    ]
    stationary_max = max(abs(float(row["bias_sigma"])) for row in stationary)
    convergence_bias = np.asarray([float(row["bias_sigma"]) for row in convergence])
    convergence_spread = float(np.max(convergence_bias) - np.min(convergence_bias))
    profile_differences = [
        abs(float(row["nonlinear_alpha_hat"]) - float(row["linear_alpha_hat"]))
        for row in profile_rows
        if int(row["nonlinear_success"]) == 1
    ]
    max_profile_difference = max(profile_differences) if profile_differences else math.nan
    max_false_positive = max(
        noise_summary, key=lambda row: float(row["false_positive_abs_z_gt_1p645"])
    )
    validation_pass = all(int(row["pass"]) == 1 for row in validation)
    lines = [
        "# Evolving-Kerr-spectrum false-EFT-bias study",
        "",
        "## Scope",
        "",
        "This calculation is a controlled adiabatic evolving-spectrum surrogate. It is not an exact time-dependent Kerr spacetime and does not use detector strain.",
        "",
        f"- Event anchor: `{config['event']['name']}`",
        f"- Final detector-frame mass: `{config['event']['mass_detector_msun']}` solar masses",
        f"- Final spin: `{config['event']['spin']}`",
        f"- Peak-time synthetic SNR: `{config['peak_snr']}`",
        f"- Representative fingerprint: `{config['representative_operator']}/{config['representative_branch']}`",
        f"- Main scan rows: `{len(main_rows)}`",
        f"- Operator-robustness rows: `{len(operator_rows)}`",
        f"- Nonlinear profile-validation rows: `{len(profile_rows)}`",
        "",
        "## Validation",
        "",
        f"- All mandatory validation checks pass: `{validation_pass}`",
        f"- Maximum stationary false bias: `{stationary_max:.3e}` sigma",
        f"- Convergence-grid bias spread: `{convergence_spread:.4f}` sigma",
        f"- Maximum absolute linear/nonlinear alpha difference in selected profile checks: `{max_profile_difference:.5g}`",
        "",
        "## Main diagnostic result",
        "",
        "The largest absolute bias in the configured scan is:",
        "",
        f"- family: `{maximum['drift_family']}`",
        f"- initial drift strength: `{100 * float(maximum['drift_strength_fractional']):.2f}%`",
        f"- relaxation pair: `{maximum['relaxation_pair']}`",
        f"- start time: `{maximum['start_time_M']} M_f`",
        f"- mode set: `{maximum['mode_set']}`",
        f"- amplitude scenario: `{maximum['amplitude_scenario']}`",
        f"- alpha_hat: `{float(maximum['alpha_hat']):.6g}`",
        f"- sigma_alpha: `{float(maximum['sigma_alpha']):.6g}`",
        f"- false bias: `{float(maximum['bias_sigma']):.3f} sigma`",
        "",
        "The largest empirical 90%-threshold false-positive fraction in the configured Gaussian-noise ensemble is:",
        "",
        f"- drift strength: `{100 * float(max_false_positive['drift_strength_fractional']):.2f}%`",
        f"- start time: `{max_false_positive['start_time_M']} M_f`",
        f"- fraction with `|z_alpha| > 1.645`: `{float(max_false_positive['false_positive_abs_z_gt_1p645']):.3f}`",
        f"- Wilson 95% interval for that fraction: `[{float(max_false_positive['false_positive_abs_z_gt_1p645_wilson95_low']):.3f}, {float(max_false_positive['false_positive_abs_z_gt_1p645_wilson95_high']):.3f}]`",
        "",
        "## Interpretation guardrails",
        "",
        "- `220_only` is retained as an expected identifiability negative control after mass-spin profiling.",
        "- A nonzero recovered alpha is a misspecification diagnostic, not evidence for modified gravity.",
        "- Drift strengths define a controlled envelope; they are not a posterior for the physical remnant evolution.",
        "- Public GW250114 products remain a calibration anchor and are not statistically combined with these synthetic results.",
        "",
        "## Generated outputs",
        "",
        "- `evolving_kerr_validation.csv`",
        "- `evolving_kerr_convergence.csv`",
        "- `evolving_kerr_main_bias_scan.csv`",
        "- `evolving_kerr_prior_informed_bias_scan.csv`",
        "- `evolving_kerr_snr_scaling.csv`",
        "- `evolving_kerr_prior_snr_report.md`",
        "- `evolving_kerr_operator_robustness.csv`",
        "- `evolving_kerr_profile_validation.csv`",
        "- `evolving_kerr_noise_ensemble.csv`",
        "- `evolving_kerr_noise_summary.csv`",
        "- five publication-facing PNG figures",
    ]
    (output / "evolving_kerr_bias_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config",
        type=Path,
        default=ROOT / "config" / "evolving_kerr_bias.json",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "results" / "evolving_kerr_bias",
    )
    args = parser.parse_args()
    config = load_config(args.config)
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    anchors = load_anchors(
        ROOT / "results" / "gw250114_kerr_qnm" / "qnm_solver_crosscheck.csv",
        float(config["event"]["mass_detector_msun"]),
        float(config["event"]["spin"]),
    )
    coefficients = load_eft_coefficients(ROOT / "data" / "beyond_kerr_qnm_selected_fits.csv")
    priors = remnant_prior_covariances(config)

    validation = validation_rows(config, anchors, coefficients)
    write_csv(output / "evolving_kerr_validation.csv", validation)
    if not all(int(row["pass"]) == 1 for row in validation):
        raise RuntimeError("Mandatory evolving-Kerr validation failed")
    print("Validation: PASS")

    convergence = run_convergence(config, anchors, coefficients)
    write_csv(output / "evolving_kerr_convergence.csv", convergence)
    print(f"Convergence rows: {len(convergence)}")

    main_rows = run_main_scan(config, anchors, coefficients)
    write_csv(output / "evolving_kerr_main_bias_scan.csv", main_rows)
    print(f"Main scan rows: {len(main_rows)}")

    prior_rows: list[dict] = []
    for prior_name, prior_covariance in priors.items():
        if prior_covariance is None:
            continue
        prior_rows.extend(
            run_main_scan(
                config,
                anchors,
                coefficients,
                prior_scenario=prior_name,
                prior_covariance=prior_covariance,
            )
        )
    write_csv(output / "evolving_kerr_prior_informed_bias_scan.csv", prior_rows)
    print(f"Prior-informed scan rows: {len(prior_rows)}")

    snr_rows = run_snr_scaling(config, anchors, coefficients, priors)
    write_csv(output / "evolving_kerr_snr_scaling.csv", snr_rows)
    print(f"SNR scaling rows: {len(snr_rows)}")

    operator_rows = run_operator_robustness(config, anchors, coefficients)
    write_csv(output / "evolving_kerr_operator_robustness.csv", operator_rows)
    print(f"Operator robustness rows: {len(operator_rows)}")

    profile_rows = run_profile_validation(config, anchors, coefficients)
    write_csv(output / "evolving_kerr_profile_validation.csv", profile_rows)
    print(f"Profile validation rows: {len(profile_rows)}")

    noise_rows, noise_summary = run_noise_ensemble(config, anchors, coefficients)
    write_csv(output / "evolving_kerr_noise_ensemble.csv", noise_rows)
    write_csv(output / "evolving_kerr_noise_summary.csv", noise_summary)
    print(f"Noise rows: {len(noise_rows)}")

    make_figures(output, main_rows, noise_summary)
    make_prior_figures(output, main_rows, prior_rows, snr_rows)
    make_prior_report(output, config, prior_rows, snr_rows)
    make_report(
        output,
        config,
        validation,
        main_rows,
        operator_rows,
        profile_rows,
        noise_summary,
        convergence,
    )
    print(f"Output: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
