"""Time-domain validation for static spherical master potentials.

The script evolves

    d^2 psi/dr_*^2 - d^2 psi/dt^2 = V(r) psi

in null coordinates using the standard second-order characteristic stencil.
It is intentionally dependency-light: only NumPy is required, with the local
`qnm` package used for the Schwarzschild reference value when available.

Current scope: Schwarzschild scalar, electromagnetic, Regge-Wheeler, and
Zerilli potentials. This is a validation layer for later static-metric
stress tests, not a precision replacement for Leaver or high-order WKB/Pade.
"""

from __future__ import annotations

import csv
import math
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / ".deps" / "python"))


try:
    import qnm  # type: ignore
except Exception:  # pragma: no cover - only used for a nicer local report
    qnm = None


@dataclass(frozen=True)
class PotentialCase:
    name: str
    sector: str
    ell: int


CASES = [
    PotentialCase("Schwarzschild_scalar_l2", "scalar", 2),
    PotentialCase("Schwarzschild_EM_l2", "electromagnetic", 2),
    PotentialCase("Schwarzschild_ReggeWheeler_l2", "gravitational_odd", 2),
    PotentialCase("Schwarzschild_Zerilli_l2", "gravitational_even", 2),
]


def schwarzschild_radius_from_tortoise(x: float) -> float:
    """Invert r_* = r + 2 log(r/2 - 1) for M=1 Schwarzschild."""

    if x < -50.0:
        return 2.0 + 2.0 * math.exp((x - 2.0) / 2.0)

    low = 2.0 + 1e-14
    high = max(4.0, x + 6.0)

    def tortoise(radius: float) -> float:
        return radius + 2.0 * math.log(radius / 2.0 - 1.0)

    while tortoise(high) < x:
        high *= 2.0

    for _ in range(90):
        mid = 0.5 * (low + high)
        if tortoise(mid) < x:
            low = mid
        else:
            high = mid
    return 0.5 * (low + high)


def schwarzschild_potential(case: PotentialCase, x_values: np.ndarray) -> np.ndarray:
    radii = np.array([schwarzschild_radius_from_tortoise(float(x)) for x in x_values])
    ell = case.ell
    f = 1.0 - 2.0 / radii

    if case.sector == "scalar":
        return f * (ell * (ell + 1) / radii**2 + 2.0 / radii**3)
    if case.sector == "electromagnetic":
        return f * ell * (ell + 1) / radii**2
    if case.sector == "gravitational_odd":
        return f * (ell * (ell + 1) / radii**2 - 6.0 / radii**3)
    if case.sector == "gravitational_even":
        lam = 0.5 * (ell - 1) * (ell + 2)
        return (
            2.0
            * f
            * (
                lam**2 * (lam + 1.0) * radii**3
                + 3.0 * lam**2 * radii**2
                + 9.0 * lam * radii
                + 9.0
            )
            / (radii**3 * (lam * radii + 3.0) ** 2)
        )
    raise ValueError(f"unknown sector: {case.sector}")


def characteristic_evolve_with_potential(
    potential_at_tortoise,
    *,
    h: float = 0.2,
    n_grid: int = 2200,
    observer_rstar: float = 20.0,
    gaussian_center: float = 0.0,
    gaussian_width: float = 3.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Return time and waveform at fixed observer tortoise radius."""

    k_obs = int(round(2.0 * observer_rstar / h))
    if abs(k_obs * h / 2.0 - observer_rstar) > 1e-10:
        raise ValueError("observer_rstar must be compatible with h/2 grid")
    if n_grid <= k_obs + 20:
        raise ValueError("n_grid too small for requested observer")

    diffs = np.arange(-n_grid, n_grid + 1, dtype=float)
    x_for_diff = 0.5 * h * diffs
    v_by_diff = potential_at_tortoise(x_for_diff)

    psi = np.zeros((n_grid + 1, n_grid + 1), dtype=np.float64)
    v_axis = np.arange(n_grid + 1, dtype=float) * h
    initial_rstar = 0.5 * v_axis
    psi[0, :] = np.exp(-((initial_rstar - gaussian_center) ** 2) / (2.0 * gaussian_width**2))
    psi[:, 0] = 0.0
    psi[0, 0] = 0.0

    coeff = (h * h) / 8.0
    offset = n_grid
    for i in range(n_grid):
        row_next = psi[i + 1]
        row_current = psi[i]
        for j in range(n_grid):
            v_here = v_by_diff[j - i + offset]
            row_next[j + 1] = (
                row_next[j]
                + row_current[j + 1]
                - row_current[j]
                - coeff * v_here * (row_next[j] + row_current[j + 1])
            )

    i_values = np.arange(0, n_grid - k_obs + 1)
    j_values = i_values + k_obs
    times = 0.5 * h * (i_values + j_values)
    waveform = psi[i_values, j_values]
    return times, waveform


def characteristic_evolve(
    case: PotentialCase,
    **kwargs,
) -> tuple[np.ndarray, np.ndarray]:
    """Return time and waveform for one of the built-in Schwarzschild cases."""

    return characteristic_evolve_with_potential(
        lambda x_values: schwarzschild_potential(case, x_values),
        **kwargs,
    )


def prony_two_mode(times: np.ndarray, waveform: np.ndarray) -> tuple[float, float, float]:
    """Estimate omega, gamma, and normalized residual from a real damped sinusoid."""

    if len(times) < 50:
        return math.nan, math.nan, math.inf

    y = waveform.astype(float).copy()
    y -= float(np.mean(y))
    scale = float(np.max(np.abs(y)))
    if scale <= 0.0 or not math.isfinite(scale):
        return math.nan, math.nan, math.inf
    y /= scale

    lhs = -y[2:]
    rhs = np.column_stack([y[1:-1], y[:-2]])
    coeffs, *_ = np.linalg.lstsq(rhs, lhs, rcond=None)
    roots = np.roots([1.0, coeffs[0], coeffs[1]])

    dt = float(np.median(np.diff(times)))
    candidates = []
    for root in roots:
        if abs(root) <= 0.0:
            continue
        exponent = np.log(root) / dt
        omega = abs(float(np.imag(exponent)))
        gamma = -float(np.real(exponent))
        if omega > 0.0 and gamma > 0.0 and math.isfinite(omega) and math.isfinite(gamma):
            candidates.append((omega, gamma))

    if not candidates:
        return math.nan, math.nan, math.inf

    omega, gamma = sorted(candidates, key=lambda item: item[0], reverse=True)[0]
    design = np.column_stack(
        [
            np.exp(-gamma * (times - times[0])) * np.cos(omega * (times - times[0])),
            np.exp(-gamma * (times - times[0])) * np.sin(omega * (times - times[0])),
            np.ones_like(times),
        ]
    )
    amps, *_ = np.linalg.lstsq(design, waveform, rcond=None)
    residual = waveform - design @ amps
    denom = np.linalg.norm(waveform - np.mean(waveform))
    rel_residual = float(np.linalg.norm(residual) / denom) if denom > 0.0 else math.inf
    return omega, gamma, rel_residual


def reference_frequency(case: PotentialCase) -> tuple[float, float, str]:
    if qnm is None:
        return math.nan, math.nan, "unavailable"

    spin_weight_by_sector = {
        "scalar": 0,
        "electromagnetic": -1,
        "gravitational_odd": -2,
        "gravitational_even": -2,
    }
    spin_weight = spin_weight_by_sector[case.sector]
    seq = qnm.modes_cache(s=spin_weight, l=case.ell, m=case.ell, n=0)
    omega, _angular_constant, _mixing = seq(a=0.0)
    return float(omega.real), abs(float(omega.imag)), "qnm"


def fit_windows(
    times: np.ndarray,
    waveform: np.ndarray,
    *,
    t_start_min: float = 40.0,
    t_start_max: float = 160.0,
    widths: tuple[float, ...] = (60.0, 70.0, 80.0, 100.0),
    step: float = 10.0,
) -> list[dict[str, float]]:
    rows = []
    for width in widths:
        start = t_start_min
        while start <= t_start_max + 1e-9:
            stop = start + width
            mask = (times >= start) & (times <= stop)
            omega, gamma, rel_residual = prony_two_mode(times[mask], waveform[mask])
            rows.append(
                {
                    "fit_t_start": start,
                    "fit_t_stop": stop,
                    "fit_width": width,
                    "td_Momega_real": omega,
                    "td_Momega_imag_abs": gamma,
                    "relative_residual": rel_residual,
                }
            )
            start += step
    return rows


def pct_delta(value: float, reference: float) -> float:
    if not math.isfinite(value) or not math.isfinite(reference) or reference == 0.0:
        return math.nan
    return 100.0 * (value - reference) / reference


def format_float(value: float, digits: int = 12) -> str:
    if not math.isfinite(value):
        return "nan"
    return f"{value:.{digits}g}"


def main() -> None:
    out_dir = ROOT / "results" / "static_master_potential_time_domain"
    out_dir.mkdir(parents=True, exist_ok=True)

    summary_rows: list[dict[str, object]] = []
    window_rows: list[dict[str, object]] = []

    for case in CASES:
        times, waveform = characteristic_evolve(case)
        np.savetxt(
            out_dir / f"{case.name}_waveform.csv",
            np.column_stack([times, waveform]),
            delimiter=",",
            header="t_over_M,psi",
            comments="",
        )

        ref_real, ref_imag, ref_source = reference_frequency(case)
        windows = fit_windows(times, waveform)
        for row in windows:
            window_rows.append({"name": case.name, "sector": case.sector, **row})

        finite_windows = [
            row
            for row in windows
            if math.isfinite(row["td_Momega_real"])
            and math.isfinite(row["td_Momega_imag_abs"])
            and math.isfinite(row["relative_residual"])
        ]
        if finite_windows:
            best = min(finite_windows, key=lambda row: row["relative_residual"])
        else:
            best = {
                "fit_t_start": math.nan,
                "fit_t_stop": math.nan,
                "fit_width": math.nan,
                "td_Momega_real": math.nan,
                "td_Momega_imag_abs": math.nan,
                "relative_residual": math.inf,
            }

        summary_rows.append(
            {
                "name": case.name,
                "sector": case.sector,
                "ell": case.ell,
                "reference_source": ref_source,
                "reference_Momega_real": ref_real,
                "reference_Momega_imag_abs": ref_imag,
                "fit_t_start": best["fit_t_start"],
                "fit_t_stop": best["fit_t_stop"],
                "fit_width": best["fit_width"],
                "td_Momega_real": best["td_Momega_real"],
                "td_Momega_imag_abs": best["td_Momega_imag_abs"],
                "delta_real_pct": pct_delta(float(best["td_Momega_real"]), ref_real),
                "delta_imag_pct": pct_delta(float(best["td_Momega_imag_abs"]), ref_imag),
                "relative_residual": best["relative_residual"],
            }
        )

    summary_csv = out_dir / "time_domain_qnm_summary.csv"
    with summary_csv.open("w", newline="", encoding="utf-8") as handle:
        fields = list(summary_rows[0].keys())
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(summary_rows)

    windows_csv = out_dir / "time_domain_qnm_fit_windows.csv"
    with windows_csv.open("w", newline="", encoding="utf-8") as handle:
        fields = list(window_rows[0].keys())
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(window_rows)

    lines = [
        "# Static Master-Potential Time-Domain Validation",
        "",
        "This report evolves Schwarzschild supplied potentials with a characteristic null-grid scheme.",
        "The extracted frequency uses a two-root Prony/linear-prediction fit over several fixed ringdown windows.",
        "",
        "Default numerical setup: `h = 0.2 M`, `r*_obs = 20 M`, Gaussian center `r* = 0 M`, width `3 M`.",
        "",
        "Scope note: this is an independent validation layer for static supplied potentials. It is not a replacement for Leaver, high-order WKB/Pade, or a full strain likelihood.",
        "",
        "| potential | sector | reference Re | TD Re | delta Re [%] | reference Abs(Im) | TD Abs(Im) | delta Abs(Im) [%] | window | residual |",
        "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in summary_rows:
        lines.append(
            "| {name} | {sector} | {reference_Momega_real:.12f} | {td_Momega_real:.12f} | {delta_real_pct:.3f} | {reference_Momega_imag_abs:.12f} | {td_Momega_imag_abs:.12f} | {delta_imag_pct:.3f} | {fit_t_start:.0f}-{fit_t_stop:.0f} | {relative_residual:.4f} |".format(
                **row
            )
        )
    lines.extend(
        [
            "",
            "Interpretation:",
            "",
            "- The Regge-Wheeler and Zerilli rows should agree with the same Schwarzschild gravitational reference.",
            "- The time-domain extraction is expected to be less precise than Leaver/qnm, but it is independent of WKB.",
            "- A case is suitable for the next static-metric stress-test layer only after it passes this benchmark and has a published master potential.",
            "",
        ]
    )
    report_path = out_dir / "time_domain_qnm_validation.md"
    report_path.write_text("\n".join(lines), encoding="utf-8")

    print(f"Summary: {summary_csv}")
    print(f"Windows: {windows_csv}")
    print(f"Report: {report_path}")
    for row in summary_rows:
        print(
            "{name}\tTD={td_Momega_real:.8f}-{td_Momega_imag_abs:.8f}i\t"
            "ref={reference_Momega_real:.8f}-{reference_Momega_imag_abs:.8f}i\t"
            "delta=({delta_real_pct:.2f}%, {delta_imag_pct:.2f}%)".format(**row)
        )


if __name__ == "__main__":
    main()
