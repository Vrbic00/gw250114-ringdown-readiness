"""Continued-fraction QNMs for the rotating hairy metric of Zhen Li.

This implements the approximate Dudley--Finley radial problem in
arXiv:2212.08112.  It is deliberately kept separate from the observational
pipeline until the published tables and the Kerr limit are reproduced.

The source leaves the deformed horizon indices ambiguous.  The default
``table_calibrated`` branch retains the Kerr horizon identity and reproduces
the published tables; ``ode_consistent`` evaluates the index from the stated
radial ODE.  Both give the standard Teukolsky/Kerr indices when
``alpha_hair = 0``.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from dataclasses import asdict, dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
LOCAL_DEPS = ROOT / ".deps" / "python"
if LOCAL_DEPS.exists():
    sys.path.insert(0, str(LOCAL_DEPS))

import numpy as np
import qnm
from qnm import angular
from scipy.optimize import root


PUBLISHED_TABLE = ROOT / "data" / "ovalle_rotating_hairy_qnm_zhen_li_2022.csv"
DEFAULT_OUTPUT = ROOT / "results" / "hairy_continued_fraction"


@dataclass(frozen=True)
class SolveResult:
    spin: float
    alpha_hair: float
    h0: float
    overtone: int
    omega_real: float
    minus_omega_imag: float
    angular_real: float
    angular_imag: float
    residual_norm: float
    root_success: bool
    iterations: int
    truncation: int
    inversion: int
    horizon_index_convention: str
    horizon_source: str

    @property
    def omega(self) -> complex:
        return complex(self.omega_real, -self.minus_omega_imag)


def _hair_exponential(radius: float, h0: float) -> float:
    denominator = 1.0 - 0.5 * h0
    if denominator <= 0.0:
        return 0.0
    return math.exp(-radius / denominator)


def approximate_horizons(spin: float, alpha_hair: float, h0: float, iterations: int = 2) -> tuple[float, float]:
    """Iterated horizons used in arXiv:2212.08112, with M=1."""
    if not 0.0 <= spin < 1.0:
        raise ValueError("spin must satisfy 0 <= a < 1")
    if alpha_hair < 0.0:
        raise ValueError("alpha_hair must be non-negative")
    if not 0.0 < h0 <= 2.0:
        raise ValueError("h0 must satisfy 0 < h0 <= 2 in units M=1")

    root_kerr = math.sqrt(1.0 - spin * spin)
    plus = 1.0 + root_kerr
    minus = 1.0 - root_kerr
    for _ in range(iterations):
        radicand_plus = 1.0 - spin * spin - alpha_hair * plus * plus * _hair_exponential(plus, h0)
        radicand_minus = 1.0 - spin * spin - alpha_hair * minus * minus * _hair_exponential(minus, h0)
        if radicand_plus <= 0.0 or radicand_minus <= 0.0:
            raise ValueError("hair parameters do not yield two real approximate horizons")
        plus = 1.0 + math.sqrt(radicand_plus)
        minus = 1.0 - math.sqrt(radicand_minus)
    if plus <= minus:
        raise ValueError("invalid approximate horizon ordering")
    return plus, minus


def separation_constant(omega: complex, spin: float, *, s: int = -2, ell: int = 2, m: int = 2, l_max: int = 24) -> complex:
    """Spin-weighted spheroidal eigenvalue in the convention of the paper."""
    spherical = complex(angular.swsphericalh_A(s, ell, m))
    return complex(angular.sep_const_closest(spherical, s, spin * omega, m, l_max))


def recurrence_constants(
    omega: complex,
    spin: float,
    alpha_hair: float,
    h0: float,
    *,
    s: int = -2,
    ell: int = 2,
    m: int = 2,
    horizon_index_convention: str = "table_calibrated",
    horizons: tuple[float, float] | None = None,
) -> tuple[complex, complex, complex, complex, complex, complex]:
    """Return the five recurrence constants and angular eigenvalue.

    These quantities depend on the trial frequency but not on the recurrence
    index.  Computing them once per continued-fraction evaluation avoids
    hundreds of redundant angular eigenvalue solves.
    """
    r_plus, r_minus = approximate_horizons(spin, alpha_hair, h0) if horizons is None else horizons
    b = r_plus - r_minus
    separation = separation_constant(omega, spin, s=s, ell=ell, m=m)

    # The preprint's compact index formulas are ambiguous.  Its tables are
    # reproduced only if the Kerr identity r_h^2+a^2=2 r_h is retained after
    # deforming the horizons (``table_calibrated``).  Applying K(r_h) from
    # the stated radial ODE gives the distinct ``ode_consistent`` convention.
    # Both reduce exactly to the standard Teukolsky recurrence for alpha=0.
    if horizon_index_convention == "table_calibrated":
        k_plus = 2.0 * r_plus * omega - spin * m
        k_minus = 2.0 * r_minus * omega - spin * m
    elif horizon_index_convention == "ode_consistent":
        k_plus = (r_plus * r_plus + spin * spin) * omega - spin * m
        k_minus = (r_minus * r_minus + spin * spin) * omega - spin * m
    else:
        raise ValueError(f"unknown horizon-index convention: {horizon_index_convention}")
    sigma_plus = k_plus / b - 1j * s
    sigma_minus = k_minus / b + 1j * s

    b1 = (-s - 1.0 - 2j * sigma_minus) * b
    b2 = 2.0 * (1j * sigma_minus - 1j * sigma_plus + s + 1.0)
    b3 = (
        2.0 * omega * omega * r_plus * r_plus
        + omega * omega * (r_plus + r_minus) ** 2
        + spin * spin * omega * omega
        - separation
        + 1j * (2.0 * s + 1.0) * (sigma_minus - sigma_plus)
        - (sigma_minus - sigma_plus) ** 2
        - 1j * s * omega * (b + 4.0 - 4.0 * r_plus)
    )
    eta = -omega * (r_plus + r_minus) - 1j * s

    c0 = b2 + b1 / b
    c1 = -2.0 * (c0 + 1.0 + 1j * (eta - omega * b))
    c2 = c0 + 2.0 * (1.0 + 1j * eta)
    c4 = (0.5 * b2 + 1j * eta) * (0.5 * b2 + 1j * eta + 1.0 + b1 / b)
    c3 = (
        -c4
        - 0.5 * b2 * (0.5 * b2 - 1.0)
        + eta * (1j - eta)
        + 1j * omega * b * c0
        + b3
    )

    return c0, c1, c2, c3, c4, separation


def recurrence_coefficients(
    index: int,
    constants: tuple[complex, complex, complex, complex, complex, complex],
) -> tuple[complex, complex, complex]:
    """Return radial ``(alpha_n, beta_n, gamma_n)`` from shared constants."""
    c0, c1, c2, c3, c4, _ = constants
    n = float(index)
    alpha_n = n * n + (c0 + 1.0) * n + c0
    beta_n = -2.0 * n * n + (c1 + 2.0) * n + c3
    gamma_n = n * n + (c2 - 3.0) * n + c4 - c2 + 2.0
    return alpha_n, beta_n, gamma_n


def continued_fraction_residual(
    omega: complex,
    spin: float,
    alpha_hair: float,
    h0: float,
    *,
    overtone: int = 0,
    truncation: int = 300,
    s: int = -2,
    ell: int = 2,
    m: int = 2,
    horizon_index_convention: str = "table_calibrated",
    horizons: tuple[float, float] | None = None,
) -> complex:
    """Evaluate the overtone-inverted radial continued fraction."""
    if overtone < 0:
        raise ValueError("overtone must be non-negative")
    if truncation <= overtone + 5:
        raise ValueError("truncation is too small for requested inversion")

    constants = recurrence_constants(
        omega,
        spin,
        alpha_hair,
        h0,
        s=s,
        ell=ell,
        m=m,
        horizon_index_convention=horizon_index_convention,
        horizons=horizons,
    )
    coeffs = [recurrence_coefficients(k, constants) for k in range(truncation + 1)]

    # Right infinite fraction, truncated from below.
    right = coeffs[truncation][1]
    for k in range(truncation - 1, overtone, -1):
        alpha_k, beta_k, _ = coeffs[k]
        gamma_next = coeffs[k + 1][2]
        if abs(right) < 1.0e-300:
            right += 1.0e-300
        right = beta_k - alpha_k * gamma_next / right

    alpha_o, beta_o, _ = coeffs[overtone]
    right_term = alpha_o * coeffs[overtone + 1][2] / right
    if overtone == 0:
        return beta_o - right_term

    left = coeffs[0][1]
    for k in range(1, overtone):
        alpha_prev = coeffs[k - 1][0]
        beta_k = coeffs[k][1]
        gamma_k = coeffs[k][2]
        if abs(left) < 1.0e-300:
            left += 1.0e-300
        left = beta_k - alpha_prev * gamma_k / left
    left_term = coeffs[overtone - 1][0] * coeffs[overtone][2] / left
    return beta_o - left_term - right_term


def kerr_seed(spin: float, overtone: int, *, s: int = -2, ell: int = 2, m: int = 2) -> complex:
    sequence = qnm.modes_cache(s=s, l=ell, m=m, n=overtone)
    omega, _, _ = sequence(a=spin)
    return complex(omega)


def solve_mode(
    spin: float,
    alpha_hair: float,
    h0: float,
    overtone: int,
    *,
    initial: complex | None = None,
    truncation: int = 300,
    s: int = -2,
    ell: int = 2,
    m: int = 2,
    tolerance: float = 1.0e-11,
    horizon_index_convention: str = "table_calibrated",
    horizons: tuple[float, float] | None = None,
    horizon_source: str = "second_order_approximation",
) -> SolveResult:
    if initial is None:
        initial = kerr_seed(spin, overtone, s=s, ell=ell, m=m)

    def real_residual(values: np.ndarray) -> np.ndarray:
        trial = complex(float(values[0]), float(values[1]))
        value = continued_fraction_residual(
            trial,
            spin,
            alpha_hair,
            h0,
            overtone=overtone,
            truncation=truncation,
            s=s,
            ell=ell,
            m=m,
            horizon_index_convention=horizon_index_convention,
            horizons=horizons,
        )
        return np.asarray([value.real, value.imag], dtype=float)

    fit = root(real_residual, np.asarray([initial.real, initial.imag], dtype=float), method="hybr", tol=tolerance)
    omega = complex(float(fit.x[0]), float(fit.x[1]))
    residual = continued_fraction_residual(
        omega,
        spin,
        alpha_hair,
        h0,
        overtone=overtone,
        truncation=truncation,
        s=s,
        ell=ell,
        m=m,
        horizon_index_convention=horizon_index_convention,
        horizons=horizons,
    )
    separation = separation_constant(omega, spin, s=s, ell=ell, m=m)
    return SolveResult(
        spin=spin,
        alpha_hair=alpha_hair,
        h0=h0,
        overtone=overtone,
        omega_real=omega.real,
        minus_omega_imag=-omega.imag,
        angular_real=separation.real,
        angular_imag=separation.imag,
        residual_norm=abs(residual),
        root_success=bool(fit.success and abs(residual) < 1.0e-7 and omega.imag < 0.0),
        iterations=int(getattr(fit, "nfev", -1)),
        truncation=truncation,
        inversion=overtone,
        horizon_index_convention=horizon_index_convention,
        horizon_source=horizon_source,
    )


def read_published_rows() -> list[dict[str, str]]:
    with PUBLISHED_TABLE.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def validation_cases() -> list[dict[str, float | int | str]]:
    """Unique cases covering every independent published table branch."""
    selected: list[dict[str, float | int | str]] = []
    seen: set[tuple[float, float, float, int]] = set()
    for row in read_published_rows():
        key = (float(row["spin_a"]), float(row["alpha"]), float(row["h0"]), int(row["n"]))
        if key in seen:
            continue
        seen.add(key)
        selected.append(
            {
                "source_table": row["source_table"],
                "spin": key[0],
                "alpha_hair": key[1],
                "h0": key[2],
                "overtone": key[3],
                "published_real": float(row["re_momega"]),
                "published_minus_imag": float(row["minus_im_momega"]),
            }
        )
    return selected


def validate(output: Path, truncation: int) -> dict[str, object]:
    output.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    by_branch: dict[tuple[int, float, float], complex] = {}

    for case in validation_cases():
        spin = float(case["spin"])
        alpha_hair = float(case["alpha_hair"])
        h0 = float(case["h0"])
        overtone = int(case["overtone"])
        published = complex(float(case["published_real"]), -float(case["published_minus_imag"]))

        # The published value is a legitimate validation seed, not an output
        # constraint.  Kerr cases additionally test independent qnm seeds.
        initial = kerr_seed(spin, overtone) if alpha_hair == 0.0 else published
        result = solve_mode(
            spin,
            alpha_hair,
            h0,
            overtone,
            initial=initial,
            truncation=truncation,
        )
        by_branch[(overtone, alpha_hair, h0)] = result.omega
        delta = result.omega - published
        row = dict(case)
        row.update(asdict(result))
        row.update(
            {
                "delta_real": delta.real,
                "delta_minus_imag": -delta.imag,
                "abs_delta": abs(delta),
                "relative_delta": abs(delta) / max(abs(published), 1.0e-30),
            }
        )
        rows.append(row)

    fields = list(rows[0])
    with (output / "published_table_reproduction.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    converged = [row for row in rows if bool(row["root_success"])]
    summary: dict[str, object] = {
        "cases": len(rows),
        "converged": len(converged),
        "max_abs_delta": max(float(row["abs_delta"]) for row in rows),
        "median_abs_delta": float(np.median([float(row["abs_delta"]) for row in rows])),
        "max_relative_delta": max(float(row["relative_delta"]) for row in rows),
        "truncation": truncation,
        "status": "PASS" if len(converged) == len(rows) else "FAIL",
    }
    (output / "validation_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--truncation", type=int, default=300)
    args = parser.parse_args()
    summary = validate(args.output, args.truncation)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
