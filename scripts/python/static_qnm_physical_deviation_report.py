"""Separate numerical validation error from physical QNM deviations."""

from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = ROOT / "results" / "static_qnm_physical_deviations"


SCHWARZSCHILD_GRAV_L2_N0 = (0.37367168441804177, 0.08896231568893546)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def as_float(value: object) -> float:
    try:
        return float(value)
    except Exception:
        return math.nan


def pct_delta(value: float, reference: float) -> float:
    return 100.0 * (value / reference - 1.0)


def add_deviation_row(
    rows: list[dict[str, object]],
    *,
    family: str,
    case: str,
    parameter_name: str,
    parameter_value: float,
    ell: int,
    n: int,
    source_frequency: str,
    omega_real: float,
    omega_imag_abs: float,
    baseline_label: str,
    baseline_real: float,
    baseline_imag_abs: float,
    validation_delta_real_pct: float,
    validation_delta_imag_pct: float,
) -> None:
    delta_re = pct_delta(omega_real, baseline_real)
    delta_im = pct_delta(omega_imag_abs, baseline_imag_abs)
    rows.append(
        {
            "family": family,
            "case": case,
            "parameter_name": parameter_name,
            "parameter_value": parameter_value,
            "ell": ell,
            "n": n,
            "source_frequency": source_frequency,
            "Momega_real": omega_real,
            "Momega_imag_abs": omega_imag_abs,
            "baseline_label": baseline_label,
            "baseline_Momega_real": baseline_real,
            "baseline_Momega_imag_abs": baseline_imag_abs,
            "physical_delta_real_pct": delta_re,
            "physical_delta_imag_pct": delta_im,
            "max_abs_physical_delta_pct": max(abs(delta_re), abs(delta_im)),
            "spectral_distance_pct": math.hypot(delta_re, delta_im),
            "validation_delta_real_pct": validation_delta_real_pct,
            "validation_delta_imag_pct": validation_delta_imag_pct,
        }
    )


def threshold_crossings(rows: list[dict[str, object]], thresholds: tuple[float, ...]) -> list[dict[str, object]]:
    grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[f"{row['family']}|n={row['n']}|{row['parameter_name']}"].append(row)

    out: list[dict[str, object]] = []
    for key, group in grouped.items():
        group = sorted(group, key=lambda row: as_float(row["parameter_value"]))
        family, n_text, parameter_name = key.split("|", 2)
        x = [as_float(row["parameter_value"]) for row in group]
        y = [as_float(row["max_abs_physical_delta_pct"]) for row in group]

        for threshold in thresholds:
            crossing = math.nan
            status = "not_reached_in_sample"
            if y and y[0] >= threshold:
                crossing = x[0]
                status = "already_above_first_sample"
            for idx in range(1, len(group)):
                y0, y1 = y[idx - 1], y[idx]
                if (y0 - threshold) * (y1 - threshold) <= 0 and y0 != y1:
                    fraction = (threshold - y0) / (y1 - y0)
                    crossing = x[idx - 1] + fraction * (x[idx] - x[idx - 1])
                    status = "linear_interpolation_sparse_grid"
                    break
            out.append(
                {
                    "family": family,
                    "mode": n_text,
                    "parameter_name": parameter_name,
                    "threshold_pct": threshold,
                    "crossing_parameter_value": crossing,
                    "status": status,
                }
            )
    return out


def svg_escape(text: object) -> str:
    return (
        str(text)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def write_deviation_svg(path: Path, rows: list[dict[str, object]]) -> None:
    grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        if as_float(row["parameter_value"]) == 0.0 and as_float(row["max_abs_physical_delta_pct"]) == 0.0:
            continue
        key = f"{row['family']} n={row['n']}"
        grouped[key].append(row)

    summary = []
    for key, group in grouped.items():
        max_row = max(group, key=lambda row: as_float(row["max_abs_physical_delta_pct"]))
        summary.append((key, max_row))
    summary.sort(key=lambda item: as_float(item[1]["max_abs_physical_delta_pct"]), reverse=True)

    width = 1040
    row_h = 48
    top = 70
    label_w = 260
    plot_w = 620
    height = top + row_h * len(summary) + 70
    max_value = max((as_float(row["max_abs_physical_delta_pct"]) for _key, row in summary), default=1.0)
    scale_max = max(10.0, math.ceil(max_value / 10.0) * 10.0)

    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<style>text{font-family:Arial,Helvetica,sans-serif;font-size:14px;fill:#1f2933}.small{font-size:12px;fill:#52606d}.title{font-size:22px;font-weight:700}.grid{stroke:#d9e2ec;stroke-width:1}.thr{stroke:#d64545;stroke-width:1.5;stroke-dasharray:4 5}</style>',
        '<text x="30" y="34" class="title">Physical static-QNM deviations from baseline</text>',
        '<text x="30" y="56" class="small">Maximum sampled absolute deviation in Re(omega) or Abs(Im omega), not validation error.</text>',
    ]
    x0 = label_w + 45
    for tick in [0, 1, 3, 5, 10, 20, 50, 100]:
        if tick > scale_max:
            continue
        x = x0 + tick / scale_max * plot_w
        cls = "thr" if tick in {1, 3, 5, 10} else "grid"
        lines.append(f'<line x1="{x:.1f}" y1="{top - 14}" x2="{x:.1f}" y2="{height - 44}" class="{cls}"/>')
        lines.append(f'<text x="{x - 8:.1f}" y="{height - 24}" class="small">{tick:g}</text>')

    for idx, (label, row) in enumerate(summary):
        y = top + idx * row_h
        value = as_float(row["max_abs_physical_delta_pct"])
        bar_w = min(value / scale_max * plot_w, plot_w)
        parameter = f"{row['parameter_name']}={as_float(row['parameter_value']):g}"
        lines.append(f'<text x="30" y="{y + 19}" font-weight="700">{svg_escape(label)}</text>')
        lines.append(f'<text x="30" y="{y + 37}" class="small">largest sample: {svg_escape(parameter)}</text>')
        lines.append(f'<rect x="{x0:.1f}" y="{y + 8}" width="{bar_w:.1f}" height="24" fill="#2f80ed"/>')
        lines.append(f'<text x="{x0 + bar_w + 7:.1f}" y="{y + 25}" class="small">{value:.2f}%</text>')
    lines.append("</svg>")
    path.write_text("\n".join(lines), encoding="utf-8")


def write_report(
    path: Path,
    rows: list[dict[str, object]],
    crossings: list[dict[str, object]],
) -> None:
    grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[f"{row['family']} n={row['n']}"].append(row)

    lines = [
        "# Static QNM Physical-Deviation Report",
        "",
        "This report separates two quantities that should not be conflated:",
        "",
        "- validation error: how accurately the local solver reproduces a published table;",
        "- physical deviation: how far a metric's QNM spectrum moves away from a Schwarzschild or zero-parameter baseline.",
        "",
        "The physical-deviation numbers below use the published/reference frequencies where available. The validation deltas are retained only as a numerical-error diagnostic.",
        "",
        "## Largest Sampled Deviations",
        "",
        "| family / mode | largest sampled parameter | delta Re [%] | delta Abs(Im) [%] | max abs delta [%] | validation-error scale [%] |",
        "| --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for key in sorted(grouped):
        group = grouped[key]
        max_row = max(group, key=lambda row: as_float(row["max_abs_physical_delta_pct"]))
        validation_scale = max(
            abs(as_float(max_row["validation_delta_real_pct"])),
            abs(as_float(max_row["validation_delta_imag_pct"])),
        )
        lines.append(
            f"| {key} | {max_row['parameter_name']}={as_float(max_row['parameter_value']):g} | "
            f"{as_float(max_row['physical_delta_real_pct']):.3f} | "
            f"{as_float(max_row['physical_delta_imag_pct']):.3f} | "
            f"{as_float(max_row['max_abs_physical_delta_pct']):.3f} | "
            f"{validation_scale:.3f} |"
        )

    lines.extend(
        [
            "",
            "## Sampled Rows",
            "",
            "| family | parameter | mode | baseline | Re shift [%] | Abs(Im) shift [%] | max abs [%] |",
            "| --- | ---: | --- | --- | ---: | ---: | ---: |",
        ]
    )
    for row in rows:
        lines.append(
            f"| {row['family']} | {row['parameter_name']}={as_float(row['parameter_value']):g} | "
            f"l={row['ell']},n={row['n']} | {row['baseline_label']} | "
            f"{as_float(row['physical_delta_real_pct']):.3f} | "
            f"{as_float(row['physical_delta_imag_pct']):.3f} | "
            f"{as_float(row['max_abs_physical_delta_pct']):.3f} |"
        )

    lines.extend(
        [
            "",
            "## Sparse Threshold Crossings",
            "",
            "These crossings are linear interpolations over sparse published tables, not final exclusion limits.",
            "",
            "| family | mode | parameter | threshold [%] | crossing | status |",
            "| --- | --- | --- | ---: | ---: | --- |",
        ]
    )
    for row in crossings:
        crossing = as_float(row["crossing_parameter_value"])
        crossing_text = "nan" if math.isnan(crossing) else f"{crossing:.4g}"
        lines.append(
            f"| {row['family']} | {row['mode']} | {row['parameter_name']} | "
            f"{as_float(row['threshold_pct']):.0f} | {crossing_text} | {row['status']} |"
        )

    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- Sub-percent validation does not mean the metric is physically close to Schwarzschild.",
            "- In the current static benchmarks, physical QNM shifts range from the percent level to tens of percent across the sampled parameter domains.",
            "- This branch can therefore be discriminating, but only after a physically meaningful parameter range and an observational or synthetic tolerance are declared.",
            "- These static deviations remain a readiness/stress-test layer; they are not direct GW250114 rotating-remnant exclusions.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []

    # Tidal charge: compare with the Schwarzschild gravitational l=2,n=0 qnm baseline.
    tidal_rows = read_csv(ROOT / "results" / "tidal_charge_time_domain" / "tidal_charge_time_domain_summary.csv")
    add_deviation_row(
        rows,
        family="Braneworld_tidal_charge",
        case="tidal_charge_q0_axial_l2_n0_baseline",
        parameter_name="q_tidal",
        parameter_value=0.0,
        ell=2,
        n=0,
        source_frequency="Schwarzschild_qnm_baseline",
        omega_real=SCHWARZSCHILD_GRAV_L2_N0[0],
        omega_imag_abs=SCHWARZSCHILD_GRAV_L2_N0[1],
        baseline_label="Schwarzschild_qnm_l2_n0",
        baseline_real=SCHWARZSCHILD_GRAV_L2_N0[0],
        baseline_imag_abs=SCHWARZSCHILD_GRAV_L2_N0[1],
        validation_delta_real_pct=0.0,
        validation_delta_imag_pct=0.0,
    )
    for row in tidal_rows:
        add_deviation_row(
            rows,
            family="Braneworld_tidal_charge",
            case=row["case"],
            parameter_name="q_tidal",
            parameter_value=as_float(row["q_tidal"]),
            ell=int(row["ell"]),
            n=int(row["n"]),
            source_frequency=row["reference_source"],
            omega_real=as_float(row["reference_Momega_real"]),
            omega_imag_abs=as_float(row["reference_Momega_imag_abs"]),
            baseline_label="Schwarzschild_qnm_l2_n0",
            baseline_real=SCHWARZSCHILD_GRAV_L2_N0[0],
            baseline_imag_abs=SCHWARZSCHILD_GRAV_L2_N0[1],
            validation_delta_real_pct=as_float(row["delta_real_pct"]),
            validation_delta_imag_pct=as_float(row["delta_imag_pct"]),
        )

    # Bardeen: alpha=0 row is the same-paper Schwarzschild/WKB baseline.
    bardeen_rows = read_csv(ROOT / "results" / "bardeen_time_domain" / "bardeen_time_domain_summary.csv")
    bardeen_baseline = next(row for row in bardeen_rows if as_float(row["alpha"]) == 0.0)
    base_re = as_float(bardeen_baseline["reference_Momega_real"])
    base_im = as_float(bardeen_baseline["reference_Momega_imag_abs"])
    for row in bardeen_rows:
        add_deviation_row(
            rows,
            family="Bardeen_NED",
            case=row["case"],
            parameter_name="alpha",
            parameter_value=as_float(row["alpha"]),
            ell=int(row["ell"]),
            n=int(row["n"]),
            source_frequency=row["reference_source"],
            omega_real=as_float(row["reference_Momega_real"]),
            omega_imag_abs=as_float(row["reference_Momega_imag_abs"]),
            baseline_label="Bardeen_alpha_0_same_table",
            baseline_real=base_re,
            baseline_imag_abs=base_im,
            validation_delta_real_pct=as_float(row["delta_real_pct"]),
            validation_delta_imag_pct=as_float(row["delta_imag_pct"]),
        )

    # Hayward fundamental: compare with Schwarzschild gravitational l=2,n=0 qnm.
    hayward_rows = read_csv(ROOT / "results" / "hayward_time_domain" / "hayward_time_domain_summary.csv")
    add_deviation_row(
        rows,
        family="Hayward_fundamental",
        case="hayward_gamma_0_axial_l2_n0_baseline",
        parameter_name="gamma",
        parameter_value=0.0,
        ell=2,
        n=0,
        source_frequency="Schwarzschild_qnm_baseline",
        omega_real=SCHWARZSCHILD_GRAV_L2_N0[0],
        omega_imag_abs=SCHWARZSCHILD_GRAV_L2_N0[1],
        baseline_label="Schwarzschild_qnm_l2_n0",
        baseline_real=SCHWARZSCHILD_GRAV_L2_N0[0],
        baseline_imag_abs=SCHWARZSCHILD_GRAV_L2_N0[1],
        validation_delta_real_pct=0.0,
        validation_delta_imag_pct=0.0,
    )
    for row in hayward_rows:
        add_deviation_row(
            rows,
            family="Hayward_fundamental",
            case=row["case"],
            parameter_name="gamma",
            parameter_value=as_float(row["gamma"]),
            ell=int(row["ell"]),
            n=int(row["n"]),
            source_frequency=row["reference_source"],
            omega_real=as_float(row["reference_Momega_real"]),
            omega_imag_abs=as_float(row["reference_Momega_imag_abs"]),
            baseline_label="Schwarzschild_qnm_l2_n0",
            baseline_real=SCHWARZSCHILD_GRAV_L2_N0[0],
            baseline_imag_abs=SCHWARZSCHILD_GRAV_L2_N0[1],
            validation_delta_real_pct=as_float(row["delta_real_pct"]),
            validation_delta_imag_pct=as_float(row["delta_imag_pct"]),
        )

    # Hayward overtone: gamma=0 row is the same-table Schwarzschild-like baseline.
    overtone_rows = read_csv(
        ROOT
        / "results"
        / "hayward_overtone_matrix_pencil"
        / "hayward_overtone_matrix_pencil_summary.csv"
    )
    overtone_baseline = next(row for row in overtone_rows if as_float(row["gamma"]) == 0.0)
    base_re = as_float(overtone_baseline["reference_Momega_real"])
    base_im = as_float(overtone_baseline["reference_Momega_imag_abs"])
    for row in overtone_rows:
        add_deviation_row(
            rows,
            family="Hayward_overtone",
            case=row["case"],
            parameter_name="gamma",
            parameter_value=as_float(row["gamma"]),
            ell=int(row["ell"]),
            n=int(row["n"]),
            source_frequency=row["reference_source"],
            omega_real=as_float(row["reference_Momega_real"]),
            omega_imag_abs=as_float(row["reference_Momega_imag_abs"]),
            baseline_label="Hayward_gamma_0_same_table",
            baseline_real=base_re,
            baseline_imag_abs=base_im,
            validation_delta_real_pct=as_float(row["delta_real_pct"]),
            validation_delta_imag_pct=as_float(row["delta_imag_pct"]),
        )

    fields = [
        "family",
        "case",
        "parameter_name",
        "parameter_value",
        "ell",
        "n",
        "source_frequency",
        "Momega_real",
        "Momega_imag_abs",
        "baseline_label",
        "baseline_Momega_real",
        "baseline_Momega_imag_abs",
        "physical_delta_real_pct",
        "physical_delta_imag_pct",
        "max_abs_physical_delta_pct",
        "spectral_distance_pct",
        "validation_delta_real_pct",
        "validation_delta_imag_pct",
    ]
    summary_csv = OUT_DIR / "static_qnm_physical_deviations.csv"
    write_csv(summary_csv, rows, fields)

    crossings = threshold_crossings(rows, thresholds=(1.0, 3.0, 5.0, 10.0))
    crossing_fields = [
        "family",
        "mode",
        "parameter_name",
        "threshold_pct",
        "crossing_parameter_value",
        "status",
    ]
    crossings_csv = OUT_DIR / "static_qnm_physical_thresholds.csv"
    write_csv(crossings_csv, crossings, crossing_fields)

    report_path = OUT_DIR / "static_qnm_physical_deviations.md"
    write_report(report_path, rows, crossings)

    figure_path = OUT_DIR / "static_qnm_physical_deviations.svg"
    write_deviation_svg(figure_path, rows)

    print(f"Physical deviations CSV: {summary_csv}")
    print(f"Thresholds CSV: {crossings_csv}")
    print(f"Report: {report_path}")
    print(f"Figure: {figure_path}")

    grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[f"{row['family']} n={row['n']}"].append(row)
    for key in sorted(grouped):
        max_row = max(grouped[key], key=lambda row: as_float(row["max_abs_physical_delta_pct"]))
        print(
            f"{key}: max={as_float(max_row['max_abs_physical_delta_pct']):.2f}% "
            f"at {max_row['parameter_name']}={as_float(max_row['parameter_value']):g}"
        )


if __name__ == "__main__":
    main()
