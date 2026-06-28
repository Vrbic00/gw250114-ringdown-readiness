"""Generate a forest plot for the public GW250114 alpha projections."""

from __future__ import annotations

import csv
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.pagesizes import landscape
from reportlab.pdfgen import canvas


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "results" / "gw250114_constraints_comparison" / "projection_constraints_long.csv"
OUT_DIR = ROOT / "results" / "gw250114_constraints_comparison"
PDF_PATH = OUT_DIR / "projection_alpha_interval_forest.pdf"
SVG_PATH = OUT_DIR / "projection_alpha_interval_forest.svg"


def read_rows() -> list[dict[str, str]]:
    with SOURCE.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    order = {"RINGDOWN": 0, "PYRING_DELTA": 1}
    op_order = {"lambda_ev": 0, "lambda_odd": 1, "epsilon1": 2, "epsilon2": 3, "epsilon3": 4}
    pol_order = {"plus": 0, "minus": 1}
    return sorted(
        rows,
        key=lambda r: (order[r["projection"]], op_order[r["operator"]], pol_order[r["polarization"]]),
    )


def as_float(row: dict[str, str], key: str) -> float:
    return float(row[key])


def branch_label(value: str) -> str:
    return {"RINGDOWN": "RINGDOWN", "PYRING_DELTA": "pyRing"}.get(value, value)


def operator_label(row: dict[str, str]) -> str:
    ops = {
        "lambda_ev": "lambda_ev",
        "lambda_odd": "lambda_odd",
        "epsilon1": "epsilon_1",
        "epsilon2": "epsilon_2",
        "epsilon3": "epsilon_3",
    }
    sign = {"plus": "+", "minus": "-"}.get(row["polarization"], row["polarization"])
    return f"{ops.get(row['operator'], row['operator'])}^{sign}"


def projection_color(value: str):
    return colors.HexColor("#2f80ed") if value == "RINGDOWN" else colors.HexColor("#c75b12")


def draw_pdf(rows: list[dict[str, str]]) -> None:
    width, height = landscape((960, 690))
    pdf = canvas.Canvas(str(PDF_PATH), pagesize=(width, height))
    dark = colors.HexColor("#1f2933")
    muted = colors.HexColor("#52606d")
    grid = colors.HexColor("#d9e2ec")

    margin_left = 170
    margin_right = 45
    margin_top = 72
    margin_bottom = 56
    plot_width = width - margin_left - margin_right
    row_gap = (height - margin_top - margin_bottom) / (len(rows) - 1)
    x_min, x_max = -4.0, 4.0

    def xcoord(x: float) -> float:
        return margin_left + (x - x_min) / (x_max - x_min) * plot_width

    pdf.setFillColor(colors.white)
    pdf.rect(0, 0, width, height, fill=True, stroke=False)

    pdf.setFillColor(dark)
    pdf.setFont("Helvetica-Bold", 20)
    pdf.drawString(32, height - 34, "Public GW250114 one-at-a-time alpha projections")
    pdf.setFillColor(muted)
    pdf.setFont("Helvetica", 11)
    pdf.drawString(32, height - 52, "Horizontal bars show 90 percent public-product intervals; the vertical line marks alpha = 0.")

    for tick in range(-4, 5):
        x = xcoord(tick)
        pdf.setStrokeColor(grid if tick != 0 else colors.HexColor("#334e68"))
        pdf.setLineWidth(0.8 if tick != 0 else 1.4)
        pdf.line(x, margin_bottom - 5, x, height - margin_top + 10)
        pdf.setFillColor(muted)
        pdf.setFont("Helvetica", 9)
        pdf.drawCentredString(x, margin_bottom - 25, str(tick))

    pdf.setFillColor(muted)
    pdf.setFont("Helvetica", 10)
    pdf.drawCentredString(margin_left + plot_width / 2, 18, "projected coupling alpha")

    for i, row in enumerate(rows):
        y = height - margin_top - i * row_gap
        color = projection_color(row["projection"])
        label = f"{branch_label(row['projection'])}: {operator_label(row)}"

        pdf.setFillColor(dark)
        pdf.setFont("Helvetica", 10)
        pdf.drawRightString(margin_left - 14, y - 3, label)

        lo = as_float(row, "alpha_q05")
        hi = as_float(row, "alpha_q95")
        best = as_float(row, "alpha_best")
        x1, x2, xb = xcoord(lo), xcoord(hi), xcoord(best)

        pdf.setStrokeColor(color)
        pdf.setLineWidth(2.1)
        pdf.line(x1, y, x2, y)
        pdf.setLineWidth(1.2)
        pdf.line(x1, y - 4, x1, y + 4)
        pdf.line(x2, y - 4, x2, y + 4)
        pdf.setFillColor(color)
        pdf.circle(xb, y, 3.4, stroke=False, fill=True)

    pdf.setFillColor(colors.HexColor("#2f80ed"))
    pdf.circle(width - 170, height - 34, 4, stroke=False, fill=True)
    pdf.setFillColor(muted)
    pdf.setFont("Helvetica", 10)
    pdf.drawString(width - 160, height - 38, "RINGDOWN")
    pdf.setFillColor(colors.HexColor("#c75b12"))
    pdf.circle(width - 86, height - 34, 4, stroke=False, fill=True)
    pdf.setFillColor(muted)
    pdf.drawString(width - 76, height - 38, "pyRing")

    pdf.showPage()
    pdf.save()


def esc(text: str) -> str:
    return (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def draw_svg(rows: list[dict[str, str]]) -> None:
    width, height = 960, 690
    margin_left = 170
    margin_right = 45
    margin_top = 72
    margin_bottom = 56
    plot_width = width - margin_left - margin_right
    row_gap = (height - margin_top - margin_bottom) / (len(rows) - 1)
    x_min, x_max = -4.0, 4.0

    def xcoord(x: float) -> float:
        return margin_left + (x - x_min) / (x_max - x_min) * plot_width

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        "<style>text{font-family:Arial,Helvetica,sans-serif}.title{font-size:20px;font-weight:700;fill:#1f2933}.small{font-size:11px;fill:#52606d}.label{font-size:10px;fill:#1f2933}.tick{font-size:9px;fill:#52606d}</style>",
        '<rect width="100%" height="100%" fill="white"/>',
        '<text x="32" y="34" class="title">Public GW250114 one-at-a-time alpha projections</text>',
        '<text x="32" y="52" class="small">Horizontal bars show 90 percent public-product intervals; the vertical line marks alpha = 0.</text>',
    ]
    for tick in range(-4, 5):
        x = xcoord(tick)
        color = "#334e68" if tick == 0 else "#d9e2ec"
        width_line = "1.4" if tick == 0 else "0.8"
        parts.append(f'<line x1="{x:.1f}" y1="{margin_bottom - 5}" x2="{x:.1f}" y2="{height - margin_top + 10}" stroke="{color}" stroke-width="{width_line}"/>')
        parts.append(f'<text x="{x:.1f}" y="{margin_bottom - 25}" class="tick" text-anchor="middle">{tick}</text>')
    parts.append(f'<text x="{margin_left + plot_width / 2:.1f}" y="18" class="small" text-anchor="middle">projected coupling alpha</text>')

    for i, row in enumerate(rows):
        y = height - margin_top - i * row_gap
        color = "#2f80ed" if row["projection"] == "RINGDOWN" else "#c75b12"
        label = f"{branch_label(row['projection'])}: {operator_label(row)}"
        lo = as_float(row, "alpha_q05")
        hi = as_float(row, "alpha_q95")
        best = as_float(row, "alpha_best")
        x1, x2, xb = xcoord(lo), xcoord(hi), xcoord(best)
        parts.append(f'<text x="{margin_left - 14}" y="{y + 3:.1f}" class="label" text-anchor="end">{esc(label)}</text>')
        parts.append(f'<line x1="{x1:.1f}" y1="{y:.1f}" x2="{x2:.1f}" y2="{y:.1f}" stroke="{color}" stroke-width="2.1"/>')
        parts.append(f'<line x1="{x1:.1f}" y1="{y - 4:.1f}" x2="{x1:.1f}" y2="{y + 4:.1f}" stroke="{color}" stroke-width="1.2"/>')
        parts.append(f'<line x1="{x2:.1f}" y1="{y - 4:.1f}" x2="{x2:.1f}" y2="{y + 4:.1f}" stroke="{color}" stroke-width="1.2"/>')
        parts.append(f'<circle cx="{xb:.1f}" cy="{y:.1f}" r="3.4" fill="{color}"/>')

    parts.extend(
        [
            '<circle cx="790" cy="34" r="4" fill="#2f80ed"/>',
            '<text x="800" y="38" class="small">RINGDOWN</text>',
            '<circle cx="874" cy="34" r="4" fill="#c75b12"/>',
            '<text x="884" y="38" class="small">pyRing</text>',
            "</svg>",
        ]
    )
    SVG_PATH.write_text("\n".join(parts), encoding="utf-8")


def main() -> None:
    rows = read_rows()
    draw_pdf(rows)
    draw_svg(rows)
    print(f"Projection forest PDF: {PDF_PATH}")
    print(f"Projection forest SVG: {SVG_PATH}")


if __name__ == "__main__":
    main()
