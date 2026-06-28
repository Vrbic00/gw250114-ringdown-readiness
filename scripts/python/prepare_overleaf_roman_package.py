"""Build an Overleaf-ready ZIP from the Roman review package."""

from __future__ import annotations

import shutil
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.pagesizes import landscape
from reportlab.pdfgen import canvas


ROOT = Path(__file__).resolve().parents[2]
PACKAGE_DIR = ROOT / "results" / "manuscript_package"
ROMAN_DIR = PACKAGE_DIR / "roman_review_package"
OUT_DIR = PACKAGE_DIR / "overleaf_roman_package_v10"
FIGURE_DIR = OUT_DIR / "figures"
TABLE_DIR = OUT_DIR / "tables"
NOTES_DIR = OUT_DIR / "notes"


def clean_out_dir() -> None:
    if OUT_DIR.exists():
        shutil.rmtree(OUT_DIR)
    FIGURE_DIR.mkdir(parents=True)
    TABLE_DIR.mkdir(parents=True)
    NOTES_DIR.mkdir(parents=True)


def draw_readiness_pdf(path: Path) -> None:
    width, height = 1100, 250
    pdf = canvas.Canvas(str(path), pagesize=landscape((width, height)))
    pdf.setFillColor(colors.white)
    pdf.rect(0, 0, width, height, fill=True, stroke=False)

    dark = colors.HexColor("#1f2933")
    small = colors.HexColor("#52606d")
    stroke = colors.HexColor("#bcccdc")
    arrow = colors.HexColor("#52606d")
    fills = [
        colors.HexColor("#f5f7fa"),
        colors.HexColor("#f5f7fa"),
        colors.HexColor("#fff8e1"),
        colors.HexColor("#e3fcef"),
        colors.HexColor("#e7f0ff"),
    ]
    strokes = [
        stroke,
        stroke,
        colors.HexColor("#f2a900"),
        colors.HexColor("#2f9e44"),
        colors.HexColor("#2f80ed"),
    ]
    boxes = [
        (28, 63, 190, 115),
        (240, 63, 190, 115),
        (452, 63, 190, 115),
        (664, 63, 190, 115),
        (876, 63, 190, 115),
    ]
    labels = [
        ("Metric/geodesic only", ["horizon, photon sphere,", "ISCO, QPO, shadow"]),
        ("Test-field QNM", ["scalar or", "electromagnetic fields", "on fixed background"]),
        ("Supplied master potential", ["axial/polar", "gravitational potential", "and boundary conditions"]),
        ("Validated static QNM", ["local solver reproduces", "published spectra below", "target error"]),
        ("Theory-backed rotating QNM", ["action, perturbations,", "Kerr limit, mode fits", "for event projection"]),
    ]

    pdf.setFillColor(dark)
    pdf.setFont("Helvetica-Bold", 22)
    pdf.drawString(28, height - 34, "Static metric ringdown-readiness ladder")
    pdf.setFillColor(small)
    pdf.setFont("Helvetica", 12)
    pdf.drawString(28, height - 56, "A line element alone is not a gravitational ringdown model.")

    for i, (x, y, w, h) in enumerate(boxes):
        pdf.setFillColor(fills[i])
        pdf.setStrokeColor(strokes[i])
        pdf.roundRect(x, y, w, h, 6, fill=True, stroke=True)
        title, body = labels[i]
        pdf.setFillColor(dark)
        pdf.setFont("Helvetica-Bold", 13)
        pdf.drawString(x + 12, y + h - 26, title)
        pdf.setFillColor(small)
        pdf.setFont("Helvetica", 12)
        for j, line in enumerate(body):
            pdf.drawString(x + 12, y + h - 52 - 18 * j, line)
        if i < len(boxes) - 1:
            x1 = x + w + 5
            x2 = boxes[i + 1][0] - 5
            yy = y + h / 2
            pdf.setStrokeColor(arrow)
            pdf.setFillColor(arrow)
            pdf.setLineWidth(1.6)
            pdf.line(x1, yy, x2 - 8, yy)
            pdf.line(x2 - 8, yy + 5, x2, yy)
            pdf.line(x2 - 8, yy - 5, x2, yy)

    pdf.showPage()
    pdf.save()


def copy_files() -> None:
    tex = (ROMAN_DIR / "manuscript_for_roman.tex").read_text(encoding="utf-8")
    tex = tex.replace(r"\bibliography{references_for_roman}", r"\bibliography{references}")
    (OUT_DIR / "main.tex").write_text(tex, encoding="utf-8")
    shutil.copyfile(ROMAN_DIR / "references_for_roman.bib", OUT_DIR / "references.bib")

    for source in (ROMAN_DIR / "figures").glob("*"):
        if source.is_file():
            shutil.copyfile(source, FIGURE_DIR / source.name)
    draw_readiness_pdf(FIGURE_DIR / "static_metric_readiness_ladder.pdf")

    for source in (ROMAN_DIR / "tables").glob("*"):
        if source.is_file():
            shutil.copyfile(source, TABLE_DIR / source.name)

    for name in [
        "README.md",
        "executive_summary_for_roman.md",
        "cover_email_to_roman.md",
        "open_points_before_submission.md",
        "manuscript_for_roman.md",
    ]:
        source = ROMAN_DIR / name
        if source.exists():
            target_name = "README_OVERLEAF.md" if name == "README.md" else name
            shutil.copyfile(source, NOTES_DIR / target_name)

    (OUT_DIR / "README_OVERLEAF.md").write_text(overleaf_readme(), encoding="utf-8")


def overleaf_readme() -> str:
    return """# Overleaf Roman Review Package

Upload this ZIP to Overleaf as a new project.

## Main file

Set `main.tex` as the main document. The recommended compiler is pdfLaTeX.

## Required structure

```text
main.tex
references.bib
figures/
  projection_alpha_interval_forest.pdf
  projection_alpha_interval_forest.svg
  projection_sigma_from_zero_comparison.png
  static_metric_readiness_ladder.pdf
  static_metric_readiness_ladder.svg
  static_qnm_physical_deviations.svg
  static_qnm_validation_scorecard.svg
tables/
  main_table1_public_projection_summary.csv
  main_table2_static_readiness_summary.csv
  projection_consistency_by_operator.csv
  projection_constraints_long.csv
  projection_constraints_summary.csv
  static_metric_readiness_audit.csv
  static_qnm_physical_deviations.csv
notes/
  README_OVERLEAF.md
  executive_summary_for_roman.md
  cover_email_to_roman.md
  open_points_before_submission.md
  manuscript_for_roman.md
```

Only `main.tex`, `references.bib`, and the files in `figures/` are needed for
compilation. The `tables/` and `notes/` folders are included for Roman's review.

The bibliography file is based on the user's existing Mark-style `references.bib`.
New entries added for this project use keys such as
`Ber-Car-Wil:2006:PRD:` and `Can-Cap-Fra:2024:ARXIV:`.

## If compilation complains

- Make sure `main.tex` is selected as the main file.
- Use pdfLaTeX.
- Recompile once after bibliography warnings.
- If REVTeX is not available, the temporary fallback is to replace the first
  line by `\\documentclass[twocolumn]{article}`, but this is only for reading,
  not for journal submission.
"""


def make_zip() -> Path:
    zip_path = Path(shutil.make_archive(str(PACKAGE_DIR / "overleaf_roman_package_v10"), "zip", root_dir=OUT_DIR))
    return zip_path


def main() -> None:
    clean_out_dir()
    copy_files()
    zip_path = make_zip()
    print(f"Overleaf package: {OUT_DIR}")
    print(f"ZIP: {zip_path}")


if __name__ == "__main__":
    main()
