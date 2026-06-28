# Overleaf Roman Review Package

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
  line by `\documentclass[twocolumn]{article}`, but this is only for reading,
  not for journal submission.
