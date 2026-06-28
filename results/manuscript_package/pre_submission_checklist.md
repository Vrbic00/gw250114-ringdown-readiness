# Pre-Submission Checklist

## Computation

- Re-run the public-data projection workflow from `REPRODUCIBILITY.md`.
- Re-run the static validation and readiness scripts.
- Re-run `pre_manuscript_audit.py`; require zero STOP items.
- Confirm all main-table numbers match regenerated CSV files.
- Confirm all figures referenced in `overleaf_roman_package_v10/main.tex` exist
  and are current.

## Manuscript

- Treat `manuscript_v1_author_style.md` as the active source draft.
- Keep `Results and Discussion` combined unless a target journal or coauthor
  asks for a different structure.
- Keep the public GW250114 branch framed as a projection, not as a full
  strain-level beyond-Kerr inference.
- Keep the static QNM branch framed as a readiness benchmark, not as a
  constraint on rotating GW250114 remnants.
- Keep the inserted per-operator projection table unless a target journal asks
  for it to be moved to a supplement.
- Add data/code availability statement with public release URLs and local
  repository/archive information.
- Replace the placeholder repository/archive statement by a real Git or Zenodo
  link before submission.

## Referee Risk Review

- Check that "projection" is not accidentally upgraded to "constraint" where a
  full likelihood is required.
- Check that RINGDOWN and pyRing are never combined as independent likelihoods.
- Check that static metrics are never described as ruled out by GW250114.
- Check that WKB/Prony reproductions are not called exact spectra.
- Check that threshold crossings are called diagnostic, not exclusion limits.
- Check that the static branch is motivated by ringdown-readiness, not as a
  detached catalog of metrics.

## Journal Strategy

- PRD route: emphasize public reproducibility, reusable tables/code, and a
  clear readiness standard.
- CQG fallback: emphasize ringdown methodology, perturbation-theory hygiene,
  and QNM validation.
- Avoid a title that sounds like direct constraints on Bardeen/Hayward from
  GW250114.

## Before Submission

- Decide author list and acknowledgements.
- Decide whether to upload code/data snapshot to Zenodo.
- Send the current Overleaf/Roman review package to Roman Konoplya only as an
  invitation to comment or join, not as a pre-filled coauthorship claim.
- Prepare arXiv version and choose arXiv categories.
- Verify journal policy for non-OA publication and arXiv posting.
- Run one final independent read as a hostile referee.
