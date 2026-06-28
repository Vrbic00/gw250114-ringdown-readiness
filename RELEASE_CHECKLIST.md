# Release checklist

Use this checklist before tagging a citable release or archiving the repository
with Zenodo.

## Before release

- Run `python scripts/python/check_release_package.py`.
- Run `python -m compileall -q scripts/python`.
- Confirm `README.md`, `DATA_SOURCES.md`, and `REPRODUCIBILITY.md` match the
  manuscript version.
- Confirm `paper/main.tex` contains the final repository or Zenodo DOI.
- Confirm `CITATION.cff` has the final title, authors, version, release date,
  and DOI.
- Confirm no raw data, local environments, ZIP exports, private PDFs, or cache
  files are tracked.

## Suggested versioning

- `v0.1.0`: pre-submission reproducibility snapshot.
- `v1.0.0`: version matching the submitted manuscript.
- `v1.1.0`: version revised after peer review.

## Zenodo

After the GitHub release is created, archive it through Zenodo and copy the DOI
back into `README.md`, `CITATION.cff`, and the manuscript data/code
availability statement.
