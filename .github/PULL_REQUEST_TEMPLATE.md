## Summary

Describe the scientific or reproducibility change.

## Checks

- [ ] `python scripts/python/check_release_package.py`
- [ ] `python -m compileall -q scripts/python`
- [ ] New or changed claims follow `results/manuscript_package/claim_guardrails.md`
- [ ] No raw GW data, local environments, PDFs, ZIPs, or cache files are committed

## Interpretation boundary

State whether the change affects public GW250114 projections, static
QNM-readiness benchmarks, manuscript text, or only repository packaging.
