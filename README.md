# GW250114 ringdown-readiness framework

This repository contains the reproducible code and compact result products for
a draft paper on public GW250114 ringdown projections and static QNM-readiness
benchmarks for phenomenological metrics.

The central claim is deliberately conservative:

> Public GW250114 ringdown/spectroscopy products, projected onto published
> higher-derivative Kerr QNM fingerprints in a one-at-a-time linearized
> Gaussian approximation, show no robust beyond-Kerr deviation. Static
> supplied-potential examples are used as QNM-readiness benchmarks, not as
> observational constraints on the rotating GW250114 remnant.

## What is included

- `scripts/`: Wolfram Language and Python scripts used for the projection,
  validation, audit, and manuscript-table layers.
- `config/`: run configurations and small registries.
- `data/`: small machine-readable theory and candidate-metric registries.
- `results/`: compact CSV/Markdown/figure outputs used in the draft.
- `paper/`: current LaTeX manuscript source, bibliography, figures, and tables.
- `notes/`: selected project notes documenting the main scientific decisions.
- `DATA_SOURCES.md` and `REPRODUCIBILITY.md`: data provenance and run order.

## What is not included

Large public data products are intentionally not committed:

- GW250114 Zenodo/GWOSC tarballs and posterior samples,
- extracted raw HDF5/DAT products,
- local Python environments and package caches,
- private reference PDFs used only to tune author style,
- Overleaf ZIP exports and other transient build files.

See `DATA_SOURCES.md` for the data provenance and download policy.

## Minimal reproduction path

The public-data projection layer requires Wolfram Language and the public
GW250114 data products described in `DATA_SOURCES.md`.

```powershell
wolframscript -file scripts/wolfram/gw250114_posterior_calibration.wl
wolframscript -file scripts/wolfram/gw250114_public_ringdown_products.wl
wolframscript -file scripts/wolfram/gw250114_ringdown_eft_projection.wl
wolframscript -file scripts/wolfram/gw250114_pyring_delta_eft_projection.wl
wolframscript -file scripts/wolfram/gw250114_constraints_comparison.wl
wolframscript -file scripts/wolfram/gw250114_pyring_filter_robustness.wl
wolframscript -file scripts/wolfram/gw250114_linearized_posterior_projection.wl
wolframscript -file scripts/wolfram/gw250114_paper_tables.wl
```

The static supplied-potential validation and readiness layer is Python based:

```powershell
python scripts/python/static_master_potential_time_domain.py
python scripts/python/tidal_charge_time_domain_benchmark.py
python scripts/python/bardeen_time_domain_benchmark.py
python scripts/python/hayward_time_domain_benchmark.py
python scripts/python/hayward_overtone_matrix_pencil.py
python scripts/python/static_qnm_scorecard.py
python scripts/python/static_qnm_readiness_audit.py
python scripts/python/static_qnm_physical_deviation_report.py
```

## Repository status

This is a pre-submission research repository. The code and result tables are
kept public-facing, but the manuscript interpretation should still follow the
guardrails in `results/manuscript_package/claim_guardrails.md`.

Remote repository: https://github.com/Vrbic00/gw250114-ringdown-readiness
