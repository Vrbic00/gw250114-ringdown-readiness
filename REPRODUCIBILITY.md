# Reproducibility guide

Run all commands from the repository root. Generated outputs are written under
`results/`, which is intentionally ignored by Git.

## 1. Prepare derived inputs

Download the public inputs described in `DATA_SOURCES.md`, then create the
compact intermediate files:

```powershell
python scripts/python/article1_independent_prior.py
wolframscript -file scripts/wolfram/gw250114_posterior_calibration.wl
python scripts/python/article1_nr_calibration.py
```

The reconstruction path needs `requirements-optional.txt`, Wolfram Language,
and the external raw files.

## 2. Static supplied-potential validation

```powershell
python scripts/python/static_master_potential_time_domain.py
python scripts/python/tidal_charge_time_domain_benchmark.py
python scripts/python/bardeen_time_domain_benchmark.py
python scripts/python/hayward_time_domain_benchmark.py
python scripts/python/hayward_overtone_matrix_pencil.py
```

These are numerical-intake benchmarks, not physical models of the rotating
GW250114 remnant.

## 3. Factorised Dudley--Finley QNM validation and grid

```powershell
python scripts/python/hairy_continued_fraction.py
python scripts/python/build_hairy_qnm_production_grid.py
python scripts/python/hairy_qnm_internal_systematics.py --grid results/hairy_qnm_production_grid/hairy_qnm_production_grid.npz
python scripts/python/effective_kerr_newman_control.py
```

The production grid uses the settings in
`config/hairy_gw250114_publication.json` and contains 35,343 direct roots.
The final command audits the exact effective Kerr--Newman-form polynomial
identity, the signed charge-squared map, and the difference between the
table-calibrated and ODE-consistent recurrences. It does not validate the
missing coupled perturbation system.

## 4. Public-posterior likelihood model

This stage requires the pyRing posterior from `DATA_SOURCES.md`:

```powershell
python scripts/python/gw250114_bayesian_spectral_eft.py
```

It creates the compact GMM bundle used by the Li-prescription event analysis.

## 5. Event-level Li-prescription calculation

```powershell
python scripts/python/gw250114_hairy_constraints.py
```

This stage also uses the public RINGDOWN deviation and time-scan files listed
in `DATA_SOURCES.md`. pyRing and RINGDOWN are separate analyses of the same
event and are not multiplied as independent likelihoods.

## 6. Controls and robustness

```powershell
python scripts/python/hairy_evolving_kerr_control.py
python scripts/python/hairy_positive_injection_recovery.py
python scripts/python/hairy_referee_robustness.py
python scripts/python/hairy_publication_computation_audit.py
```

## Expected numerical checkpoints

Principal checkpoints reported for the publication snapshot include:

- 81 published table cases reproduced;
- 35,343 direct production-grid modes;
- maximum depth change below `1e-8`;
- maximum interpolation relative error below `5e-4`;
- effective Kerr--Newman-form polynomial identity at machine precision;
- successful stationary-Kerr and evolving-Kerr controls.

The exact public file inventory is recorded in `MANIFEST.csv`.

## What is deliberately absent

This repository does not contain the article source, bibliography, Overleaf
project, submission package, publication figures, large raw data, virtual
environments, caches, or the full generated-result tree.
