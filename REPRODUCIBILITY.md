# Reproducibility notes

The workflow has two independent layers.

## 1. Public GW250114 projection

This layer projects public RINGDOWN/pyRing posterior products onto published
higher-derivative Kerr QNM fingerprints. It is a linearized public-data
projection, not a full strain-level EFT likelihood.

Main outputs:

```text
results/gw250114_constraints_comparison/projection_constraints_long.csv
results/gw250114_constraints_comparison/projection_consistency_by_operator.csv
results/gw250114_paper_tables/table1_main_projected_constraints.csv
```

Main interpretation:

- `alpha = 0` remains inside every one-at-a-time 90 percent projected interval.
- RINGDOWN and pyRing branches are used as consistency checks and are not
  statistically combined as independent likelihoods.
- Gaussian intervals may be narrower than a full strain-level Bayesian
  analysis because nonlinear mass-spin-coupling correlations and prior-volume
  effects are not modeled.

## 2. Static QNM-readiness benchmarks

This layer tests whether supplied master-potential examples can reproduce
published QNM frequencies. It is not an observational constraint on GW250114.

Main outputs:

```text
results/static_qnm_scorecard/static_qnm_validation_summary.csv
results/static_qnm_readiness_audit/static_metric_readiness_audit.csv
results/static_qnm_physical_deviations/static_qnm_physical_deviations.csv
```

Main interpretation:

- The current supplied-potential validations are sub-percent checks against
  published references.
- Physical deviations from Schwarzschild or zero-parameter baselines are model
  comparisons, not detector residuals or exclusion statistics.
- Metric-only, QPO-only, or shadow-only models are not gravitational-ringdown
  ready unless they also provide perturbation equations, boundary conditions,
  and reproducible QNM spectra.
