# Case Configuration

Edit `spherical_cases.wl` to add candidate metrics or supplied master
potentials.

For exploratory runs, copy `custom_cases_template.wl` to `custom_cases.wl` and
run:

```powershell
wolframscript -file scripts/wolfram/deviation_detectability_report.wl config/custom_cases.wl
```

Custom configs write to `results/deviation_detectability_<config-name>/` by
default. You can override that with a second argument:

```powershell
wolframscript -file scripts/wolfram/deviation_detectability_report.wl config/custom_cases.wl results/my_custom_report
```

Every custom file should include its own reference rows, for example
`Schwarzschild` for metric-only comparisons or `Schwarzschild_ReggeWheeler_l2`
for Regge-Wheeler-like supplied-potential comparisons.

## Metric-Only Case

Use this when only the background metric is known:

```wolfram
<|
  "name" -> "MyMetric",
  "kind" -> "candidate metric",
  "parameters" -> "M=1, alpha=0.1",
  "reference" -> "Schwarzschild",
  "horizonRole" -> "LargestPositive",
  "f" -> 1 - 2/r + 0.1/r^3
|>
```

This enables background checks and the eikonal QNM proxy only.

Use `"horizonRole" -> "SmallestPositive"` for Schwarzschild-de Sitter-like
cases where the largest positive root is the cosmological horizon.

## Supplied-Potential Case

Use this when a perturbation potential is known or proposed:

```wolfram
<|
  "name" -> "MyPotential_l2",
  "sector" -> "gravitational_candidate",
  "level" -> "Level 4 supplied potential",
  "description" -> "Candidate master potential",
  "reference" -> "Schwarzschild_ReggeWheeler_l2",
  "f" -> 1 - 2/r,
  "V" -> (1 - 2/r) (ell (ell + 1)/r^2 - 6/r^3) (1 + 0.1/r^2)
|>
```

This enables the WKB deviation/detectability report. A full Level 5 claim
requires theory-backed perturbation equations or a trusted master potential.

## Parameter Scans

Edit `parameter_scans.wl` to scan parameterized metric or potential families.

Run:

```powershell
wolframscript -file scripts/wolfram/parameter_threshold_scan.wl
```

The scanner writes a full sampled CSV, threshold-crossing CSV, plots, and a
Markdown report under:

```text
results/parameter_threshold_scan/
```

Tolerance profiles such as `optimistic_1pct`, `near_future_3pct`,
`screening_5pct`, and `loose_10pct` are illustrative precision scenarios.
They are not observational posterior widths unless explicitly calibrated later.

## Kerr Spin Scan

Edit `kerr_spin_scan.wl` to change the Kerr modes or spin grid used by:

```powershell
wolframscript -file scripts/wolfram/kerr_spin_contamination_scan.wl
```

This scan estimates how much ordinary Kerr spin moves QNM frequencies away from
the `a=0` Schwarzschild-limit reference.

## Event Registry

Edit `event_registry.csv` to add or update public events without changing code.
Then run:

```powershell
wolframscript -file scripts/wolfram/event_suitability_classifier.wl
```

The classifier is only a bookkeeping layer. It does not replace LVK parameter
estimation or a dedicated literature review for a newly added event.

To combine the registry with the Kerr spin-contamination scan, run:

```powershell
wolframscript -file scripts/wolfram/event_spin_floor_report.wl
```

This estimates whether each event is sensible for static-remnant, Kerr-ringdown,
or inspiral/progenitor use.

For pre-merger/progenitor diagnostics, run:

```powershell
wolframscript -file scripts/wolfram/inspiral_progenitor_diagnostics.wl
```

## Metric Intake Registry

Edit `metric_intake_registry.csv` when a new metric appears in the literature.
Then run:

```powershell
wolframscript -file scripts/wolfram/metric_intake_classifier.wl
```

The intake classifier decides whether the metric currently supports only
geodesic/QPO diagnostics, test-field QNMs, supplied gravitational QNMs, or a
full theory-backed ringdown analysis.

## Rotating Geodesic Cases

Edit `rotating_geodesic_cases.wl` to add stationary axisymmetric metrics in
Boyer-Lindquist-like coordinates. The current script expects metric components
`gtt`, `gtp`, `gpp`, `grr`, and `gthth` as functions of `{r, theta, a,
deformation}`.

Run:

```powershell
wolframscript -file scripts/wolfram/rotating_geodesic_diagnostics.wl
```

This computes equatorial horizon, ergosphere, prograde photon orbit, ISCO, and
ISCO orbital/vertical frequencies, then compares every case with Kerr at the
same spin.

To estimate grid-level threshold crossings from those results, run:

```powershell
wolframscript -file scripts/wolfram/rotating_geodesic_thresholds.wl
```

For the denser Johannsen scan at `a=0.68`, run:

```powershell
wolframscript -file scripts/wolfram/rotating_geodesic_diagnostics.wl config/rotating_geodesic_johannsen_dense.wl results/rotating_geodesic_johannsen_dense
wolframscript -file scripts/wolfram/rotating_geodesic_thresholds.wl results/rotating_geodesic_johannsen_dense/rotating_geodesic_deviations.csv results/rotating_geodesic_johannsen_dense_thresholds
```

For the rotating Bardeen geodesic and validity scan, run:

```powershell
wolframscript -file scripts/wolfram/rotating_geodesic_diagnostics.wl config/rotating_geodesic_bardeen.wl results/rotating_geodesic_bardeen
```

For the two-dimensional spin-deformation degeneracy grid, edit
`johannsen_spin_degeneracy.wl` and run:

```powershell
wolframscript -file scripts/wolfram/johannsen_spin_degeneracy.wl
```

For the rotating Bardeen spin-deformation degeneracy map, run:

```powershell
wolframscript -file scripts/wolfram/johannsen_spin_degeneracy.wl config/bardeen_spin_degeneracy.wl results/bardeen_spin_degeneracy
```

For the QPO observable layer, edit `qpo_observable_scan.wl` and run:

```powershell
wolframscript -file scripts/wolfram/qpo_observable_layer.wl
```

For the rotating Bardeen QPO observable layer, run:

```powershell
wolframscript -file scripts/wolfram/qpo_observable_layer.wl config/qpo_observable_bardeen.wl results/qpo_observables_bardeen
```

For the GRO J1655-40 observational RPM benchmark, edit
`qpo_observational_fit.wl` and run:

```powershell
wolframscript -file scripts/wolfram/qpo_observational_fit.wl
```

The measured triplet is stored externally in `data/qpo_observations.csv`.

For the XTE J1859+226 triplet, run:

```powershell
wolframscript -file scripts/wolfram/qpo_observational_fit.wl config/qpo_observational_fit_xte_j1859.wl results/qpo_observational_fit_xte_j1859
```

For the GRO J1655-40 profile-chi-square scan, edit
`qpo_profile_likelihood.wl` and run:

```powershell
wolframscript -file scripts/wolfram/qpo_profile_likelihood.wl
```

For the XTE J1859+226 statistical profile, run:

```powershell
wolframscript -file scripts/wolfram/qpo_profile_likelihood.wl config/qpo_profile_likelihood_xte_j1859.wl results/qpo_profile_likelihood_xte_j1859
```

For the horizon-constrained rotating Bardeen profiles, run:

```powershell
wolframscript -file scripts/wolfram/qpo_profile_likelihood.wl config/qpo_profile_likelihood_bardeen_gro.wl results/qpo_profile_likelihood_bardeen_gro
wolframscript -file scripts/wolfram/qpo_profile_likelihood.wl config/qpo_profile_likelihood_bardeen_xte_j1859.wl results/qpo_profile_likelihood_bardeen_xte_j1859
```

The XTE J1859+226 positive-`alpha13` extension and final merge use:

```powershell
wolframscript -file scripts/wolfram/qpo_profile_likelihood.wl config/qpo_profile_likelihood_xte_j1859_extension.wl results/qpo_profile_likelihood_xte_j1859_extension
wolframscript -file scripts/wolfram/qpo_profile_likelihood.wl config/qpo_profile_likelihood_xte_j1859_finalize.wl results/qpo_profile_likelihood_xte_j1859
```

The configured fractional QPO systematics are sensitivity scenarios added in
quadrature to the statistical frequency errors. They are not calibrated
accretion-flow uncertainties.

`qpo_profile_likelihood_validation.wl` contains a small multi-start regression
grid for checking the main continuation scan at its minima, Kerr points, and
profile edges.

For correlated QPO nuisance modes and a finite-width emitting annulus, edit
`qpo_robustness_profile.wl` and run:

```powershell
wolframscript -file scripts/wolfram/qpo_robustness_profile.wl
```

For the rotating Bardeen robustness profile of GRO J1655-40, run:

```powershell
wolframscript -file scripts/wolfram/qpo_robustness_profile.wl config/qpo_robustness_bardeen_gro.wl results/qpo_robustness_bardeen_gro
```

`qpo_robustness_profile_finalize.wl` reloads the main profile and merges
targeted validation rows before regenerating the final CSV, report, and plot.

For the conditional two-source joint deformation profile, edit
`qpo_joint_profile.wl` and run:

```powershell
wolframscript -file scripts/wolfram/qpo_joint_profile.wl
```

For the conditional common-`g/M` rotating Bardeen profile, run:

```powershell
wolframscript -file scripts/wolfram/qpo_joint_profile.wl config/qpo_joint_profile_bardeen.wl results/qpo_joint_profile_bardeen
```

## Higher-Derivative Ringdown Bridge

The initial Cano et al. (2023) 220/330 polynomial reproduction is:

```powershell
wolframscript -file scripts/wolfram/higher_derivative_qnm_bridge.wl
```

The complete selected-mode spectrum from the public BeyondKerrQNM fits is:

```powershell
wolframscript -file scripts/wolfram/higher_derivative_qnm_complete.wl
```

This second calculation includes `220`, `221`, `222`, `330`, and `440`, both
polarizations, and all five EFT operators. It uses the numerical Python `qnm`
Kerr sequences generated by `scripts/python/qnm_solver_crosscheck.py`.

The fit import can be regenerated from a checked-out BeyondKerrQNM repository:

```powershell
python scripts/python/import_beyond_kerr_qnm_fits.py <Fits-directory> data/beyond_kerr_qnm_selected_fits.csv --commit <commit-hash>
```

These scripts produce theory-backed spectral predictions, not a strain
likelihood or an observational exclusion.

For the first synthetic identifiability check using the complete `220+221`
spectra, run:

```powershell
wolframscript -file scripts/wolfram/ringdown_synthetic_likelihood.wl
```

This profiles each one-at-a-time EFT fingerprint over remnant mass and spin
using illustrative `{log f, log tau}` widths. It is a bridge toward an
event-level likelihood, not yet a GW250114 strain analysis.

For a stricter synthetic damped-sinusoid likelihood with linear
amplitude/phase profiling, start-time profiling, and a simple colored-noise
inner product, run:

```powershell
wolframscript -file scripts/wolfram/ringdown_toy_time_domain_likelihood.wl
```

After downloading the public preferred NRSur7dq4 posterior sample file
`data/raw/posterior_samples_NRSur7dq4.h5`, calibrate the remnant posterior and
Kerr QNM pushforward with:

```powershell
wolframscript -file scripts/wolfram/gw250114_posterior_calibration.wl
```

This summarizes the full IMR posterior and its induced Kerr QNM covariance; it
is not a ringdown-only EFT likelihood.

To compare free remnant profiling with a Gaussian remnant prior inferred from
that public NRSur7dq4 posterior, run:

```powershell
wolframscript -file scripts/wolfram/ringdown_posterior_informed_projection.wl
```

This remains synthetic because the `{log f, log tau}` measurement widths are
scenario choices. The posterior is used only as a GR-informed prior on
`M_f, chi`.

After extracting `data/raw/GW250114_data_release.tar.gz` into
`data/raw/GW250114_data_release/`, summarize the public ringdown-specific
products with:

```powershell
wolframscript -file scripts/wolfram/gw250114_public_ringdown_products.wl
```

This reads the public RINGDOWN and pyRing `220+221` deviation products and
creates a first one-dimensional `df_221 -> EFT alpha` proxy. That proxy is not
yet a full EFT likelihood.

For an approximate event-level Gaussian projection in the public RINGDOWN
variables `{log f_220, log f_221, df_221}`, run:

```powershell
wolframscript -file scripts/wolfram/gw250114_ringdown_eft_projection.wl
```

This profiles remnant mass and spin linearly and maps one EFT coupling at a
time into the public ringdown posterior covariance. It is still not a
strain-level EFT likelihood.

For a complementary pyRing projection that includes the 221 damping-time
deviation posterior, run:

```powershell
wolframscript -file scripts/wolfram/gw250114_pyring_delta_eft_projection.wl
```

This maps `{log(1 + domega_221), log(1 + dtau_221)}` to the imported
higher-derivative QNM frequency and damping-time sensitivities. It is a
phenomenological deviation-posterior projection, not a full EFT waveform
likelihood.

To build the side-by-side public constraints table and consistency check, run:

```powershell
wolframscript -file scripts/wolfram/gw250114_constraints_comparison.wl
```

This compares the RINGDOWN and pyRing projections without statistically
combining them, since they are analyses of the same event and are not
independent likelihoods.

To test robustness of the pyRing frequency-plus-damping projection against the
lower-tail `domega_221` filter choice, run:

```powershell
wolframscript -file scripts/wolfram/gw250114_pyring_filter_robustness.wl
```

The sweep includes stricter, public, looser, and positive-domain-only filter
choices.

To map the actual public posterior samples into linearized `alpha` posterior
samples, run:

```powershell
wolframscript -file scripts/wolfram/gw250114_linearized_posterior_projection.wl
```

This is a non-Gaussian sanity check of the projection layer. It still assumes
the same linearized EFT response and is not a full strain-level likelihood.

To build manuscript/supplement-ready CSV tables from the generated outputs,
run:

```powershell
wolframscript -file scripts/wolfram/gw250114_paper_tables.wl
```
