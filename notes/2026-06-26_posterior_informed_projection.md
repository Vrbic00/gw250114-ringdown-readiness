# Posterior-Informed Spectral Projection

This note records the first use of the public GW250114 NRSur7dq4 posterior as
a remnant prior inside the synthetic EFT spectral projection.

## Implemented

- Added `config/ringdown_posterior_informed_projection.wl`.
- Added `scripts/wolfram/ringdown_posterior_informed_projection.wl`.
- The script reads `nrSur7dq4_selected_posterior_samples.csv`, estimates the
  covariance of `{delta ln M_f, delta chi}`, and compares:
  - free remnant profiling;
  - the direct NRSur7dq4 IMR prior;
  - a loose 3x-wider prior.

## Scientific Caveat

The NRSur7dq4 posterior is a GR full-IMR product. Using it as a prior helps
calibrate scales and break degeneracies, but it is not an independent
ringdown-only no-hair test.

## Next Step

Search for a public ringdown-specific posterior or mode-frequency posterior
for GW250114. If none is available in compact form, keep the current
observational use explicitly IMR-prior-informed.
