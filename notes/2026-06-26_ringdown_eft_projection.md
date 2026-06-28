# Public Ringdown EFT Projection

This note records the first approximate event-level EFT projection using
public GW250114 ringdown variables.

## Implemented

- Added `config/gw250114_ringdown_eft_projection.wl`.
- Added `scripts/wolfram/gw250114_ringdown_eft_projection.wl`.
- The observable vector is `{log f_220, log f_221, df_221}` from the public
  RINGDOWN HDF5 posterior.
- The likelihood is approximated as Gaussian around the posterior mean.
- Remnant mass and spin are profiled as linear nuisance directions.
- One higher-derivative EFT coupling is enabled at a time.

## Caveat

This is not a full EFT likelihood. It uses frequency-sector information only
and does not yet include damping-time deviations, amplitude correlations,
non-Gaussian structure, or strain-level likelihood evaluation.

## Next Step

Add the damping-time/deviation sector when the correct public variable mapping
is selected, then compare the RINGDOWN and pyRing projections under the same
observable model.
