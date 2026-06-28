# GW250114 Linearized Posterior Projection

This note records the first non-Gaussian sanity check of the projected EFT
constraints.

## Implemented

- Added `config/gw250114_linearized_posterior_projection.wl`.
- Added `scripts/wolfram/gw250114_linearized_posterior_projection.wl`.
- Mapped every public RINGDOWN posterior sample to a linearized `alpha`
  estimate with mass/spin nuisance profiling.
- Mapped every filtered pyRing 221 deviation posterior sample to a linearized
  `alpha` estimate in `{log(1 + domega_221), log(1 + dtau_221)}`.

## Main Result

- No empirical 90 percent interval excludes `alpha = 0`.
- Maximum empirical `abs(median)/sd`:
  - RINGDOWN: `0.485`
  - pyRing: `1.286`
- The no-deviation conclusion is therefore not an artifact of using symmetric
  Gaussian intervals.

## Caveat

This is still a linearized projection and not a full strain-level EFT
likelihood.
