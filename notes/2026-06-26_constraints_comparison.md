# GW250114 Constraints Comparison

This note records the first side-by-side public constraints table.

## Implemented

- Added `config/gw250114_constraints_comparison.wl`.
- Added `scripts/wolfram/gw250114_constraints_comparison.wl`.
- Compared the public RINGDOWN projection
  `{log f_220, log f_221, df_221}` with the pyRing projection
  `{log(1 + domega_221), log(1 + dtau_221)}`.
- Exported long-form constraints, operator-by-operator consistency checks,
  projection summaries, observable summaries, and a comparison plot.

## Main Result

- No one-at-a-time EFT row has `alpha = 0` outside the 90 percent interval.
- No RINGDOWN/pyRing operator pair is flagged as a 90 percent interval tension.
- The largest normalized projection difference is `0.535 sigma`.

## Publication Role

This is the current best draft constraints table. It should be presented as an
approximate public-data projection and pipeline consistency check, not as a
final LVK-style EFT bound.
