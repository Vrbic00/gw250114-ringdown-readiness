# GW250114 pyRing Filter Robustness

This note records the first robustness sweep for the pyRing
frequency-plus-damping projection.

## Implemented

- Added `config/gw250114_pyring_filter_robustness.wl`.
- Added `scripts/wolfram/gw250114_pyring_filter_robustness.wl`.
- Swept the lower-tail `domega_221` filter through stricter, public, looser,
  and positive-domain-only choices.

## Main Result

- No scenario pushes `alpha = 0` outside a 90 percent interval.
- The maximum nominal shift over all tested scenarios is `1.29 sigma`.
- For every operator branch, the common 90 percent interval over all tested
  filters still contains zero.

## Publication Role

This supports the statement that the current no-deviation conclusion is not an
artifact of the public pyRing lower-tail filter choice. It does not replace a
full non-Gaussian posterior treatment or strain-level likelihood.
