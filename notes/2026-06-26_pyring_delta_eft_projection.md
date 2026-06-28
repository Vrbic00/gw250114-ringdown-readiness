# GW250114 pyRing Delta EFT Projection

This note records the first frequency-plus-damping observational projection in
the project.

## Implemented

- Added `config/gw250114_pyring_delta_eft_projection.wl`.
- Added `scripts/wolfram/gw250114_pyring_delta_eft_projection.wl`.
- Used the public `posterior_with_qnm_frequencies.dat` pyRing 220+221
  deviation posterior from the GW250114 spectroscopy release.
- Applied the same lower-tail `domega_221` filter used in the public Figure 4
  script and required positive `1 + domega_221` and `1 + dtau_221`.
- Built a Gaussian projection in
  `{log(1 + domega_221), log(1 + dtau_221)}`.

## Main Result

No one-at-a-time higher-derivative EFT operator gives a robust nonzero
deviation. The largest nominal shift is about `1.2 sigma`, and zero remains
inside every 90 percent interval.

Best constrained examples:

```text
epsilon1 minus: alpha = -0.02726 +/- 0.02381
epsilon2 minus: alpha = -0.02890 +/- 0.02386
epsilon2 plus:  alpha = -0.07656 +/- 0.06406
```

## Caveat

This is a Gaussian projection of phenomenological pyRing deviation samples.
It is useful for a constraints-table prototype, but it is not a full
strain-level EFT likelihood.
