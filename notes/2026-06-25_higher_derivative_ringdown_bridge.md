# Higher-Derivative Ringdown Bridge

## Scientific Input

The project now uses a genuine perturbation-theory framework rather than
promoting a metric-only deformation to gravitational ringdown.

Primary sources:

- Cano et al. (2023), arXiv:2307.07431
- Cano et al. (2024), arXiv:2409.04517
- public data: https://github.com/pacmn91/BeyondKerrQNM

The theory supplies:

- an EFT action through eight derivatives;
- perturbative rotating black-hole backgrounds;
- modified gravitational Teukolsky equations;
- complex QNM frequency shifts at first order in the EFT couplings.

The dimensionless couplings are

```text
alpha_ev  = (ell/M)^4 lambda_ev
alpha_odd = (ell/M)^4 lambda_odd
alpha_i   = (ell/M)^6 epsilon_i
```

and the spectral model is

```text
M omega_lmn = M omega_lmn^Kerr + alpha_q delta omega_lmn^q
              + O(alpha_q^2).
```

## Reproduction Layer

`scripts/wolfram/higher_derivative_qnm_bridge.wl` implements the 2023
polynomial tables for 220 and 330.

At `chi = 0.7`, all 16 stored polynomials reproduce the rounded values in
Tables V and VI. The largest complex difference is `0.001187`, consistent
with the three-decimal published values.

The Berti Kerr fits differ from the numerical Python `qnm` sequences at
`chi = 0.68` by about:

- `0.51%` in the real 220 frequency and `0.79%` in its imaginary part;
- `0.53%` in the real 330 frequency and `0.74%` in its imaginary part.

The complete bridge therefore uses numerical `qnm` Kerr values.

## Complete Selected Spectrum

The public 2024 fit repository was inspected on 2026-06-25:

```text
commit 0afe6281bec6a6224bfd55fe4600d5966c6a7135
license GPL-3.0
```

The repository reports about 27.3 MB in total, but its full `Fits` directory
contains only about 64 kB. The local generated CSV selects:

```text
220, 221, 222, 330, 440
```

for all five operators and both polarizations. It contains 570 polynomial
coefficient rows and records source provenance per row.

`scripts/wolfram/higher_derivative_qnm_complete.wl` evaluates these fits at
the GW250114 central remnant values:

```text
M_f,det = 68.1 Msun
chi     = 0.68
```

## Numerical Findings

The maximum complex relative sensitivity per unit coupling is:

| mode | max sensitivity | operator / polarization |
| --- | ---: | --- |
| 220 | 3.478 | epsilon2 minus |
| 221 | 7.639 | epsilon1 minus |
| 222 | 16.144 | epsilon1 minus |
| 330 | 5.824 | epsilon2 minus |
| 440 | 9.001 | epsilon1 plus |

Overtone sensitivity is operator and polarization dependent. Relative to 220:

- the 221 amplification ranges from about `0.54` to `3.30`;
- the 222 amplification ranges from about `0.75` to `6.97`.

Thus, "overtones are more sensitive" is a useful tendency, not a universal
statement for every operator and polarization.

## Interpretation

Reproduction:

- published 2023 polynomial tables and 2024 public fits.

Numerical checks:

- paper-table regression at `chi = 0.7`;
- numerical Kerr control with Python `qnm`;
- parity-breaking branches in the public fits are opposite as expected.

Physical interpretation:

- the project now has theory-backed complex frequency shifts for both 220 and
  221, matching the principal mode content used in the GW250114 spectroscopy
  analysis.

Not yet established:

- no EFT coupling has been constrained by GW250114;
- no mode amplitudes or excitation coefficients are modeled;
- no detector likelihood, start-time marginalization, or noise treatment is
  included;
- first-order validity becomes more delicate for sensitive overtones;
- the public repository commit is newer than the paper and is therefore
  recorded exactly.

## Next Block

Construct a synthetic 220+221 injection/recovery likelihood:

1. generate a Kerr injection at the GW250114 central remnant parameters;
2. fit mass, spin, amplitudes, phases, and one EFT coupling at a time;
3. test identifiability and mass-spin-coupling degeneracies;
4. vary ringdown start time and SNR;
5. proceed to real GW250114 data only if the synthetic recovery is stable.
