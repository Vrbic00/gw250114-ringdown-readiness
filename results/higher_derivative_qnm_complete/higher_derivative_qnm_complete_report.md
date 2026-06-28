# Complete Higher-Derivative QNM Spectrum for GW250114

This report evaluates the public polynomial fits associated with Cano et al. (2024), arXiv:2409.04517.

## Data Provenance

- Repository: https://github.com/pacmn91/BeyondKerrQNM
- Imported commit: `0afe6281bec6a6224bfd55fe4600d5966c6a7135`
- License: `GPL-3.0`
- Selected modes: `220, 221, 222, 330, 440`
- The repository commit is newer than the paper release; the exact commit is recorded to keep the calculation reproducible.

## Physical Model

- The source supplies a higher-derivative EFT action, perturbative rotating backgrounds, modified gravitational Teukolsky equations, and numerical QNM fits.
- The fits contain both polarizations for all five EFT operators and include the 221 and 222 overtones.
- Frequencies obey `M omega = M omega_Kerr + alpha_q delta omega_q + O(alpha_q^2)`.
- The Kerr baseline below uses the numerical Python `qnm` sequences, not the approximate Berti fit.

## GW250114 Kerr Baseline

| mode | Re(Momega) | Im(Momega) | f [Hz] | tau [ms] | Q |
| --- | ---: | ---: | ---: | ---: | ---: |
| 220 | 0.523975 | -0.081513 | 248.619 | 4.115 | 3.2141 |
| 221 | 0.511941 | -0.246529 | 242.909 | 1.361 | 1.0383 |
| 222 | 0.489747 | -0.416451 | 232.378 | 0.805 | 0.5880 |
| 330 | 0.830913 | -0.083698 | 394.257 | 4.008 | 4.9637 |
| 440 | 1.125432 | -0.085118 | 534.002 | 3.941 | 6.6110 |

## Mode Sensitivity Summary

The final two columns use the illustrative coupling `alpha_q = 0.001` for the most sensitive operator/polarization of each mode.

| mode | largest operator | polarization | max |delta omega/omega| per alpha | Delta f/f [%] | Delta tau/tau [%] |
| --- | --- | --- | ---: | ---: | ---: |
| 220 | epsilon2 | minus | 3.478 | -0.315 | -1.000 |
| 221 | epsilon1 | minus | 7.639 | -0.734 | -0.874 |
| 222 | epsilon1 | minus | 16.144 | -2.093 | 0.392 |
| 330 | epsilon2 | minus | 5.824 | -0.569 | -1.348 |
| 440 | epsilon1 | plus | 9.001 | -0.901 | -0.804 |

## Overtone Amplification

Ratios below compare the magnitude `|delta omega/omega|` with the 220 fundamental at the same spin.

| operator | polarization | 221 / 220 | 222 / 220 |
| --- | --- | ---: | ---: |
| lambda_ev | plus | 1.228 | 1.138 |
| lambda_ev | minus | 2.366 | 4.499 |
| lambda_odd | plus | 1.522 | 1.732 |
| lambda_odd | minus | 1.522 | 1.732 |
| epsilon1 | plus | 0.535 | 0.746 |
| epsilon1 | minus | 3.298 | 6.970 |
| epsilon2 | plus | 1.098 | 1.754 |
| epsilon2 | minus | 2.157 | 3.633 |
| epsilon3 | plus | 0.837 | 1.998 |
| epsilon3 | minus | 0.837 | 1.998 |

## Scientific Reading

- This removes the main spectral mismatch of the earlier bridge: the theory now supplies both 220 and 221, the two modes central to the current GW250114 spectroscopy result.
- Overtones are often substantially more sensitive than the fundamental, but this also means their first-order EFT regime can be narrower.
- At `chi = 0.68`, the relevant corotating modes lie in the paper's intended accuracy range near `chi approximately 0.7`.
- A QNM spectrum is not yet a ringdown likelihood. Mode amplitudes, excitation, start time, detector noise, remnant priors, and polarization content still have to be modeled.
- Couplings must be tested one at a time first. Mixed parity-preserving and parity-breaking terms require the full polarization eigenvalue rule rather than naive addition.

## Next Defensible Step

Build a GW250114 spectral likelihood for the 220+221 content with Kerr remnant mass and spin as nuisance parameters, then map each one-at-a-time EFT coupling into the complex mode frequencies. Start with a synthetic injection/recovery and only then connect to strain or published posterior products.

## Generated Files

- `gw250114_complete_kerr_baseline.csv`
- `gw250114_complete_eft_sensitivities.csv`
- `gw250114_overtone_amplification.csv`
- `gw250114_mode_sensitivity_summary.csv`
- `gw250114_complete_eft_sensitivity_heatmap.png`