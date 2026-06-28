# Static QNM Physical-Deviation Report

This report separates two quantities that should not be conflated:

- validation error: how accurately the local solver reproduces a published table;
- physical deviation: how far a metric's QNM spectrum moves away from a Schwarzschild or zero-parameter baseline.

The physical-deviation numbers below use the published/reference frequencies where available. The validation deltas are retained only as a numerical-error diagnostic.

## Largest Sampled Deviations

| family / mode | largest sampled parameter | delta Re [%] | delta Abs(Im) [%] | max abs delta [%] | validation-error scale [%] |
| --- | ---: | ---: | ---: | ---: | ---: |
| Bardeen_NED n=0 | alpha=0.7698 | 109.800 | -14.565 | 109.800 | 0.683 |
| Braneworld_tidal_charge n=0 | q_tidal=2 | -21.214 | -11.423 | 21.214 | 0.206 |
| Hayward_fundamental n=0 | gamma=1.18 | 6.035 | -16.596 | 16.596 | 0.092 |
| Hayward_overtone n=1 | gamma=1.18 | 6.981 | -17.407 | 17.407 | 0.136 |

## Sampled Rows

| family | parameter | mode | baseline | Re shift [%] | Abs(Im) shift [%] | max abs [%] |
| --- | ---: | --- | --- | ---: | ---: | ---: |
| Braneworld_tidal_charge | q_tidal=0 | l=2,n=0 | Schwarzschild_qnm_l2_n0 | 0.000 | 0.000 | 0.000 |
| Braneworld_tidal_charge | q_tidal=0.1 | l=2,n=0 | Schwarzschild_qnm_l2_n0 | -1.705 | -0.744 | 1.705 |
| Braneworld_tidal_charge | q_tidal=0.4 | l=2,n=0 | Schwarzschild_qnm_l2_n0 | -6.121 | -2.543 | 6.121 |
| Braneworld_tidal_charge | q_tidal=0.7 | l=2,n=0 | Schwarzschild_qnm_l2_n0 | -9.814 | -4.454 | 9.814 |
| Braneworld_tidal_charge | q_tidal=1 | l=2,n=0 | Schwarzschild_qnm_l2_n0 | -13.025 | -6.140 | 13.025 |
| Braneworld_tidal_charge | q_tidal=2 | l=2,n=0 | Schwarzschild_qnm_l2_n0 | -21.214 | -11.423 | 21.214 |
| Bardeen_NED | alpha=0 | l=2,n=0 | Bardeen_alpha_0_same_table | 0.000 | 0.000 | 0.000 |
| Bardeen_NED | alpha=0.3 | l=2,n=0 | Bardeen_alpha_0_same_table | 8.971 | -1.387 | 8.971 |
| Bardeen_NED | alpha=0.7698 | l=2,n=0 | Bardeen_alpha_0_same_table | 109.800 | -14.565 | 109.800 |
| Hayward_fundamental | gamma=0 | l=2,n=0 | Schwarzschild_qnm_l2_n0 | 0.000 | 0.000 | 0.000 |
| Hayward_fundamental | gamma=0.1 | l=2,n=0 | Schwarzschild_qnm_l2_n0 | 0.440 | -0.831 | 0.831 |
| Hayward_fundamental | gamma=0.35 | l=2,n=0 | Schwarzschild_qnm_l2_n0 | 1.599 | -3.201 | 3.201 |
| Hayward_fundamental | gamma=0.6 | l=2,n=0 | Schwarzschild_qnm_l2_n0 | 2.849 | -6.121 | 6.121 |
| Hayward_fundamental | gamma=0.85 | l=2,n=0 | Schwarzschild_qnm_l2_n0 | 4.192 | -9.834 | 9.834 |
| Hayward_fundamental | gamma=1.1 | l=2,n=0 | Schwarzschild_qnm_l2_n0 | 5.592 | -14.718 | 14.718 |
| Hayward_fundamental | gamma=1.18 | l=2,n=0 | Schwarzschild_qnm_l2_n0 | 6.035 | -16.596 | 16.596 |
| Hayward_overtone | gamma=0 | l=2,n=1 | Hayward_gamma_0_same_table | 0.000 | 0.000 | 0.000 |
| Hayward_overtone | gamma=0.02 | l=2,n=1 | Hayward_gamma_0_same_table | 0.148 | -0.212 | 0.212 |
| Hayward_overtone | gamma=0.04 | l=2,n=1 | Hayward_gamma_0_same_table | 0.295 | -0.459 | 0.459 |
| Hayward_overtone | gamma=0.06 | l=2,n=1 | Hayward_gamma_0_same_table | 0.464 | -0.786 | 0.786 |
| Hayward_overtone | gamma=0.08 | l=2,n=1 | Hayward_gamma_0_same_table | 0.834 | -1.223 | 1.223 |
| Hayward_overtone | gamma=0.1 | l=2,n=1 | Hayward_gamma_0_same_table | 1.368 | -1.015 | 1.368 |
| Hayward_overtone | gamma=1.18 | l=2,n=1 | Hayward_gamma_0_same_table | 6.981 | -17.407 | 17.407 |

## Sparse Threshold Crossings

These crossings are linear interpolations over sparse published tables, not final exclusion limits.

| family | mode | parameter | threshold [%] | crossing | status |
| --- | --- | --- | ---: | ---: | --- |
| Braneworld_tidal_charge | n=0 | q_tidal | 1 | 0.05865 | linear_interpolation_sparse_grid |
| Braneworld_tidal_charge | n=0 | q_tidal | 3 | 0.188 | linear_interpolation_sparse_grid |
| Braneworld_tidal_charge | n=0 | q_tidal | 5 | 0.3239 | linear_interpolation_sparse_grid |
| Braneworld_tidal_charge | n=0 | q_tidal | 10 | 0.7174 | linear_interpolation_sparse_grid |
| Bardeen_NED | n=0 | alpha | 1 | 0.03344 | linear_interpolation_sparse_grid |
| Bardeen_NED | n=0 | alpha | 3 | 0.1003 | linear_interpolation_sparse_grid |
| Bardeen_NED | n=0 | alpha | 5 | 0.1672 | linear_interpolation_sparse_grid |
| Bardeen_NED | n=0 | alpha | 10 | 0.3048 | linear_interpolation_sparse_grid |
| Hayward_fundamental | n=0 | gamma | 1 | 0.1178 | linear_interpolation_sparse_grid |
| Hayward_fundamental | n=0 | gamma | 3 | 0.3288 | linear_interpolation_sparse_grid |
| Hayward_fundamental | n=0 | gamma | 5 | 0.504 | linear_interpolation_sparse_grid |
| Hayward_fundamental | n=0 | gamma | 10 | 0.8585 | linear_interpolation_sparse_grid |
| Hayward_overtone | n=1 | gamma | 1 | 0.06979 | linear_interpolation_sparse_grid |
| Hayward_overtone | n=1 | gamma | 3 | 0.2099 | linear_interpolation_sparse_grid |
| Hayward_overtone | n=1 | gamma | 5 | 0.3446 | linear_interpolation_sparse_grid |
| Hayward_overtone | n=1 | gamma | 10 | 0.6813 | linear_interpolation_sparse_grid |

## Interpretation

- Sub-percent validation does not mean the metric is physically close to Schwarzschild.
- In the current static benchmarks, physical QNM shifts range from the percent level to tens of percent across the sampled parameter domains.
- This branch can therefore be discriminating, but only after a physically meaningful parameter range and an observational or synthetic tolerance are declared.
- These static deviations remain a readiness/stress-test layer; they are not direct GW250114 rotating-remnant exclusions.
