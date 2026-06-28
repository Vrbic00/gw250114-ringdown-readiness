# Static QNM Readiness Scorecard

This table joins the local time-domain validations for supplied static master potentials.

| family | rows | max abs delta Re [%] | max abs delta Abs(Im) [%] | verdicts |
| --- | ---: | ---: | ---: | --- |
| Bardeen_NED | 3 | 0.174 | 0.873 | PASS:sub_percent_validation |
| Braneworld_tidal_charge | 5 | 0.206 | 0.454 | PASS:sub_percent_validation |
| Hayward | 6 | 0.057 | 0.718 | PASS:sub_percent_validation |
| Hayward_overtone | 7 | 0.925 | 0.673 | PASS:sub_percent_validation |
| Schwarzschild | 4 | 0.144 | 0.432 | PASS:sub_percent_validation |

Interpretation:

- `PASS:sub_percent_validation` means the local time-domain extraction reproduces the reference table to within 1 percent in both real frequency and damping rate.
- These are static supplied-potential validations. They do not replace rotating-remnant GW250114 constraints.
- The Hayward overtone rows use a multi-mode matrix-pencil extraction and should be read as an overtone-sensitivity check, not as a final precision spectroscopy solver.
- The scorecard is intended as the paper-facing readiness layer for metrics often used in geodesic, QPO, and shadow phenomenology.
