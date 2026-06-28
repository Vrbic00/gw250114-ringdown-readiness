# Static Master-Potential Time-Domain Validation

This report evolves Schwarzschild supplied potentials with a characteristic null-grid scheme.
The extracted frequency uses a two-root Prony/linear-prediction fit over several fixed ringdown windows.

Default numerical setup: `h = 0.2 M`, `r*_obs = 20 M`, Gaussian center `r* = 0 M`, width `3 M`.

Scope note: this is an independent validation layer for static supplied potentials. It is not a replacement for Leaver, high-order WKB/Pade, or a full strain likelihood.

| potential | sector | reference Re | TD Re | delta Re [%] | reference Abs(Im) | TD Abs(Im) | delta Abs(Im) [%] | window | residual |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Schwarzschild_scalar_l2 | scalar | 0.483643872211 | 0.483667272272 | 0.005 | 0.096758775978 | 0.096646873538 | -0.116 | 80-180 | 0.0006 |
| Schwarzschild_EM_l2 | electromagnetic | 0.457595511630 | 0.457680605079 | 0.019 | 0.095004425819 | 0.095067222485 | 0.066 | 100-200 | 0.0006 |
| Schwarzschild_ReggeWheeler_l2 | gravitational_odd | 0.373671684418 | 0.373132623641 | -0.144 | 0.088962315689 | 0.089346713235 | 0.432 | 50-150 | 0.0026 |
| Schwarzschild_Zerilli_l2 | gravitational_even | 0.373671684418 | 0.373593229317 | -0.021 | 0.088962315689 | 0.088994025235 | 0.036 | 160-220 | 0.0014 |

Interpretation:

- The Regge-Wheeler and Zerilli rows should agree with the same Schwarzschild gravitational reference.
- The time-domain extraction is expected to be less precise than Leaver/qnm, but it is independent of WKB.
- A case is suitable for the next static-metric stress-test layer only after it passes this benchmark and has a published master potential.
