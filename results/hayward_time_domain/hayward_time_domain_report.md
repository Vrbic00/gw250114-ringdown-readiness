# Hayward Axial Gravitational Time-Domain Benchmark

Reference: Bolokhov & Skvortsova 2025, arXiv:2508.19989, Table V.

Metric and potential:

```text
f(r) = 1 - 2 r^2 / (r^3 + gamma)
V = f [2 f/r^2 - f'/r + ((ell+2)(ell-1))/r^2]
```

| gamma | r_h | r_peak | reference Re | TD Re | delta Re [%] | reference Abs(Im) | TD Abs(Im) | delta Abs(Im) [%] | window | residual |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.10 | 1.974346 | 3.255110 | 0.375316 | 0.375531 | 0.057 | 0.088223 | 0.087590 | -0.718 | 120-260 | 0.0033 |
| 0.35 | 1.903392 | 3.186710 | 0.379646 | 0.379605 | -0.011 | 0.086115 | 0.086259 | 0.167 | 110-250 | 0.0009 |
| 0.60 | 1.818579 | 3.110894 | 0.384318 | 0.384306 | -0.003 | 0.083517 | 0.083406 | -0.132 | 130-270 | 0.0010 |
| 0.85 | 1.708958 | 3.025137 | 0.389335 | 0.389307 | -0.007 | 0.080214 | 0.080311 | 0.120 | 140-280 | 0.0008 |
| 1.10 | 1.530247 | 2.925109 | 0.394567 | 0.394516 | -0.013 | 0.075869 | 0.076113 | 0.322 | 170-310 | 0.0017 |
| 1.18 | 1.383623 | 2.888891 | 0.396224 | 0.396280 | 0.014 | 0.074198 | 0.074130 | -0.092 | 130-270 | 0.0004 |

Interpretation:

- This benchmark reproduces the published fundamental-mode time-domain/Prony values for the Hayward axial gravitational potential.
- The first overtone is intentionally left for a separate multi-exponential extraction layer.
- The result is a static readiness/stress test, not a rotating GW250114 remnant constraint.
