# Bardeen Axial Gravitational Time-Domain Benchmark

Reference: Ulhoa 2013, arXiv:1303.3143, Tables I-III.

Metric and potential:

```text
f(r) = 1 - 2 r^2 / (r^2 + alpha^2)^(3/2)
V = f [(l(l+1)+2(f-1))/r^2 + f'/r + f'' + 2 k L]
k = 8 pi,  L = 3 alpha^2 / (r^2 + alpha^2)^(5/2)
```

| alpha | r_h | r_peak | reference Re | TD Re | delta Re [%] | reference Abs(Im) | TD Abs(Im) | delta Abs(Im) [%] | window | residual |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.000000000 | 2.000000 | 3.280776 | 0.373162 | 0.373812 | 0.174 | 0.089217 | 0.088964 | -0.284 | 120-200 | 0.0009 |
| 0.300000000 | 1.929617 | 2.942084 | 0.406640 | 0.406196 | -0.109 | 0.087980 | 0.087212 | -0.873 | 70-190 | 0.0004 |
| 0.769800359 | 1.088662 | 1.761271 | 0.782892 | 0.783066 | 0.022 | 0.076223 | 0.076743 | 0.683 | 70-190 | 0.0039 |

Interpretation:

- This is a source-dependent nonlinear-electrodynamics Bardeen axial gravitational potential, not a metric-only proxy.
- The comparison is against the paper's third-order WKB tables, so sub-percent agreement is not expected in every case.
- The extremal case is numerically delicate because the horizon is degenerate and the tortoise coordinate is more singular.
