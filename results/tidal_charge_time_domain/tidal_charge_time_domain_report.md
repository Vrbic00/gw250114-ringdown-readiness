# Braneworld Tidal-Charge Time-Domain Benchmark

Reference: Toshmatov et al. 2016, arXiv:1605.02058, Table II.

Metric and potential:

```text
f(r) = 1 - 2/r - q/r^2
V_axial = f [l(l+1)/r^2 - 2(3r + 2q)/r^4]
```

The tortoise coordinate is shifted so that the axial potential peak is near `r* = 0` for each `q`.

| q=Q*/M^2 | r_peak | reference Re | TD Re | delta Re [%] | reference Abs(Im) | TD Abs(Im) | delta Abs(Im) [%] | window | residual |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.1 | 3.354611 | 0.367300 | 0.367295 | -0.001 | 0.088300 | 0.088701 | 0.454 | 130-250 | 0.0017 |
| 0.4 | 3.559609 | 0.350800 | 0.350575 | -0.064 | 0.086700 | 0.086867 | 0.192 | 60-180 | 0.0016 |
| 0.7 | 3.745222 | 0.337000 | 0.337155 | 0.046 | 0.085000 | 0.085048 | 0.057 | 140-200 | 0.0011 |
| 1.0 | 3.916084 | 0.325000 | 0.325368 | 0.113 | 0.083500 | 0.083595 | 0.114 | 150-210 | 0.0013 |
| 2.0 | 4.411898 | 0.294400 | 0.295006 | 0.206 | 0.078800 | 0.078757 | -0.054 | 110-230 | 0.0007 |

Interpretation:

- This is the first non-Schwarzschild supplied gravitational-potential benchmark in the static branch.
- The comparison is against the paper's sixth-order WKB table, not an exact Leaver spectrum.
- Agreement at the percent level is sufficient for adopting the case as a validation/stress-test target before moving to Bardeen and Hayward.
