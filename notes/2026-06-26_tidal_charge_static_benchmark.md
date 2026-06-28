# Tidal-Charge Static QNM Benchmark

This note records the first non-Schwarzschild supplied-potential benchmark in
the static master-potential branch.

Reference:

```text
Toshmatov et al., "Quasinormal frequencies of black hole in the braneworld",
arXiv:1605.02058, Table II.
```

Local run:

```powershell
& 'C:\Users\vrb0015\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' scripts/python/tidal_charge_time_domain_benchmark.py
```

Output:

```text
results/tidal_charge_time_domain/
```

## Model

The static braneworld metric is

```text
f(r) = 1 - 2/r - q/r^2
```

where `q = Q*/M^2 > 0` is the positive tidal-charge parameter used in the
paper. The axial gravitational potential used locally is

```text
V_axial = f [l(l+1)/r^2 - 2(3r + 2q)/r^4]
```

for `M = 1`, `l = 2`.

## Validation Result

The time-domain/Prony extraction reproduces the paper's sixth-order WKB table
very closely:

| `Q*/M^2` | delta Re [%] | delta Abs(Im) [%] |
| ---: | ---: | ---: |
| `0.1` | `-0.001` | `+0.454` |
| `0.4` | `-0.064` | `+0.193` |
| `0.7` | `+0.046` | `+0.056` |
| `1.0` | `+0.113` | `+0.114` |
| `2.0` | `+0.206` | `-0.055` |

This establishes the static branch as more than a Schwarzschild-only
validation. It can reproduce a published non-Schwarzschild gravitational-QNM
benchmark with sub-percent accuracy.

## Interpretation

This is still a static benchmark and not a rotating GW250114 remnant
constraint. Its value for the article is methodological: it shows that the
pipeline can take a published gravitational master potential, evolve it
independently in time domain, and compare it against the literature. That is
exactly the readiness test we want to apply to popular static metrics.
