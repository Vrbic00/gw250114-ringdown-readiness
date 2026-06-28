# Hayward Static QNM Benchmark

This note records the second regular-black-hole target in the static
master-potential branch.

Reference:

```text
Bolokhov & Skvortsova, "Gravitational quasinormal modes of the Hayward
spacetime", arXiv:2508.19989.
```

Local run:

```powershell
& 'C:\Users\vrb0015\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' scripts/python/hayward_time_domain_benchmark.py
```

Output:

```text
results/hayward_time_domain/
```

## Model

The Hayward/asymptotic-safety-inspired metric is

```text
f(r) = 1 - 2 r^2 / (r^3 + gamma)
```

for `M = 1`. The axial gravitational potential implemented locally is

```text
V = f [2 f/r^2 - f'/r + ((ell+2)(ell-1))/r^2]
```

This is the potential given explicitly in the source paper under the
anisotropic-fluid axial-sector assumption.

## Validation Result

The local time-domain/Prony extraction reproduces the paper's Table V Prony
fits for `ell=2,n=0`:

| `gamma` | delta Re [%] | delta Abs(Im) [%] |
| ---: | ---: | ---: |
| `0.10` | `+0.057` | `-0.718` |
| `0.35` | `-0.011` | `+0.168` |
| `0.60` | `-0.003` | `-0.133` |
| `0.85` | `-0.007` | `+0.121` |
| `1.10` | `-0.013` | `+0.321` |
| `1.18` | `+0.014` | `-0.092` |

All tested deviations are below about `0.8%`.

## Interpretation

Hayward is now a successful second regular-black-hole stress test. It is
particularly useful for the article because the source paper emphasizes that
the first overtone is more sensitive to the deformation parameter. The current
local implementation validates the fundamental mode only; a separate
multi-exponential extraction layer is needed before using the overtone table
quantitatively.
