# Bardeen Static QNM Benchmark

This note records the first popular regular-black-hole target in the static
master-potential branch.

Reference:

```text
Ulhoa, "On Quasinormal Modes for Gravitational Perturbations of Bardeen
Black Hole", arXiv:1303.3143.
```

Local run:

```powershell
& 'C:\Users\vrb0015\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' scripts/python/bardeen_time_domain_benchmark.py
```

Output:

```text
results/bardeen_time_domain/
```

## Model

The Bardeen metric is

```text
f(r) = 1 - 2 r^2 / (r^2 + alpha^2)^(3/2)
```

for `m = 1`. The axial gravitational potential implemented locally is the
source-dependent nonlinear-electrodynamics potential from the paper:

```text
V = f [(l(l+1)+2(f-1))/r^2 + f'/r + f'' + 2 k L]
k = 8 pi
L = 3 alpha^2 / (r^2 + alpha^2)^(5/2)
```

This is important: it is not a metric-only Regge-Wheeler proxy. It uses a
specific matter/source interpretation of the Bardeen spacetime.

## Validation Result

The local time-domain/Prony extraction reproduces Ulhoa's third-order WKB
tables for `l=2,n=0`:

| `alpha` | delta Re [%] | delta Abs(Im) [%] |
| ---: | ---: | ---: |
| `0` | `+0.174` | `-0.284` |
| `0.3` | `-0.109` | `-0.873` |
| `4/sqrt(27)` | `+0.022` | `+0.682` |

The agreement is comfortably within the expected accuracy of comparing a
time-domain extraction against a third-order WKB table.

## Interpretation

Bardeen is now a successful high-value static regular-black-hole stress test.
It supports the central methodological claim:

```text
A popular metric can be admitted into the static gravitational-ringdown
branch only after a source-dependent master potential is supplied and
validated against published QNM data.
```

This is not a rotating GW250114 constraint, but it is a strong readiness audit
for static metrics used in QPO/shadow phenomenology.
