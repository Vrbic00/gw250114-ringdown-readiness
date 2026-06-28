# Static Time-Domain QNM Validation

The static master-potential branch now has an independent time-domain
validation script:

```powershell
& 'C:\Users\vrb0015\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' scripts/python/static_master_potential_time_domain.py
```

Output directory:

```text
results/static_master_potential_time_domain/
```

## Method

The script evolves supplied Schwarzschild master potentials on a null grid,
using the characteristic finite-difference stencil for

```text
d^2 psi/dr_*^2 - d^2 psi/dt^2 = V(r) psi
```

The QNM frequency is extracted with a two-root Prony/linear-prediction fit
over multiple fixed ringdown windows. The reference values come from the
local Python `qnm` package.

Default numerical setup:

- null-grid spacing `h = 0.2 M`;
- observer position `r*_obs = 20 M`;
- initial Gaussian center `r* = 0 M`;
- initial Gaussian width `3 M`.

## Validation Result

The clean validation run gives:

| potential | delta Re [%] | delta Abs(Im) [%] |
| --- | ---: | ---: |
| Schwarzschild scalar `l=2,n=0` | `+0.005` | `-0.116` |
| Schwarzschild EM `l=2,n=0` | `+0.019` | `+0.066` |
| Schwarzschild Regge-Wheeler `l=2,n=0` | `-0.144` | `+0.432` |
| Schwarzschild Zerilli `l=2,n=0` | `-0.021` | `+0.036` |

This is good enough for the next static-QNM stress-test layer. It is still
not a precision Leaver or high-order WKB/Pade solver; its role is an
independent validation and reproducible cross-check for supplied potentials.

## Next Use

The next safe implementation step is a simple non-Schwarzschild supplied
potential with published QNM tables:

1. Reissner-Nordstrom or braneworld tidal charge as the first charged/RN-like
   benchmark.
2. Bardeen axial gravitational perturbations as the first high-value
   regular-black-hole metric.
3. Hayward axial gravitational perturbations as the overtone-sensitive regular
   black-hole target.
