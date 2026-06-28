# Hayward First-Overtone Matrix-Pencil Check

Reference: Bolokhov & Skvortsova (2025), arXiv:2508.19989, Table III WKB-8
entries for the Hayward axial gravitational `ell=2,n=1` mode.

## Purpose

The fundamental Hayward time-domain benchmark validates the supplied axial
gravitational potential at the `n=0` level. The overtone table is a stronger
stress test because the extracted waveform is multi-mode and a single Prony fit
is not stable enough to isolate `n=1`.

This block therefore adds a matrix-pencil extraction layer. It is intended as a
readiness and sensitivity check, not as a final precision QNM solver.

## Run Command

```powershell
& 'C:\Users\vrb0015\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' scripts/python/hayward_overtone_matrix_pencil.py
```

Output files:

```text
results/hayward_overtone_matrix_pencil/hayward_overtone_matrix_pencil_summary.csv
results/hayward_overtone_matrix_pencil/hayward_overtone_matrix_pencil_modes.csv
results/hayward_overtone_matrix_pencil/hayward_overtone_matrix_pencil_report.md
```

## Extraction Settings

- Same Hayward axial gravitational potential and null-grid time-domain
  waveform as the fundamental-mode benchmark.
- Rank-8 matrix-pencil fit.
- Fixed fitting window: `60 <= t/M <= 120`.
- Mode selection uses the expected hierarchy relative to the already extracted
  fundamental mode: lower real frequency and substantially stronger damping.
  The overtone reference value is not used as the selection target.

## Results

| gamma | reference Re | TD Re | delta Re [%] | reference Abs(Im) | TD Abs(Im) | delta Abs(Im) [%] |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.00 | 0.346003 | 0.344213 | -0.517 | 0.273555 | 0.273679 | 0.045 |
| 0.02 | 0.346515 | 0.344699 | -0.524 | 0.272975 | 0.273281 | 0.112 |
| 0.04 | 0.347024 | 0.345213 | -0.522 | 0.272300 | 0.272894 | 0.218 |
| 0.06 | 0.347609 | 0.345732 | -0.540 | 0.271405 | 0.272481 | 0.397 |
| 0.08 | 0.348888 | 0.346260 | -0.753 | 0.270209 | 0.272028 | 0.673 |
| 0.10 | 0.350735 | 0.347491 | -0.925 | 0.270779 | 0.270829 | 0.018 |
| 1.18 | 0.370157 | 0.369653 | -0.136 | 0.225937 | 0.225743 | -0.086 |

All tested entries remain within the configured sub-percent criterion.

## Interpretation

This result makes Hayward more valuable than a fundamental-only static example:
it demonstrates that the local pipeline can recover an overtone-sensitive
spectral feature from the supplied master potential.

The limitation is important. Matrix-pencil extraction is a controlled
time-domain diagnostic, but it is not a replacement for a final Leaver solver,
high-order WKB/Pade calculation, or a carefully benchmarked multi-mode
spectral method. In the paper this should be phrased as an overtone readiness
check, not as final high-precision spectroscopy.
