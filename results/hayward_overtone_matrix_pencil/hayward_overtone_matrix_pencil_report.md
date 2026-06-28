# Hayward First-Overtone Matrix-Pencil Check

Reference: Bolokhov & Skvortsova 2025, arXiv:2508.19989, Table III WKB-8 entries.

Extraction settings: rank-8 matrix pencil, fixed window `60 <= t/M <= 120`, `h=0.2 M`.

| gamma | reference Re | TD Re | delta Re [%] | reference Abs(Im) | TD Abs(Im) | delta Abs(Im) [%] | status |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 0.00 | 0.346003 | 0.344213 | -0.517 | 0.273555 | 0.273679 | 0.045 | PASS:sub_percent_overtone |
| 0.02 | 0.346515 | 0.344699 | -0.524 | 0.272975 | 0.273281 | 0.112 | PASS:sub_percent_overtone |
| 0.04 | 0.347024 | 0.345213 | -0.522 | 0.272300 | 0.272894 | 0.218 | PASS:sub_percent_overtone |
| 0.06 | 0.347609 | 0.345732 | -0.540 | 0.271405 | 0.272481 | 0.397 | PASS:sub_percent_overtone |
| 0.08 | 0.348888 | 0.346260 | -0.753 | 0.270209 | 0.272028 | 0.673 | PASS:sub_percent_overtone |
| 0.10 | 0.350735 | 0.347491 | -0.925 | 0.270779 | 0.270829 | 0.018 | PASS:sub_percent_overtone |
| 1.18 | 0.370157 | 0.369653 | -0.136 | 0.225937 | 0.225743 | -0.086 | PASS:sub_percent_overtone |

Interpretation:

- The first overtone is extracted from the same time-domain waveform as the fundamental mode, but requires a multi-mode matrix-pencil fit.
- The selection rule uses the fundamental-mode damping hierarchy, not the reference overtone frequency itself.
- This is a useful validation of overtone sensitivity, but a final precision overtone solver should still use Leaver, high-order WKB/Pade, or a more controlled multi-mode fit.
