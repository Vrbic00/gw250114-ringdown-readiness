# Static Metric QNM Readiness Audit

This report converts the candidate registry and validation scorecard into a paper-facing referee filter.

Core rule: a static line element, geodesic observable, QPO fit, or shadow calculation is not enough for a gravitational ringdown claim. A candidate needs at least a supplied gravitational master potential plus a reproduced QNM spectrum.

## Status Counts

| audit status | count |
| --- | ---: |
| future_theory_specific_project | 3 |
| negative_control_metric_only | 1 |
| ready_for_next_implementation | 1 |
| source_model_caution | 1 |
| validated_static_gravitational_qnm | 4 |

## Validated Families

| family | rows | max abs delta Re [%] | max abs delta Abs(Im) [%] | modes |
| --- | ---: | ---: | ---: | --- |
| Bardeen_NED | 3 | 0.174 | 0.873 | l=2,n=0 |
| Braneworld_tidal_charge | 5 | 0.206 | 0.454 | l=2,n=0 |
| Hayward | 6 | 0.057 | 0.718 | l=2,n=0 |
| Hayward_overtone | 7 | 0.925 | 0.673 | l=2,n=1 |
| Schwarzschild | 4 | 0.144 | 0.432 | l=2,n=0 |

## Candidate Audit

| candidate | readiness | audit status | validation families | referee action |
| --- | --- | --- | --- | --- |
| schwarzschild_rw_zerilli | A+ | validated_static_gravitational_qnm | Schwarzschild | May be used as a static supplied-potential benchmark, with the non-rotating limitation explicit. |
| reissner_nordstrom_gravito_em | A- | ready_for_next_implementation | none | Good candidate for the next reproducible benchmark before any observational language. |
| braneworld_tidal_charge_rn_like | A- | validated_static_gravitational_qnm | Braneworld_tidal_charge | May be used as a static supplied-potential benchmark, with the non-rotating limitation explicit. |
| bardeen_ned_axial | A- | validated_static_gravitational_qnm | Bardeen_NED | May be used as a static supplied-potential benchmark, with the non-rotating limitation explicit. |
| hayward_axial_gravitational | A- | validated_static_gravitational_qnm | Hayward; Hayward_overtone | May be used as a static supplied-potential benchmark, with the non-rotating limitation explicit. |
| simpson_visser_black_bounce_source_ambiguity | B- | source_model_caution | none | Useful as a conceptual caution: the same metric can imply different perturbation physics. |
| scalarized_einstein_gauss_bonnet_static | C+ | future_theory_specific_project | none | Do not fold into the current paper unless the full perturbation setup is reproduced. |
| generalized_einstein_maxwell_scalar | C+ | future_theory_specific_project | none | Do not fold into the current paper unless the full perturbation setup is reproduced. |
| abg_stvg_odd_gravitational | C | future_theory_specific_project | none | Do not fold into the current paper unless the full perturbation setup is reproduced. |
| regular_metric_scalar_only_cases | D | negative_control_metric_only | none | Do not promote to gravitational ringdown constraints; cite as a readiness failure or outlook item. |

## Manuscript Use

- Main text: use the validated families as positive examples of a static supplied-potential readiness audit.
- Appendix or outlook: list ready-but-not-yet-implemented candidates such as RN/gravito-electromagnetic perturbations.
- Negative-control paragraph: metric-only regular black holes, shadow/QPO-only models, and source-ambiguous examples cannot be promoted to gravitational ringdown constraints without perturbation physics.
- GW250114 connection: keep this branch separate from rotating-remnant constraints; its role is community-facing model triage.

Figures generated with this audit:

- `static_qnm_validation_scorecard.svg`
- `static_metric_readiness_ladder.svg`
