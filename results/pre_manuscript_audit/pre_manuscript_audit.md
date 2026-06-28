# Pre-Manuscript Audit

Overall readiness: `GO_TO_METHODS_DRAFT`.

Audit counts: `7` GO, `0` CAUTION, `0` STOP.

## Computation And Result Checks

| block | status | evidence | caveat | recommended action |
| --- | --- | --- | --- | --- |
| gw250114_public_projection | GO | zero-outside-90 count=0; max sigma from zero=1.222. | Linearized one-at-a-time projection of public posterior products, not a strain-level EFT likelihood. | Safe to state no robust projected beyond-Kerr deviation, with public-product caveats. |
| ringdown_pyring_consistency | GO | max normalized projection difference=0.535 (epsilon1/plus); non-overlap rows=0; zero-fail rows=0. | RINGDOWN and pyRing products should be compared, not statistically combined. | Keep comparison table in main text or early supplement; avoid combined intervals. |
| robustness_checks | GO | filter zero-outside=0, max filter sigma=1.290; empirical zero-outside=0, max empirical nominal sigma=1.286. | Filter and empirical checks are robustness diagnostics, not independent detections. | Keep as support for null result; place detailed rows in supplement. |
| static_physical_deviation_layer | GO | max physical shift=109.80%; max validation-error scale=0.925%; sparse threshold rows=16, interpolated=16. | Threshold crossings are sparse-grid diagnostics and require physical parameter priors before constraint language. | Use to answer whether the static test is discriminating; avoid observational exclusion wording. |
| static_qnm_validation | GO | 25 rows, families=Bardeen_NED, Braneworld_tidal_charge, Hayward, Hayward_overtone, Schwarzschild, max validation delta=0.925%. | Several references are WKB/Prony tables rather than exact Leaver spectra. | Use as readiness/stress-test validation; do not call it final precision spectroscopy. |
| static_readiness_audit | GO | future_theory_specific_project=3, negative_control_metric_only=1, ready_for_next_implementation=1, source_model_caution=1, validated_static_gravitational_qnm=4 | Audit is a classification standard, not a population-complete literature review. | Frame as a reusable referee filter and explicitly invite extension. |
| table_figure_manifest | GO | 12 referenced tables/figures exist; 2 marked diagnostic. | Diagnostic items must be labeled as such in captions. | Use manifest as the draft table/figure control list. |

## Referee Risk Register

| risk | severity | affected section | mitigation | writing rule |
| --- | --- | --- | --- | --- |
| Linearized one-at-a-time EFT couplings miss nonlinear and multi-parameter degeneracies. | high | Methods and limitations | State one-at-a-time linear response explicitly; reserve multi-parameter EFT for outlook. | Do not claim a full theory-space exclusion. |
| Public posterior products are not independent likelihoods. | high | GW250114 projection | Compare RINGDOWN and pyRing side by side; do not combine intervals. | Use 'projection' and 'consistency check', not 'combined constraint'. |
| Static spherical QNM branch is not a rotating GW250114 remnant model. | high | Static readiness audit | Keep static branch as readiness/stress-test audit and negative-control standard. | Never say GW250114 rules out Bardeen/Hayward/tidal-charge static metrics. |
| Manuscript may look like two papers joined together. | medium | Framing | Make the unifying theme 'ringdown-readiness and public reproducibility'. | Use static branch as a methodological audit supporting the main standard. |
| Some static references are WKB/Prony tables, not exact spectra. | medium | Static validation | Report validation against the same published approximation level and label precision limits. | Use 'reproduces published table' rather than 'exact QNM'. |
| Sparse physical-deviation threshold crossings can look like constraints. | medium | Static physical deviations | Call crossings sparse-grid diagnostics and require parameter priors for constraint language. | Use 'threshold crossing' not 'excluded parameter'. |
| Target journal may demand a stronger original physics result. | medium | Submission strategy | Aim PRD only if Methods and audit standard are sharp; keep CQG fallback framing ready. | Emphasize reusable benchmark tables, code, and clear limitations. |

## Go/No-Go Decision

There is no current computational STOP item. The next step can be Methods drafting, provided the draft follows the claim guardrails:

- no combined RINGDOWN+pyRing constraint;
- no strain-level likelihood claim;
- no static-metric observational exclusion;
- sparse static thresholds remain diagnostic;
- static branch is presented as a ringdown-readiness audit.
