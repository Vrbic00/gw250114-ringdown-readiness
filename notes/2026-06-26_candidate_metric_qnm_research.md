# Candidate Metric/QNM Research

This note records a targeted search for additional metrics or theories that
could be run through the GW250114 public ringdown projection.

## Referee Filter

A candidate is useful for the main ringdown paper only if it supplies:

- gravitational QNM frequencies, not only geodesic, scalar, or shadow proxies;
- rotating black-hole coverage near the GW250114 remnant spin `chi ~ 0.68`;
- complex frequencies or shifts, including damping;
- ideally both `220` and `221`, because the current public-data leverage comes
  from the 220+221 ringdown/spectroscopy products;
- tables, fits, or public code that can be imported reproducibly.

## Main Finding

The filter is severe. Most interesting metrics do not have the perturbation
data needed for a defensible GW250114 ringdown constraint. The currently
implemented Cano et al. higher-derivative Kerr QNM data remain the only
clearly GW250114-ready source found: they provide rotating Kerr QNM correction
fits for `l=2,3,4` and overtones `n=0,1,2`, including the crucial `221` mode.

## Best New Candidate

The best external candidate is dynamical Chern-Simons gravity through the 2025
METRICS calculation. It reaches spins `a <= 0.75`, close to the GW250114
remnant, and provides fitting polynomials for leading-order frequency shifts.
The immediate limitation is mode coverage: the abstract lists `022`, `033`,
and `032`, but not `221`. Therefore it is promising for a future/fundamental
mode extension, but it does not drop cleanly into the present 220+221
pipeline.

## Other Promising but Incomplete Candidates

- EdGB second-order slow-rotation QNMs are physically attractive and the paper
  explicitly targets spins around `~0.7`, but the available results appear
  fundamental-mode focused rather than 221-overtone ready.
- Rapidly rotating EGBd and shift-symmetric EsGB non-perturbative spectra are
  important recent developments, but the public material found is not yet in a
  simple machine-readable 220+221 fit form.
- A follow-up EGBd spectrum-method paper gives valuable details on the
  perturbation PDEs, boundary conditions, and spectral eigenvalue problem. This
  strengthens EGBd as a future-contact/extraction target, but still does not
  supply the current 220+221 overtone table.
- Nonrotating cubic-EFT overtone/excitation work and exotic static regular
  black-hole spectra are useful methodological references for why overtones are
  sensitive to beyond-GR physics. They are not direct GW250114 additions
  because the remnant is rotating.
- Kerr-Newman has analytic overtone models and exact numerical literature, but
  the most overtone-friendly model relies on a modified geodesic
  correspondence. It is better as an appendix/toy comparison than as a main
  Q1-level constraint.

## Rejected for Main Ringdown Constraints

Parametric metrics such as Johannsen/KRZ and rotating regular Bardeen/Hayward
metrics are useful for geodesic/QPO/shadow diagnostics. They should not be used
as ringdown constraints unless the associated gravitational perturbation
equations and QNM spectra are supplied. This is exactly the trap a referee
would object to: changing the metric is not enough to define the ringdown
theory.

## Recommendation

Do not add a weak metric case to the main results just to broaden the paper.
Instead:

1. Keep the main manuscript centered on higher-derivative Kerr QNM
   fingerprints.
2. Add a short "candidate-theory intake" paragraph/table explaining why generic
   metrics are excluded.
3. Optionally add dynamical Chern-Simons as a future-work or appendix target
   after extracting the METRICS fit coefficients.

The most useful next technical step, if we want one more theory, is to extract
the dCS METRICS fit coefficients and test whether a defensible 220-only or
future 330-compatible projection can be included without weakening the
GW250114 220+221 story.
