# Static Master-Potential/QNM Feasibility Audit

This note answers whether the project can safely add a static spherical
master-potential/QNM module without discovering halfway through that no
defensible examples exist.

## Bottom Line

Yes, the direction is feasible if the scope is kept static, axial-first, and
validation-first.

The project should not try to infer gravitational ringdown from an arbitrary
metric alone. It should accept a case into the stronger branch only when the
literature supplies a one-dimensional master equation, a clear source/theory
interpretation, and published QNM values for validation.

## What We Can Validate

The current local code already has first-order and Iyer-Will third-order WKB
for supplied one-dimensional potentials. It is already validated at the
Schwarzschild level: WKB3 reproduces the real part of the known `l=2,n=0`
Schwarzschild scalar, electromagnetic, and Regge-Wheeler frequencies at roughly
the sub-percent level, while WKB1 is much worse.

For a publishable static branch, the next module should add:

1. a higher-order WKB/Pade or time-domain solver for supplied `V_l(r)`;
2. automated horizon and potential-peak checks;
3. a reference-table comparator;
4. a readiness verdict that separates metric-only, test-field, supplied
   gravitational potential, and full-theory cases.

## Green Candidates

### Schwarzschild Regge-Wheeler/Zerilli

This is the non-negotiable anchor. It is already partly implemented locally and
has abundant reference frequencies. It validates conventions, signs, tortoise
coordinate, boundary conditions, and WKB/time-domain extraction.

### Reissner-Nordstrom / gravito-electromagnetic modes

This is the clean charged GR benchmark. The main technical caution is that the
physical perturbations are coupled gravito-electromagnetic modes, so we should
not pretend it is a single pure-gravity deformation of the Regge-Wheeler
potential. Still, it is a very good validation case before moving to regular
black holes.

### Braneworld tidal-charge RN-like metric

This may be the best first non-Schwarzschild target. The metric is simple,
static, and spherically symmetric, and the literature gives scalar,
electromagnetic, axial gravitational, and polar gravitational perturbations
with published WKB results. It is also close to the kind of phenomenological
static metric literature the project wants to audit.

The caveat is interpretive: the paper uses a simplifying assumption about the
projected bulk Weyl perturbation. That is acceptable if clearly labeled.

## Strong Regular-BH Targets

### Bardeen

Bardeen is strategically important because it is heavily reused in
QPO/shadow/geodesic papers. There is an asymptotically flat axial
gravitational QNM paper with a master equation and WKB tables, plus a more
systematic nonlinear-electrodynamics perturbation framework and later
Bardeen-(A)dS axial/polar work.

This is probably the first "popular metric" we should attack after the
benchmark cases.

The caveat is source dependence: Bardeen can be interpreted through nonlinear
electrodynamics or through an effective anisotropic source, and the
perturbation equations depend on that choice. This is not a problem; it is
actually part of the scientific point. A metric is not enough.

### Hayward

Hayward is also feasible. A recent gravitational-QNM paper reports axial
perturbations, WKB-Pade, time-domain evolution, and first-overtone behavior.
This makes it attractive because overtone sensitivity connects naturally to
the project's ringdown theme.

The caveat is that it is recent and must be extracted carefully. It should be
second after Bardeen, not the first implementation target.

## Useful But Not First

### Simpson-Visser / black-bounce

This is conceptually excellent because it shows that the same metric can lead
to different ringdown behavior depending on the matter-source interpretation.
That supports the project's critique of metric-only phenomenology. It is not
the cleanest first numerical target.

### Scalarized Einstein-Gauss-Bonnet

This is physically stronger than regular-metric phenomenology, but it requires
numerical background functions and careful hyperbolicity/stability domains. It
is a follow-up project, not the first static module.

### Generalized Einstein-Maxwell-scalar frameworks

These are valuable as theory umbrellas and may provide notebooks or general
master equations, but they are not a single ready benchmark. They become useful
once a specific model is selected.

## No-Go For Gravitational Claims

Cases with only scalar-field QNMs, electromagnetic test-field QNMs, shadows,
photon spheres, or QPO/geodesic observables should remain metric-only or
test-field examples. They are useful as negative controls in the readiness
audit, but they should not be presented as gravitational ringdown constraints.

## Recommended Implementation Order

1. Strengthen the supplied-potential solver:
   - keep WKB3;
   - add time-domain evolution and Prony/ringdown extraction;
   - optionally add higher-order WKB/Pade later.
2. Validate on Schwarzschild Regge-Wheeler/Zerilli.
3. Add Reissner-Nordstrom or braneworld tidal charge as the first
   non-Schwarzschild benchmark.
4. Add Bardeen axial gravitational perturbations as the first high-value
   regular-BH target.
5. Add Hayward axial gravitational perturbations as the overtone-sensitive
   exotic target.
6. Turn scalar-only/QPO/shadow-only metrics into a public "not ringdown-ready"
   scorecard.

## Referee-Safe Framing

The strong claim is not:

```text
We rule out metric X.
```

The defensible claim is:

```text
We provide a reproducible readiness and static-QNM stress test. Metrics with
published gravitational master potentials can be quantitatively benchmarked;
metrics with only geodesic, shadow, QPO, scalar-field, or electromagnetic
test-field calculations are not yet gravitational-ringdown-ready.
```

That is sharp enough to be useful to referees, while still being technically
defensible.
