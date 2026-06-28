<!-- Near-final manuscript for expert feedback before journal submission. -->

Working author list: J. Vrba
# A Ringdown-Readiness Framework for Beyond-Kerr Claims: Public GW250114 Projections and Static QNM Benchmarks

## Abstract

In the present work we formulate a reproducible ringdown-readiness framework
for beyond-Kerr claims. The motivation is simple: a phenomenological black-hole
metric is not ready for gravitational-ringdown use only because it has a
horizon, photon ring, ISCO, QPO fit, or shadow. A ringdown test also requires
the perturbation problem, boundary conditions, and reproducible quasinormal
mode spectra. We illustrate the framework in two complementary ways. First, we
use public GW250114 spectroscopy products and published higher-derivative Kerr
QNM fingerprints to make a transparent public-data projection. In this worked
example the Kerr value `alpha = 0` remains inside every one-at-a-time 90
percent interval, with the largest nominal displacement equal to `0.497 sigma`
for the RINGDOWN branch and `1.222 sigma` for the pyRing branch. The two public
branches are mutually consistent, with maximum normalized projection
difference `0.535 sigma`. This null result is not advertised as an evidence
ratio or as a full beyond-Kerr reanalysis; it is a reproducible benchmark for
what can be inferred from public ringdown products alone. Second, we validate
`25` static supplied-potential rows for Schwarzschild, braneworld tidal charge,
Bardeen, and Hayward examples as static QNM benchmarks. These static benchmarks
are not GW250114 constraints on rotating remnants. They demonstrate, however,
where metric-only phenomenology ends and where gravitational-ringdown readiness
begins.

## 1. Introduction

Black-hole ringdown spectroscopy is becoming an important test of the Kerr
description of compact-object remnants \cite{Berti2009QNMReview,Isi2019NoHair,Gho-Bri-Buo:2021:PRD:}.
In general relativity, the late post-merger signal is described by a discrete
set of quasinormal modes. Their frequencies and damping times are fixed by the
remnant mass and spin. When the signal-to-noise ratio is high enough, measured
ringdown quantities can therefore be compared with the Kerr pattern and used
as a test of the strong-field regime.

This program is powerful, but it also imposes a strict theoretical standard.
A metric alone is not a gravitational ringdown model. A line element may be
useful for geodesic phenomenology, QPO modelling, ray tracing, or shadow
studies \cite{Johannsen2015Metric,KRZ2016Parametrization,Pedrotti2024Eikonal},
but gravitational ringdown additionally requires a perturbation problem, the
relevant propagating degrees of freedom, boundary conditions, and reproducible
QNM spectra. Without these ingredients it is not clear what quantity should be
compared with the observed ringdown frequencies and damping times.

The aim of this paper is to make this distinction operational. We propose a
ringdown-readiness framework for beyond-Kerr claims and illustrate it with two
controlled examples. The first example is observationally motivated: we use
public GW250114 RINGDOWN and pyRing products
\cite{LVK_GW250114_Spectroscopy,GW250114_ZenodoRelease} to project published
higher-derivative Kerr QNM fingerprints \cite{Cano2024BeyondKerrQNM,Cano2023RotatingQNM,BeyondKerrQNMRepo}
onto one-at-a-time effective couplings. This is not a full strain-level
beyond-Kerr inference. It is a reproducible public-data projection showing how
theory-backed rotating QNM fingerprints can be confronted with public
spectroscopy products.

The second example is theory-readiness motivated. We validate supplied
gravitational master-potential examples for static spherical backgrounds,
including Schwarzschild, braneworld tidal charge, Bardeen, and Hayward cases
\cite{Berti2009QNMReview,Toshmatov2016TidalCharge,Ulhoa2013Bardeen,Bolokhov2025Hayward}.
This part is not a direct constraint on the GW250114 remnant, which is a
rotating object. Its role is different: it demonstrates a reproducible
benchmark level between metric-only phenomenology and a theory-backed
rotating-ringdown calculation.

The central message is therefore methodological. Public GW250114 products can
be used for transparent sanity checks when a model supplies rotating QNM
fingerprints. Static supplied-potential benchmarks can be used to test whether
phenomenological metrics have reached at least a reproducible perturbative
level. Metric-only models remain below the threshold required for
observational gravitational-ringdown claims.

This framing is meant to be useful in practice. It does not forbid the use of
phenomenological metrics in QPO, shadow, plasma, or geodesic calculations. It
only says that, when the same metric is promoted to a ringdown statement, the
paper must also supply the perturbation system or an equivalent validated QNM
fingerprint. The readiness ladder is therefore a referee-facing diagnostic:
it tells the reader what has actually been tested and what has merely been
assumed.

## 2. Methods

### 2.1 Readiness Logic

We use the term ringdown-ready in a narrow sense. A model is not ringdown-ready
only because it specifies a metric. It must also specify enough perturbation
physics to compute gravitational QNM spectra. In the present paper we use the
following ladder:

```text
metric/geodesic only
test-field QNM
supplied gravitational master potential
validated static gravitational QNM
theory-backed rotating gravitational QNM.
```

The public GW250114 projection uses the highest level available in this work:
published rotating Kerr-like QNM fingerprints for higher-derivative
corrections. The static branch reaches the validated static gravitational-QNM
level for selected examples. It does not claim to replace a theory-backed
rotating calculation.

### 2.2 Public GW250114 Inputs

We use public GW250114 ringdown/spectroscopy products in two separate
projection branches. The RINGDOWN branch uses the posterior samples for

```text
y_R = {log f_220, log f_221, df_221}.
```

The pyRing branch uses the 221 fractional deviation posterior

```text
y_P = {log(1 + domega_221), log(1 + dtau_221)}.
```

These two branches are treated as complementary public-product projections of
the same event. We compare them for consistency, but we do not combine them
statistically because they are not independent likelihoods.

### 2.3 Kerr Baseline

For each mode, the Kerr reference is represented by the complex dimensionless
frequency \cite{Berti2006Spectroscopy,Stein2019QNM}

```text
M omega_lmn(chi) = M omega_R,lmn(chi) - i M/tau_lmn(chi).
```

The physical frequency and damping time scale as

```text
f_lmn(M, chi) propto Re[M omega_lmn(chi)] / M,
tau_lmn(M, chi) propto M / Abs(Im[M omega_lmn(chi)]).
```

The RINGDOWN projection uses the Kerr baseline for the event remnant mass and
spin, and in this branch the linearized remnant-mass and remnant-spin
directions are profiled over. The pyRing variables are already fractional
deviations from the Kerr prediction, and therefore the baseline in the pyRing
branch is zero.

### 2.4 Higher-Derivative QNM Fingerprints

For each one-at-a-time EFT operator and polarization branch, the imported QNM
fingerprints give linear responses such as

```text
d log f_lmn / d alpha,
d log tau_lmn / d alpha.
```

The projected coupling `alpha` is therefore a local linear coordinate along a
published higher-derivative QNM fingerprint. We do not interpret it here as a
full nonlinear or multi-parameter exclusion of theory space.

### 2.5 Gaussian Projection

For a public posterior sample vector `y`, we calculate the sample mean `mu` and
covariance `C`. The local Gaussian chi-square is

```text
chi2(p) = (mu - y_model(p))^T C^-1 (mu - y_model(p)).
```

For the RINGDOWN branch, the linearized model residual is

```text
r = mu - y_Kerr = N theta + s alpha,
```

where `theta = {delta log M_f, delta chi_f}` contains nuisance directions and
`s` is the EFT fingerprint column. The mass nuisance direction contributes
`{-1, -1, 0}` to `{log f_220, log f_221, df_221}`. The spin nuisance direction
is evaluated numerically from the Kerr QNM sequence. For each EFT fingerprint
we solve the weighted least-squares problem and profile over the nuisance
directions.

For the pyRing branch, the model is

```text
y_P = s alpha,
s = {d log f_221 / d alpha, d log tau_221 / d alpha}.
```

No remnant-mass or remnant-spin nuisance profiling is used in this branch,
because the public variables are already given as fractional deviations
relative to the Kerr prediction.

For both branches, the reported Gaussian 90 percent interval is

```text
alpha_hat +/- 1.6448536269514722 sigma_alpha.
```

The result should be read as a public-data projection. It is not a Bayes
factor, not a model-selection statistic, and not a substitute for a full
strain-level likelihood with priors, waveform systematics, and nonlinear
parameter correlations.

The RINGDOWN construction profiles over the local mass and spin directions,
but this is still only a Gaussian projection of public posterior products. It
does not marginalize over the full nonlinear degeneracy between the remnant
parameters and the beyond-Kerr coupling, and it does not include alternative
ringdown start times, waveform systematics, calibration uncertainty, or a
penalty for the additional coupling. We keep these limitations explicit
because they define the boundary between the present public-data exercise and
a true event-level beyond-Kerr inference. In particular, the quoted intervals
should not be read as guaranteed conservative bounds. They can be narrower
than intervals from a full strain-level Bayesian analysis if nonlinear
mass-spin-coupling correlations or prior-volume effects are important.

### 2.6 Static Supplied-Potential QNM Validation

The static branch is separated from the GW250114 rotating-remnant projection.
For a supplied master potential, the local code evolves

```text
d_t^2 Psi - d_x^2 Psi + V(r(x)) Psi = 0,
```

where `x` is the tortoise coordinate. The extracted waveform is fitted with
Prony, linear-prediction, or matrix-pencil methods according to the mode
content. The validation target is not nature. It is the corresponding
published QNM table or the local Schwarzschild `qnm` reference.

This step tests reproducibility of a perturbative calculation. It does not
decide whether the corresponding metric is realized astrophysically. It also
does not make a static model directly comparable to the rotating GW250114
remnant.

## 3. Results and Discussion

### 3.1 Public GW250114 Projection

The public RINGDOWN and pyRing branches give a consistent null result in the
linearized one-at-a-time higher-derivative projection. In the RINGDOWN branch,
all ten operator/polarization rows retain `alpha = 0` inside the configured
90 percent Gaussian interval. The largest nominal displacement from zero is
`0.497 sigma`, in the `lambda_odd plus` row.

The identical `0.497 sigma` value appearing across the RINGDOWN rows is a
consequence of the local profiling geometry rather than an independent
coincidence in every operator. The RINGDOWN observable vector has three
components, and the projection profiles over two nuisance directions
associated with remnant mass and spin. The remaining projected residual is
therefore effectively one-dimensional. Different one-at-a-time fingerprint
columns mainly change the scale and sign of the inferred `alpha`, while the
normalized distance of the public posterior mean from the profiled Kerr
surface remains the same at this order.

In the pyRing deviation branch, all ten rows again retain `alpha = 0` inside
the 90 percent interval. The largest nominal displacement is `1.222 sigma`, in
the `epsilon1 plus` row.

The two branches should be interpreted as complementary projections of
overlapping public analyses, not as independent constraints. Their
operator-by-operator comparison gives a maximum normalized projection
difference of `0.535 sigma` for the `epsilon1 plus` row. No row shows a
90 percent interval tension between the RINGDOWN and pyRing projections.

This is the sense in which the worked GW250114 example is a null result. It
does not prove that higher-derivative corrections are absent. It shows that,
within the public products and linearized one-at-a-time treatment used here,
there is no robust projected displacement from the Kerr value. The full
per-operator values of `alpha_hat`, `sigma_alpha`, the 90 percent intervals,
and the branch-to-branch comparison are kept in the numerical tables so that
the summary statement can be checked without rerunning the code.

### 3.2 Robustness Checks

The pyRing lower-tail filter sweep supports the same null statement. Across
the stricter, public, looser, and positive-domain-only filter choices,
`alpha = 0` remains inside every configured 90 percent interval. The largest
nominal filter-sweep displacement is `1.290 sigma`.

The empirical linearized posterior-sample projection gives the same qualitative
result without using only symmetric Gaussian intervals. In this check, no
empirical 90 percent interval excludes `alpha = 0`. The largest nominal
`abs(median)/sd` is `1.286`.

These checks support the stability of the public projection, but they are not
independent evidence and they do not replace a full Bayesian event-level
analysis. They address two limited failure modes of the Gaussian summary: the
choice of pyRing lower-tail filter and the use of a symmetric covariance
approximation. They do not test start-time dependence or waveform
model-dependence, because those choices are already built into the public
products used here.

### 3.3 Static QNM Benchmarks

The static branch is included for one reason: it turns the readiness ladder
from a slogan into a reproducible test. We take examples that are not directly
usable for GW250114, but that do supply a gravitational master potential or a
published QNM reference, and ask whether the spectrum can be reproduced before
any phenomenological interpretation is attempted. The present validation
scorecard contains `25` rows across Schwarzschild, braneworld tidal charge,
Bardeen NED, Hayward fundamental, and Hayward first-overtone examples. All
validation families pass the configured sub-percent criterion. The largest
validation discrepancy is `0.925%`, for the Hayward first-overtone
matrix-pencil extraction. This number measures reproduction of a published
reference table, not agreement with observation.

In the static summary, `max physical shift` means the maximum sampled absolute
relative change, in percent, of either `Re(M omega)` or `Abs(Im(M omega))`
with respect to the stated Schwarzschild or zero-parameter baseline. It is not
a complex-norm posterior, not a detector-space residual, and not an
observational exclusion statistic.

After validation, the same static examples can be used to measure physical QNM
shifts relative to a Schwarzschild or zero-parameter baseline. These shifts are
larger than the numerical validation errors:

```text
Bardeen NED:               max sampled shift 109.8% at alpha=0.7698
braneworld tidal charge:   max sampled shift 21.214% at q_tidal=2
Hayward fundamental:       max sampled shift 16.596% at gamma=1.18
Hayward first overtone:    max sampled shift 17.407% at gamma=1.18
```

This shows that the static test is not only a numerical reproduction exercise.
The threshold-crossing table gives sparse-grid 1, 3, 5, and 10 percent
crossings. However, these are diagnostic interpolation results, not
observational exclusion intervals.

### 3.4 Interpretation of the Static Branch

The static benchmarks should be read as readiness tests, not as GW250114
constraints. They answer a limited but useful question: if a phenomenological
metric is supplied together with a gravitational master potential, can the QNM
spectrum be reproduced and used in a controlled stress test? For the examples
considered here the answer is yes at the static supplied-potential level.

This is different from saying that Bardeen, Hayward, tidal-charge, or other
static backgrounds are ruled out, preferred, or disfavoured by GW250114. Such a
claim would require a rotating solution, the correct perturbation system, mode
identification for the remnant spin, priors on the extra parameters, and an
observational likelihood at the same approximation level.

The usefulness of the static branch is therefore diagnostic. It identifies
where a model sits on the ladder from metric-only phenomenology to
ringdown-ready theory. A model used in QPO or shadow studies may be useful at
the geodesic level, but that does not make it a gravitational-ringdown model.
The readiness standard makes this distinction explicit.

### 3.5 Limitations and Next Steps

The public GW250114 projection is intentionally conservative. It is linearized,
tests one effective coupling at a time, and uses public posterior products
rather than a full strain-level beyond-Kerr likelihood. It does not compute
Bayes factors or information criteria, and it does not explore multi-parameter
EFT directions. It also does not scan ringdown start time, because the goal is
to project the released spectroscopy products rather than to redo the event
analysis. These are limitations of scope rather than hidden assumptions.

The most direct extension of the observational branch is a full event-level
likelihood using the strain data, detector noise model, ringdown start-time
treatment, and beyond-Kerr waveform response in one inference problem. A second
extension is to include multi-parameter EFT directions and nonlinear responses.

The most direct extension of the readiness branch is to add more
theory-complete or supplied-potential cases. Reissner-Nordstrom or
gravito-electromagnetic perturbations provide a clean charged benchmark.
Slow-rotation supplied-potential examples would be especially valuable,
because they would reduce the conceptual distance between static benchmarks
and high-spin merger remnants without requiring a full generic rotating
solution. Dynamical Chern-Simons or other theory-backed rotating calculations
\cite{Chung2025DCSMetrics,Li2025DCSSlowRotation,Pierini2022EdGB} could be
more directly connected to high-spin remnant spectroscopy if the required mode
set is available. The criterion should remain the same: the model must provide
enough perturbation information to compute or reproduce gravitational QNM
spectra.

The main value of the present work is therefore not a claim of new physics. It
is a reproducible public-data example plus a ringdown-readiness standard. The
first part keeps observational statements controlled, while the second part
keeps theoretical phenomenology controlled. Together they define a practical
standard for future high-SNR ringdown events and for the alternative metrics
proposed to interpret them.

## Data and Code Availability

All scripts, configuration files, intermediate tables, validation summaries,
and manuscript-supporting outputs are included in the accompanying repository.
The public GW250114 inputs are drawn from the GWOSC/Zenodo event release. The
higher-derivative QNM fingerprints are imported from the public
BeyondKerrQNM/Cano et al. fit data. Static benchmark inputs are derived from
the cited published master potentials and QNM tables. The working repository
placeholder is `https://github.com/jvrba/gw250114-ringdown-readiness`; this
placeholder should be replaced by the final public archive before journal
submission.
