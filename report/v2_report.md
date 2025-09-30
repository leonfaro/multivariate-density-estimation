---
bibliography:
- biblio.bib
---

:::: titlepage
::: center
**Multivariate Density Estimation\
Comparing Transformation Random Forests, Triangular Transport Maps, and
Copulas**

Master's Thesis in Biostatistics (STA495)

by

Léon Kia Faro\

13-795-026

supervised by

Prof. Dr. Torsten Hothorn

Zurich, September 2025
:::
::::

# Abstract {#abstract .unnumbered}

This thesis evaluates three approaches to multivariate density
estimation for tabular data within a single, consistent pipeline:
separable triangular transport maps (TTM-Sep), Transformation Random
Forests (TRTF), and copulas (used only for
low-dimensional diagnostics, $K\!\le\!3$). All methods use standardized
inputs and a common evaluation protocol so that likelihoods,
diagnostics, and compute are directly comparable. In the configuration
studied (monotone CDF smoothing with the default forest aggregation), TRTF and
TTM-Sep yield the same triangular-likelihood form, which enables
like-for-like evaluation.

On Half-Moon ($n=250$), mean joint negative log-likelihoods (NLL; lower
is better) were $1.71$ (TRTF), $1.93$ (TTM-Sep), and $1.54$ (copula). On
a four-dimensional autoregressive generator they were $4.53$, $5.66$,
and $5.45$, respectively; permutation averages confirm order sensitivity
for triangular maps. On MiniBooNE ($K=43$; sum test log-likelihood),
TRTF reached $-30.01$ under the standard preprocessing and training
budget used here; published flow models report values around $-12$ to
$-16$ under their settings. These numbers are not strictly comparable
but indicate the relative accuracy of this configuration.

Overall, TRTF tends to lead within the separable family at low
dimension, while higher-dimensional datasets expose the limits of
separable structure. We report robustness checks (ordering), calibration
diagnostics, and the numerical safeguards used, and we outline
directions toward richer parameterizations within the same evaluation
frame.

# Introduction {#ch:intro}

Multivariate density estimation underpins likelihood-based modeling,
simulation, and uncertainty quantification for tabular data. Yet
striking a balance between flexibility, computation, and
interpretability remains difficult in practice. Three complementary
model families are widely used: (i) normalizing flows that learn
expressive change-of-variables maps
[@rezende2015variational; @papamakarios2017masked; @dinh2017real; @durkan2019neural; @kingma2018glow; @papamakarios2021normalizing],
(ii) transformation models and forests that target conditional
distributions through flexible parametric building blocks
[@hothorn2017transformation; @hothorn2018conditional; @hothorn2021transformation],
and (iii) copulas that decouple marginals from dependence
[@sklar1959fonctions; @nelsen2006introduction; @joe2014dependence; @nagler2017kdecopula].

We compare these approaches within a common **transport** frame that
puts them on the same footing
[@rosenblatt1952remarks; @knothe1957contributions; @bogachev2005triangular; @ramgraber2025friendly].
Concretely, we work in standardized evaluation coordinates $u$
constructed from training-split statistics only (to avoid leakage and
keep results comparable across estimators). Every model then couples $u$
to a simple Gaussian reference and reports densities via an exact change
of variables. The high-level frame is summarized in
Section [1.2](#sec:ch1-frame){reference-type="ref"
reference="sec:ch1-frame"} and detailed in
Chapter [2](#ch:background){reference-type="ref"
reference="ch:background"};
Figure [5.1](#fig:transport-schematic){reference-type="ref"
reference="fig:transport-schematic"} in
Appendix [5](#ch:appendix){reference-type="ref" reference="ch:appendix"}
shows the pipeline at a glance.

For clarity, efficiency, and interpretability we emphasize separable
triangular transport maps as a running theme. The lower-triangular
Jacobian yields $\mathcal{O}(K)$ evaluation, back-substitution enables
one-pass sampling, and the structure stabilizes conditional shapes
across contexts. These properties make comparisons transparent while
preserving explicit marginals and an interpretable dependence mechanism.

*Notation.* We use $\varphi$ and $\Phi$ for the univariate standard
normal density and CDF, respectively; later sections also write $\eta$
for the $K$-variate standard normal density. All references to the
evaluation and reporting protocol (datasets, preprocessing, metrics, and
timing) are consolidated in
Section [3.3](#sec:evaluation-protocol){reference-type="ref"
reference="sec:evaluation-protocol"}.

## Thesis and Problem Statement {#sec:ch1-problem}

We standardize elementwise with training-split statistics only;
$\mu,\sigma\in\mathbb{R}^K$. Equivalently,
$u=\operatorname{diag}(\sigma)^{-1}(x-\mu)$.

This thesis investigates tabular multivariate density estimation within
a unified transport-based evaluation frame. We compare separable
triangular transport maps (TTM-Sep), Transformation Random Forests
(TRTF), and copula baselines.

A transport perspective couples standardized data to a Gaussian
reference through a monotone lower-triangular map. This structure yields
exact likelihoods, transparent conditionals, exact inversion by back
substitution, and linear per-sample evaluation. The Rosenblatt and
Knothe rearrangements justify the triangular coupling for any variable
order [@rosenblatt1952remarks; @knothe1957contributions].

We standardize features using training statistics only.
Equation [\[eq:transport-standardise\]](#eq:transport-standardise){reference-type="eqref"
reference="eq:transport-standardise"} defines the standardized
coordinates $u$ used for evaluation. All derivatives and Jacobians are
computed in $u$. The diagonal affine correction in
Equation [\[eq:transport-affine\]](#eq:transport-affine){reference-type="eqref"
reference="eq:transport-affine"} reports log densities on the original
scale $x$. This convention keeps objectives, diagnostics, and
comparisons interoperable across estimators and datasets. All log
quantities are reported in nats.

We denote the $K$-variate standard normal density by $\eta$, and the
univariate density and CDF by $\varphi$ and $\Phi$. Abbreviations for
models and references appear in
Table [1.1](#tab:model-abbrev){reference-type="ref"
reference="tab:model-abbrev"}. Copulas decouple marginals from
dependence and serve as interpretable baselines.

Separable triangular maps decompose each component into a context shift
and a univariate monotone shape as in
Equation [\[eq:transport-separable\]](#eq:transport-separable){reference-type="eqref"
reference="eq:transport-separable"}. The decomposition fixes conditional
shape across contexts and stabilizes the triangular determinant. Under
strictly increasing conditional CDFs after standard monotone smoothing
and with the forest aggregation, TRTF implements the same separable
triangular likelihood via the probability integral transform. Copulas
preserve explicit marginals and model dependence on the unit hypercube;
in this thesis they serve strictly as low-dimensional ($K\!\le\!3$)
diagnostic baselines and are not evaluated on high-$K$ datasets.

We evaluate all estimators under the single protocol referenced above,
with matched preprocessing and reporting to keep results comparable.
Section [3.3](#sec:evaluation-protocol){reference-type="ref"
reference="sec:evaluation-protocol"} defines metrics and timing
conventions.

Figure [5.1](#fig:transport-schematic){reference-type="ref"
reference="fig:transport-schematic"} in
Appendix [5](#ch:appendix){reference-type="ref" reference="ch:appendix"}
visualizes the pipeline by showing standardization
$u = T_{\text{std}}(x)$, the triangular transport branch containing
TTM-Sep and TRTF, and the copula branch. Both branches feed the reported
outputs, namely log density, conditionals, sampling, calibration, and
compute, under the shared frame.

The central problem is to determine when separability is appropriate for
tabular data. We study how TRTF and copulas position themselves against
direct triangular transports inside the same reporting convention
(reported log densities on $x$ apply the affine correction in
Eq. [\[eq:transport-affine\]](#eq:transport-affine){reference-type="eqref"
reference="eq:transport-affine"}). Ordering effects, conditional
calibration, and computational trade-offs address this question.

On synthetic data, TRTF tends to outperform separable TTM variants yet
shares their separability limits; on the MiniBooNE benchmark it improves
on Gaussian references but trails published flow baselines.
Chapter [3](#ch:dataanalysis){reference-type="ref"
reference="ch:dataanalysis"} presents the evidence and discusses these
comparisons.

## The Transport Frame on One Page {#sec:ch1-frame}

To avoid duplication, the canonical derivations and notation live in
Chapter [2](#ch:background){reference-type="ref"
reference="ch:background"}. This section serves only as a map: we
standardize with train-only statistics
(Eq. [\[eq:transport-standardise\]](#eq:transport-standardise){reference-type="eqref"
reference="eq:transport-standardise"}), evaluate likelihoods via the
pullback
(Eq. [\[eq:transport-pullback\]](#eq:transport-pullback){reference-type="eqref"
reference="eq:transport-pullback"}), exploit the triangular determinant
factorization
(Eq. [\[eq:transport-det\]](#eq:transport-det){reference-type="eqref"
reference="eq:transport-det"}), and apply the affine correction for
reporting
(Eq. [\[eq:transport-affine\]](#eq:transport-affine){reference-type="eqref"
reference="eq:transport-affine"}). Separable components are defined in
Eq. [\[eq:transport-separable\]](#eq:transport-separable){reference-type="eqref"
reference="eq:transport-separable"}.
Figure [5.1](#fig:transport-schematic){reference-type="ref"
reference="fig:transport-schematic"} in
Appendix [5](#ch:appendix){reference-type="ref" reference="ch:appendix"}
illustrates the pipeline.

Notation remains consistent. We write $\eta$ for the $K$-variate
standard normal density, and $\varphi$ and $\Phi$ for the univariate
standard normal density and CDF. We reserve $u$ for standardized
coordinates and $x$ for original coordinates, and we compute all
derivatives with respect to $u$. These choices align symbols across
Chapters [1](#ch:intro){reference-type="ref"
reference="ch:intro"}--[3](#ch:dataanalysis){reference-type="ref"
reference="ch:dataanalysis"} and prevent ambiguity in later diagnostics
and tables.

This one-page frame removes duplicated exposition from
Chapter [2](#ch:background){reference-type="ref"
reference="ch:background"}. It establishes where logs and Jacobians live
and makes complexity, inversion, and units explicit before the
comparisons that follow.
Section [1.1](#sec:ch1-problem){reference-type="ref"
reference="sec:ch1-problem"} documented the motivation, and
Section [1.3](#sec:ch1-contributions){reference-type="ref"
reference="sec:ch1-contributions"} states the resulting contributions
and research questions.

## Contributions and Research Questions {#sec:ch1-contributions}

This section states the contributions and the research questions, and
maps them to the chapters and figures that deliver the evidence. We
adopt the shared transport frame summarized in
Section [1.2](#sec:ch1-frame){reference-type="ref"
reference="sec:ch1-frame"};
Chapter [2](#ch:background){reference-type="ref"
reference="ch:background"} records notation and assumptions, and
Figure [5.1](#fig:transport-schematic){reference-type="ref"
reference="fig:transport-schematic"} anchors the comparisons.

The first contribution formalizes a unified likelihood view for
separable triangular transport maps, Transformation Random Forests, and
copula baselines. Where we claim an equivalence between TRTF and
separable triangular transports, it holds under the conditions made
explicit in Section [2.2.1](#sec:transport-trtf){reference-type="ref"
reference="sec:transport-trtf"} (strictly increasing conditional CDFs
after monotone smoothing with the forest aggregation).
Chapter [2](#ch:background){reference-type="ref"
reference="ch:background"} and
Figure [5.1](#fig:transport-schematic){reference-type="ref"
reference="fig:transport-schematic"} establish the conventions and
remove duplication in later chapters.

The second contribution provides empirical benchmarks under a single
protocol with matched preprocessing and reporting. We evaluate TTM-Sep,
TTM-Marg, and TRTF on synthetic generators and real tabular data;
copulas are included only as low-dimensional ($K\!\le\!3$) diagnostic
baselines (e.g., Half-Moon, 4D) and are not used in high-$K$ studies.
The protocol records three families of measurements: average test
log-likelihoods, conditional diagnostics based on probability integral
transforms, and compute indicators for training and per-sample
evaluation. Section [3.3](#sec:evaluation-protocol){reference-type="ref"
reference="sec:evaluation-protocol"} defines the protocol,
Section [3.4](#sec:synthetic-results){reference-type="ref"
reference="sec:synthetic-results"} presents the synthetic and
autoregressive results, and
Section [3.5](#sec:realdata){reference-type="ref"
reference="sec:realdata"} positions our measurements against published
normalizing-flow baselines where appropriate.

The third contribution distills practical guidance from the unified
frame and the benchmarks. We state operational choices that preserve
comparability, highlight ordering sensitivity and separability limits,
and summarize when copulas serve as informative baselines.
Chapter [4](#ch:conclusion){reference-type="ref"
reference="ch:conclusion"} consolidates these points as actionable
recommendations and records limitations that motivate richer
parameterizations or alternative predictors.

Two questions drive the empirical study and bind the contributions to
specific measurements. The first question asks how TRTF compares with
TTM-Sep and copula baselines on synthetic data. All estimators share the
transport frame in this comparison. We answer by reporting average test
negative log-likelihoods, conditional negative log-likelihood
decompositions, and probability integral transform diagnostics, with
timing summaries that quantify practical cost.
Section [3.4](#sec:synthetic-results){reference-type="ref"
reference="sec:synthetic-results"} provides the corresponding tables and
figures.

The second question asks how closely our TRTF results on real benchmarks
approach the published performance of modern normalizing flows under the
standard preprocessing. We answer by placing our test log-likelihoods
beside reported numbers from the literature. The gaps are interpreted
through the separable Jacobian constraint and compute profiles.
Section [3.5](#sec:realdata){reference-type="ref"
reference="sec:realdata"} reports these comparisons, and
Chapter [4](#ch:conclusion){reference-type="ref"
reference="ch:conclusion"} interprets their implications for model
choice.

Taken together, these commitments make the comparisons interpretable,
keep units and complexity explicit, and prepare the reader for the
empirical evidence that answers the two questions under a single,
transparent evaluation frame.

::: {#tab:model-abbrev}
  Label        Meaning
  ------------ -----------------------------------------------------------------------------
  TTM-Marg     Marginal triangular transport map (per-dimension; no context)
  TTM-Sep      Separable triangular transport map (additive: $g_k$ shift + monotone $h_k$)
  TRTF         Transformation Random Forests (axis-parallel splits)
  True-Marg    Oracle marginal density
  True-Joint   Oracle conditional joint density
  Copula       Copula baseline (Gaussian or nonparametric)

  : Model abbreviations used throughout the thesis.
:::

# Methodological Background {#ch:background}

## Transport Frame and Notation {#sec:transport-frame}

This section fixes the standardized coordinate system, notation, and
algebraic identities used throughout the thesis. The motivation and
schematic in Figure [5.1](#fig:transport-schematic){reference-type="ref"
reference="fig:transport-schematic"} housed in
Appendix [5](#ch:appendix){reference-type="ref" reference="ch:appendix"}
remain valid; here we strip the exposition down to the formulas needed
in later chapters. We summarize the standardized pullback likelihood,
state the triangularity assumption, and record the Jacobian
factorization that drives evaluation and inversion.

We work with observations on the original scale $x \in \mathbb{R}^K$.
Training-split statistics define a fixed standardization map
$$u \;=\; T_{\mathrm{std}}(x) \;=\; (x-\mu)\oslash\sigma,\qquad \sigma_k>0,\label{eq:transport-standardise}$$
where $\mu$ and $\sigma$ denote the empirical mean and standard
deviation estimated on the training split and $\oslash$ denotes
elementwise division. In words, we shift and rescale features once,
using training data only, and keep all derivatives and Jacobians in
$u$-space to avoid leakage and to ensure comparability across
estimators.

The standardized density $\pi_U$ is coupled to a simple reference
through a monotone triangular map $S:u\mapsto z$. Throughout the thesis
the reference is the $K$-variate standard normal density $\eta(z)$. The
pullback identity then reads
$$\pi_U(u) \;=\; \eta\!\left(S(u)\right)\,\left|\det\nabla_u S(u)\right|,\label{eq:transport-pullback}$$
which evaluates the reference at $S(u)$ and applies the exact volume
correction given by the Jacobian determinant. Reporting log densities on
the original scale requires only the diagonal affine correction implied
by standardization,
$$\log \pi_X(x) \;=\; \log \pi_U\!\left(T_{\mathrm{std}}(x)\right) - \sum_{k=1}^{K}\log\sigma_k.\label{eq:transport-affine}$$
We therefore differentiate with respect to $u$, and we convert to
$x$-scale only at reporting time.

The transport is assumed to be lower triangular and componentwise
monotone,
$$S(u) \;=\; \big(S_1(u_1), S_2(u_{1:2}), \ldots, S_K(u_{1:K})\big), \qquad \partial_{u_k}S_k(u_{1:k})>0,\label{eq:transport-triangular}$$
so the Jacobian $\nabla_u S(u)$ is lower triangular. Its determinant
factorizes into a sum of one-dimensional log derivatives,
$$\log \big|\det \nabla_u S(u)\big| \;=\; \sum_{k=1}^{K}\log \partial_{u_k}S_k(u_{1:k}).\label{eq:transport-det}$$
The factorization yields $\mathcal{O}(K)$ evaluation cost per-sample,
improves numerical stability, and guarantees global invertibility:
strictly monotone diagonal derivatives let us recover $x$ by solving $K$
one-dimensional monotone equations in sequence, mirroring the Rosenblatt
and Knothe rearrangements
[@rosenblatt1952remarks; @knothe1957contributions].

Table [2.1](#tab:transport-notation){reference-type="ref"
reference="tab:transport-notation"} consolidates the notation used in
this transport frame. All derivatives and Jacobians act on $u$; the
affine correction
[\[eq:transport-affine\]](#eq:transport-affine){reference-type="eqref"
reference="eq:transport-affine"} converts log densities back to $x$ for
reporting. The remainder of this chapter adopts this frame.
Section [2.2](#sec:transport-separable){reference-type="ref"
reference="sec:transport-separable"} details the separable triangular
parameterization used for direct transports.
Section [2.2.1](#sec:transport-trtf){reference-type="ref"
reference="sec:transport-trtf"} shows how Transformation Random Forests
induce the same triangular likelihood via the probability integral
transform under the conditions stated there (strictly increasing
conditional CDFs after monotone smoothing with the forest aggregation).
Section [2.3](#sec:transport-copula){reference-type="ref"
reference="sec:transport-copula"} places copulas in the same reporting
convention.

::: {#tab:transport-notation}
  Symbol                      Meaning
  --------------------------- ---------------------------------------------------
  $x \in \mathbb{R}^K$        Original features on the data scale
  $T_{\mathrm{std}}$          Standardization map using training $(\mu,\sigma)$
  $u = T_{\mathrm{std}}(x)$   Standardized evaluation coordinates
  $z \in \mathbb{R}^K$        Reference coordinates after transport
  $S:u\mapsto z$              Monotone lower-triangular transport map
  $\nabla_u S(u)$             Jacobian of $S$ with respect to $u$
  $\eta(z)$                   $K$-variate standard normal density
  $\varphi(t)$, $\Phi(t)$     Univariate standard normal density and CDF
  $\pi_U$, $\pi_X$            Densities on $u$- and $x$-space, respectively
  $\mu$, $\sigma$             Training mean vector and positive scales
  $K$                         Dimension of the feature vector

  : Notation for the transport frame used in
  Chapters [2](#ch:background){reference-type="ref"
  reference="ch:background"}
  and [3](#ch:dataanalysis){reference-type="ref"
  reference="ch:dataanalysis"}. All derivatives and Jacobians are taken
  with respect to $u$; log densities on $x$-space apply the affine
  correction in
  Equation [\[eq:transport-affine\]](#eq:transport-affine){reference-type="eqref"
  reference="eq:transport-affine"}.
:::

## Separable Triangular Maps and Transformation Random Forests as Transport {#sec:transport-separable}

This section unifies separable triangular maps and Transformation Random
Forests (TRTF) within the transport frame fixed in
Section [2.1](#sec:transport-frame){reference-type="ref"
reference="sec:transport-frame"}. Both estimators realize a monotone
lower-triangular map $S:u\mapsto z$ that couples the standardized target
to the Gaussian reference $\eta$. The use of triangular transports
builds on modern measure-transport literature; see, for instance,
triangular transformations and their properties in
@bogachev2005triangular.
Figure [5.1](#fig:transport-schematic){reference-type="ref"
reference="fig:transport-schematic"} in
Appendix [5](#ch:appendix){reference-type="ref" reference="ch:appendix"}
illustrates the shared backbone and locates [TTM-Sep]{.smallcaps} and
[TRTF]{.smallcaps} on the transport branch introduced in
Chapter [1](#ch:intro){reference-type="ref" reference="ch:intro"}. We
focus on shared likelihood identities, modeling assumptions, and limits
of separability, and defer implementation details to
Chapter [3](#ch:dataanalysis){reference-type="ref"
reference="ch:dataanalysis"} and
Appendix [5](#ch:appendix){reference-type="ref"
reference="ch:appendix"}.

The goal is to state a single likelihood for both constructions, clarify
what separability permits, and identify failure modes that motivate
richer parameterizations. We do not pursue non-additive TRTF predictors,
cross-term triangular maps, or ordering heuristics in this section;
Chapter [3](#ch:dataanalysis){reference-type="ref"
reference="ch:dataanalysis"} evaluates those choices empirically and
Appendix [5](#ch:appendix){reference-type="ref" reference="ch:appendix"}
documents routines and defaults.

We adopt the notation introduced in
Section [2.1](#sec:transport-frame){reference-type="ref"
reference="sec:transport-frame"}. Coordinates satisfy
$u=T_{\mathrm{std}}(x)$, the reference density is $\eta(z)$, and the
pullback identity
[\[eq:transport-pullback\]](#eq:transport-pullback){reference-type="eqref"
reference="eq:transport-pullback"} gives
$\pi_U(u)=\eta(S(u))\,|\det\nabla_u S(u)|$. The map is lower-triangular
with strictly positive diagonal partial derivatives, which yields the
sum decomposition in
Equation [\[eq:transport-det\]](#eq:transport-det){reference-type="eqref"
reference="eq:transport-det"}. These conventions keep derivatives in
$u$-space and apply the affine correction
[\[eq:transport-affine\]](#eq:transport-affine){reference-type="eqref"
reference="eq:transport-affine"} only when reporting $\log \pi_X(x)$.

We restrict attention to separable triangular maps. Component $k$
decomposes into a context shift and a univariate monotone shape,
$$S_k(u_{1:k}) \;=\; g_k(u_{1:k-1}) + h_k(u_k),\qquad \log \partial_{u_k}S_k(u_{1:k}) \;=\; \log h_k'(u_k),\label{eq:transport-separable}$$
which fixes context effects in $g_k$ and reserves $h_k$ for the
one-dimensional marginal shape. Intuitively, earlier coordinates
translate the location, while the conditional shape along $u_k$ remains
fixed across contexts. The Jacobian contribution depends only on $u_k$,
which reduces per-sample evaluation cost and simplifies inversion.

::: shaded
**Assumptions.** Unless stated otherwise, we assume:

- *Lower-triangularity:* $S$ has the structure in
  Eq. [\[eq:transport-triangular\]](#eq:transport-triangular){reference-type="eqref"
  reference="eq:transport-triangular"}.

- *Strict monotone coordinates:* $\partial_{u_k} S_k(u_{1:k}) > 0$ for
  all $k$ and all arguments.

- *Separable component:*
  Eq. [\[eq:transport-separable\]](#eq:transport-separable){reference-type="eqref"
  reference="eq:transport-separable"} holds, so conditional shape along
  $u_k$ is fixed across contexts.
:::

Substituting the standard normal reference into
[\[eq:transport-pullback\]](#eq:transport-pullback){reference-type="eqref"
reference="eq:transport-pullback"} produces a separable objective,
$$\log \pi_U(u) \;=\; \sum_{k=1}^{K}\Big[\log \varphi\!\big(S_k(u_{1:k})\big) + \log h_k'(u_k)\Big],\label{eq:transport-likelihood}$$
where $\varphi$ denotes the univariate standard normal density.
Equation [\[eq:transport-likelihood\]](#eq:transport-likelihood){reference-type="eqref"
reference="eq:transport-likelihood"} splits the log density into a
reference fit and an exact volume correction. In plain language, the
model evaluates how Gaussian each transformed coordinate appears, then
corrects for the local stretch induced by $h_k$. The same decomposition
produces linear per-sample time in $K$ and stable accumulation of log
derivatives.

The negative log-likelihood per-sample takes the quadratic-plus-barrier
form
$$\mathcal{L}(u) \;=\; \sum_{k=1}^{K}\Big[\tfrac{1}{2}\,S_k(u_{1:k})^2 - \log h_k'(u_k)\Big],\label{eq:transport-loss}$$
which follows because
$\log \varphi(t) = -\tfrac{1}{2}t^2 - \tfrac{1}{2}\log(2\pi)$ and
constants independent of the parameters drop out.
Equation [\[eq:transport-loss\]](#eq:transport-loss){reference-type="eqref"
reference="eq:transport-loss"} pulls each component toward the reference
while preventing degenerate derivatives through the log barrier. In
practice we enforce $h_k'(u_k)>0$ by construction and control tails with
mild regularization; implementation choices appear in
Chapter [3](#ch:dataanalysis){reference-type="ref"
reference="ch:dataanalysis"} and
Appendix [5](#ch:appendix){reference-type="ref"
reference="ch:appendix"}.

Separable structure encodes clear modeling assumptions. Conditional
variance, skewness, and modality do not change with the preceding
coordinates once $g_k$ shifts location. Consequently, separable maps can
underfit heteroskedastic or multimodal conditionals, which manifests as
U-shaped or inverted-U probability integral transform (PIT) diagnostics.
Variable ordering also matters for finite bases because triangular
transports are anisotropic, even though a Knothe--Rosenblatt
rearrangement exists for any ordering
[@rosenblatt1952remarks; @knothe1957contributions]. These caveats guide
the robustness checks in
Chapter [3](#ch:dataanalysis){reference-type="ref"
reference="ch:dataanalysis"}.

### Transformation Random Forests within the Transport Frame {#sec:transport-trtf}

Transformation Random Forests
[@hothorn2017transformation; @hothorn2018conditional; @hothorn2021transformation]
fit into the same transport frame through the probability integral
transform. Let $\widehat F_k(\cdot \mid u_{1:k-1})$ denote the strictly
increasing conditional CDF returned by a TRTF for coordinate $k$. (In
practice, forest CDFs can be stepwise; we assume a measurable, strictly
increasing version after standard monotone smoothing so that inversion
and derivatives are well-defined.) The induced triangular component is
$$S_k(u_{1:k}) \;=\; \Phi^{-1}\!\Big(\widehat F_k(u_k \mid u_{1:k-1})\Big),\label{eq:transport-trtf-map}$$
which maps conditionals to standard normal margins. In plain language,
TRTF predicts a conditional CDF, then the probit transform places the
result on the Gaussian reference scale. Differentiating
$\Phi\!\big(S_k(u_{1:k})\big)=\widehat F_k(u_k \mid u_{1:k-1})$ with
respect to $u_k$ yields
$$\widehat \pi_k(u_k \mid u_{1:k-1}) \;=\; \varphi\!\big(S_k(u_{1:k})\big)\,\partial_{u_k}S_k(u_{1:k}),\label{eq:transport-trtf-likelihood}$$
which is exactly the pullback factor in
Equation [\[eq:transport-likelihood\]](#eq:transport-likelihood){reference-type="eqref"
reference="eq:transport-likelihood"}. Summing over $k$ recovers
Equation [\[eq:transport-likelihood\]](#eq:transport-likelihood){reference-type="eqref"
reference="eq:transport-likelihood"} in standardized coordinates.

The additive-predictor TRTF used in this thesis yields a separable
transport. Under the model
$$\widehat F_k(u_k \mid u_{1:k-1}) \;=\; \Phi\!\big(h_k(u_k) + g_k(u_{1:k-1})\big),\label{eq:transport-trtf-additive}$$
we obtain
$$S_k(u_{1:k}) \;=\; h_k(u_k) + g_k(u_{1:k-1}),\qquad \partial_{u_k}S_k(u_{1:k}) \;=\; h_k'(u_k),\label{eq:transport-trtf-separable}$$
so TRTF implements the same separable triangular likelihood as the
direct parameterization in
Equation [\[eq:transport-separable\]](#eq:transport-separable){reference-type="eqref"
reference="eq:transport-separable"}. The map is monotone in $u_k$ by
construction, the Jacobian depends only on $u_k$, and inversion proceeds
by back-substitution identical to the separable map. This equivalence
underpins the empirical comparisons in
Chapter [3](#ch:dataanalysis){reference-type="ref"
reference="ch:dataanalysis"}.

The equivalence also clarifies limits. Additive TRTF predictors shift
location but cannot alter conditional shape with context, which mirrors
the separable constraint. Axis-aligned partitions stabilize estimation,
yet they do not remove residual multimodality when the conditional shape
varies with $u_{1:k-1}$. These limits are visible in PIT diagnostics and
conditional negative log-likelihood decompositions on synthetic studies.

We emphasize operational scope and supporting references. All
derivatives and Jacobians are computed in standardized coordinates,
evaluation uses the triangular pullback, and reported log densities on
the original scale include the affine correction
[\[eq:transport-affine\]](#eq:transport-affine){reference-type="eqref"
reference="eq:transport-affine"}. Implementation details on basis
choices for $h_k$, feature construction for $g_k$, regularization,
derivative clipping, timing, and memory footprints appear in
Chapter [3](#ch:dataanalysis){reference-type="ref"
reference="ch:dataanalysis"} and
Appendix [5](#ch:appendix){reference-type="ref"
reference="ch:appendix"}, which also provides pseudo-code for both
estimators. Figure [5.1](#fig:transport-schematic){reference-type="ref"
reference="fig:transport-schematic"} in
Appendix [5](#ch:appendix){reference-type="ref" reference="ch:appendix"}
visualizes how the [TTM-Sep]{.smallcaps} and [TRTF]{.smallcaps} branches
share the same computational path from standardized data to reported
likelihoods.

In summary, separable triangular maps and additive-predictor TRTF
realize the same lower-triangular likelihood once the data are
standardized and the conditional CDFs are strictly increasing after
monotone smoothing. The shared structure yields exact likelihoods, exact
inversion, transparent conditionals, and linear per-sample complexity,
but it restricts context-dependent shape.
Section [2.3](#sec:transport-copula){reference-type="ref"
reference="sec:transport-copula"} positions copulas within the same
reporting convention to decouple marginals from dependence.

## Copula Baselines {#sec:transport-copula}

This section positions copulas within the unified transport frame and
links their reported likelihoods to the evaluation conventions used for
triangular maps and Transformation Random Forests. Copulas decouple
marginal modeling from dependence modeling by pairing univariate
marginals with a separate dependence density on the unit hypercube
[@nelsen2006introduction; @joe2014dependence].
Figure [5.1](#fig:transport-schematic){reference-type="ref"
reference="fig:transport-schematic"} in
Appendix [5](#ch:appendix){reference-type="ref" reference="ch:appendix"}
displays the copula branch beside the triangular branch and highlights
how both yield comparable reported log densities under the shared
evaluation pipeline.

We begin with pseudo-observations built from training-split marginals.
Let $\widehat F_k$ denote the strictly increasing empirical or smoothed
CDF of $X_k$ estimated on the training split. Define the
pseudo-observations and their probit transform as
$$v_k \;=\; \widehat F_k(x_k),\qquad z_k \;=\; \Phi^{-1}(v_k),\label{eq:copula-probit}$$
which map each coordinate to $(0,1)$ and then to $\mathbb{R}$ through
the probit function. In plain language, the marginals become uniform
scores, and $z$ records those scores on a Gaussian scale. Mid-ranks and
clamping near $(0,1)$ stabilize the transformation in finite samples.

The copula representation combines marginal densities with a dependence
factor. The joint log density on the original scale satisfies
$$\log \widehat \pi_X(x) \;=\; \sum_{k=1}^{K} \log \widehat f_k(x_k) + \log c\!\left(v_1,\ldots,v_K\right),\label{eq:copula-logdensity}$$
where $c$ denotes the copula density on $(0,1)^K$.
Equation [\[eq:copula-logdensity\]](#eq:copula-logdensity){reference-type="eqref"
reference="eq:copula-logdensity"} separates the task into two parts: fit
interpretable marginals and correct for dependence through $\log c$. The
reported quantity already lives on the original scale, so the affine
correction in
Equation [\[eq:transport-affine\]](#eq:transport-affine){reference-type="eqref"
reference="eq:transport-affine"} is unnecessary.
Figure [5.1](#fig:transport-schematic){reference-type="ref"
reference="fig:transport-schematic"} in
Appendix [5](#ch:appendix){reference-type="ref" reference="ch:appendix"}
uses the equivalent shorthand $\log c(z)$ because dependence is
evaluated through $z=\Phi^{-1}(v)$.

The independence baseline fixes a lower bound for dependence modeling.
Setting $c \equiv 1$ yields
$$\log \widehat \pi_X^{\mathrm{ind}}(x) \;=\; \sum_{k=1}^{K} \log \widehat f_k(x_k),\label{eq:copula-independence}$$
which treats coordinates as independent after marginal fitting.
Chapter [3](#ch:dataanalysis){reference-type="ref"
reference="ch:dataanalysis"} uses this baseline as a reference point in
evaluation tables and figures.

The Gaussian copula specifies elliptical dependence through a
correlation matrix $\Sigma$. With $z=\Phi^{-1}(v)$, the copula density
admits the closed form
$$c_{\Sigma}(v) \;=\; |\Sigma|^{-1/2}\,\exp\!\Big(-\tfrac{1}{2}\,z^{\top}(\Sigma^{-1} - I)z\Big),\label{eq:copula-gaussian}$$
which reduces dependence estimation to fitting $\Sigma$ on the
transformed scores. In plain language, the Gaussian copula bends the
joint shape away from independence according to $\Sigma$ while
preserving the learned marginals.

A low-dimensional nonparametric variant avoids elliptical assumptions at
small $K$. We fit a kernel density $\widehat f_Z$ on $z=\Phi^{-1}(v)$
and recover the copula density by
$$c(v) \;=\; \frac{\widehat f_Z\!\big(\Phi^{-1}(v)\big)}{\prod_{k=1}^{K} \varphi\!\big(\Phi^{-1}(v_k)\big)},\label{eq:copula-kde}$$
which applies the change of variables from $z$ back to $v$ and yields a
proper copula density. In words, the kernel density models the joint
shape of the probit scores, and division by the product of standard
normal densities restores the unit-cube scale. This approach is viable
only for small $K$, where kernel density estimation remains accurate and
stable. We implement this baseline via the `kdecopula` package
[@nagler2017kdecopula].
Chapter [3](#ch:dataanalysis){reference-type="ref"
reference="ch:dataanalysis"} employs it strictly as a diagnostic
baseline.

The transport frame keeps reporting interoperable across modeling
branches despite distinct parameterizations. Triangular maps and TRTF
evaluate the pullback likelihood in standardized coordinates and apply
the fixed affine correction
[\[eq:transport-affine\]](#eq:transport-affine){reference-type="eqref"
reference="eq:transport-affine"} when mapping back to $x$. Copulas
operate on $x$ directly through
Equation [\[eq:copula-logdensity\]](#eq:copula-logdensity){reference-type="eqref"
reference="eq:copula-logdensity"}, yet
Figure [5.1](#fig:transport-schematic){reference-type="ref"
reference="fig:transport-schematic"} in
Appendix [5](#ch:appendix){reference-type="ref" reference="ch:appendix"}
shows how the probit scores $z$ maintain comparability with the Gaussian
reference used above. This alignment keeps objectives and diagnostics
consistent across Chapters [2](#ch:background){reference-type="ref"
reference="ch:background"}
and [3](#ch:dataanalysis){reference-type="ref"
reference="ch:dataanalysis"}.

Modeling choices and limits follow from the chosen copula. The Gaussian
copula imposes elliptical dependence and may misrepresent tail behavior
or localized asymmetry. The nonparametric variant mitigates these issues
only at small dimension and sufficient sample size. The independence
baseline provides a transparent reference when dependence is weak or
data are scarce. These caveats motivate treating copulas as
interpretable baselines rather than definitive high-dimensional models
in the empirical study of
Chapter [3](#ch:dataanalysis){reference-type="ref"
reference="ch:dataanalysis"}. Sklar's theorem underlies all
constructions above and formalizes the decoupling of marginals from
dependence [@sklar1959fonctions].

# Data Analysis and Validation {#ch:dataanalysis}

This chapter turns the commitments of
Chapters [1](#ch:intro){reference-type="ref" reference="ch:intro"} and
[2](#ch:background){reference-type="ref" reference="ch:background"} into
a practical modeling program. Our aim is to express three model
families---triangular transport maps (TTM), transformation random
forests (TRTF), and copulas---within a common transport framework so
that likelihoods, calibration, and computational cost are directly
comparable. Every method we study standardizes the data, learns a
monotone triangular map to a simple reference, and evaluates Jacobians
in the standardized space. That alignment keeps objectives, diagnostics,
and reported log-densities interoperable.

**Model abbreviations.** We follow
Table [1.1](#tab:model-abbrev){reference-type="ref"
reference="tab:model-abbrev"} for concise labels (TTM-Marg, TTM-Sep,
TRTF, True-Marg/True-Joint, Copula) across tables and figures.

## Datasets and Preprocessing {#sec:datasets-preprocessing}

This section fixes data sources, generators, and preprocessing so
likelihoods, calibration, and compute remain comparable across models.
All estimators operate in standardized coordinates, evaluate Jacobians
in that space, and report log densities on the original scale using the
common affine correction. We keep a single triangular-map direction
$S:u \rightarrow z$ across methods to avoid mixed objectives. To avoid
symbol collisions, $\sigma$ denotes standardization scales only;
logistic gates are written $\operatorname{logistic}(\cdot)$ throughout.

We standardize features with training-split statistics only.
Equation [\[eq:transport-standardise\]](#eq:transport-standardise){reference-type="eqref"
reference="eq:transport-standardise"} defines
$u = T_{\mathrm{std}}(x) = (x-\mu)\oslash \sigma$ with $\sigma_k > 0$.
We evaluate $\log \pi_U(u)$ through the pullback identity in
Equation [\[eq:transport-pullback\]](#eq:transport-pullback){reference-type="eqref"
reference="eq:transport-pullback"}, apply the triangular factorization
from
Equation [\[eq:transport-det\]](#eq:transport-det){reference-type="eqref"
reference="eq:transport-det"}, then convert to $\log \pi_X(x)$ using the
diagonal correction in
Equation [\[eq:transport-affine\]](#eq:transport-affine){reference-type="eqref"
reference="eq:transport-affine"}. We report average test negative
log-likelihoods (NLL) in nats. Negative per-dimension NLL values can
occur because valid densities may exceed one on subdomains.
Figure [5.1](#fig:transport-schematic){reference-type="ref"
reference="fig:transport-schematic"} in
Appendix [5](#ch:appendix){reference-type="ref" reference="ch:appendix"}
shows the standardized pipeline shared by transport maps, Transformation
Random Forests, and copulas.

We use fixed train, validation, and test splits with proportions
$0.60/0.20/0.20$ unless a benchmark provides official splits. Synthetic
studies report results for $n \in \{25, 50, 100, 250\}$ and use $n=250$
for headline tables and figures; for real-data benchmarks we use $N$ to
denote the training budget. The canonical four-dimensional ordering is
$(1,2,3,4)$. Robustness to ordering is assessed by averaging over all
$4! = 24$ permutations. We adopt the natural column order for real
datasets. We fix seeds $\{11, 13, 17, 19, 23\}$ for data generation and
model fitting, and we average repeated runs with standard errors to
quantify stochastic variability.
Figure [3.2](#fig:autoregressive-joint-calibration){reference-type="ref"
reference="fig:autoregressive-joint-calibration"} displays the $20\%$
test split for the synthetic calibration study.

The Half-Moon dataset provides a curved, bimodal joint in $K=2$. We draw
a class $Y \sim \mathrm{Bernoulli}(0.5)$, an angle
$\Theta \sim \mathrm{Unif}[0,\pi]$, and additive noise
$\varepsilon \sim \mathcal{N}(0, \sigma^2 I_2)$ with $\sigma = 0.10$.
For $Y = 0$ we set $m(\Theta) = (\cos \Theta, \sin \Theta)$. For $Y = 1$
we set $m(\Theta) = (1 - \cos \Theta, -\sin \Theta + 0.5)$. The observed
$X = m(\Theta) + \varepsilon$. The "True joint" oracle evaluates the
mixture density by numerical quadrature over $\Theta$ with the known
Gaussian noise, and the "True conditional" oracle conditions on $Y$.
Figure [3.1](#fig:halfmoon-panels){reference-type="ref"
reference="fig:halfmoon-panels"} (p. ) shows representative contour
plots at $n=250$, which align with this generator.
Table [3.2](#tab:halfmoon-nll){reference-type="ref"
reference="tab:halfmoon-nll"} (p. ) reports the corresponding NLLs.

The four-dimensional autoregressive generator combines Gaussian,
exponential, beta, and gamma components to induce heteroskedasticity,
skew, and conditional multimodality. The first coordinate is
$X_1 \sim \mathcal{N}(0,1)$. The second coordinate is independent
$X_2 \sim \mathrm{Exp}(\lambda_0)$ with rate $\lambda_0 = 1$. The third
coordinate lies on $(0,1)$ and is a context-gated mixture of two beta
laws,
$X_3 \mid X_{1:2} \sim w\,\mathrm{Beta}(\alpha_1, \beta_1) + (1 - w)\,\mathrm{Beta}(\alpha_2, \beta_2)$.
We set $(\alpha_1, \beta_1) = (2.5, 5.0)$ and
$(\alpha_2, \beta_2) = (5.0, 2.5)$. The mixing weight is
$w = \operatorname{logistic}(\gamma_0 + \gamma_1 X_1 + \gamma_2(X_2 - 1))$
with $\operatorname{logistic}(\cdot)$ the logistic function and
$(\gamma_0, \gamma_1, \gamma_2) = (0, 1.5, 1.0)$. The fourth coordinate
is positive and conditionally heteroskedastic,
$X_4 \mid X_{1:3} \sim \tilde{w}\,\mathrm{Gamma}(k_1, r_1(X_2)) + (1 - \tilde{w})\,\mathrm{Gamma}(k_2, r_2(X_2))$.
We set shapes $(k_1, k_2) = (3, 6)$, rates $r_1(X_2) = 1 + 0.5 X_2$ and
$r_2(X_2) = 0.75 + 0.25 X_2$, and gate
$\tilde{w} = \operatorname{logistic}(\delta_0 + \delta_1 X_1 + \delta_3(X_3 - 0.5))$
with $(\delta_0, \delta_1, \delta_3) = (0, 1.0, 3.0)$. The "True joint"
baseline uses these closed-form conditionals to evaluate the exact joint
density, while the "True marginal" baseline uses the corresponding
univariate marginals and ignores dependence.

::: {#tab:autoregressive-config}
  Coordinate           Distribution    Parameters / gate
  -------------------- --------------- -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  $X_1$                Normal          $\mathcal{N}(0,1)$
  $X_2$                Exponential     $\mathrm{rate} = \lambda_0 = 1$
  $X_3 \mid X_{1:2}$   Mixture Beta    $\mathrm{Beta}(2.5,5.0)$ / $\mathrm{Beta}(5.0,2.5)$; $w = \operatorname{logistic}(\gamma_0 + \gamma_1 X_1 + \gamma_2(X_2-1))$, $\gamma=(0,1.5,1.0)$
  $X_4 \mid X_{1:3}$   Mixture Gamma   $\mathrm{Gamma}(k_1{=}3, r_1{=}1+0.5X_2)$ / $\mathrm{Gamma}(k_2{=}6, r_2{=}0.75+0.25X_2)$; $\tilde{w} = \operatorname{logistic}(\delta_0 + \delta_1 X_1 + \delta_3(X_3-0.5))$, $\delta=(0,1.0,3.0)$

  : Configuration for the four-dimensional autoregressive generator used
  in the synthetic study. The beta and gamma coordinates are
  two-component mixtures with logistic gates; fixed parameters and gates
  match the prose above.
:::

Mixture weights use the logistic gate
$\operatorname{logistic}(a)=\exp(a)/(1+\exp(a))$, which coincides with
the two-component softmax and therefore keeps probabilities in $(0,1)$
that sum to one. For completeness, the general softmax takes a vector
$a$ and returns $\mathrm{softmax}(a)_i = \exp(a_i) / \sum_j \exp(a_j)$.
This normalization is essential for the beta and gamma mixtures because
it translates linear predictors into valid probability weights while
preserving differentiability.

Table [3.1](#tab:autoregressive-config){reference-type="ref"
reference="tab:autoregressive-config"} (p. ) summarizes the mixture
families, gates, and fixed parameters by dimension.
Tables [3.4](#tab:autoregressive-perm){reference-type="ref"
reference="tab:autoregressive-perm"} (p. )
and [3.6](#tab:autoregressive-perm-avg){reference-type="ref"
reference="tab:autoregressive-perm-avg"} (p. ) then summarize the
permutation and sample-size studies used later in this chapter.

The MiniBooNE benchmark follows the published preprocessing to ensure
comparability with flow-based baselines. We remove 11 outliers with
value $-1000$, drop seven features with extreme mass at a single value,
and retain $K = 43$ attributes. We use the fixed train, validation, and
test splits from the benchmark, apply train-only standardization, and
avoid any extra pruning of correlated features. We report all
log-likelihoods in nats and retain the published naming for flow
comparators in later tables.
Section [3.5](#sec:realdata){reference-type="ref"
reference="sec:realdata"} records these steps and provides the dataset
context. Table [3.8](#tab:uci-loglik){reference-type="ref"
reference="tab:uci-loglik"} reproduces the flow baselines that motivate
our TRTF runs.

Additional UCI datasets appear only when we retain them for real-data
context. POWER keeps household electricity attributes after jittering
the minute-of-day encoding, dropping the calendar date and
reactive-power column, and adding uniform noise to break ties. GAS keeps
the `ethylene_CO` subset, treats the series as i.i.d., removes strongly
correlated attributes, and retains an eight-dimensional representation.
HEPMASS keeps only the positive class from the "1000" split and discards
five variables with repeated values to avoid density spikes. These
preprocessing steps follow the same train-only standardization and
reporting conventions described above.
Section [3.5](#sec:realdata){reference-type="ref"
reference="sec:realdata"} provides the corresponding background and
positions these datasets within our evaluation.

All models use the same standardized frame and direction for evaluation,
which keeps objectives, diagnostics, and reported quantities
interoperable across triangular transports, TRTF, and copula baselines.
This alignment is necessary for the conditional decompositions,
probability integral transform (PIT) calibration checks, and compute
summaries presented later in this chapter.

## Models and Implementation {#sec:models-implementation}

This section specifies the estimators and implementation details that
keep likelihoods, calibration, and compute directly comparable across
models. All estimators share the transport direction
$S:u \rightarrow z$, operate in standardized coordinates, and report log
densities on the original scale using the affine correction from
Chapter [2](#ch:background){reference-type="ref"
reference="ch:background"}.
Figure [5.1](#fig:transport-schematic){reference-type="ref"
reference="fig:transport-schematic"} in
Appendix [5](#ch:appendix){reference-type="ref" reference="ch:appendix"}
and Table [2.1](#tab:transport-notation){reference-type="ref"
reference="tab:transport-notation"} summarize the shared pipeline and
notation.

We implement separable lower-triangular transport maps denoted TTM-Sep.
Component $k$ decomposes into a context shift and a univariate monotone
shape,
$$S_k(u_{1:k}) = g_k(u_{1:k-1}) + h_k(u_k), \qquad \partial_{u_k} S_k(u_{1:k}) = h_k'(u_k) > 0,
  \label{eq:ttm-separable-def}$$ so the Jacobian contribution depends
only on $u_k$. The structure yields linear per-sample complexity in $K$
and exact inversion by back-substitution.

We minimize the Gaussian pullback objective induced by the shared
reference,
$$\mathcal{L}(u) = \sum_{k=1}^K \Big[ \tfrac{1}{2} S_k(u_{1:k})^2 - \log h_k'(u_k) \Big],
  \label{eq:ttm-separable-loss}$$ which follows from the
change-of-variables identity in
Equation [\[eq:transport-pullback\]](#eq:transport-pullback){reference-type="eqref"
reference="eq:transport-pullback"} combined with the triangular
determinant factorization in
Equation [\[eq:transport-det\]](#eq:transport-det){reference-type="eqref"
reference="eq:transport-det"}. The quadratic term pulls the transformed
coordinates toward the reference, and the log-derivative term prevents
degenerate solutions. We solve the regularized problem with
bound-constrained optimization and enforce monotone structure by
construction.

We construct $h_k$ with monotone one-dimensional bases that combine
identity, integrated sigmoids, softplus-like edge terms, and integrated
radial basis functions. Nonnegativity constraints on the derivative
coefficients guarantee $h_k'(u_k) \ge 0$. We linearize tails to
stabilize likelihoods as $|u_k|$ grows. Ridge penalties apply to all
basis coefficients, and optional sparsity penalties shrink context
shifts when multicollinearity inflates variance. During training and
evaluation we clip log-derivatives to $[-H, H]$ to avoid numerical
overflow in the Jacobian sum; the bound $H$ is tuned on the validation
split.

We build $g_k$ from low-degree polynomial features of $u_{1:k-1}$ and
drop predecessors whose inclusion does not improve validation
likelihood. This pruning keeps $\nabla_u S(u)$ sparse and improves
stability in small-sample regimes. Ordering matters for finite bases, so
headline results use the natural variable order while robustness studies
vary the order as described in
Section [3.1](#sec:datasets-preprocessing){reference-type="ref"
reference="sec:datasets-preprocessing"}. When heuristics are applied, we
evaluate two candidates on the validation split---identity and a
Cholesky-pivoted ordering with optional Gaussianization---and select the
ordering with the lower validation NLL.
Appendix [5](#ch:appendix){reference-type="ref" reference="ch:appendix"}
records how the ordering is stored and reapplied at prediction time.

We reference a cross-term variant, denoted TTM-X, only to delimit scope.
The variant augments the separable component with low-rank interactions,
$$S_k(u_{1:k}) = g_k(u_{1:k-1}) + h_k(u_k) + \sum_{j<k} \alpha_{kj}\,q_j(u_j)\,r_k(u_k),
  \label{eq:ttm-cross}$$ where $q_j$ and $r_k$ are monotone features and
constraints ensure $\partial_{u_k} S_k(u_{1:k}) > 0$. We exclude TTM-X
from headline tables because the interactions alter identifiability and
complicate calibration. The definition clarifies the naming used in the
synthetic analyses.

We implement Transformation Random Forests with the forest aggregation and
denote the model TRTF. Implementations rely on the `partykit` toolkit
for recursive partitioning [@partykit2015] together with the
transformation-model framework
[@hothorn2017transformation; @hothorn2018conditional]. Let
$\widehat{F}_k(\cdot \mid u_{1:k-1})$ denote the strictly increasing
conditional CDF returned by the forest. The induced triangular component
is
$$S_k(u_{1:k}) = \Phi^{-1}\!\big(\widehat{F}_k(u_k \mid u_{1:k-1})\big),
  \label{eq:trtf-transport}$$ and differentiation yields
$\varphi\big(S_k(u_{1:k})\big)\,\partial_{u_k} S_k(u_{1:k}) = \widehat{\pi}_k(u_k \mid u_{1:k-1})$.
Under the forest aggregation
$\widehat{F}_k(u_k \mid u_{1:k-1}) = \Phi\big(h_k(u_k) + g_k(u_{1:k-1})\big)$
we obtain $S_k = h_k + g_k$ and $\partial_{u_k} S_k = h_k'(u_k)$, which
matches
Equation [\[eq:ttm-separable-def\]](#eq:ttm-separable-def){reference-type="eqref"
reference="eq:ttm-separable-def"} exactly in the transport frame.
Consequently TRTF shares the likelihood in
Equation [\[eq:ttm-separable-loss\]](#eq:ttm-separable-loss){reference-type="eqref"
reference="eq:ttm-separable-loss"}, inherits exact inversion, and
differs operationally through forest training and aggregation.

We keep copulas as dependence baselines with explicit scope. We fit only
low-dimensional nonparametric copulas for $K \le 3$, using
probit-transformed pseudo-observations and kernel density estimation on
the Gaussian scale before mapping back to the unit cube with the
appropriate Jacobian. The independence baseline evaluates the product of
fitted marginals. We omit a Gaussian copula from the main experiments to
preserve consistency with the low-$K$ nonparametric dependence analyzed
in Section [3.5](#sec:realdata){reference-type="ref"
reference="sec:realdata"}.

We adopt a single inversion and evaluation convention across estimators.
Training, Jacobians, and conditional evaluations occur in standardized
coordinates. Sampling draws $z \sim \mathcal{N}(0, I)$, applies $S^{-1}$
by back-substitution, and converts to the original scale with the stored
affine parameters. This convention prevents mixed objectives and keeps
all reported quantities interoperable.

We ensure reproducibility and comparability with fixed seeds, cached
standardization parameters, and shared reporting utilities.
Appendix [5](#ch:appendix){reference-type="ref" reference="ch:appendix"}
provides pseudo-code for TRTF fitting and prediction, nonparametric
copulas, marginal and separable triangular maps, and the shared
transport core that implements ordering, bases, derivatives, and
evaluation. The appendix also records optimizer choices, timing hooks,
and object layouts used in the experiments.

## Evaluation Metrics and Protocol {#sec:evaluation-protocol}

This section defines the metrics and procedures applied across all
models so likelihoods, calibration, and compute remain directly
comparable. We evaluate every estimator in standardized coordinates,
apply the triangular determinant, and report log densities on the
original scale using the affine correction from
Chapter [2](#ch:background){reference-type="ref"
reference="ch:background"}.
Figure [5.1](#fig:transport-schematic){reference-type="ref"
reference="fig:transport-schematic"} in
Appendix [5](#ch:appendix){reference-type="ref" reference="ch:appendix"}
and Table [2.1](#tab:transport-notation){reference-type="ref"
reference="tab:transport-notation"} summarize the shared pipeline and
notation.

We distinguish pointwise log density from dataset averages. The test log
likelihood (LL) is the mean of $\log \hat{\pi}_X(x)$ over the test
split, and the test negative log likelihood (NLL) is its negative.
Tables note "LL (higher is better)" or "NLL (lower is better)" to avoid
ambiguity. Reported log densities on the original scale equal the
standardized quantity minus $\sum_k \log \sigma_k$ as given by
Equation [\[eq:transport-affine\]](#eq:transport-affine){reference-type="eqref"
reference="eq:transport-affine"}. Consequently, datasets with larger
training scales introduce large constant offsets that the affine
correction removes.

Triangular models exploit the separable pullback in standardized
coordinates. With $u = T_{\mathrm{std}}(x)$, the log density decomposes
as
$$\log \hat{\pi}_U(u) = \sum_{k=1}^{K} \Big[ \log \varphi\big(S_k(u_{1:k})\big) + \log \partial_{u_k} S_k(u_{1:k}) \Big],
  \label{eq:evaluation-triangular}$$ so the determinant factorization in
Equation [\[eq:transport-det\]](#eq:transport-det){reference-type="eqref"
reference="eq:transport-det"} yields linear per-sample cost in $K$. In
plain language, the model checks how Gaussian each transformed
coordinate looks, then adds the exact log-Jacobian contribution from its
one-dimensional derivative. The affine correction in
Equation [\[eq:transport-affine\]](#eq:transport-affine){reference-type="eqref"
reference="eq:transport-affine"} converts $\log \hat{\pi}_U$ to
$\log \hat{\pi}_X$ for reporting.

We report per-dimension conditional NLLs for triangular models to
localize error. For each coordinate,
$$\mathrm{NLL}_k = -\frac{1}{N_{\mathrm{test}}} \sum_{i=1}^{N_{\mathrm{test}}} \log \hat{\pi}\big(x_{ik} \mid x_{i,1:k-1}\big),
  \label{eq:evaluation-conditional-nll}$$ and the joint NLL equals
$\sum_{k=1}^{K} \mathrm{NLL}_k$ by construction. Copulas lack a unique
triangular factorization, so we report only their joint NLL. Negative
per-dimension NLL values can occur because valid densities may exceed
one on subdomains. These conventions align with
Equations [\[eq:transport-trtf-likelihood\]](#eq:transport-trtf-likelihood){reference-type="eqref"
reference="eq:transport-trtf-likelihood"}
and [\[eq:transport-trtf-separable\]](#eq:transport-trtf-separable){reference-type="eqref"
reference="eq:transport-trtf-separable"}, which link separable
transports and Transformation Random Forests inside the common frame
under the additive-predictor and monotone-smoothing assumptions stated
in Section [2.2.1](#sec:transport-trtf){reference-type="ref"
reference="sec:transport-trtf"}.

Calibration assesses whether predictive probabilities align with
empirical frequencies. For triangular models we form conditional
probability integral transform (PIT) values
$V_{ik} = \widehat{F}_k(u_{ik} \mid u_{i,1:k-1})$ on the test split and
expect independent $\mathrm{Unif}(0,1)$ draws under correct calibration
[@gneiting2007probabilistic]. We summarize departures from uniformity
with the Kolmogorov--Smirnov statistic
$D_n = \sup_t \lvert \widehat{F}_n(t) - t \rvert$ and report the
associated $p$-value [@massey1951kolmogorov]. We complement the scalar
summary with brief PIT distribution descriptions when patterns recur
across seeds. For copulas we assess marginal PITs and low-dimensional
slices where dependence is transparent. Systematic U-shaped or
inverted-U PIT indicates under- or over-dispersion and motivates richer
parameterizations.

Compute metrics quantify practical cost alongside fit. We record
wall-clock training time on the training split and per-sample evaluation
time on the test split. Triangular transports scale linearly in $K$ and
approximately linearly in the number of basis functions. Transformation
Random Forests scale with the number and depth of trees per conditional
during training, while prediction remains linear after aggregation.
Copula training is dominated by correlation estimation or kernel density
fitting, followed by fast evaluation. We also track peak resident memory
when caching affects runtime. All timings use the deterministic pipeline
defined in Chapter [3](#ch:dataanalysis){reference-type="ref"
reference="ch:dataanalysis"} and are averaged over seeds with standard
errors.

Protocol choices keep comparisons stable and reproducible:

1.  Standardize with training-split statistics, fit a single map
    $S:u \rightarrow z$, and evaluate Jacobians in standardized space.

2.  Compute LL, NLL, and conditional decompositions in standardized
    coordinates, then apply the affine correction once for reporting.

3.  Evaluate PIT diagnostics, Kolmogorov--Smirnov statistics, and
    compute metrics on the fixed test split with seeds
    $\{11, 13, 17, 19, 23\}$ and quote means with $\pm$ two standard
    errors across seeds (SE $= s/\sqrt{m}$ over $m$ seeds).

Appendix [5](#ch:appendix){reference-type="ref" reference="ch:appendix"}
lists routine interfaces that support exact re-execution; figure
captions and table notes repeat units, splits, and the "higher/lower is
better" convention for clarity.

Numerical safeguards prevent unstable likelihoods from dominating
summaries. We enforce strictly monotone derivatives by construction and
clip log-derivative contributions to $[-H, H]$ during training and
evaluation. We tune $H$ and regularization on the validation split and
reuse the selected bound on the test split. We report clipping status
inline per study and keep the exact bound values with the corresponding
experiment logs to avoid duplication. This practice controls overflow in
the Jacobian sum without masking systematic misfit that PIT diagnostics
would reveal.

Scope limits clarify non-goals. We do not report AIC or BIC because
effective parameter counts are not comparable across estimators in this
frame. We also do not adjust $p$-values for multiple PIT checks;
instead, we treat Kolmogorov--Smirnov results as diagnostics and
corroborate them with effect sizes and plots. These limits keep the
evaluation focused on likelihood, calibration, and compute under a
single, transparent protocol.

## Synthetic Results and Diagnostics {#sec:synthetic-results}

This section reports synthetic results for the Half-Moon and
four-dimensional generators under the protocol in
Section [3.3](#sec:evaluation-protocol){reference-type="ref"
reference="sec:evaluation-protocol"}. We summarize mean test negative
log likelihoods, per-dimension conditional NLLs, calibration evidence,
and ordering robustness, referencing the corresponding tables and
figures.

The Half-Moon generator stresses conditional shape in two dimensions.
Table [3.2](#tab:halfmoon-nll){reference-type="ref"
reference="tab:halfmoon-nll"} lists mean joint NLLs with $\pm$ two
standard errors: TRTF achieved $1.71 \pm 0.09$ nats, TTM-Sep achieved
$1.93 \pm 0.08$ nats, and TTM-Marg achieved $2.02 \pm 0.07$ nats. The
copula baseline reached $1.54 \pm 0.09$ nats and bracketed the
triangular transports. The oracle references set $0.78 \pm 0.10$ nats
for the true marginal density and $0.70 \pm 0.12$ nats for the true
joint. Per-dimension NLLs confirm that the first coordinate is harder:
TRTF reported $(1.23, 0.47)$, while TTM-Sep reported $(1.28, 0.65)$.
Figure [3.1](#fig:halfmoon-panels){reference-type="ref"
reference="fig:halfmoon-panels"} shows contours consistent with these
rankings and with the standardized pipeline in
Figure [5.1](#fig:transport-schematic){reference-type="ref"
reference="fig:transport-schematic"} of
Appendix [5](#ch:appendix){reference-type="ref"
reference="ch:appendix"}.Clipping status: not triggered in these runs
(no log-derivative terms reached the bound).

*(mean NLL in nats).*

::: {#tab:halfmoon-nll}
  Model         Mean joint NLL    Conditional NLL 1   Conditional NLL 2
  ------------ ----------------- ------------------- -------------------
  True-Marg     $0.78 \pm 0.10$        $0.39$              $0.39$
  True-Joint    $0.70 \pm 0.12$        $0.35$              $0.35$
  TRTF          $1.71 \pm 0.09$        $1.23$              $0.47$
  TTM-Marg      $2.02 \pm 0.07$        $1.28$              $0.74$
  TTM-Sep       $1.93 \pm 0.08$        $1.28$              $0.65$
  Copula        $1.54 \pm 0.09$        $0.77$              $0.77$

  : Half-Moon ($n=250$): mean test negative log-likelihood (NLL; nats;
  lower is better). Values are means $\pm$ 2SE.
:::

![Half-Moon ($n=250$) log-density contours for the true joint, TRTF, TTM
variants, and the copula mixture. Each panel overlays the train/test
samples; contour levels correspond to the highest density regions at
$50\%$, $70\%$, and
$90\%$.](figure/halfmoon_panels_seed007_N250.png){#fig:halfmoon-panels
width="85%"}

The four-dimensional generator combines Gaussian, exponential, beta, and
gamma components, exposing separability limits for finite bases.
Table [3.3](#tab:autoregressive-nll){reference-type="ref"
reference="tab:autoregressive-nll"} (p. ) reports the canonical ordering
$(1,2,3,4)$. TRTF aligned closely with the exponential coordinate,
recording $1.51$ nats compared with $1.49$ for the true joint reference.
TTM-Sep over-penalized that coordinate at $1.88$ nats, and TTM-Marg
overfit at $2.57$ nats. The beta coordinate yielded negative NLLs for
the oracles because valid densities can exceed one on $(0,1)$; values
were $-0.79$ for the true joint and $-0.48$ for the true marginal. TRTF
reached $-0.25$, while TTM-Sep and the copula baseline reported $0.07$
and $0.05$ nats, respectively. The gamma coordinate remained most
challenging, with $1.99$ nats for TRTF and $2.41$ nats for TTM-Sep.
Joint sums were $4.53$ nats for TRTF, $5.66$ nats for TTM-Sep, $6.83$
nats for TTM-Marg, and $5.45$ nats for the copula, compared with $3.80$
nats for the true joint oracle.
Figure [3.2](#fig:autoregressive-joint-calibration){reference-type="ref"
reference="fig:autoregressive-joint-calibration"} (p. ) compares
predicted and true joint log densities, highlighting calibration gaps
relative to the identity line.Clipping status: not triggered at $n=250$
under the selected configuration (see Appendix
Table [5.1](#tab:ttmsep-n25-overflow){reference-type="ref"
reference="tab:ttmsep-n25-overflow"} for the small-sample $n=25$ edge
case).

![Four-dimensional autoregressive generator ($n=250$): joint log-density
calibration for each estimator (axes in nats). Panels are ordered
left-to-right, top-to-bottom as True-Joint, True-Marg, TRTF, TTM-Marg,
TTM-Sep, and Copula. Gray dots mark the $20\%$ test split (50 samples).
The dotted red line denotes perfect calibration and the blue line is a
LOWESS
smoother.](figure/logdensity_joint_N250.png){#fig:autoregressive-joint-calibration
width="85%"}

*(mean NLL in nats).*

::: {#tab:autoregressive-nll}
  Dim   Distribution     True-Marg   True-Joint      TRTF   TTM-Marg   TTM-Sep   Copula
  ----- -------------- ----------- ------------ --------- ---------- --------- --------
  1     Normal              $1.29$       $1.28$    $1.28$     $1.29$    $1.29$   $1.30$
  2     Exponential         $1.75$       $1.49$    $1.51$     $2.57$    $1.88$   $1.87$
  3     Beta               $-0.48$      $-0.79$   $-0.25$     $0.28$    $0.07$   $0.05$
  4     Gamma               $2.05$       $1.83$    $1.99$     $2.69$    $2.41$   $2.22$
  $K$   Sum (joint)         $4.61$       $3.80$    $4.53$     $6.83$    $5.66$   $5.45$

  : Four-dimensional autoregressive generator ($n=250$, permutation
  $1,2,3,4$): mean conditional and joint NLL (nats; lower is better).
  Values are means over test samples (no SE shown).
:::

Ordering affected finite-basis triangular maps, and permutation averages
quantify that sensitivity.
Table [3.4](#tab:autoregressive-perm){reference-type="ref"
reference="tab:autoregressive-perm"} (p. ) summarizes test NLLs over all
$4! = 24$ permutations: TRTF averaged $4.65$ nats, TTM-Sep averaged
$5.62$ nats, TTM-Marg averaged $6.83$ nats, and the copula baseline
averaged $5.45$ nats. The joint and marginal oracles remained stable at
$3.80$ and $4.61$ nats, respectively. These effects confirm anisotropy
and motivate the ordering heuristics described in
Section [3.2](#sec:models-implementation){reference-type="ref"
reference="sec:models-implementation"} when bases are finite. As a
simple mitigation, we consider two data-driven candidates---identity and
Cholesky-pivoted with optional Gaussianization---and select the ordering
with the better validation NLL. Appendix
Figure [5.2](#fig:ordering-heuristics-4d){reference-type="ref"
reference="fig:ordering-heuristics-4d"} visualizes the potential
improvement window by marking the canonical, median, and
best-over-permutations joint NLLs for TRTF and TTM-Sep at $n=250$.

*(mean NLL in nats).*

::: {#tab:autoregressive-perm}
  Model          Dim 1   Dim 2   Dim 3   Dim 4    Sum
  ------------ ------- ------- ------- ------- ------
  True-Marg       1.22    1.13    1.15    1.11   4.61
  True-Joint      1.03    0.93    0.94    0.91   3.80
  TRTF            1.33    1.19    1.09    1.04   4.65
  TTM-Marg        1.77    1.67    1.73    1.66   6.83
  TTM-Sep         1.59    1.38    1.36    1.29   5.62
  Copula          1.42    1.34    1.36    1.32   5.45

  : Four-dimensional autoregressive generator ($n=250$): mean test NLL
  (nats; lower is better) averaged over all $24$ permutations of
  $(1,2,3,4)$.
:::

*(mean NLL in nats).*

::: {#tab:autoregressive-perm-spread}
  Model           Min   Median    Max
  ------------ ------ -------- ------
  True-Marg      4.61     4.61   4.61
  True-Joint     3.80     3.80   3.80
  TRTF           4.46     4.59   5.23
  TTM-Marg       6.83     6.83   6.83
  TTM-Sep        5.48     5.60   5.78
  Copula         5.45     5.45   5.45

  : Permutation spread of joint NLLs (nats) over all $24$ permutations
  for $n=250$. Values report $\min/\mathrm{median}/\max$ across
  orderings (lower is better).
:::

Sample size influenced stability and ranking, especially in the sparse
regime. Table [3.6](#tab:autoregressive-perm-avg){reference-type="ref"
reference="tab:autoregressive-perm-avg"} (p. ) aggregates joint NLLs
across permutations for $n \in \{25, 50, 100, 250\}$. TRTF decreased
from $38.18$ to $4.64$ nats as $n$ increased, while TTM-Sep decreased
from $6.35$ to $5.61$ nats across the stable regimes. The TTM-Sep result
at $n=25$ exhibited numerical overflow and is reported in Appendix
Table [5.1](#tab:ttmsep-n25-overflow){reference-type="ref"
reference="tab:ttmsep-n25-overflow"} marked with an asterisk ($^{\ast}$)
as out of scope; it is excluded from main-text comparisons. The copula
decreased from $9.02$ to $5.45$ nats and tracked TTM-Sep once
$n \ge 100$.

*(mean NLL in nats).*

::: {#tab:autoregressive-perm-avg}
  Model          $n=25$   $n=50$   $n=100$   $n=250$
  ------------ -------- -------- --------- ---------
  True-Marg       10.50     4.75      4.91      4.61
  True-Joint       4.35     4.23      3.55      3.80
  TRTF            38.18     6.10      4.59      4.64
  TTM-Marg        49.36     7.43      7.72      6.83
  TTM-Sep            --     6.35      6.08      5.61
  Copula           9.02     6.66      6.02      5.45

  : Note: The TTM-Sep entry at $n=25$ is omitted from the main table due
  to numerical overflow; see Appendix
  Table [5.1](#tab:ttmsep-n25-overflow){reference-type="ref"
  reference="tab:ttmsep-n25-overflow"}, where it is marked with an
  asterisk ($^{\ast}$) as out of scope.
:::

Calibration assessments align with the likelihood evidence.
Figure [3.2](#fig:autoregressive-joint-calibration){reference-type="ref"
reference="fig:autoregressive-joint-calibration"} (p. ) shows joint
log-density calibration against the oracle, with residual structure
visible for triangular transports in the canonical ordering. Conditional
PIT diagnostics and Kolmogorov--Smirnov distances, computed as in
Section [3.3](#sec:evaluation-protocol){reference-type="ref"
reference="sec:evaluation-protocol"}, exhibited the same qualitative
patterns across seeds, so we omit redundant tables.

These studies indicate that TRTF closes part of the gap to oracle
likelihoods while preserving the triangular evaluation frame. Separable
maps remain competitive at moderate sample sizes but exhibit ordering
sensitivity and sparse-regime fragility, and copulas provide competitive
baselines in low dimensions.
Section [3.5](#sec:realdata){reference-type="ref"
reference="sec:realdata"} turns to real-data benchmarks and compute
summaries under the same protocol.

##### Calibration numbers.

To complement the visual diagnostics,
Table [3.7](#tab:ks-synth){reference-type="ref"
reference="tab:ks-synth"} reports median Kolmogorov--Smirnov (KS)
distances of probability-integral-transform (PIT) values per coordinate,
aggregated across seeds ($\pm$ 2SE). Lower is better. We include entries
for methods with an accessible marginal CDF in our implementation.

*(median KS of PIT per coordinate; lower is better).* 

::: {#tab:ks-synth}
  -------------- -------------- ------------------- --------------- ------------------- ------------------- -------------------
    Dataset        True (Joint)     True (marginal)   Random Forest        Marginal Map       Separable Map           Copula NP
   Half-Moon                 --   0.079 $\pm$ 0.015     NA $\pm$ NA   0.078 $\pm$ 0.015   0.090 $\pm$ 0.018   0.067 $\pm$ 0.001
  4D generator               --                  --              --                  --                  --                  --
                                                                                                            
  -------------- -------------- ------------------- --------------- ------------------- ------------------- -------------------

  : Calibration via PIT--KS on synthetic datasets: median KS distance
  per coordinate (mean $\pm$ 2SE across seeds). Entries marked '--'
  indicate that the CDF was not available in the corresponding backend.
:::

## Real-Data Benchmarks and Compute {#sec:realdata}

This section presents real-data evidence on MiniBooNE and the UCI
tabular benchmarks under the transport frame introduced in
Chapters [1](#ch:intro){reference-type="ref" reference="ch:intro"}
and [2](#ch:background){reference-type="ref" reference="ch:background"}.
We keep preprocessing identical to the published flow literature where
applicable, align likelihood reporting through standardized coordinates
and the affine correction in
Equation [\[eq:transport-affine\]](#eq:transport-affine){reference-type="eqref"
reference="eq:transport-affine"}, and pair test log likelihoods with
compute summaries so that score differences reflect modeling assumptions
rather than inconsistent units.

##### Preprocessing.

We treat dataset-specific preprocessing as part of each estimator to
preserve comparability. MiniBooNE follows @papamakarios2017masked: we
remove $11$ outliers at $-1000$, drop $7$ near-constant attributes,
retain $K=43$ variables, and rely on the official train, validation, and
test splits. We standardize with training statistics only, evaluate
Jacobians in standardized coordinates, and apply the diagonal affine
correction once at reporting time. The UCI datasets follow the same
rule. POWER receives jitter on the minute-of-day encoding, removal of
the calendar-date and reactive-power attributes, and a small uniform
perturbation to break ties. GAS keeps the `ethylene_CO` subset and
removes strongly correlated attributes to yield an eight-dimensional
representation. HEPMASS keeps the positive class from the "1000" split
and discards five repeated-value variables to avoid density spikes.
These steps match the literature conventions and keep the reported
likelihoods interpretable.

##### Flow baselines.

Published normalizing flows compose invertible layers with permutations
or autoregressive sublayers and report strong test log likelihoods on
the UCI suite and MiniBooNE
[@rezende2015variational; @dinh2017real; @kingma2018glow; @durkan2019neural; @papamakarios2021normalizing].
Table [3.8](#tab:uci-loglik){reference-type="ref"
reference="tab:uci-loglik"} reproduces the published average test
log-likelihoods per example together with $\pm$ two standard errors
reported by @papamakarios2017masked and appends our TRTF measurements
trained with $N=2500$ observations. Higher values indicate better fits.
We report TRTF as means $\pm$ 2SE under the same evaluation pipeline.

*(average LL; nats per example).*

::: {#tab:uci-loglik}
  Model                        POWER                GAS             HEPMASS           MiniBooNE
  --------------- ------------------ ------------------ ------------------- -------------------
  Gaussian          $-7.74 \pm 0.02$   $-3.58 \pm 0.75$   $-27.93 \pm 0.02$   $-37.24 \pm 1.07$
  MADE              $-3.08 \pm 0.03$    $3.56 \pm 0.04$   $-20.98 \pm 0.02$   $-15.59 \pm 0.50$
  MADE MoG           $0.40 \pm 0.01$    $8.47 \pm 0.02$   $-15.15 \pm 0.02$   $-12.27 \pm 0.47$
  Real NVP (5)      $-0.02 \pm 0.01$    $4.78 \pm 1.80$   $-19.62 \pm 0.02$   $-13.55 \pm 0.49$
  Real NVP (10)      $0.17 \pm 0.01$    $8.33 \pm 0.14$   $-18.71 \pm 0.02$   $-13.84 \pm 0.52$
  MAF (5)            $0.14 \pm 0.01$    $9.07 \pm 0.02$   $-17.70 \pm 0.02$   $-11.75 \pm 0.44$
  MAF MoG (5)        $0.30 \pm 0.01$    $9.59 \pm 0.02$   $-17.39 \pm 0.02$   $-11.68 \pm 0.44$
  TRTF (ours)       $-7.17 \pm 0.39$   $-2.41 \pm 0.37$   $-25.47 \pm 0.37$   $-30.01 \pm 1.26$

  : UCI: average test log-likelihood per example (nats; higher is
  better). Baselines (first seven rows): means $\pm$ 2SE as reported by
  @papamakarios2017masked. TRTF (ours): single-seed measurements at
  $N=2500$ (no SE). Entries marked "--" indicate configurations not
  executed in this draft.
:::

##### MiniBooNE.

Table [3.8](#tab:uci-loglik){reference-type="ref"
reference="tab:uci-loglik"} shows that the Gaussian reference yields
$-37.24 \pm 1.07$ nats, providing a weak baseline. MADE reaches
$-15.59 \pm 0.50$ nats, the Real NVP variants lie near $-13.7$ nats, and
MAF MoG improves to $-11.68 \pm 0.44$ nats. Our TRTF result attains
$-30.01 \pm 1.26$ nats at $N=2500$, improving over the Gaussian baseline
yet trailing the flow families by a wide margin. This ranking is
consistent with the separable Jacobian and the forest aggregation discussed
in Section [3.2](#sec:models-implementation){reference-type="ref"
reference="sec:models-implementation"}. The high dimensionality of
MiniBooNE amplifies residual misfit through the triangular
determinant.Clipping: validation-tuned bound $H$ applied; the exact
value is recorded with the experiment logs.

##### POWER.

POWER offers a milder conditional structure and lower dimensionality.
Table [3.8](#tab:uci-loglik){reference-type="ref"
reference="tab:uci-loglik"} reports that TRTF records $-7.17 \pm 0.39$
nats at $N=2500$, which falls short of the flow baselines. Real NVP with
ten steps reaches $0.17 \pm 0.01$ nats, while MAF MoG attains
$0.30 \pm 0.01$ nats. The gap indicates that the current TRTF
configuration underutilizes structure in this benchmark; additional
seeds or hyperparameter tuning may recover the performance previously
observed at smaller sample sizes.Clipping: validation-tuned bound $H$
applied; the exact value is recorded with the experiment logs.

##### GAS and HEPMASS.

The TRTF results on GAS and HEPMASS yield $-2.41 \pm 0.37$ and
$-25.47 \pm 0.37$ nats, respectively. Both scores remain below the flow
baselines, emphasizing that the present configuration sacrifices
likelihood accuracy for interpretability. Additional seeds and tuning
remain planned, yet we retain the current numbers to document the
outcome of the standardized pipeline at $N=2500$.Clipping:
validation-tuned bound $H$ applied; the exact values are recorded with
the experiment logs.

##### Sample size sensitivity.

Figure [3.3](#fig:n-sensitivity){reference-type="ref"
reference="fig:n-sensitivity"} plots test negative log likelihood versus
sample size $N$ for the UCI benchmarks, aggregating seeds at each
budget. The new $N=2500$ runs extend the trajectories: GAS continues the
mild decreasing trend, HEPMASS and MiniBooNE remain sensitive to
additional data, and POWER shows a deterioration relative to the
mid-range budgets. The figure reports one standard error bars (zero when
only a single seed is available), restates that lower curves indicate
better fits because the vertical axis plots NLL, and mirrors the
diagnostic procedures in
Section [3.3](#sec:evaluation-protocol){reference-type="ref"
reference="sec:evaluation-protocol"}.

![Test negative log-likelihood (NLL; nats; lower is better) versus
sample size $N$ on the UCI benchmarks. Points denote averages across
seeds; vertical bars show one standard error
(1SE).](figure/N_sensitivity_all.png){#fig:n-sensitivity width="85%"}

##### Compute metrics.

Likelihood comparisons require compute summaries because similar
accuracy at very different costs leads to different recommendations.
Training time is wall-clock time to fit the model on the training split
with fixed seeds and deterministic preprocessing. Evaluation time is the
wall-clock time per $10^5$ joint log-density evaluations on the test
split, averaged over seeds. These definitions mirror the compute
discussion in
Section [3.3](#sec:evaluation-protocol){reference-type="ref"
reference="sec:evaluation-protocol"}, use the same standardized inputs
across datasets, and yield the budget-specific totals collected in
Table [3.9](#tab:real-compute){reference-type="ref"
reference="tab:real-compute"}.

::: {#tab:real-compute}
  Dataset       $N=25$   $N=50$   $N=100$   $N=250$   $N=500$   $N=1000$   $N=2500$
  ----------- -------- -------- --------- --------- --------- ---------- ----------
  POWER            $1$      $1$       $2$       $6$      $39$      $115$      $130$
  GAS              $1$      $1$       $2$       $5$      $39$      $138$      $600$
  HEPMASS          $1$      $2$       $4$       $9$      $12$      $153$      $721$
  MiniBooNE        $3$      $4$       $8$      $20$      $27$      $202$     $2007$

  : TRTF wall-clock training plus evaluation time (seconds) as a
  function of the training budget $N$. Runs use the standardized inputs,
  seeds, and transport direction shared across datasets. Dashes denote
  configurations that were not executed in the current draft.
:::

##### Interpretation.

The real-data evidence aligns with the synthetic diagnostics in
Section [3.4](#sec:synthetic-results){reference-type="ref"
reference="sec:synthetic-results"}. MiniBooNE exposes the limits of
separable structure in high dimensions, and the updated POWER value
shows that the present TRTF configuration no longer matches flow
baselines once the training budget increases to $N=2500$. GAS and
HEPMASS also trail the published flows, illustrating that
interpretability and exact inversion come at a likelihood cost under the
current hyperparameters.
Table [3.9](#tab:real-compute){reference-type="ref"
reference="tab:real-compute"} documents the corresponding compute
budgets and confirms the anticipated near-linear growth in wall-clock
time.

## Reproducibility {#sec:reproducibility}

We avoid AIC or BIC because effective parameter counts differ across
estimators, and we do not treat small likelihood differences as
practically significant when $\pm 2$ SE intervals overlap. This
subsection consolidates the settings needed to reproduce the reported
numbers.

\- Data splits and direction - Synthetic: fixed train/validation/test
proportions $0.60/0.20/0.20$; evaluations use the shared direction
$S:u\to z$ in standardized coordinates and apply the diagonal affine
correction once for reporting. - Real data: use official splits where
provided (MiniBooNE) and the same standardized evaluation pipeline;
otherwise adopt the same $0.60/0.20/0.20$ convention.

\- Seeds - Synthetic generators and model fits: seeds
$\{11,13,17,19,23\}$ across repeats; permutation studies average over
all $4!=24$ orderings in the 4D case. - Real data (UCI + MiniBooNE):
single-seed runs with seed $42$ for training/evaluation in this draft.

\- Standardization and evaluation - Standardize features with
training-split $(\mu,\sigma)$ only; compute all derivatives/Jacobians in
$u$; report on $x$ via the affine correction in
Eq. [\[eq:transport-affine\]](#eq:transport-affine){reference-type="eqref"
reference="eq:transport-affine"}. - TRTF uses forest aggregation and
monotone CDF smoothing so that the induced likelihood matches the
separable triangular form
(Sec. [2.2.1](#sec:transport-trtf){reference-type="ref"
reference="sec:transport-trtf"}).

\- Hyperparameters and tuning - TTM-Sep: monotone one-dimensional bases
for $h_k$ (identity, integrated sigmoids, softplus-like edge terms,
integrated RBFs); low-degree polynomial features for $g_k$; ridge
regularization on all coefficients; log-derivative clipping to $[-H,H]$
(bound $H$ tuned on validation). Degree and penalty strengths are
selected by validation; ordering is fixed to the natural order in
headline tables and varied in robustness checks. - TTM-Sep: monotone
one-dimensional bases for $h_k$ (identity, integrated sigmoids,
softplus-like edge terms, integrated RBFs); low-degree polynomial
features for $g_k$; ridge regularization on all coefficients;
log-derivative clipping to $[-H,H]$ (bound $H$ tuned on validation).
Degree and penalty strengths are selected by validation; ordering is
fixed to the natural order in headline tables and, when heuristics are
enabled, chosen as the better of identity vs. Cholesky-pivoted (with
optional Gaussianization) according to validation NLL. - TRTF: additive
predictor with forest aggregation; strictly increasing conditional CDFs
after standard monotone smoothing; remaining fit options follow package
defaults unless stated; we record the number of trees, depth and split
rules in the experiment logs. - Copulas (diagnostics only for $K\le 3$):
probit pseudo-observations and kernel density copula via `kdecopula`
with default bandwidth selection; independence and Gaussian baselines
are used only for reference in text where noted. - Exact choices (e.g.,
basis sizes, ridge penalties, selected $H$) are captured alongside each
run in the experiment logs and summarized inline where relevant; we
avoid duplicate tables in the PDF.

Final safeguard settings used for the reported results. For Half-Moon
($n=250$) and 4D ($n=250$), TTM-Sep used degree$_g=2$, ridge
$\lambda=0$, and no log-derivative clipping was activated (no terms hit
the bound). The $n=25$ 4D case overflowed under $\lambda=0$; reruns with
$\lambda>0$ and tighter $H$ removed the failure but are omitted as out
of scope. Real-data tables report TRTF only, so derivative clipping does
not apply there. Exact package versions and per-run settings (including
any tuned $H$) are recorded with the experiment logs.

\- Software and hardware - R with packages: `tram`, `trtf`, `partykit`,
`mlt`, `dplyr`, `parallel`, and `knitr`/LaTeX for the report. We record
package versions via `sessionInfo()` in run logs. - Single-threaded BLAS
by default; optional parallel training for TRTF via
`options(trtf.train_cores=4)` when available. - CPU-only runs on a
laptop-class machine; logs include hardware notes (CPU model, RAM) and
wall-clock timings (Table [3.9](#tab:real-compute){reference-type="ref"
reference="tab:real-compute"}).

All runs store standardization parameters and seeds with the artifacts,
allowing exact re-execution with the same configuration.
Appendix [5](#ch:appendix){reference-type="ref" reference="ch:appendix"}
provides routine interfaces and object layouts to support this.

##### Bridge to Chapter [4](#ch:conclusion){reference-type="ref" reference="ch:conclusion"}. {#bridge-to-chapter-chconclusion.}

The real-data study closes
Chapter [3](#ch:dataanalysis){reference-type="ref"
reference="ch:dataanalysis"} by positioning separable triangular
transports and TRTF within the UCI and MiniBooNE landscape. TRTF offers
exact inversion, linear evaluation, and transparent conditional
structure, yet trails modern flows on MiniBooNE.
Chapter [4](#ch:conclusion){reference-type="ref"
reference="ch:conclusion"} interprets these trade-offs and distills
guidance for practitioners choosing between separable transports,
transformation forests, and copula baselines on tabular data.

# Interpretation and Conclusion {#ch:conclusion}

This chapter synthesizes the empirical evidence gathered in
Chapter [3](#ch:dataanalysis){reference-type="ref"
reference="ch:dataanalysis"}, interprets the behavior of the estimators
within the unified transport frame, and prepares the concluding guidance
that follows. We retain the shared preprocessing, likelihood
conventions, and diagnostic procedures so that numerical comparisons
remain meaningful across synthetic and real datasets. Copulas enter our
study only as low-dimensional ($K\!\le\!3$) diagnostic baselines (e.g.,
Half-Moon, 4D) and are not evaluated on high-$K$ datasets.

## Interpretation of Results {#sec:interpretation-results}

This section interprets the empirical evidence under the unified
transport frame. We focus on TRTF, TTM-Sep, and,
where applicable, copula baselines (only for $K\!\le\!3$) evaluated with
matched preprocessing, metrics, and units. Synthetic studies report NLL,
real datasets report LL, and we apply the shared affine correction.
These commitments keep objectives, diagnostics, and compute
interoperable across estimators.

TRTF often leads within the separable family because the forest aggregation
shifts conditional location while the underlying monotone shapes remain
stable. The likelihood identities equate TRTF with separable triangular
maps, so observed gaps arise from how each estimator realizes context
shifts and stabilizes derivatives. On Half-Moon ($K=2$), TRTF achieved
an NLL of $1.71$ while TTM-Sep reached $1.93$, and the first coordinate
remained the main source of residual error.
Table [3.2](#tab:halfmoon-nll){reference-type="ref"
reference="tab:halfmoon-nll"} records the per-dimension decomposition
and associated uncertainty bands, showing that location adjustments
dominate the remaining discrepancies when separability holds
approximately in low dimensions.

The four-dimensional generator sharpens this interpretation by isolating
coordinates with different conditional structure. TRTF matched the
exponential coordinate with an NLL of $1.51$ compared with $1.49$ for
the oracle, whereas TTM-Sep over-penalized that coordinate. The beta
coordinate produced negative NLLs for the oracles because valid
densities can exceed one on $(0,1)$; TRTF approached those values at
$-0.25$. The gamma coordinate remained the most challenging, with TRTF
at $1.99$ and TTM-Sep at $2.41$. Joint sums favored TRTF at $4.53$
versus $5.66$, consistent with concentrated gains on location-dominated
coordinates. Table [3.3](#tab:autoregressive-nll){reference-type="ref"
reference="tab:autoregressive-nll"} lists these values, and
Figure [3.2](#fig:autoregressive-joint-calibration){reference-type="ref"
reference="fig:autoregressive-joint-calibration"} visualizes the
residual curvature relative to the identity line.

These comparisons reveal where separability fails to adapt to
context-dependent shape. Under a separable map, conditional variance,
skewness, and modality remain fixed after the location shift.
Probability-integral-transform diagnostics display U-shaped or
inverted-U patterns when dispersion misaligns, indicating under- or
over-dispersion rather than pure location error. The calibration plots
corroborate the per-dimension NLLs and localize remaining structure to
the beta and gamma coordinates, where separability is least appropriate.
Figure [3.2](#fig:autoregressive-joint-calibration){reference-type="ref"
reference="fig:autoregressive-joint-calibration"} summarizes these
deviations under the canonical ordering.

Ordering sensitivity stems from finite parameterizations, not from the
triangular theory itself. A Knothe--Rosenblatt rearrangement exists for
any order, yet limited bases introduce anisotropy that affects fit.
Averaging over all $24$ permutations yielded joint NLLs of $4.65$ for
TRTF and $5.62$ for TTM-Sep, leaving a $0.97$ nat gap that persisted
despite order changes, while the copula baseline averaged $5.45$.
Table [3.4](#tab:autoregressive-perm){reference-type="ref"
reference="tab:autoregressive-perm"} consolidates these
permutation-averaged results and underlines the value of data-driven
orderings when available.

Small-sample regimes amplified numerical fragility through the
log-Jacobian accumulation. TRTF decreased from $38.18$ to $4.64$ joint
NLL as $n$ grew from $25$ to $250$, reflecting stabilization with
additional data. TTM-Sep spiked to $6{,}829.45$ at $n=25$ and dropped to
$5.61$ at $n=250$, indicating overflow rather than intrinsic misfit.
Table [3.6](#tab:autoregressive-perm-avg){reference-type="ref"
reference="tab:autoregressive-perm-avg"} reports these trajectories, and
Section [3.2](#sec:models-implementation){reference-type="ref"
reference="sec:models-implementation"} documents the derivative clipping
and ridge penalties that mitigate this failure mode when samples are
scarce.

High dimensionality converts small calibration errors into large
likelihood gaps because the triangular determinant accumulates
coordinate-wise discrepancies. MINIBOONE with $K=43$ illustrates this
accumulation: published flows achieved LL values between $-15.59$ and
$-11.68$, whereas TRTF reached $-30.01$ under the shared preprocessing.
Table [3.8](#tab:uci-loglik){reference-type="ref"
reference="tab:uci-loglik"} positions TRTF beside the flow baselines and
shows that the improvement over the Gaussian reference remains clear
even though an approximately $18$ nat gap persists to the strongest
flow.

Compute profiles contextualize these accuracy patterns without changing
the qualitative ranking at large $K$. At $N=1000$, TRTF required $115$ s
on POWER, $138$ s on GAS, $153$ s on HEPMASS, and $202$ s on MINIBOONE,
matching the near-linear growth in the training budget and
$\mathcal{O}(K)$ evaluation cost.
Table [3.9](#tab:real-compute){reference-type="ref"
reference="tab:real-compute"} summarizes these wall-clock measurements
and highlights that separable estimators remain practical in moderate
dimensions, yet accuracy dominates the choice once $K \approx 40$.

Taken together, the transport frame delineates when separability
suffices and when richer models become necessary. TRTF leads within the
separable family when location shifts capture most structure, exhibits
ordering sensitivity only through finite bases, and stabilizes with
modest sample sizes under the safeguards of
Section [3.2](#sec:models-implementation){reference-type="ref"
reference="sec:models-implementation"}. Performance degrades in high
dimensions where shape changes and interactions matter, at which point
non-separable models offer clear likelihood gains. These conclusions
motivate the guidance that will follow in the concluding subsection of
this chapter.

## Conclusions, Limitations, and Outlook {#sec:conclusion-outlook}

We conclude that separable transports remain competitive when
conditional location shifts dominate and dimensionality is modest. TRTF
led TTM-Sep on Half-Moon ($1.71$ versus $1.93$ NLL) and matched the
exponential coordinate in the four-dimensional generator, supporting
this interpretation. Conditional decompositions and calibration plots
indicate that residual error concentrates in context-dependent shapes,
particularly on the beta and gamma components. These findings align with
permutation averages that favor TRTF and quantify finite-basis
anisotropy. Tables [3.2](#tab:halfmoon-nll){reference-type="ref"
reference="tab:halfmoon-nll"}--[3.4](#tab:autoregressive-perm){reference-type="ref"
reference="tab:autoregressive-perm"} together with
Figure [3.2](#fig:autoregressive-joint-calibration){reference-type="ref"
reference="fig:autoregressive-joint-calibration"} document this evidence
under the shared protocol.

Performance on MINIBOONE reveals the cost of separability at higher
dimension. TRTF improved the Gaussian reference yet remained about
$18$ nats behind the best published flow, consistent with accumulated
Jacobian error across $43$ coordinates. POWER exhibited the opposite
regime: under identical preprocessing, the reported flows outperformed
TRTF (Table [3.8](#tab:uci-loglik){reference-type="ref"
reference="tab:uci-loglik"} lists TRTF at $-7.17$ versus flow baselines
near $0.17$ to $0.30$). These contrasts suggest that conditional shape
and dimension jointly determine whether separable structure suffices.
Table [3.8](#tab:uci-loglik){reference-type="ref"
reference="tab:uci-loglik"} reports these comparisons in a common unit.

Compute profiles remained practical and scaled near-linearly with the
training budget. Training plus evaluation required $115$ s at $N=1000$
on POWER and $202$ s on MINIBOONE, with longer totals at $N=2500$ that
preserved the same trend. These measurements keep separable transports
viable for exploratory analysis and model diagnostics.
Table [3.9](#tab:real-compute){reference-type="ref"
reference="tab:real-compute"} records the budgeted timings and the
shared pipeline settings.

Several limitations qualify these conclusions. Separable maps fix
conditional shape and therefore cannot resolve heteroskedasticity or
conditional multimodality. Ordering remained a material source of
variance under finite bases, as shown by the $0.97$ nat permutation gap
despite stable rankings at moderate sample sizes. In our $n=250$
synthetic runs, the TRTF versus TTM-Sep ranking did not change across
the $24$ permutations
(Table [3.4](#tab:autoregressive-perm){reference-type="ref"
reference="tab:autoregressive-perm"}); ordering affected magnitudes
rather than the lead. Simple ordering heuristics (identity or
Cholesky-pivoted with optional Gaussianization; see
Section [3.2](#sec:models-implementation){reference-type="ref"
reference="sec:models-implementation"}) reduced variance but did not
alter this pattern. Small-sample regimes created numerical fragility
through steep log-Jacobian terms, which clipping and ridge
regularization mitigate but do not eliminate. Real-data tables still
contain missing GAS and HEPMASS entries, and single-seed settings
persist for some runs, limiting external comparability.
Tables [3.4](#tab:autoregressive-perm){reference-type="ref"
reference="tab:autoregressive-perm"}--[3.9](#tab:real-compute){reference-type="ref"
reference="tab:real-compute"} catalog these caveats within the
standardized protocol.

The outlook follows directly from the evidence. Data-driven orderings
are likely to reduce anisotropy without abandoning the lower-triangular
map. Low-rank cross-terms in triangular transports and non-additive
predictors in TRTF may adapt conditional shapes while preserving
monotone structure, exact inversion, and linear per-sample evaluation.
We excluded these richer variants by design in
Chapter [1](#ch:intro){reference-type="ref" reference="ch:intro"}
(Non-goals) due to compute and calibration overhead; they remain
promising future work once resources permit. Expanded calibration
reporting, including probability-integral-transform summaries and
Kolmogorov--Smirnov distances, should remain part of any
deployment-grade evaluation. Completing GAS and HEPMASS under the same
protocol will improve generality and sharpen the accuracy-versus-compute
trade-off. These steps target smaller likelihood gaps on high-$K$
datasets while retaining the interpretability and reproducibility
provided by the transport frame.

# Appendix {#ch:appendix}

## Unified Transport Schematic {#app:transport-schematic-figure}

Figure [5.1](#fig:transport-schematic){reference-type="ref"
reference="fig:transport-schematic"} provides the full schematic of the
unified transport pipeline referenced throughout the thesis. The
landscape layout preserves readability for the granular annotations on
each modeling branch.

::: landscape
<figure id="fig:transport-schematic">

<figcaption>Unified evaluation pipeline shared by triangular transport
maps, Transformation Random Forests, and copulas. The diagram shows how
standardized features flow through the triangular pullback or copula
dependence block before reporting log densities, conditional
diagnostics, samples, calibration metrics, and compute summaries. The
shared preprocessing, Jacobian accumulation, and reporting path define
the evaluation protocol used across the thesis.</figcaption>
</figure>
:::

## Pseudo-code Summaries for Model Routines

This appendix records consolidated pseudo-code for the core R
implementations used in the experiments. Each summary captures inputs,
main processing stages, and outputs so the execution flow is transparent
without consulting the source code files.

### Transformation Random Forest (TRTF) {#app:trtf}

**Routine:**`fit_TRTF(S, config, seed, cores)` (calls `mytrtf`).

1.  Validate that the training matrix is numeric, set the RNG seed, and
    label columns as $X_1,\ldots,X_K$.

2.  Fit an intercept-only transformation model `BoxCox` for each $X_k$
    to provide baseline monotone transformations.

3.  For $k = 2,\ldots,K$:

    1.  Build the formula $X_k \sim X_1 + \cdots + X_{k-1}$.

    2.  Choose `mtry = max(1, floor((k-1)/2))` and standard `ctree`
        controls (`minsplit`, `minbucket`, `maxdepth`).

    3.  Fit a transformation forest with `traforest` and store the
        conditional model (one forest per $k$).

4.  Return a `mytrtf` object containing baseline transformations,
    conditional forests, variable-importance scores, and the seed.

5.  **Prediction (`predict.mytrtf`):**

    1.  Convert new data to the same column naming scheme and evaluate
        $X_1$ through its baseline transformation model to obtain
        marginal log densities.

    2.  For each conditional forest ($k\geq 2$) evaluate the log density
        of $X_k$ given $X_{1:(k-1)}$, extracting the diagonal when the
        forest returns a log density matrix.

    3.  Stack the per-dimension log densities (`logdensity_by_dim`) or
        sum them to obtain the joint log likelihood (`logdensity`).

### Nonparametric Copula Baseline {#app:copula}

**Routine:**`fit_copula_np(S, seed)`.

1.  Inspect the training matrix and optional class labels; detect
    whether the dedicated copula packages are available.

2.  If prerequisites fail (dimension $K \neq 2$ or labels missing), fall
    back to independent univariate kernel density estimates per
    dimension and store them for later interpolation.

3.  Otherwise, for each class label:

    1.  Fit one-dimensional `kde1d` models to each marginal $X_1$ and
        $X_2$.

    2.  Convert training samples to pseudo-observations using mid-ranks
        scaled by $(n+1)^{-1}$ and clamp to
        $(\varepsilon, 1-\varepsilon)$.

    3.  Fit a two-dimensional kernel copula with `kdecopula::kdecop`
        (method `TLL2`).

    4.  Store marginals, copula fit, and effective sample size for the
        class.

4.  Record class priors and return a `copula_np` object.

5.  **Prediction (`predict.copula_np`):**

    1.  In fallback mode evaluate each univariate KDE at the requested
        points and sum log densities.

    2.  In copula mode compute marginal log densities and CDF values,
        evaluate the copula density, and either:

        1.  Average over class-specific log densities weighted by priors
            (mixture prediction), or

        2.  Use the class labels supplied at prediction time.

    3.  Return per-dimension log densities or their sum depending on the
        requested type.

### Triangular Transport Core Utilities {#app:ttm-core}

**Module:**`ttm_core.R` (shared by marginal and separable TTM fits).

1.  Provide train-only standardization helpers that cache feature means
    and standard deviations and reapply them to new data.

2.  Define basis builders: polynomial features for predecessor
    coordinates $g_k$, monotone basis functions $f_k$ for the current
    coordinate, and their derivatives.

3.  Implement optional ordering heuristics (identity or Cholesky
    pivoting with optional Gaussianization) and persist selected
    permutations.

4.  Expose a dispatcher `ttm_forward(model, X)` that:

    1.  Standardizes inputs using stored parameters.

    2.  For marginal maps apply affine transformations $a_k + b_k x_k$
        with precomputed coefficients.

    3.  For separable maps constructs $g_k$ and $f_k$, computes
        $S_k = g_k + f_k$, and records the Jacobian diagonal
        $\partial_{x_k} S_k$.

5.  Provide `ttm_ld_by_dim` to combine the forward map with the Gaussian
    reference, yielding per-dimension log densities used by all TTM
    variants.

### Marginal Triangular Transport Map {#app:ttm-marg}

**Routine:**`fit_ttm_marginal(data, seed)`.

1.  Split data into train/test subsets if only a matrix is provided;
    otherwise accept a prepared list.

2.  Standardize training features and, for each dimension $k$, compute
    closed-form coefficients $(a_k, b_k)$ that minimize the Gaussian
    pullback objective subject to $b_k > 0$.

3.  Store model parameters (standardization, per-dimension coefficients,
    ordering) and time measurements.

4.  During prediction call `ttm_forward` with the marginal coefficients
    and convert Jacobian diagonals to log densities via `ttm_ld_by_dim`;
    aggregate per-dimension contributions when the joint log density is
    requested.

### Separable Triangular Transport Map {#app:ttm-sep}

**Routine:**`fit_ttm_separable(data, degree_g, lambda, seed)`.

1.  Prepare train/test splits and standardize training features as in
    the marginal case.

2.  For each coordinate $k$:

    1.  Build polynomial features $g_k$ on previous coordinates (degree
        set by `degree_g`).

    2.  Build monotone basis functions $f_k$ on the current coordinate
        and their derivatives.

    3.  If `degree_g = 0`, use the marginal closed-form solution to
        recover affine parameters.

    4.  Otherwise solve the regularized optimization problem
        $\min_c \frac{1}{2}\lVert (I - \Phi_{\text{non}} M)c \rVert^2 - \sum \log (B c) + \lambda\,\text{penalty}(c)$
        using `optim` with L-BFGS-B while enforcing positivity of the
        derivative.

    5.  Store coefficients $c_{\text{non}}$ and $c_{\text{mon}}$ for the
        coordinate.

3.  Assemble the model list with standardization parameters,
    coefficients, and metadata; record training/prediction timings.

4.  At prediction time re-use `ttm_forward` and `ttm_ld_by_dim` to
    obtain per-dimension and joint log densities.

### Evaluation Utilities {#app:evaluation}

**Module:**`evaluation.R` (experiment orchestration).

1.  Define convenience helpers such as `stderr(x)` and `add_sum_row` for
    table post-processing.

2.  `prepare_data(n, config, seed)` samples from the configured
    data-generating process, splits the sample into
    train/validation/test sets, and returns both the matrix of draws and
    the split structure.

3.  `fit_models(S, config)` fits the oracle TRUE density and the TRTF
    baseline on a split, times their evaluations, and returns the fitted
    objects together with per-dimension log-likelihood arrays.

4.  `calc_loglik_tables(models, config, X_te, ...)` aggregates negative
    log-likelihoods (nats) for TRUE (marginal and joint), TRTF, TTM, and
    separable TTM, formats the results with standard-error bands,
    appends a summary row, and renames columns for presentation.

5.  `eval_halfmoon(mods, S, out_csv)` ensures all requisite models are
    available (TRTF, TTM variants, copula baseline), evaluates them on
    the half-moon test split, computes joint and per-dimension negative
    log-likelihoods, and optionally persists the metrics as CSV
    artifacts.

These structured summaries allow reproducing the algorithmic flow of
each model without navigating the full R implementation.

### Supplementary Results {#app:supplementary}

*(mean NLL in nats).*

::: {#tab:ttmsep-n25-overflow}
  Model                   $n=25$
  ------------------ -----------
  TTM-Sep$^{\ast}$     $6829.45$

  : Note: $^{\ast}$ out-of-scope setting. This value reflects numerical
  overflow of the separable map in the sparse regime ($n=25$). Stronger
  derivative clipping and ridge regularization
  (Section [3.2](#sec:models-implementation){reference-type="ref"
  reference="sec:models-implementation"}) remove this failure in reruns.
  We mark this configuration as out of scope and exclude it from
  main-text comparisons; the table remains for transparency.
:::

<figure id="fig:ordering-heuristics-4d">

<figcaption>Ordering sensitivity and mitigation window on the
four-dimensional generator at <span
class="math inline"><em>n</em> = 250</span> (joint NLL; nats; lower is
better). Markers show the best (min over <span
class="math inline">24</span> permutations), canonical ordering, and
permutation median for each method (values from Chapter <a
href="#ch:dataanalysis" data-reference-type="ref"
data-reference="ch:dataanalysis">3</a>: Tables <a
href="#tab:autoregressive-nll" data-reference-type="ref"
data-reference="tab:autoregressive-nll">3.3</a> and <a
href="#tab:autoregressive-perm-spread" data-reference-type="ref"
data-reference="tab:autoregressive-perm-spread">3.5</a>). A simple
heuristic selects the better of two candidates—identity and
Cholesky-pivoted with optional Gaussianization—aiming to move toward the
“min” marker while keeping evaluation linear.</figcaption>
</figure>
