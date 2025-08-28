---

# A FRIENDLY INTRODUCTION TO TRIANGULAR TRANSPORT

**Maximilian Ramgraber**
Delft University of Technology
Delft, Netherlands
m.ramgraber@tudelft.nl

**Daniel Sharp**
Massachusetts Institute of Technology
Cambridge, USA
dannys4@mit.edu

**Mathieu Le Provost**
Massachusetts Institute of Technology
Cambridge, USA
mleprovo@mit.edu

**Youssef Marzouk**
Massachusetts Institute of Technology
Cambridge, USA
ymarz@mit.edu

March 28, 2025

## ABSTRACT

Decision making under uncertainty is a cross-cutting challenge in science and engineering. Most approaches to this challenge employ probabilistic representations of uncertainty. In complicated systems accessible only via data or black-box models, however, these representations are rarely known. We discuss how to characterize and manipulate such representations using **triangular transport maps**, which approximate any complex probability distribution as a transformation of a simple, well-understood distribution. The particular structure of triangular transport guarantees many desirable mathematical and computational properties that translate well into solving practical problems. Triangular maps are actively used for density estimation, (conditional) generative modelling, Bayesian inference, data assimilation, optimal experimental design, and related tasks. While there is ample literature on the development and theory of triangular transport methods, this manuscript provides a detailed introduction for scientists interested in employing measure transport without assuming a formal mathematical background. We build intuition for the key foundations of triangular transport, discuss many aspects of its practical implementation, and outline the frontiers of this field.

## 1 Motivation

**Who is this tutorial for?** This manuscript is an accessible introduction to **triangular transport**, a powerful and versatile method for generative modelling and Bayesian inference. In particular, triangular transport underpins effective algorithms for data assimilation, solving inverse problems, and performing simulation-based inference, with applications across myriad scientific disciplines.

This tutorial targets researchers with an interest in applied statistical methods but without a formal background in mathematics. Consequently, we will focus more on intuition, general concepts, and implementation, referring the reader to other relevant articles for more formal exposition and theory.

**How does triangular transport work?** Like other measure transport methods, triangular (measure) transport is a framework to transform one probability distribution into another. This operation is highly useful, as it allows us to characterize a complex **target distribution** π by transforming a simpler, known **reference distribution** η. Throughout this manuscript, we use orange to denote variables associated with the (problem-specific) target distribution and green to denote variables associated with the (user-defined) reference distribution, e.g., a standard Gaussian. The idea of coupling two distributions is used in a wide range of applications. In **generative modelling**, for example, we are usually interested in creating samples of a target distribution π, such as the distribution of 400 × 400 pixel images of cats. Measure transport methods approach this challenge by first learning a transport map S that transforms η to π, and then use S to convert reference samples z ∼ η (here: 400 × 400 pixel images of Gaussian white noise) into samples from the target distribution x̄ = S(z) ∼ π (here: 400 × 400 pixel images of cats).

Some of these methods - among them triangular transport - can also characterize conditionals π(a|b∗) of the joint target distribution π(a, b) of two random variables a and b. Here blue denotes the fact that b∗ is a deterministic, fixed value. Conditioning operations often arise as stochastic generalizations of evaluating deterministic processes (see Figure 1). As we will describe in Section 2, conditioning is also central to the Bayesian approach to statistical inference, which is an important tool across many scientific disciplines. We distinguish here between generating from π(a|b∗) (“generate an image of a grey cat”) and π(a,b) (“generate a color and a cat of that color”) by assuming that b∗ is determined outside of our control, either by user or application.

[Image of Figure 1: Progression from a fully deterministic to a fully stochastic system. (A) Numerical models are usually represented as deterministic functions. (B) In the presence of input uncertainty, deterministic functions encode a deterministic coupling which yields uncertain output (see Section 2.2). (C) If the function itself is uncertain, this coupling “blurs” into a joint probability distribution. Function evaluation now corresponds to characterizing a conditional distribution. Mind that (A) and (B) can also be parsed as degenerate joint probability distributions.]

**What makes triangular transport special?** At the heart of triangular transport is their eponymous triangular structure. This structure sets them apart from other measure transport methods such as normalizing flows (e.g., Rezende and Mohamed, 2015; Kobyzev et al., 2021), which compose together many simpler but somewhat ad hoc transformations, often interleaved with permutations, and even from conditional normalizing flows (e.g., Van Den Oord et al., 2016), which parameterize normalizing flows in order to represent block triangular (rather than strictly triangular) maps. Triangular structure – discussed in greater detail in subsequent sections – has many important practical properties:

*   **Parsimony**: The parameterization of the map function is at the user’s discretion (see Section 3.1). This means we can adjust the map’s overall complexity, down to the complexity with which it resolves individual variables and variable dependencies. This allows us to implement nonlinear maps that are as complex as necessary, and yet as simple as possible (see Section 4.2).
*   **Sparsity**: Triangular maps have a natural ability to exploit conditional independence. This improves their computational efficiency, which enables these maps to scale to high-dimensional settings. Further, such structure makes them highly robust to spurious correlations and smaller sample sizes (see Section 2.3.3).
*   **Numerical convenience**: Constructing triangular maps boils down to parameterizing simple one-dimensional monotone functions, a task with a rich body of supporting literature. Because of this, these maps are easy to optimize and invert, which we investigate in detail.
*   **Explainability**: Triangular maps have a clear correspondence between their constituent elements and the statistical features they represent (see Section 2.3). We can then readily describe different factorizations of the target distribution using the elements of a triangular map, with a particular focus on combinations of various conditional and marginals of the target.

**In what applications has triangular transport been successful?** Triangular transport has been applied to a wide range of statistical problems in many different disciplines. Examples of such applications include:

*   **Bayesian inference**: Triangular transport lends itself exceptionally well to the sampling of conditional distributions. As such, it has found application in both variational (El Moselhy and Marzouk, 2012) and simulation-based inference (Marzouk et al., 2017; Rubio et al., 2023; Baptista et al., 2024a), for large-scale inverse problems (Brennan et al., 2020) and in applications with multiscale structure (Parno et al., 2016).
*   **Data assimilation**: Triangular transport provides true nonlinear generalizations of popular filtering (Spantini et al., 2022) and smoothing (Ramgraber et al., 2023a,b) algorithms such as the ensemble Kalman filter and smoother, and their many variants (Grange et al., 2024).
*   **Density estimation and generative modelling**: The coupling learned by triangular transport is highly useful for the estimation (Wang and Marzouk, 2022; Martinez-Sanchez et al., 2024; López-Marrero et al., 2024) and sampling (Irons et al., 2022) of non-Gaussian probability distributions, even in high dimensions (Katzfuss and Schäfer, 2024).
*   **Optimal experimental design**: Due to the close connections between conditional densities and expected information gain or mutual information, triangular transport maps are useful for estimating common objectives in Bayesian optimal experimental design (Huan et al., 2024; Koval et al., 2024; Li et al., 2024).

Further applications of triangular transport include methods for joint state-parameter inference in state-space models (Spantini et al., 2018; Grashorn et al., 2024; Zhao and Cui, 2024), solving Fokker–Planck equations (Zeng et al., 2023), stochastic programming (Backhoff et al., 2017), and even the discovery of causal models from data (Akbari et al., 2023; Xi et al., 2023).

**How is this tutorial structured?** In the following, we will provide an intuition-focused introduction to the theoretical basics of triangular transport (Section 2), discuss practical aspects related to their implementation in code (Section 3), and conclude with a brief overview of interesting research directions (Section 5). For researchers interested chiefly in practical implementation, a flow chart of the most important steps in the construction and application of triangular transport is provided in Figure 2. First, we define the target π and reference η (Section 2.3). Then, we structure, parameterize (Section 3.1), and optimize (Section 3.2) the triangular map. Finally, we can deploy the map in the application of our choice, with some practical heuristics listed in Section 4. The tutorial will involve several recurring variables, which are summarized in Table 1.

[Image of Figure 2: Flowchart of the main steps required to define, build, optimize, and apply triangular maps, with links to the relevant sections.]

## 2 Theory

### 2.1 Bayesian inference

To begin, let us briefly revisit some basic concepts of Bayesian inference which will serve to motivate the operations explored in the following sections. In short, Bayesian inference is based on *Bayes’ theorem*: given two random variables (RVs) **a**, **b** with joint probability density function (pdf) p(a,b), we see

`p(a|b*) = (p(a)p(b*|a)) / p(b*)` (1)

Equation (1) subsumes three sequential operations (see Figure 3; e.g., Gelman et al., 2013):

1.  First, the marginal prior p(a) is combined with a conditional observation model (sometimes also called likelihood model) p(b|a), yielding a joint probability distribution p(a,b) = p(a)p(b|a) over all possible combinations of a and b. This joint distribution describes how the variable of interest a and the predicted observations b relate to each other.
2.  Next, this joint distribution is conditioned on a specific observation b∗. In practice, this means evaluating p(a,b) for all possible values a while keeping b fixed at the value of b∗. This extracts a slice p(a,b∗) of this joint distribution at h∗ along different values of a.
3.  Since the probability densities along this slice do not generally integrate to 1, this slice does not constitute a valid pdf. The final step thus normalizes the probability densities against the slice’s probability mass p(b∗) = ∫ p(t,b∗) dt, yielding the posterior pdf p(a|b∗).

In summary, Bayes’ theorem first constructs a joint distribution p(a,b) from a prior p(a) and an observation model p(b|a), then conditions it on a specific observation value b∗ and re-normalizes. In consequence, one could reformulate Equation (1) equivalently as:

`p(a|b*) = p(a,b*) / ∫ p(t,b*)dt` (2)

where p(a,b∗) evaluates the joint pdf p(a,b) for all possible a while keeping b fixed at b∗, and the denominator acts as a normalizing constant. This equation, or reformulation thereof, lie at the heart of all Bayesian inference algorithms. Unfortunately, it is generally impossible to formulate p(a,b) in closed form, which in turn makes it difficult to evaluate the posterior p(a|b∗). To overcome this challenge, different Bayesian inference methods use different strategies. As we shall see in the following, triangular transport solves this challenge by first approximating an almost arbitrary joint pdf p(a,b) by using the concept of measure transport (Section 2.2.1). Crucially, this transformation then allows us to evaluate any of its conditionals p(a|b∗) (Section 2.3.2).

**Table 1: Recurring notation and variables.**

| | | | |
| :--- | :--- | :--- | :--- |
| **bold font** | vector-valued variable or function | Roman font | scalar-valued variable or function |
| **S** | Target-to-reference map | S_k | k-th map component function |
| **R** | Reference-to-target map (see Section 3.2.2) | T | twice Archimedes’ constant, i.e., 6.283185... |
| π | target distribution of interest | η | reference distribution, often standard Gaussian |
| orange | variable associated with π | green | variable associated with η |
| x∼π | target random variable | z∼η | reference random variable |
| X^i | ith realization of x∼π | Z^i | ithe realization of z∼η |
| S^#_η | pullback distribution | S_±π | pushforward distribution |
| K | number of target dimensions | N | ensemble size |
| p | generic probability density function (pdf) | a, b | generic random variables |
| y∗ | conditioning variable | x◦ | conditioned variable x−∼p(x|y−) |
| c | basis function coefficient | r | rectifier (r:R→R+) |
| f | monotone function | g | nonmonotone function |

[Image of Figure 3: Schematic illustration of Bayes theorem, using a Beta mixture prior p(a) and an observation model p(b|a) = Ν (μ = 2a^3 - a, σ = 0.075). Left: We can create a joint distribution p(a,b) from a prior p(a) and the observation model p(b|a). Right: conditioning this joint distribution on a specific value b* retrieves a posterior p(a|b*).]

The remainder of this tutorial drops this generic notation for two RVs a and b jointly distributed as p(a,b) in favour of a single RV x distributed according to a distribution π. You can think of x as an augmented RV x = [b,a], and of the joint density as π(x) = p(a,b). Transport methods generally operate on this joint distribution p(a,b); we focus primarily on the ones that will allow us to characterize the desired conditionals p(a|b).

### 2.2 The change-of-variables formula

The key to understand transport methods is the **change-of-variables formula**. This formula allows us to relate a RV x associated with a complicated target pdf π, known only to proportionality or through samples, to a second RV z, associated with a much simpler, user-specified reference pdf η through an invertible, differentiable transformation S. In essence, the change-of-variables formula describes what happens to pdfs when they are subjected to specific transformations. For scalar-valued RVs x and z, the change-of-variables formula is defined as

`π(x) = S*η(x) = η(S(x)) |∂S(x)/∂x|`
`η(z) = S*π(z) = π(S⁻¹(z)) |∂S⁻¹(z)/∂z|` (3)

where z = S(x) and x = S⁻¹(z), and the pullback density¹ S\*η(x) obtains by applying the inverse map S⁻¹ to the reference distribution η. The alternate form in the second line of Equation (3) reflects the fact that since S is invertible, each distribution π and η can be expressed in terms of the other through either the (forward) map S or its inverse S⁻¹. Consequently, the **pushforward** pdf² S±π(z) likewise approximates η by applying the forward map S to the target³ π. For multivariate RVs **x** and **z**, the same principle applies. Given a multivariate monotone function **S**, Equation (3) generalizes to:

`π(x) = S*η(x) = η(S(x))|det ∇ₓS(x)|`
`η(z) = S*π(z) = π(S⁻¹(z))|det ∇₂S⁻¹(z)|` (4)

The change-of-variables formula in Equation (3) has a surprisingly intuitive interpretation. It states that the probability density η(z) at a point of interest z after the transformation equals the original probability density π(x) at the pre-transformation point x = S⁻¹(z), adjusted for any deformation the transformation might have induced at this location |∂S(x)/∂x|. Equation (4) extends this notion to multi-dimensional systems. Intuitively, the absolute value of the determinant of the Jacobian of **S**, namely |det ∇ₓS(x)|, measures the inflation/deflation of an infinitesimal volume centered about **x** by the map **S**. If |det ∇ₓS(x)| = 1, the map **S** preserves infinitesimal volumes about **x**. The compensation term, similar to ‘u-substitution’ in calculus, is necessary because transformations **S** “stretch” or “squeeze” the spaces they are applied to. Accounting for this spatial distortion ensures that the probability mass is preserved and thus the transformed distribution η(z) remains a valid pdf.

The change-of-variables formula allows us to describe how a probability distribution π changes when subjected to a specific transformation S. An example is provided in Figure 4, which shows how a non-Gaussian target pdf π can be related to a Gaussian reference pdf η through an invertible transformation; this invertibility is equivalent to monotonicity in the scalar case.

[Image of Figure 4: Illustration of the change-of-variables formula. A monotone function S allows us to relate a RV x associated with a pdf π to a RV z associated with a pdf η. We can evaluate η(z) by evaluating the target π at the pre-image of z, that is to say π(S⁻¹(z)), then multiplying it with the absolute inverse map’s gradient |∂zS⁻¹(z)|.]

#### 2.2.1 Connection to transport methods

The change-of-variables formula has three knobs to tune: the original distribution π, the map S, and its transformed output η. When considering the change-of-variables in elementary calculus courses, π and S are often assumed known, and we seek its transformed output η. Transport methods choose a slightly different approach. We assume π and η to be known (at least partially), and instead seek the specific map S which relates the two distributions to each other.

As discussed in Section 1, we generally do not know the target distribution π in closed form. Often, the target density π is known only partially, either through samples X ∼ π or up to proportionality (π̄ = mπ, where m > 0 is an unknown constant). On the other hand, the reference distribution η is defined as a simple, well-known distribution, often a standard Gaussian pdf N(0,I) (Figure 5). The map S is identified by minimizing an objective function over a specified class of functions, see the discussion in Section 3.2. By finding this map, we learn how to construct the unknown target distribution π by transforming a well-defined reference distribution η.

[Image of Figure 5: A transport map S relates a RV x associated with the target π to a RV z associated with the reference η. If the map is monotone, it can be applied both ways.]

Among its other uses, learning the map S allows us to cheaply draw new samples from the target distribution π. This is achieved by first sampling the reference η, then applying the inverse map S⁻¹ to the resulting reference samples z. This is especially useful in applications where sampling the target conventionally would involve computationally expensive simulations of, e.g., partial differential equations (emulation), or systems in which the sample-generating process is not known exactly (generative modelling; e.g., Baptista et al., 2024c).

### 2.3 Triangular maps and their uses

Many classes of functions are viable choices for the transport map S. Common examples include normalizing flows and GANs. However, an especially useful class among these are **triangular transport maps**. In addition to sampling the target π, triangular maps are unique in allowing us to sample conditionals of π, which makes them a flexible and versatile tool for Bayesian inference. Triangular maps are structured as follows:

`S(x) = [S₁(x₁), S₂(x₁, x₂), ..., Sₖ(x₁, ..., xₖ)]ᵀ = [z₁, z₂, ..., zₖ]ᵀ = z` (5)

where the full map S: ℝᴷ → ℝᴷ is comprised of map components Sₖ: ℝᵏ → ℝ, k=1,...,K, each of which depends only on the first k entries of the target RV x = [x₁, ..., xₖ]ᵀ and we enforce that ∂ₓₖSₖ(x₁, ..., xₖ) > 0 for any feasible choice of x₁, ..., xₖ. When all of S₁, ..., Sₖ satisfy this, we call the map S “monotone”. The eponymous triangular nature of S refers to the fact that the map’s partial derivatives with regards to x are lower-triangular; that is to say, the Jacobian matrix VS has all zeros above its diagonal. This structure—also known as a Knothe–Rosenblatt rearrangement (Knothe, 1957; Rosenblatt, 1952) or KR map—has a number of highly desirable properties:

1.  Over all functions satisfying this structure, there is one that uniquely couples π and η under mild conditions (e.g., Marzouk et al., 2017).
2.  This triangular structure allows us to evaluate the **determinant of the map’s Jacobian** det ∇S(x) efficiently as the product of its diagonal entries (Marzouk et al., 2017), which proves highly useful for the map’s optimization (see Section 3.2), not to mention when performing the density estimation itself for the changed variables:

    `det ∇S(x) = Πₖ(∂Sₖ(x₁, ..., xₖ) / ∂xₖ)` (6)
3.  Triangular maps are **easily invertible**. In particular, we pick a class of functions that are monotone in their last input xₖ (with no requirements on the other inputs x₁, ..., xₖ₋₁), i.e. ∂ₓ Sₖ > 0 everywhere; this guarantees the monotonicity and thus eases our inversion computation. We will discuss ways to guarantee this property in Section 3.1.1.
4.  Perhaps most importantly, triangular maps naturally **factorize the target distribution** into a product of marginal conditional pdfs. We will investigate this property in greater detail in the following sections.

#### 2.3.1 Map inversion

In the forward map evaluation (Equation (5)), each of the map’s constituent map component functions Sₖ can be evaluated independently, even in parallel, and then assembled into the full reference vector z. However, the same does not hold for the inverse map:

`S⁻¹(z) = [S₁⁻¹(z₁), S₂⁻¹(z₂; x₁), ..., Sₖ⁻¹(zₖ; x₁, ..., xₖ₋₁)]ᵀ = [x₁, x₂, ..., xₖ]ᵀ = x` (7)

Here, the inverse map’s component functions Sₖ⁻¹ must be evaluated in sequence and cannot be evaluated independently. This process begins by inverting the first map component S₁⁻¹(z₁), a trivial one-dimensional root finding problem, which yields x₁. This output serves as auxiliary input for the second map component’s inversion S₂⁻¹(z₂;x₁), yielding another one-dimensional root-finding problem, which provides x₂. All subsequent map component inversions are similar one-dimensional root-finding problems that likewise depend on the outcomes of previous inversions. This dependence of each inversion Sₖ⁻¹ on each of x₁, ..., xₖ₋₁, the outcomes of previous inversions, effectively factorizes the target distribution as a product of marginal conditionals (Villani, 2007):

`π(x) = π(x₁)π(x₂|x₁) ... π(xₖ|x₁, ..., xₖ₋₁)` (8)
`S₁⁻¹(z₁) S₂⁻¹(z₂;x₁) ... Sₖ⁻¹(zₖ;x₁, ..., xₖ₋₁)`

where each term corresponds to, and is in turn sampled by, one of the inverse map components indicated in the underbraces. In other words, for a particular sample i, each row of Equation (7) can be used to generate a sample Xₗⁱ from a particular marginal distribution conditioned on X₁, ..., Xₖ₋₁:

`Xₖⁱ = Sₖ⁻¹(Zₖⁱ; X₁ⁱ, ..., Xₖ₋₁ⁱ) ~ Sₖⁱⁱηₖ = π(xₖ|x₁, ..., xₖ₋₁)` (9)

where Sₖⁱⁱηₖ is the pullback of the one-dimensional marginal reference ηₖ.

#### 2.3.2 Sampling conditionals

As it turns out, the factorization of the target distribution π in Equations (8) and (9) also allows us to sample conditionals of π, including the Bayesian posterior p(a|b) (assuming p := π and \[b, a] := \[x₁:ₖ, xₖ₊₁:ₖ]). This can be achieved by manipulating the inversion process. First, observe that the factorization in Equation (8) can be aggregated into two blocks:

`π(x) = π(x₁:ₖ)π(xₖ₊₁:ₖ|x₁:ₖ)` (10)
`S⁻¹₁:ₖ(z₁:ₖ) S⁻¹ₖ₊₁:ₖ(zₖ₊₁:ₖ; x₁:ₖ)`

Similarly, we can aggregate the map component functions into two blocks:

`S⁻¹(z) = [S⁻¹₁:ₖ(z₁:ₖ), S⁻¹ₖ₊₁:ₖ(zₖ₊₁:ₖ; x₁:ₖ)]ᵀ = [x₁:ₖ, xₖ₊₁:ₖ]ᵀ = x` (11)

Instead of evaluating of Equation (7) from top to bottom (S⁻¹₁ to S⁻¹ₖ), if we are interested in sampling conditionals, we may skip the upper map block S⁻¹₁:ₖ and replace its corresponding output x₁:ₖ = \[x₁, . . . , xₖ] with arbitrary, user-specified values x\*₁:ₖ = \[x\*₁, . . . , x\*ₖ]. This results in the following truncated inversion for the lower map block S⁻¹ₖ₊₁:ₖ,

`S⁻¹ₖ₊₁:ₖ(zₖ₊₁:ₖ; x\*₁:ₖ) = [S⁻¹ₖ₊₁(zₖ₊₁; x\*₁:ₖ), ..., S⁻¹ₖ(zₖ; x\*₁:ₖ, ..., x\*ₖ₋₁)]ᵀ = [x\*ₖ₊₁, ..., x\*ₖ]ᵀ = x\*ₖ₊₁:ₖ` (12)

Resuming the inversion starting with the map component inverse S⁻¹ₖ₊₁ thus yields samples x\*ₖ₊₁:ₖ from the conditional π(xₖ₊₁:ₖ|x\*₁:ₖ) instead of the target π(x). Equivalent to Equation (8), the corresponding conditional distribution now factorizes as:

`π(xₖ₊₁:ₖ|x\*₁:ₖ) = π(xₖ₊₁|x\*₁:ₖ) ... π(xₖ|x\*₁:ₖ, xₖ₊₁:ₖ₋₁)` (13)
`S⁻¹ₖ₊₁(zₖ₊₁; x\*₁:ₖ) ... S⁻¹ₖ(zₖ; x\*₁:ₖ, ..., x\*ₖ₋₁)`

This means that the manipulated triangular map inversion in Equation (7) allows us to sample conditionals of the target distribution π for arbitrary x\*₁:ₖ, or to estimate the density of this conditional distribution. Recalling Section 2.1, it is plain to see why this operation proves extremely useful in Bayesian inference: If we define our target pdf π as the joint distribution p(a, b) (using the notation of Section 2.1) between the RV of interest a and the observation predictions b, and consider the manipulated samples x\*₁:ₖ to be the observations b\* of b, then the manipulated inversion in Equation (12) samples the posterior p(a|b\*). An illustration of the forward, inverse, and conditional mapping operations is provided in Figure 6.

[Image of Figure 6: (A) The forward map (top row, from right) and its inverse (top row, from left) operate via implicit intermediate distributions (center), as the map components transform the distributions one marginal at a time. (B) Supplying the inverse map with a manipulated intermediate distribution (bottom row, center) as a starting point will instead sample conditionals of the target distribution.]

#### 2.3.3 Conditional independence

A related, very useful property of triangular transport maps is that they naturally allow for the exploitation of conditional independence by construction. Recalling the map’s generic factorization of the conditionals of target π in Equation (8), we might ask ourselves what happens in systems in which we can exploit conditional independence. For example, if we have a target distribution π(x₁:₄) and conditional independence properties x₃ ⊥⊥ x₁|x₂ and x₄ ⊥⊥ x₁, x₂|x₃ (corresponding to Markov structure), the map’s factorization could be reduced as follows:

`π(x₁:₄) = π(x₁)π(x₂|x₁)π(x₃|x₂, x₁)π(x₄|x₃, x₁, x₂) = π(x₁)π(x₂|x₁)π(x₃|x₂)π(x₄|x₃)` (14)

We refer to this reduction as **sparsification**. Triangular transport maps allows us to leverage these conditional independence properties by simply dropping the corresponding arguments from the map components Sₖ:

`S(x₁:₄) = [S₁(x₁), S₂(x₁, x₂), S₃(x₂, x₃), S₄(x₃, x₄)]ᵀ = [z₁, z₂, z₃, z₄]ᵀ = z₁:₄` (15)

Equivalently, its inverse map would be:

`S⁻¹(z₁:₄) = [S⁻¹₁(z₁), S⁻¹₂(z₂; x₁), S⁻¹₃(z₃; x₂), S⁻¹₄(z₄; x₃)]ᵀ = [x₁, x₂, x₃, x₄]ᵀ = x₁:₄` (16)

Making use of conditional independence properties in this way is useful for two reasons:

1.  **Robustness**: Any conditional independence we can enforce by construction is statistical information the map does not have to pry from the samples or the model, improving the overall fidelity and robustness of the approximation to π in settings with finite ensemble size.
2.  **Efficiency**: The removal of superfluous dependencies reduces the number of input arguments to many of the map components Sₖ, decreasing the evaluation and inversion complexity (e.g., see Section 3.1.1 for evaluation complexity). This property is often called sparsification as it turns the Jacobian ∇S into a sparse matrix.

This second property is the key to applying transport methods in high-dimensional systems. As each map component function Sₖ generically depends on all previous arguments x₁:ₖ₋₁, the computational demand explodes with the dimension of the target π when optimizing or evaluating the map. With sufficient conditional independence, however, sparse maps can overcome this dramatic increase in complexity. For instance, the Markov structure in Equations (15) and (16) results in a sparse map S with numerical complexity scaling linearly in the target dimension K.

A comprehensive account of the link between conditional independence and the sparsity of triangular maps is given in Spantini et al. (2018), through two main lines of results. First, given a sparse undirected probabilistic graphical model that encodes Markov properties of the target distribution π, it is shown how to predict the sparsity pattern of the triangular map S. This process relies on an ordered graph elimination algorithm, and can thus be performed before learning the map itself. But the resulting sparsity pattern depends on the chosen ordering of the random variables, which underscores the fact that triangular maps are intrinsically anisotropic objects: a good ordering is necessary in order to maximize sparsity. We comment further on methods for finding such orderings in Section 5. Second, Spantini et al. (2018) show that a property somehow dual to sparsity is decomposability: given the Markov structure of some distribution π on ℝᴷ, the inverse of the map S defined by S♯π = η can be represented as the composition of finitely many low-dimensional triangular transport maps, where low dimensionality here means that each component map is a function only of a small number of inputs. The exact structure of this decomposition follows from, again, an ordered decomposition of the original graph. These results allow very high-dimensional problems to be broken into many smaller and more manageable parts, given some conditional independence.

#### 2.3.4 Block-triangular maps

Recalling the block structure in Section 2.3.2, we have so far assumed that the map component “blocks” in Equation (11) are internally triangular; i.e., S₁:ₖ has output \[S₁(x₁), . . . , Sₖ(x₁:ₖ)] and similar for Sₖ₊₁:ₖ. We may instead use any invertible functions S₁:ₖ : ℝᵏ → ℝᵏ and Sₖ₊₁:ₖ : ℝᴷ → ℝᴷ⁻ᵏ, for example certain neural networks (e.g., Baptista et al., 2024c). This still permits sampling conditionals and can be numerically advantageous in high-dimensional systems, but can compromise the ability to exploit conditional independence within the map blocks. In the following, we will assume that the transport map is fully triangular, even where we adopt the block-triangular structure for ease of notation.

## 3 Implementation

In the preceding section, we have discussed the theoretical foundations of triangular transport, but stopped short of defining the precise method to evaluate Sₖ, or the optimization problem which identifies the specific map S from π to η. In this section, we will fill these gaps, and provide guidance on the practical implementation of triangular transport methods in code. First, though, a slight interlude about the “work” a map performs to transform π and η into one another. Without a loss of generality, suppose π and η have the same mean and covariance. In the case that both are Gaussian, this implies that they are the same distribution, and the transport map S is trivially the identity. Once π becomes non-Gaussian, the transport map must be strictly nonlinear, for any linear map merely transforms these first two moments! The nonlinearity of a map is therefore directly correlated with how “different” our reference and target distributions are, measured in one way or another (e.g., the Kullback–Leibler divergence, discussed below). This theoretical idea underpins many practical considerations; for example, we will assume that π and η approximated have the same tails (i.e., π(x) ≈ Cη(x) for some C ∈ ℝ when x is sufficiently large), implying the map must become close to linear sufficiently far from the origin.

### 3.1 Structuring the map components Sₖ

A core component of any triangular transport map S is the definition of the constituent map components Sₖ. In this subsection, we will introduce and discuss important features in the definition of these map components Sₖ.

#### 3.1.1 Monotonicity

Let us recall from Section 2.2 that each map component function must be monotone in its last argument xₖ, that is to say, ∂ₓₖ Sₖ > 0 on the entire domain (i.e., for any choice x₁:ₖ₋₁). Together with the triangular structure and linearity in the tails, this property guarantees that the resulting map S is bijective.

This may not be immediately intuitive, so an illustration of this effect is provided in Figure 7, which shows an example in ℝ². Here, invertibility ensures that every tuple (X₁, X₂) maps to a unique tuple (Z₁, Z₂). Graphically, this means that if we overlay Figure 7C and Figure 7D on top of each other, any pair of contour lines between the subplots must not intersect more than once. This is achieved due to triangularity and monotonicity:

*   The triangular structure ensures that the contours of every map component S₁:ₖ₋₁(x₁:ₖ₋₁) independent of xₖ and are thus constant along xₖ. In consequence, the contours of S₁ in Figure 7C are constant along x₂, which aligns them parallel to x₂ (see Figure 7C).
*   The monotonicity of S₁(x₁) relates each X₁ to a unique Z₁. We then obtain a unique tuple (Z₁, Z₂) when we ensure that the second set of contours from S₂(X₁, x₂) increases monotonously along x₂ (see Figure 7D). This is achieved through the monotonicity requirement ∂ₓₖ Sₖ > 0.

[Image of Figure 7: (A) Monotonicity guarantees that a function remains invertible. (B) Non-monotone map component functions have non-unique inverses. How do multivariate maps S remain monotone? (C) The map component function S₁(x₁) = z₁ is monotone in x₁ but constant in x₂, as it does not depend on this input. (D) The map component function S₂(x₁, x₂) = z₂ is monotone in x₂ but can be nonmonotone in x₁. Color indicates the magnitude of the map component’s output in all subplots.]

The same principle extends to higher dimensions: the triangular structure means that the (k − 1)-dimensional contours of S₁:ₖ₋₁(x₁:ₖ₋₁) are aligned along xₖ. Monotonicity of Sₖ along xₖ then ensures that for (X₁, . . . , Xₖ₋₁) and any xₖ, we obtain a unique tuple (Z₁, . . . , Zₖ). There are several strategies that can ensure each map component satisfies these monotonicity requirements. We will explore these strategies in the following.

**Monotonicity through integration** A general, powerful, but also computationally demanding method to enforce monotonicity makes use of a rectifier and integration (e.g., Baptista et al., 2023). In this formulation, our starting point is a smooth, differentiable, but nonmonotone function Sₖⁿᵒⁿ : ℝᵏ → ℝ. An example of such a nonmonotone function is:

`S₃ⁿᵒⁿ(x₁, x₂, x₃) = c₀ + c₁x₁ + c₂x₂ + c₃x₁x₂ + c₄x₁²` (nonmonotone part g(x₁,x₂)) `+ c₅x₃ + c₆x₁x₃² + c₇x₂x₃` (pre-monotone part ĝ(x₁,x₂,x₃)) (17)

where the coefficients {cᵢ} parameterize Sₖⁿᵒⁿ. In Equation (17), we have conceptually separated the terms of S₃ⁿᵒⁿ into a nonmonotone part g, which does not depend on x₃, and a “pre-monotone” part ĝ, which does. In general, ĝ will not be monotone in x₃, and consequently neither will be S₃ⁿᵒⁿ. However, we can monotonize it by first applying a rectifier r : ℝ → ℝ⁺ to ĝ, then integrating the rectified output over xₖ:

`Sₖ(x₁, . . . , xₖ) = g(x₁, . . . , xₖ₋₁) + f(x₁, . . . , xₖ)`
`where f(x₁, . . . , xₖ) = ∫₀ˣᵏ r(ĝ(x₁, . . . , xₖ₋₁, t)) dt` (18)

The two steps combined in Equation (18) are as follows:

1.  **Rectification**: First, ĝ is sent through the rectifier r, a function which maps its output to strictly positive values. Examples of useful rectifiers include exponential functions, the softplus function (i.e., log(1 + exp(x))) or Exponential Linear Units (Clevert et al., 2016).
2.  **Integration**: Numerically or analytically integrating the output of the resulting strictly positive function over xₖ yields a function f monotone in xₖ, which in turn also guarantees the monotonicity of Sₖ in xₖ.

This process is visualized in Figure 8, and has a number of important advantages. First off, the use of an integrated rectifier places very few limitations on the structure of Sₖⁿᵒⁿ. Furthermore, it also permits the introduction of cross-terms such as x₁xₖ between the last argument xₖ and previous dimensions x₁:ₖ₋₁, which are required for many of the more complex mapping operations (see Section 3.1.2). The drawback of this approach is that in practice evaluating Equation (18) often demands one-dimensional numerical integration via, e.g., quadrature, increasing the computational demand of the map’s evaluation, optimization, and inversion.

An important practical consideration when using the integrated map formulation is to ensure that the pre-monotone term ĝ reverts to a constant in the tails of xₖ. This can be ensured by defining ĝ as a combination of a constant term and Hermite functions (see Section 3.1.3), which revert to zero in the tails. To understand why this is important, consider the process illustrated in Figure 8. Through rectification and integration, increasingly negative values of ĝ become flat stretches in Sₖ. If one or both tails of the monotone term f become near-flat, it is possible for the root finding during the map’s inversion (see Section 2.3.1) to fail for outlying values. However, if ĝ instead reverts to a positive constant in the tails, the resulting integrated monotone term f will extrapolate linearly and thus generally remain more robust to outliers during the inversion. This is an artifact from the discussion above where, while π and η may “look” very different for the bulk of the pdf, we expect the tails of both of them to behave similarly (Baptista et al., 2023).

[Image of Figure 8: Integrated rectifiers create monotonicity. (A) We start with an arbitrary smooth, non-monotone function ĝ(x). (B) Applying a monotone rectifier, for example exp(ĝ(x)), yields strictly positive function output. (C) Integrating a strictly positive function yields a monotone function f. Color indicates the magnitude of the nonmonotone function in (A).]

**Monotonicity through variable separation** A more computationally efficient way to ensure monotonicity is to formulate a map component function Sₖ which is both separable in xₖ and linear in the coefficients, though not generally linear in the inputs x. An example for such a map component function is provided below:

`S₃(x₁, x₂, x₃) = c₀ + c₁x₁ + c₂x₂ + c₃x₁x₂ + c₄x₁²` (nonmonotone part g(x₁,x₂)) `+ c₅x₃ + c₆erf(x₃)` (monotone part f(x₃)) (19)

where cᵢ once more are the map component’s coefficients, which are optimized to find the map from π to η, and c₅, c₆ > 0. The key aspects in Equation (19) are a clear separation into a nonmonotone term g, which may depend on all arguments except the last x₁:ₖ₋₁, and a monotone term f, which may only depend on the last argument xₖ. To ensure that Sₖ remains monotone in xₖ, all terms in f must likewise be monotone basis functions (e.g., xₖ and erf(xₖ)).

The advantages of this separable formulation are two-fold: First, Equation (19) has no need for numerical integration, which reduces computational demand substantially. Second, as we shall see in Section 3.2, this separable formulation allows for extremely efficient map optimization. The price for this efficiency is that variable separation does not allow for cross-terms between the last argument and previous arguments (e.g., x₁xₖ), which limits the map’s expressivity. We will explore the consequences of this limitation in the next section.

#### 3.1.2 Parameterization

Closely related to the concern of monotonicity is how each map component Sₖ relates xₖ, its last argument, to x₁, . . . , xₖ₋₁. This decision leads to three different kinds of map parameterizations, which permit maps of different complexity. In ascending order of complexity, these parameterizations are:

**Marginal maps** If both the nonmonotone part g and the monotone part f depend exclusively on input xₖ, such map component functions Sₖ result in a marginal (or diagonal) map. As the name implies, such maps only transform the marginals of x without capturing any dependencies between its dimensions. An example of a marginal map would be:

`S₃(x₃) = c₀` (nonmonotone part g(-)) `+ c₁x₃ + c₂x₃³` (monotone part f(x₃)) (20)

Recalling the implications of removing dependencies on earlier arguments x₁:ₖ₋₁ from the map component functions (see Section 2.3.3), marginal maps are of no direct use for the conditioning operations in Section 2.3.2. This is because marginal maps implicitly assume all entries of x are marginally independent, which in turn implies that π(xₖ₊₁:ₖ|x₁:ₖ) = π(xₖ₊₁:ₖ), that is to say, nothing can be learned from x₁:ₖ about xₖ₊₁:ₖ. As marginal Gaussianization schemes, however, they can find limited use in Gaussian anamorphosis (Schöniger et al., 2012; Zhou et al., 2011) or as preconditioning tools in certain graph detection methods (Liu et al., 2009). Marginal maps are only included here for completeness’ sake, and we use g(-) to denote that the nonmonotone part is constant in x, connecting better to the more complex parameterizations.

**Separable maps** Closely related to the monotonicity scheme of the same name in the previous section, separable maps separate the map component Sₖ into the sum of a function g : ℝᵏ⁻¹ → ℝ that takes in x₁:ₖ₋₁ and a univariate monotone function f : ℝ → ℝ that just takes in xₖ:

`S₂(x₁, x₂, x₃) = c₀ + c₁x₁ + c₂x₁² + c₃x₁x₂` (nonmonotone part g(x₁,x₂)) `+ c₄x₃ + c₅x₃³` (monotone part f(x₃)) (21)

Separable maps are one level of complexity above marginal maps. Such maps can consider nonlinear dependencies between the different dimensions, but do not permit cross-terms with the last argument xₖ, which in turn limits the complexity of the target distributions π they can recover⁹. This limitation is subtle, and will be discussed shortly. An important advantage of separable maps is that, when linear in the coefficients, their optimization can be partially solved in closed-form, which can improve computational efficiency substantially. More detail on this is provided in Appendix A.

**Cross-term maps** The most versatile map parameterization are cross-term maps, which permit the greatest expressiveness at the cost of increased computational demand. As mentioned previously, cross-terms are basis functions which depend both on the last argument xₖ and preceding arguments x₁:ₖ₋₁. The presence of these terms requires that such maps must use integrated rectifiers (Section 3.1.1) to ensure the monotonicity of Sₖ:

`S₂(x₁, x₂) = c₀ + c₁x₁ + c₂x₁²` (nonmonotone part g(x₁)) `+ ∫₀ˣ² exp(c₃x₁²t + c₄t + c₅x₁t²) dt` (monotone part f(x₁,x₂)) (22)

The presence of cross-terms permits increased control over the local details in the transformations from π to η. In principle, it is also possible to make Equation (22) separable in xₖ by removing the dependencies of f on x₁:ₖ₋₁. Since the resulting Sₖ would not be linear in the coefficients c, however, we cannot make use of a more efficient optimization scheme (Appendix A). In the following subsection, we will develop some intuition about the strengths and limitations of different map parameterizations.

[Image of Figure 9: Pullback approximations S♯η to different target distributions π using nonlinear map components Sₖ for marginal maps, separable maps, and cross-term maps. Some distributions require cross-term maps, for others simpler parameterizations may suffice. The variable ordering is \[x₁, x₂], with x₁ plotted on the horizontal axis and x₂ on the vertical axis. Color represents the position of the reference samples z relative to the mean of η. The number of coefficients (i.e., parameters) for each map is plotted in the bottom-right corner.]

**Choosing a parameterization** To provide more insight into the effects and limitations of each parameterization choice, we have illustrated the pullback S♯η of a transport map approximated using these parameterizations for each of three different target pdfs π in Figure 9. In every case, we begin by sampling the target distributions, proceed to optimize the different maps, then apply the inverse map (Section 2.3.1) to transform reference samples z ∼ η into samples from the pullback x ∼ S♯η. We may make the following observations:

*   **Marginal maps** only reproduce the marginal densities of all three target distributions. As their structure implies independence between the marginals (Section 2.3.3), however, the pullback approximation does not approximate the true target pdfs well. This effect can be subtle: in the corner-multimodal distribution, it is only noticeable through the emergence of a fourth phantom mode.
*   **Separable maps** provide much better approximations to the target pdfs, but have some interesting caveats: while they successfully recover the wavy target distribution, they struggle with the distribution once it is rotated by 90°, and they yield a slight – if less pronounced – fourth phantom mode for the corner-multimodal target. We will discuss the reason for this below.
*   **Cross-term maps** prove most versatile, recovering all three target distributions well at the cost of a larger parameter space.

In practice, cross-terms provide greater control over local features. To understand why, let us take a closer look at a single map component inversion S⁻¹ₖ, which forms the basis of the map inversion (see Section 2.3.1) that ultimately yields the pushforward samples in Figure 9. Consider a generic map component function comprised of a nonmonotone part g and a monotone part f:

`zₖ = Sₖ(x₁:ₖ₋₁, xₖ) = g(x₁:ₖ₋₁) + f(x₁:ₖ₋₁, xₖ)` (23)

As discussed in Section 2.3.1, the inversion of Sₖ relies on one-dimensional root finding. Reformulating this expression yields the root finding objective for S⁻¹ₖ, conditioned on the outcomes x₁:ₖ₋₁ of previous map component inversions:

`zₖ − g(x₁:ₖ₋₁) = f(x₁:ₖ₋₁, xₖ)` (24)

If we now apply the simplifications of the three map parameterizations above, we obtain three different objectives for the root finding problem:

`zₖ − g(-) = f(xₖ)` (marginal maps)
`zₖ − g(x₁:ₖ₋₁) = f(xₖ)` (separable maps)
`zₖ − g(x₁:ₖ₋₁) = f(x₁:ₖ₋₁, xₖ)` (cross-term maps) (25)

These equations provide some insight into the differences between the pullback densities S♯η we have observed in Figure 9. Each term above takes on a different role during the conditional inversion: the reference samples zₖ are an independent input, unaffected by the map parameterization choice, and define a random initial offset for the root finding. The nonmonotone term g acts as an additional dynamic offset for the inversion based on previous values, and the monotone term f defines the shape of the inverse function for the one-dimensional root finding over the unknown xₖ. From this perspective, the differences between the map parameterizations emerge from how each handles the dependence on previous values x₁:ₖ₋₁:

1.  **Marginal maps** (Figure 10, top row) permit only a constant nonmonotone term g(-), which results in no dynamic offset. Their monotone term f(xₖ) likewise only depends on xₖ. In consequence, marginal maps extract the same xₖ from a specific zₖ for any given value of x₁:ₖ₋₁.
2.  **Separable maps** (Figure 10, center row) also feature monotone term f(xₖ) that depends only on xₖ, which keeps the inverse function’s shape constant for all values of x₁:ₖ₋₁. However, their off-diagonal term g(x₁:ₖ₋₁) can induce varying offsets for different x₁:ₖ₋₁, which introduces a simple dependence on previous values.
3.  **Cross-term maps** (Figure 10, bottom row) can vary both the dynamic offset (nonmonotone term g(x₁:ₖ₋₁)) and the inverse function’s shape (monotone term f(x₁:ₖ₋₁, xₖ)) with different x₁:ₖ₋₁, and thus provide the most expressive maps, but often at higher computational cost.

These insights reveal why separable maps succeeded to recover the wavy distribution, but failed to capture its features once it is rotated (Figure 9). For the original wavy distribution, its (vertical) conditionals at different horizontal positions are mostly of the same shape, just shifted vertically, which makes them a perfect fit for recovery via separable maps. Once the distribution is rotated, however, the true vertical conditionals become more challenging to recover. Note that if we were to change the order of the target vector to \[x₂, x₁] instead of \[x₁, x₂], separable maps would succeed in recovering the rotated wavy distribution and instead fail for the standard one. This illustrates the subtle influence of input variable ordering on a map’s approximation ability.

[Image of Figure 10: Conditional inversion x₂ = S⁻¹₂(z₂; x₁) for marginal (top row), separable (center row), and cross-term (bottom row) maps. Each column corresponds to a conditional distribution for a different fixed x\*₁. Each term in Equation (25) has a different effect: z₂ defines the horizontal starting position for the inversion, f(x₁, x₂) defines the shape of coupling, and g(x₁) defines its horizontal offset. Marginal maps keep f and g fixed, and thus yield the same inverse x₂ for any previous x₁. Separable maps also keep f constant, but can adjust the offset g. Finally, cross-term maps can adjust both f and g for different x₁, providing the greatest degree of flexibility.]

#### 3.1.3 Basis functions

A very important part of the parameterization of the map component function Sₖ is its construction: what basis functions should be used to represent Sₖ and how much nonlinearity should they permit? In a nutshell, simpler – perhaps separable – parameterizations are more computationally efficient, but can struggle to recover features of the target pdf π that do not resemble the reference η, as discussed at the beginning of this section. By contrast, more complex parameterizations add the flexibility required to capture “localized” features of π (e.g., skewness, multi-modality, . . . ), but risk unfavourable bias-variance trade-offs. Depending on the properties of the problem of interest, different parameterization choices promise different advantages:

**Polynomial basis functions** In function approximation, polynomials are a canonical choice for forming an approximation; in particular, harkening to polynomial chaos expansion and other traditional stochastic function approximation literature (Le Maître and Knio, 2010; Ernst et al., 2012), we can choose orthogonal polynomials for our approximation function class. Since this work focuses on unbounded x, we will consider the probabilists’ **Hermite polynomial** family {Heⱼ}∞ⱼ=₁, which has the orthogonality property

`∫₋∞∞ Heᵢ(x)Heⱼ(x) (1/√τ)exp(-x²/2) dx = j!δᵢⱼ` (26)

where τ = 6.283185 . . . is twice Archimedes’ constant. Broadly, Equation (26) shows what it means for polynomials of different orders (e.g., linear, quadratic, cubic, etc.) to be orthogonal under the Gaussian weight. Similar to orthogonal vectors in linear algebra, orthogonal polynomials (more generally, orthogonal functions) are often chosen as a basis for function approximation, as their elements contain no “overlapping information” with respect to a particular weight (which may or may not match η). As these polynomials are dependent on the pdf, a few other orthogonal polynomial families with their respective weights and support include Legendre polynomials (respecting a constant weight on (−1, 1)) and Laguerre polynomials (respecting weight e⁻ˣ on (0, ∞)) (Szegő, 1939). These polynomial families may be suited for different problems depending on characteristics of η and π (Wang and Marzouk, 2022).

Generally, if a polynomial family {ψⱼ}∞ⱼ=₁ is orthogonal with respect to a weight ρ (e.g., the probabilist Hermite polynomials take ρ as the Gaussian weight), one can approximate a given well-behaved function f(x) with f̂(x) = c₁ψ₁(x) + c₂ψ₂(x) + · · · + cₚψₚ(x) using a finite number P of these polynomials. Guarantees for this approximation are often given when measuring the error according to our weight ρ via ∫(f(x) − f̂(x))² dρ(x) (Ernst et al., 2012). If we believe that η and π are very similar, choosing polynomials orthogonal with respect to our reference η (which we know the orthogonal polynomial family of) will also perform well when the input is distributed according to π (which we generally won’t know the orthogonal polynomial family of). On the other hand, if η and π are not very similar, such guarantees are not particularly helpful; practitioners see this manifest as a need for more complex map parameterizations (see Section 3.1.2).

As they have unbounded support, Hermite polynomials exert global influence and often require high-order terms to capture fine distributional features. While versatile, the use of polynomial basis functions has a few pitfalls which demand caution:

1.  **Importance of standardization**: In response to the above, we would like to make π resemble η more closely by pre-processing the data we have from π. When working with a Gaussian reference η, it is thus generally advised to standardize x by trying to match the mean and variance of η. We do this by (i) subtracting the empirical mean X̄ and (ii) scaling each marginal to unit variance by dividing it through the unbiased estimator of the empirical standard deviation σₖ before formulating and applying a transport map S. All subsequent map operations are then implemented in this standardized space. Finally, the map’s outputs are transformed back into target space by reversing the standardization. In effect, standardization wraps a separate linear transformation around the transport map, and thereby does not affect the validity of the operations inside this wrapper. In the subsequent sections, we will assume by default that the samples have been standardized. As pointed out in Morrison et al. (2022), invertible diagonal transformations retain the conditional independence structure of the target distribution. In other words, there is no loss of information while working in the marginally rescaled space.
2.  **Edge control**: A practical challenge of polynomials is their growth in the tails for higher-order terms. This often leads to volatility if the map is evaluated or inverted far from zero. To address this issue, one useful option is the use of edge-controlled terms such as **Hermite functions** Hⱼ. This variant of basis functions are defined as probabilist’s Hermite polynomials Heⱼ of order j multiplied with a Gaussian weight, which reverts the polynomial’s output to zero far from zero:

    `Hⱼ(x) = Heⱼ(x)exp(-x²/4)` (27)

    A nice feature of this form is that Hᵢ(x)Hⱼ(x) = Heᵢ(x)Heⱼ(x)exp(-x²/2), and so {Hⱼ}∞ⱼ=₁ will inherit the same error properties as {Heⱼ} when measuring error without any weight ρ, i.e., ∫(f(x) − f̂(x))² dx.

Illustrations of the first few orders of Hermite functions are provided in Figure 11A. Note that expressing the map exclusively in terms of Hermite functions causes the map component Sₖ to revert to zero in the tails, which can be undesirable whenever extrapolation may be required. In practice, limiting the Gaussian weight term to polynomials of order two or larger is a good compromise, thus we recommend retaining the linear terms without a weight.

A practical issue with Hermite functions can emerge for target distributions π which include sharp features. The largest maximizer x\*ⱼ of Heⱼ grows at a rate of x\*ⱼ = √j + O(j⁻¹ᐟ⁶), which means that if Hⱼ has nontrivial values on (−r, r), we will have that H₄ⱼ has nontrivial values on (−2r, 2r); e.g., it is reasonable to estimate from Figure 11(A) that r > 4 for j = 2 and so H₈ has nontrivial values on (−8, 8). Due to this “global” nature of polynomials influencing the edges of the Hermite functions, maps for such distributions often involve high-order terms with high-magnitude coefficients of opposing signs partially cancelling each other. While this is often unproblematic within the support of the training samples, the resulting high-magnitude coefficients can extend the influence of some Hermite functions farther into the tails, resulting in undesirable tail effects. As an alternative, we might prefer to employ a weight term that reverts the weights to zero at a finite distance. For instance, using a cubic spline weight yields a different edge-controlled Hermite polynomial

`HECⱼ(x) = Heⱼ(x)(2min(1, |x|/r)³ - 3min(1, |x|/r)² + 1)` (28)

[Image of Figure 11: Two types of edge-controlled basis functions based on Hermite polynomials. (A) Hermite functions are multiplied with a Gaussian weight, causing them to approach zero in the tails. (B) Replacing the Gaussian weight with a scaled cubic spline leads to a different edge-controlled basis function. This enforces limited support, ensuring the function reverts to zero at a finite distance r.]

**Radial basis functions** A second useful class of basis functions is radial basis functions (RBFs). RBFs exert more local influence than combinations of polynomial basis functions, and generally require two additional parameters: a position parameter µ and a scale (or bandwidth) parameter σ. To avoid the need to optimize these parameters along with the map’s coefficients cᵢ, one can estimate µ and σ from the empirical quantiles along the marginals of the training samples {Xᵢ}ᴺᵢ=₁ (Spantini et al., 2022). For the positions µᵢ, they propose to place each RBF at the empirical quantiles qᵢ/(j+1)(xₖ) for i = 1, . . . , j, where j is the total number of RBFs along the k-th marginal. Based on the resulting µᵢ, the corresponding scales σᵢ are determined by averaging the distances to neighbouring RBFs. To simplify notation, it can be useful to define local coordinates for each RBF:

`xˡᵒᶜ,ⁱₖ = (xₖ - µᵢ) / √(τσᵢ)` (29)

where the τ is again twice the Archimedes’ constant. Radial basis functions are particularly useful for multimodal distributions because they render the map selectively expressive wherever the target distribution’s samples are located. Where a superposition of many high-order polynomials may be required to recover isolated modes, often only a few RBFs may suffice. Besides conventional RBFs, three useful related basis functions are left edge terms (LET), integrated RBFs (iRBF), and right edge terms (RET):

`RBF(xˡᵒᶜ,ⁱₖ) = (1/√τ)exp(-(xˡᵒᶜ,ⁱₖ)²)`
`iRBF(xˡᵒᶜ,ⁱₖ) = (1/2)(1 + erf(xˡᵒᶜ,ⁱₖ))`
`LET(xˡᵒᶜ,ⁱₖ; σᵢ) = (1/2)(σᵢ/√(τxˡᵒᶜ,ⁱₖ))(1 - erf(xˡᵒᶜ,ⁱₖ)) - (4σᵢ/√τ)exp(-(xˡᵒᶜ,ⁱₖ)²) `
`RET(xˡᵒᶜ,ⁱₖ; σᵢ) = (1/2)(σᵢ/√(τxˡᵒᶜ,ⁱₖ))(1 + erf(xˡᵒᶜ,ⁱₖ)) + (4σᵢ/√τ)exp(-(xˡᵒᶜ,ⁱₖ)²) `(30)

Illustrations of these functions are provided in Figure 12. Note that iRBFs, LETs, and RETs are monotone functions, and thus an excellent choice for a linear separable map (see Section 3.1.2). Other sigmoid-like functions and monotone “edge terms” that revert to linear in the tails include logistic and softplus functions, which may be computationally faster.

[Image of Figure 12: Illustration of different useful RBF basis functions.]

### 3.2 Optimization

With the structure and parameterization of the map component functions Sₖ defined, we may now address the question of how to identify the specific map S that relates the target distribution π to the reference distribution η. In general, we have two different ways to estimate the map, depending on the type of information available: **maps from samples**, and **maps from densities**. Both follow a very similar approach, but differ fundamentally in the direction in which they define the transport map:

**Maps from samples**
*   require target samples x ∼ π
*   require evaluations of the reference η
*   seek a map S such that S♯π = η

**Maps from densities**
*   require reference samples z ∼ η
*   require evaluations of the unnormalized target π̃
*   seek a map R such that R♯η = π

#### 3.2.1 Maps from samples

In the preceding sections, we have often implicitly assumed a **maps from samples** scenario. In other words, the target pdf π is only known through samples, i.e., we have a collection of samples Xᵢ ∼ π. In this case, we seek to minimize the “distance” between the target pdf π and the pullback pdf S♯η. While several notions of distance between pdfs are possible, the **Kullback–Leibler divergence** is commonly used for its computational tractability. For two arbitrary pdfs p, q and a generic RV a, we define their Kullback–Leibler divergence D(p∥q) as:

`D(p∥q) = Eₚ[log(p(a)/q(a))] = ∫p(a)log(p(a)/q(a))da` (31)

Thus, we can formulate our optimization problem as finding the map Sₒₚₜ that minimizes the difference, measured in terms of the Kullback–Leibler divergence, between the target π and its approximation S♯η ≈ π:

`Sₒₚₜ = arg min_{S∈F} D(π∥S♯η)` (32)

where F is some appropriate class of functions (see Section 3.1). This is often referred to as the minimization of the **forward Kullback–Leibler divergence**. We note that D(π∥S♯η) can be expanded as:

`D(π∥S♯η) = Eπ[log(π(x)/S♯η(x))] = Eπ[logπ(x)] - Eπ[logS♯η(x)]` (33)

where the first term Eπ[logπ(x)] does not depend on the map S. Thus, we can alternatively view Equation (32) as an attempt to maximize the log-likelihood (or minimize the negative log-likelihood) of the map’s pullback density S♯η ≈ π over the target density π:

`Sₒₚₜ = arg min_{S∈F} Eπ[-logS♯η(x)]` (34)

[Image of Figure 13: The optimization objective of maps from samples maximizes the log-likelihood of a set of fixed samples x ∼ π over the map’s pullback density S♯η ≈ π. In practical terms, the optimization seeks maps S which mold the reference pdf η to the target samples x.]

Substituting in Equation (4) for S♯η(x) into the cost function Equation (34) and expanding the logarithms yields:

`Eπ[-logS♯η(x)] = Eπ[-logη(S(x)) - log|det∇ₓS(x)|]` (35)

Due to the lower triangular structure of S, we know that ∇ₓS is a lower triangular matrix. As established in Equation (6), the determinant of a lower triangular matrix is given by the product of its diagonal entries:

`log det∇ₓS(x) = log[Πₖ(∂Sₖ(x)/∂xₖ)] = Σₖ log(∂Sₖ(x)/∂xₖ)` (36)

Mind that due to the monotonicity of Sₖ in xₖ (see Section 3.1.1), the derivative ∂ₓₖSₖ will always be positive, and the logarithm in the right-hand side sum of Equation (36) will be defined. Plugging this identity into the cost function, we only now enforce that our reference η is multivariate Gaussian by substituting its pdf into the expectation.

`Eπ[-logS♯η(x)] = Eπ[-log(1/τ^(D/2)) + (1/2)Σₖ Sₖ(x)² - Σₖ log(∂Sₖ(x)/∂xₖ)]` (37)

where τ is again twice the Archimedes’ constant. Recognizing that the first term is constant with respect to S, it can be discarded from the objective function. Using N samples from the target distribution Xᵢ ∼ π and merging the two sums over K, we obtain a Monte Carlo approximation of the loss function:

`J(S) = Σᵢ Σₖ [(1/2)Sₖ(Xᵢ)² - log(∂Sₖ(Xᵢ)/∂xₖ)]` (38)

From this expression, we can recognize that the summands of the full objective function J(S) are independent of each other. As a consequence, we can define separate objective functions Jₖ(Sₖ) for each map component Sₖ:

`Jₖ(Sₖ) = Σᵢ [(1/2)Sₖ(Xᵢ)² - log(∂Sₖ(Xᵢ)/∂xₖ)]` (39)

Equation (39) may then be minimized with off-the-shelf software developed for other optimization and machine learning tasks. If the linear separable formulation in Section 3.1.1 is chosen, it is only necessary to optimize the monotone part f of each map component Sₖ numerically, as the coefficients for the nonmonotone terms g can be derived in closed form. More detail on this is provided in Appendix A. Note that since the objective functions in Equation (39) are independent of each other, the map component functions Sₖ can be optimized in arbitrary order, even in parallel. Each objective function Jₖ balances a first quadratic term that drives samples Sₖ(Xᵢ)² to zero with a second term that acts as a log-barrier function to prevent null diagonal entries in the gradient of the transport map (Le Provost et al., 2021).

#### 3.2.2 Maps from densities

Alternatively, we choose a different optimization strategy in the case that we can evaluate the target pdf π, but only up to a constant of proportionality; i.e., we have access to evaluating density π̃ = mπ for some unknown m > 0. A key difference when learning a map from density is that we choose to define the triangular map in the opposite direction as before: where we previously considered a forward map S that maps target samples X to reference samples Z, we now define a forward map R mapping reference samples Z to target samples X. Correspondingly, the pushforward distribution R♯η now approximates the target distribution π via the forward map R, and the pullback distribution R♯π now approximates the reference distribution η via the inverse map R⁻¹. To find the optimal Rₒₚₜ, we seek to minimize the Kullback–Leibler divergence between the reference distribution η and the pullback distribution R♯π (Marzouk et al., 2017):

`Rₒₚₜ ∈ arg min_{R∈F} D(η∥R♯π)` (40)

where F is some appropriate class of functions (see Section 3.1.3) and ∈ is used to suggest multiple minimizers (see Baptista et al. (2023) for a discussion of why this could happen). This is often referred to as the minimization of the **reverse Kullback–Leibler divergence**. An illustration of the minimization in Equation (40) is provided in Figure 14.

[Image of Figure 14: The optimization objective of maps from densities maximizes the log-likelihood of a set of fixed samples Z ∼ η over the map’s pullback density R♯π̃ ≈ η. In practical terms, the optimization seeks maps R which mold the unnormalized target pdf π̃ to the reference samples Z.]

Similar to Equation (33) in the maps from samples scenario, we expand D(η∥R♯π) as:

`D(η∥R♯π) = Eη[log(η(z)/R♯π(z))] = Eη[logη(z)] - Eη[logR♯π(z)]` (41)

where the first term Eη[logη(z)] does not depend on the map R. Thus, Equation (40) can be viewed as an attempt to maximize the log-likelihood (or minimize the negative log-likelihood) of the map’s pullback density R♯π ≈ η over the reference density η:

`Rₒₚₜ ∈ arg min_{R∈F} Eη[-logR♯π(z)]` (42)

We may now substitute the pullback density from the change-of-variables formula (Equation (4)):

`η(z) ≈ R♯π(z) = π(R(z))|det∇₂R(z)|` (43)

Equivalent to Equation (36), we express the log determinant of the forward map R as a product of its diagonal entries:

`log det∇₂R(z) = log[Πₖ(∂Rₖ(z)/∂zₖ)] = Σₖ log(∂Rₖ(z)/∂zₖ)` (44)

Substituting the two equations above into the cost function of Equation (42), we obtain:

`Eη[-logR♯π(z)] = -Eη[logπ(R(z))] - Eη[Σₖ log(∂Rₖ(z)/∂zₖ)]` (45)

Since we assume we only have access to an unnormalized target density π̃, we can expand the first term Eη[logπ(R(z))] using the identity π = π̃/m, where m is the unknown normalization factor:

`Eη[-logR♯π(z)] = -Eη[logπ̃(R(z))] + Eη[logm] - Eη[Σₖ log(∂Rₖ(z)/∂zₖ)]` (46)

Recognizing that the normalization factor m does not depend on the map R, we can discard it from the objective function:

`Rₒₚₜ ∈ arg min_{R∈F} -Eη[logπ̃(R(z))] - Eη[Σₖ log(∂Rₖ(z)/∂zₖ)]` (47)

Using N samples from the reference distribution Zᵢ ∼ η, we obtain a Monte Carlo approximation of the loss function:

`J(R) = Σᵢ [-logπ̃(R(Zᵢ)) - Σₖ log(∂Rₖ(Zᵢ)/∂zₖ)]` (48)

This yields the objective function for optimizing maps from densities. Two comments are in order regarding this optimization problem. First, observe that opposed to maps from samples, we have to first generate an ensemble of N i.i.d. samples Zᵢ ∼ η. Second, since we now find a map R to pulls the reference back to the target, we still avoid evaluating the inverse map R⁻¹ when calculating the loss J(R) for some candidate map R ∈ F. However, opposed to optimizing a map from samples (Section 3.2.1), this objective function cannot generally be subdivided into individual objectives for each constituent map component Rₖ. As such, all map components are usually found via a single optimization problem.

We note that there are other loss functions that we could consider beyond the Kullback–Leibler divergence, for instance the variance diagnostic (El Moselhy and Marzouk, 2012; Richter et al., 2020). We do not consider these options here in greater detail, as the Kullback–Leibler objective is widely used and quite tractable. Moreover, in the maps-from-samples setting, a minimizer of the forward Kullback–Leibler divergence corresponds to a maximum likelihood estimator of the transport map within the chosen function class; this link is useful for both interpretation and theory (Wang and Marzouk, 2022).

#### 3.2.3 Regularization

When the number of samples N is small relative to either the target’s dimensionality K or the number of parameters for the map components, using more complex maps (i.e., increasing the number of the parameters) risks overfitting. In turn, this might lead to numerical instabilities. To prevent these issues, we can introduce an L₁ or L₂ regularization penalty to the objective functions introduced in the preceding section. In the maps from samples (Section 3.2.1) case, we take map component Sₖ for some fixed k, and let cᵢ parameterize this component via, e.g., polynomial coefficients. The regularized objective function would be:

`Jᴸ¹ₖ(Sₖ) = Σᵢ [(1/2)Sₖ(Xᵢ)² - log(∂Sₖ(Xᵢ)/∂xₖ)] + Σᵢ λᵢ|cᵢ|` L₁ regularization
`Jᴸ²ₖ(Sₖ) = Σᵢ [(1/2)Sₖ(Xᵢ)² - log(∂Sₖ(Xᵢ)/∂xₖ)] + Σᵢ λᵢcᵢ²` L₂ regularization (49)

where we add a λ-weighted penalty term to the objective functions, which penalizes non-zero parameters cᵢ for each of the i = 1, . . . , Pₖ basis functions of map component Sₖ (see Section 3.1.3). This creates an optimum closer to where the parameters are zero, which can be interpreted as penalizing map complexity. The regularization factors λᵢ are user-specified hyperparameters that should be tuned to the problem at hand, ensuring that the penalization is large enough to have an effect and small enough to avoid excessive biasing. For maps from densities, the optimization objective follows an equivalent structure to Equation (49). Similar to other applications with regularization parameters, a bit of experimentation is required to find a suitable balance between regularization and bias.

In practice, the use of L₁ regularization is a common strategy to discover sparsity in the triangular map, which we recall corresponds to conditional independence (see Section 2.3.3). In many cases, it is more efficient to impose sparsity and parsimony by construction. We will discuss practical strategies to learn a parsimonious degree of map complexity in Section 4.2.

## 4 Practical heuristics

With the theoretical foundations and implementation-related details established, we now discuss a few important features and tricks that help to further enhance the potential of triangular transport maps.

### 4.1 Composite maps for conditional sampling

For complex target distributions π, it is often infeasible to define and optimize a map of sufficient complexity to capture all features of π. In such cases, the map S will fail to completely normalize the target π, making certain features persist in the pushforward S♯π. At the same time, the pullback S♯η may not be a good approximation of π. Figure 15 illustrates this for two different levels of map complexity.

[Image of Figure 15: Simple maps capture only simple features. (A) For a non-Gaussian target π, a simple linear map does not yield a Gaussian pushforward S♯π ≠ η. (B) Likewise, its pullback S♯η ≠ π does not provide a good approximation to the target. (C & D) More complex nonlinear maps yield better approximations. Color in (A) and (C) indicates the sample’s original cluster membership, whereas in (B) and (D) color indicates a sample position’s angle w.r.t. η’s mean.]

This imperfection has important consequences for Bayesian inference with triangular transport maps. Since the (conventional) map’s conditional inverse only samples a conditional of the pullback distribution S♯η, not of the real target π, the quality of the resulting posterior samples will depend on the mismatch between S♯η and π. Yet, we are often given a sample X = (X₁:ₖ, Xₖ₊₁:ₖ) of a joint distribution π(x₁:ₖ, xₖ₊₁:ₖ) and would like to create a high-quality sample X\*ₖ₊₁:ₖ conditioned on some fixed x\*₁:ₖ. To do this, we apply a composite map to each sample x, defined as

`x\*ₖ₊₁:ₖ = S⁻¹ₖ₊₁:ₖ(· ; x\*₁:ₖ) ∘ Sₖ₊₁:ₖ(xₖ₊₁:ₖ; x₁:ₖ)` (50)

Note that x\*₁:ₖ and x₁:ₖ are two different objects: The former comes from the real world (i.e., it is a fixed value, independent from the random sample x). For instance, x\*₁:ₖ may be specific measurement values on which we want to condition π. By contrast, the latter are jointly sampled with xₖ₊₁:ₖ from π. While this may initially seem complicated, the key idea of composite maps is as follows:

1.  We have an approximate map S. Its forward evaluation thus does not transform samples from the true target x ∼ π into samples from the reference z ∼ η, but instead yields samples from the pushforward z̃ ∼ S♯π ≠ η. These pushforward samples z̃ preserve features of the target distribution π the map did not capture. For instance, a linear map S can only shift, scale, and rotate the target distribution π, which is insufficient if π is non-Gaussian (Figure 15A).
2.  Likewise, an approximate inverse map does not transform samples from the reference z ∼ η into samples from the true target x ∼ π, instead yielding samples from the pullback x̃ ∼ S♯η ≠ π which lack any features of π that the map has not captured. For instance, the pullback S♯η of a linear map S only samples a Gaussian approximation of the target distribution (Figure 15B).
3.  Applying an invertible approximate forward map, directly followed by its inverse, cancels the map’s approximation error and always restores the original target samples. By construction, we know that the composition of S with its inverse S⁻¹ maps a sample to itself, i.e., if z̃ = S(x) ∼ S♯π, then x = S⁻¹(z̃).
4.  This restoration of unresolved features persists in part during the conditioning operation. Instead of reference samples z ∼ η (Figure 16B), we might thus use pushforward samples z̃ ∼ S♯π to preserve some features of π. This generally reduces the map’s approximation error (Figure 16C); see Proposition 11 in Baptista (2022).

[Image of Figure 16: Composite maps yield better conditionals for imperfect maps. (A) The true solution for the conditioning operation. For the simple linear map, (B) using reference samples z ∼ η directly samples a conditional of the pullback S♯η ≠ π, which here is unrepresentative of the true conditional. (C) Using pushforward samples z̃ ∼ S♯η instead as part of a composite inversion yields a better approximation to the true conditional. The black line in the x\*₂ histograms in B, C, and D shows the true conditional distribution from subplot A for comparison.]

Figure 16 illustrates the composite map’s preservation of uncaptured features during conditional inversion for a Gaussian mixture target π (Figure 16A). The pullback approximation S♯η from a linear transport map is Gaussian, and thus clearly inadequate for a multimodal Gaussian mixture target π. Inverting the map with reference samples z ∼ η only samples conditionals of the (multivariate Gaussian) pullback (Figure 16B). However, using reference samples from the pushforward z ∼ S♯π preserves complex features of the target π, yielding better conditional samples x\*₂ (Figure 16C).

Keep in mind that if we instead use a sufficiently complex, nonlinear map formulation (e.g., Figure 15D), the map approximation retrieves the target distribution sufficiently well. In this case, there is little practical difference between the use of reference z ∼ η or pushforward samples z̃ ∼ S♯π, as both distributions will be quite similar (S♯π ≈ η). In the general case, however, we may not be certain about the quality of our map approximation. In consequence, we recommend the use of composite maps as a safer strategy for conditional sampling.

### 4.2 Map adaptation

Identifying a suitable degree of map complexity can be a challenging task. For one, conditional independence structures (see Section 2.3.3) are not always known a-priori. Similarly, finding a suitable degree of complexity for each map component to achieve an optimal bias-variance trade-off is not an easy task. To address these issues, we may draw on map adaptation algorithms. An example of such an algorithm is provided by Baptista et al. (2023).

In brief, the method starts with a triangular map that is viable and particularly simple; for example, one may use the identity map Sₖ(xₖ) = xₖ. The algorithm then proposes a list of candidate basis functions for each Sₖ (Figure 17). These candidate terms either (i) extend the dependence of Sₖ to a previous state dimension x₁:ₖ₋₁, or (ii) incrementally increase the complexity of an existing dependence by adding higher degrees of nonlinearity.

To determine which candidate basis functions are added to the map, the algorithm calculates the gradient of the optimization objective function (see Section 3.2) with respect to the candidates’ coefficients. It then adds the candidate corresponding to the steepest derivative. The idea behind this approach is that candidate terms with steep optimization objective gradients promise to improve the map’s approximation significantly, whereas those with flat gradients add little beyond complexity. The algorithm then proposes new candidates to the expanded map further and repeats the procedure until a user defined stopping criterion is met.

This process gradually expands the complexity and expressiveness of each map component. Since this process is iterative, this adaptation algorithm can be expensive if the number of candidate basis functions is large, as is the case for Sₖ with many variable dependencies, particularly if the map parameterization permits cross-terms (Figure 17). The use of maps without cross-terms can drastically reduce this computational demand, as it would restrict the expansion and exploration of candidates to basis functions corresponding to the outer edges in Figure 17, making a two-dimensional result look like the letter “L”.

[Image of Figure 17: Adaptive transport maps gradually expand the complexity and expressiveness of the map components, adding terms corresponding to the largest absolute gradient of the objective function (Equation (39)) among a set of candidate terms. This identifies the necessary degree of map complexity gradually. The lower-row grids illustrate the order and composition of the map components at each iteration. Each filled cell is an accepted candidate, with coordinates (3, 5) in the S₂ block corresponding to a term H₃(x₁)H₅(x₂).]

## 5 Summary & Outlook

### 5.1 State of the art

Measure transport methods comprise an active and quickly evolving field in statistics and machine learning. Within this broader field, triangular maps occupy an important niche: transparent and versatile, they provide a powerful set of tools for conditional sampling from limited information. This, in turn, makes them highly useful for Bayesian inference.

To summarize, triangular transport methods have a number of important advantages:

1.  **Parsimony**: The ability to fine-tune the parameterization of a triangular map, coupled with the sparse variable dependence that a triangular map inherits from conditional independence, allows triangular transport to be adapted to the computational demands, data availability, and distributional complexity of a given problem – increasing efficiency and accuracy.
2.  **Numerical convenience**: Triangular maps are simple to learn and invert, and many common parameterizations provide optimization guarantees for the problem of learning maps from samples (see, e.g., Appendix A and Baptista et al. (2023)).
3.  **Transparency**: Triangular maps provide a close correspondence between the map component functions Sₖ and a specific factorization of the target distribution π. This link makes it easy to predict the impact of changes to the map on the resulting statistical model, and often allows us to decompose an inference or sampling task into smaller problems that are themselves interpretable.

At the same time, triangular maps face some important challenges. Ensuring the scalability of triangular transport with the dimension of the target distribution demands exploiting conditional independence to produce a sparse map, which in turn requires an appropriate variable ordering. We reiterate that, under mild assumptions, a triangular map exists between any target and reference (Santambrogio, 2015, Section 2.3), which includes any permutation of our target random variables. In practice, however, these maps may vary in sparsity and other notions of complexity: see Section 2.3.3 and Figure 9, respectively, and note that the rotated “wavy” distribution in Figure 9 actually corresponds to a different variable ordering. There are many approaches to finding orderings that maximize sparsity: for example, minimum degree algorithms (e.g., Amestoy et al., 1996; Cuthill and McKee, 1969) from graph theory, or other algorithms inspired by sparse Cholesky factorizations (Baptista et al., 2024d; Schäfer et al., 2021). Setting aside sparsity, however, the impact of variable ordering on the difficulty of approximating a given map component function Sₖ is generally more difficult to predict, as it depends on finer properties of the target distribution at hand.

### 5.2 Broader context within transport

There are many ways to parameterize a triangular transport map that are not discussed here – for example, performing a closed-form integration over a parameterization of the target pdf. This can be done using squared polynomials (Zanger et al., 2024), parameterizing the square-root of the pdf (Dolgov et al., 2020), or even a composition of polynomials parameterizing this square-root pdf (Cui and Dolgov, 2022), often in the context of numerical tensor methods that permit efficient integration. Regarding the approximation of triangular transport maps, there are several works elucidating approximation rates by, e.g., polynomials or neural networks (Zech and Marzouk, 2022a,b; Baptista et al., 2024b; Westermann and Zech, 2023); other works analyse the statistical consistency and convergence of triangular maps learned from finite samples (Wang and Marzouk, 2022; Irons et al., 2022).

More generally, triangular maps sit within the broader field of measure transport, which has seen remarkable computational advances in recent years. Many other approaches to constructing transport maps exchange triangular structure for different assumptions or desiderata. For example, **optimal transport** methods seek a mapping between two distributions that is as “close” to the identity function as possible, where closeness is encoded by a particular integrated transport cost (for instance, the distance between the input and output of the map) (Peyré et al., 2019). Many works seek to approximate optimal transport maps by searching over an appropriate function class – e.g., by parameterizing classes of convex functions (rather than our monotone triangular formulation) and writing the map as the gradient of such a function (Makkuva et al., 2020; Wang et al., 2023). There are some connections between triangular maps and optimal transport maps. First, each component Sₖ of a monotone triangular map is optimal in its last (scalar) input xₖ, conditioned on the first k − 1 inputs x₁:ₖ₋₁. Conditional optimal transport maps (Tabak et al., 2021; Carlier et al., 2016; Baptista et al., 2024c; Pooladian et al., 2025) generalize this idea from strictly triangular to block triangular structure. As described in Carlier et al. (2009), triangular maps also arise as the limit of optimal transport maps obtained with increasingly anisotropic quadratic transport cost.

**Normalizing flows** (Papamakarios et al., 2021) are another widely used class of transport methods, which typically do not seek to approximate some canonical (i.e., optimal or triangular) map, but rather use a composition of many simpler invertible transformations to construct a transport map, increasing the number of functions (layers) in the composition until the desired expressivity is reached. Autoregressive normalizing flows in fact use triangular maps as a building block, with parameterizations ranging in complexity from relatively simple to the complex rectified formulations discussed in Section 3.1.1 (Wehenkel and Louppe, 2019; Jaini et al., 2019). These flows interleave such triangular layers with permutations of the variables, which mitigates ordering issues described above at the cost of sacrificing the conditional independence and sparsity properties of triangular maps. Such constructions are partly inspired by earlier autoregressive models that directly parameterize marginal conditional densities (Bond-Taylor et al., 2021). Finally, many contemporary flow-based generative models (Song et al., 2020; Albergo et al., 2023; Lipman et al., 2023; Liu et al., 2023) can be understood as dynamical or continuous-time representations of transport maps; here the transport is represented by the flow map of a system of ordinary differential equations (ODEs) or by the corresponding evolution of densities under a system of stochastic differential equations. Some models of this kind explicitly employ triangular structure (Heng et al., 2021), but triangular maps have also proven useful in the analysis of more general neural ODE models (Marzouk et al., 2024).

Within this broader landscape, what is the role of triangular transport? We believe that the greatest strength of triangular transport lies in its extraordinary capacity for **parsimony**: Triangular maps may not only naturally exploit conditional independence – a special property among transport-based samplers – but they allow us to fine-tune the level of nonlinearity with which we resolve statistical dependencies. In consequence, these methods allow us to tailor the maps precisely to the system’s demands by building some structure of the target distribution into the map itself. In the small sample size regime, we believe this property holds the key to challenging the prevalence of linear methods, which are still the state of the art for practical settings from meteorological data assimilation to subsurface parameter inference.

### 5.3 Where to go from here

As a relatively new method, much of the potential of triangular maps has yet to be realized. This includes a number of promising research directions in applications, methodology, and theory. To unlock the full potential of their parsimony, a key effort for the operationalization of triangular transport is the development of efficient adaptation strategies. One of the greatest challenges of the parameterizations discussed in this paper is their scaling with the dimension of the target distribution. One way to mitigate this “curse of dimensionality” is to only include terms one knows will help, which is the problem that the adaptive algorithms above attempt to tackle. As briefly discussed in Section 4.2, these algorithms automatically seek to identify parsimonious maps for general inference problems. Increasing the efficiency of these methods will ease their applicability to higher-dimensional systems. Promising research avenues include exploring connections to graph structure learning algorithms that estimate conditional independence structure (Baptista et al., 2024d; Liaw et al., 2025; Drton and Maathuis, 2017) or developing information criteria (Lunde, 2025; Konishi and Kitagawa, 1996) to help identify suitable levels of map complexity.

A closely-related subject is research into the properties, advantages, and drawbacks of different parameterizations of triangular transport. While we have explored many practical details in this tutorial, we have still only scratched the surface of this topic. A promising avenue is to focus on known computational bottlenecks of existing parameterizations. As a simple example, it is known the erf function (see Equation (30)) has unstable evaluations near the extremes of its domain and is often computationally expensive; exploring more efficient alternatives is of cross-cutting utility. Similarly, some parameterization heuristics such as centering RBFs at marginal quantiles of the target distribution remain largely ad hoc, especially given strong dependencies in π; what can we do to improve these heuristics? What kinds of parameterizations should be employed when the distributions at hand are heavy-tailed: do the “edge terms” of maps need to increase sub- or super-linearly (depending on the direction of the map), and does the use of composite maps relieve the need to model tails when performing conditional sampling? More broadly, is there a practical role for nonparametric representations of transport maps (Pooladian and Niles-Weed, 2021) in high-dimensional problems? And whether in parametric or nonparametric settings, how can regularization or penalization schemes for learning maps from limited information be more rigorously designed and scaled? Relatedly, it is common nowadays to have multi-fidelity data sources (e.g., models that can be resolved at different discretization levels); how can we efficiently learn maps via samples and model evaluations of differing fidelities?

Furthermore, there are many open questions at the intersection of triangular transport with numerical analysis, relevant to probabilistic modelling and uncertainty quantification. For example, there has been recent work marrying deterministic quadrature or quasi-Monte Carlo methods with transport (Cui et al., 2025a; Klebanov and Sullivan, 2023; Liu, 2024), promising higher-order convergence rates for the approximation of expectations weighted by complex distributions. It is interesting to consider how triangular structure could be leveraged in this setting and to understand broader links to cubature (Cools, 1997) or other high-dimensional integration schemes.

Finally, as we have argued that sparsity provides a route to scalability, it is important to understand the impact of approximate sparsity in transport, and hence approximate conditional independence in continuous non-Gaussian distributions – e.g., the error incurred by discarding weak conditional dependencies. (See Johnson and Willsky (2007); Jog and Loh (2015) for analyses in the Gaussian case.) This line of work has direct links to localization methods widely used in data assimilation, whose analysis is a topic of much current research (Al-Ghattas and Sanz-Alonso, 2024; Tong and Morzfeld, 2023; Gottwald and Reich, 2024). Moreover, sparsity is only one form of low dimensionality. Other work seeks scalable inference and probabilistic modelling via explicit dimension reduction – for instance, by searching for low-dimensional subspaces that capture important interactions between parameters and data, or by identifying low-dimensional “updates” from prior to posterior in the Bayesian setting (Baptista et al., 2022; Zahm et al., 2022; Constantine, 2015; Fukumizu et al., 2007; Härdle and Stoker, 1989; Cui et al., 2014). These notions correspond naturally to specific transport representations (Brennan et al., 2020; Cao et al., 2024; Cui et al., 2025b) and enable the solution of problems that would otherwise – without the exploitation of structure – be intractable.

### 5.4 Toolboxes

We hope that this tutorial provided an accessible introduction to the theory and implementation of triangular transport. As we conclude, we want to leave you, the reader, with some concrete numerical tools to explore your own applications of triangular transport – if you do not want to code your own. The code to reproduce the figures and examples in this tutorial is available on GitHub and has been implemented using the Triangular Transport Toolbox, written in Python. Another implementation is available via the Monotone Parameterization Toolbox (MParT) (Parno et al., 2022), which is written in C++ but features bindings to Python, Julia, and MATLAB, with a focus on core functionality and efficiency. A package often used is TransportMaps, a Python toolbox with many features appearing in modern papers. TransportBasedInference.jl has similarly been used in previous work on triangular transport and nonlinear data assimilation. Regarding other parameterizations, there is code to reproduce sum-of-squares polynomial transport in Julia using SequentialMeasureTransport.jl and square-root pdf tensor-based methods using the Deep Tensor Toolbox in MATLAB.

## 6 Acknowledgements

We would like to acknowledge and express our deep gratitude to Mathieu Le Provost, whose invaluable contributions and insights helped shape this work. Mathieu sadly passed away before the completion of this article, and his presence is greatly missed.

The research of MR leading to these results has received funding from the Swiss National Science Foundation under the Early PostDoc Mobility grant P2NEP2 191663 and the Dutch Research Council under the Veni grant VI.Veni.232.140. MR, MLP, and YM also acknowledge support from the Office of Naval Research Multidisciplinary University Research Initiative on Integrated Foundations of Sensing, Modeling, and Data Assimilation for Sea Ice Prediction under award number N00014-20-1-2595. MLP and YM acknowledge support from the National Science Foundation (award PHY-2028125). DS and YM acknowledge support from the US Department of Energy (DOE), Office of Advanced Scientific Computing Research, under grants DE-SC0021226 (FASTMath SciDAC Institute) and DE-SC0023188.

## References

*(The list of references is extensive and has been omitted for brevity in this response, but is present in the original document from page 38 to 43.)*

---

## A Appendix 1: Optimizing linear separable maps

As established in Section 3.1.2, choosing a linear separable map parameterization permits more efficient map optimization. First, for the purpose of this material we denote gradient of some scalar-valued function r(z) with respect to z evaluated at z\* as ∇₂r|z\* and the derivative of a univariate function s(z) with respect to scalar z evaluated at z\* as ∂₂s|z\*. Then, recall for a separable map that is linear in the coefficients, each map component function Sₖ is defined as:

`Sₖ(x₁, . . . , xₖ₋₁, xₖ) = g(x₁:ₖ₋₁) + f(xₖ)`
`= cⁿᵒⁿₖ,₁ψⁿᵒⁿₖ,₁(x₁:ₖ₋₁) + · · · + cⁿᵒⁿₖ,ₘψⁿᵒⁿₖ,ₘ(x₁:ₖ₋₁) + cᵐᵒⁿₖ,₁ψᵐᵒⁿₖ,₁(xₖ) + · · · + cᵐᵒⁿₖ,ₙψᵐᵒⁿₖ,ₙ(xₖ)`
`= Ψⁿᵒⁿₖ(x₁:ₖ₋₁)cⁿᵒⁿₖ + Ψᵐᵒⁿₖ(xₖ)cᵐᵒⁿₖ` (51)

where ψⁿᵒⁿₖ,ⱼ : ℝᵏ⁻¹ → ℝ and ψᵐᵒⁿₖ,ⱼ : ℝ → ℝ are the j-th basis functions of map component Sₖ, associated with the nonmonotone and monotone terms, respectively. Then, Ψⁿᵒⁿₖ : ℝᵏ⁻¹ → ℝ¹ˣᵐ and Ψᵐᵒⁿₖ : ℝ → ℝ¹ˣⁿ are vectors of basis function evaluations for g and f, and cⁿᵒⁿₖ ∈ ℝᵐˣ¹ and cᵐᵒⁿₖ ∈ ℝⁿˣ¹ are the corresponding column vectors of coefficients. Therefore, the expressions Ψⁿᵒⁿₖ(x₁:ₖ₋₁)cⁿᵒⁿₖ and Ψᵐᵒⁿₖ(xₖ)cᵐᵒⁿₖ are both inner product functions of x, i.e. scalar-valued. For maps from samples (Section 3.2.1), we consider samples X₁, . . . , Xₙ ∼ π and recall the following optimization objective for Sₖ:

`Jₖ(Sₖ) = Σᵢ [(1/2)Sₖ(Xᵢ)² - log(∂ₓₖSₖ|Xᵢ)]` (52)

Plugging in Equation 51, we obtain

`Jₖ(Sₖ) = Σᵢ (1/2)(Ψⁿᵒⁿₖ(Xᵢ₁:ₖ₋₁)cⁿᵒⁿₖ + Ψᵐᵒⁿₖ(Xᵢₖ)cᵐᵒⁿₖ)² - log(∂ₓₖΨᵐᵒⁿₖ|Xᵢₖ cᵐᵒⁿₖ)` (53)

Note that the samples Xᵢ are defined by the target distribution π, and are thus fixed during optimization. In consequence, the basis function evaluation vectors Ψⁿᵒⁿₖ and Ψᵐᵒⁿₖ are also fixed for a given map parameterization (see Section 3.1.2), and are thus independent of the coefficients cⁿᵒⁿₖ and cᵐᵒⁿₖ. We can simplify Equation 53 further by absorbing the sum over the first term, and defining a new variable for the partial derivative in the second term. To this end, we form matrices Pⁿᵒⁿₖ ∈ ℝᴺˣᵐ and Pᵐᵒⁿₖ ∈ ℝᴺˣⁿ, and vectors bᵢₖ ∈ ℝ¹ˣⁿ, which are each defined element-wise as

`[Pⁿᵒⁿₖ]ᵢⱼ = ψⁿᵒⁿₖ,ⱼ(Xᵢ₁:ₖ₋₁), [Pᵐᵒⁿₖ]ᵢⱼ = ψᵐᵒⁿₖ,ⱼ(Xᵢₖ), [bᵢₖ]ⱼ = ∂ₓₖψᵐᵒⁿₖ,ⱼ|Xᵢₖ`

where the ij entry is the jth indexed basis function evaluated at sample index i for both matrices Pⁿᵒⁿₖ and Pᵐᵒⁿₖ, and the log of the derivative of Ψᵐᵒⁿₖ with respect to Xₖ is evaluated at Xᵢₖ for the vectors bᵢₖ. As above, these matrices and vectors can be pre-computed. The optimization objective now simplifies further to

`Jₖ(Sₖ) = (1/2)∥Pⁿᵒⁿₖcⁿᵒⁿₖ + Pᵐᵒⁿₖcᵐᵒⁿₖ∥² - Σᵢ log(bᵢₖcᵐᵒⁿₖ)`

This function is similar to that of an objective for an interior point method and, remarkably, becomes quadratic in cⁿᵒⁿₖ. At this point, we add L₂ (i.e., Tikhonov) regularization on both cᵐᵒⁿₖ and cⁿᵒⁿₖ according to the guidance in Section 3.

`Jₖ(Sₖ; λ) = (1/2)∥Pⁿᵒⁿₖcⁿᵒⁿₖ + Pᵐᵒⁿₖcᵐᵒⁿₖ∥² - Σᵢ log(bᵢₖcᵐᵒⁿₖ) + (λ/2)(∥cⁿᵒⁿₖ∥² + ∥cᵐᵒⁿₖ∥²)` (54)

In practice, this means the optimal coefficients ĉⁿᵒⁿₖ for the nonmonotone basis function evaluations minimizing (54) must satisfy

`0 ≡ ∇cⁿᵒⁿₖ Jₖ(Sₖ; λ)|ĉⁿᵒⁿₖ = (PⁿᵒⁿₖᵀPⁿᵒⁿₖ + λI)ĉⁿᵒⁿₖ + PⁿᵒⁿₖᵀPᵐᵒⁿₖcᵐᵒⁿₖ` (55)

In consequence, for a given choice of coefficients parameterizing the monotone functions, cᵐᵒⁿₖ, we can find the optimal choice of ĉⁿᵒⁿₖ by solving the normal equations; for this scenario, we assume that m < N, i.e. the number of samples surpasses the number of basis functions (and thus PⁿᵒⁿₖᵀPⁿᵒⁿₖ is full rank). The normal equations are a well-studied class of problems in numerical linear algebra (Trefethen and Bau, 2022). Assuming linearly independent basis functions, this affords a solution for ĉⁿᵒⁿₖ as a function of cᵐᵒⁿₖ given as

`ĉⁿᵒⁿₖ = -(PⁿᵒⁿₖᵀPⁿᵒⁿₖ + λI)⁻¹PⁿᵒⁿₖᵀPᵐᵒⁿₖcᵐᵒⁿₖ = -Mₖ,λPᵐᵒⁿₖcᵐᵒⁿₖ` (56)

Following best practices from numerical linear algebra, this inversion should not be computed explicitly; rather, one should use a numerical solver for the systems induced. Then, substituting Expression (56) into the optimization objective, we obtain a new objective for cᵐᵒⁿₖ (i.e., entirely independent from the nonmonotone coefficients),

`Jᵐᵒⁿₖ(cᵐᵒⁿₖ; λ) = (1/2)∥Aₖ,λcᵐᵒⁿₖ∥² - Σᵢ log(bᵢₖcᵐᵒⁿₖ) + (λ/2)(∥Dₖ,λcᵐᵒⁿₖ∥² + ∥cᵐᵒⁿₖ∥²)` (57)

where Aₖ,λ, Dₖ,λ and bᵢₖ can be precomputed prior to the optimization routine via evaluation of the basis functions. Remarkably, this method translates the original loss function into a very simple constrained convex optimization problem

`ĉᵐᵒⁿₖ = arg min_{cᵐᵒⁿₖ≥0} (1/2)∥Aₖ,λcᵐᵒⁿₖ∥² - Σᵢ log(bᵢₖcᵐᵒⁿₖ) + (λ/2)(∥Dₖ,λcᵐᵒⁿₖ∥² + ∥cᵐᵒⁿₖ∥²)`
`ĉⁿᵒⁿₖ = -Mₖ,λPᵐᵒⁿₖĉᵐᵒⁿₖ` (58)

Note that the (element-wise) constraint of cᵐᵒⁿₖ ≥ 0 is vital to maintain monotonicity; this can be enforced explicitly during optimization using particular optimization algorithms, e.g., L-BFGS-B (i.e. Low-memory BFGS with box constraints), implicitly by constructing an objective with a log-barrier term, or employing a convex transformation of the optimization objective (e.g. optimize over pᵐᵒⁿₖ := log cᵐᵒⁿₖ). Since we often will have coefficients with zero values, the first methodology might be preferable (though no empirical results appear here). This optimization objective, notably, requires no evaluations of the map during optimization, that is to say, the optimization is as fast as the combination of the implementation of linear algebra algorithms called and the optimization routine used.

In the case where we have more parameters than samples, it is worth noting that using this formulation is remarkably sensitive to λ because PⁿᵒⁿₖᵀPⁿᵒⁿₖ is no longer full rank. Thus the calculation of Mₖ,λ solves a system with possibly poor numerical properties (the industry of methods for ill-conditioned systems is dedicated to such problems). For vanishingly small λ, this corresponds to the problem of overfitting and the fact that we have infinite choices of ĉⁿᵒⁿₖ for any given choice of cᵐᵒⁿₖ.