# Triangular Transport Methods (TTM) - Core Concepts and Implementation
Read `Theory.md` and 'roadmap.md' before starting any task.
## I. Introduction and Goal

**Primary Goal:** To characterize, approximate, and manipulate complex, high-dimensional target probability distributions (`pi`) by learning an invertible transformation (`S` or `S^-1`) that couples `pi` to a simple, user-defined reference distribution (`eta`), typically a standard K-dimensional Gaussian.

**Core Idea:**
1.  Relate a random variable `x ~ pi` (target) to `z ~ eta` (reference) via `z = S(x)` and `x = S^-1(z)`.
2.  The map `S` is learned/optimized based on information about `pi` and the chosen `eta`.
3.  `S` has a specific lower-triangular structure (Knothe-Rosenblatt rearrangement), providing mathematical and computational benefits.

**Key Notation:**
| Symbol         | Description                                       |
|----------------|---------------------------------------------------|
| `pi`           | Target probability distribution                   |
| `eta`          | Reference probability distribution (e.g., `N(0,I)`)|
| `x`            | K-dim random variable from `pi`, `x = [x_1,...,x_K]^T` |
| `z`            | K-dim random variable from `eta`, `z = [z_1,...,z_K]^T` |
| `S`            | Target-to-reference map: `R^K -> R^K`, `z = S(x)`   |
| `S_k`          | k-th component of `S`: `R^k -> R`, `z_k = S_k(x_1,...,x_k)` |
| `S^-1`         | Reference-to-target map (inverse of `S`)          |
| `S_k^-1`       | k-th component of `S^-1`                          |
| `K`            | Dimensionality of `x` and `z`                     |
| `grad_x S(x)`  | Jacobian matrix of `S` w.r.t. `x`                 |
| `det(...)`     | Determinant of a matrix                           |
| `S_#eta`       | Pullback of `eta` by `S^-1` (target approx.)      |
| `S_#pi`        | Pushforward of `pi` by `S` (reference approx.)    |
| `g(...)`       | Non-monotone part of `S_k`                        |
| `f(...)`       | Monotone part of `S_k` (in its last argument)     |
| `g_hat(...)`   | Pre-monotone part of `S_k` before integration     |
| `r(...)`       | Rectifier function `R -> R+`                      |

## II. Theoretical Foundations

### A. Change-of-Variables Formula
Relates densities `pi(x)` and `eta(z)`:
*   `pi(x) = eta(S(x)) * |det(grad_x S(x))|`  (1)
*   `eta(z) = pi(S^-1(z)) * |det(grad_z S^-1(z))|` (2)
The Jacobian determinant `|det(grad_x S(x))|` accounts for local volume change.
In TTM, `pi` is often partially known (e.g., via samples, or `pi(x) propto \tilde{pi}(x)`), `eta` is known. The goal is to find `S`.

### B. Triangular Map Structure
The map `S: R^K -> R^K` is defined as:
`S(x) = [ S_1(x_1) ; S_2(x_1, x_2) ; ... ; S_K(x_1, ..., x_K) ]^T = z` (3)
Each `S_k(x_1, ..., x_k)` maps to a scalar `z_k`.

**Monotonicity Requirement:** Each `S_k` must be strictly monotone in its last argument `x_k`:
`dS_k(x_1, ..., x_k) / dx_k > 0` for all feasible `x_1, ..., x_k`.

**Properties:**
1.  **Efficient Jacobian Determinant:** `grad_x S(x)` is lower-triangular.
    `det(grad_x S(x)) = product_{k=1 to K} (dS_k(x_1, ..., x_k) / dx_k)` (4)
2.  **Invertibility:** `S` is invertible component-wise.
3.  **Distribution Factorization:** Corresponds to:
    `pi(x) = pi(x_1) * pi(x_2|x_1) * ... * pi(x_K|x_1,...,x_{K-1})`

### C. Map Inversion (`S^-1`)
Computed sequentially: `x = S^-1(z)`
*   `x_1 = S_1^-1(z_1)`
*   `x_2 = S_2^-1(z_2; x_1)` (solve `S_2(x_1, x_2_unknown) = z_2` for `x_2_unknown`)
*   ...
*   `x_K = S_K^-1(z_K; x_1, ..., x_{K-1})` (5)
Each step `S_k^-1(...)` is a 1D root-finding problem.

**Sampling from `pi` (approximate):**
1.  Draw `Z_sample ~ eta`.
2.  Compute `X_sample = S^-1(Z_sample)` using (5). `X_sample` is a sample from `S_#eta approx pi`.
   `X_k_sample = S_k^-1(Z_k_sample; X_1_sample, ..., X_{k-1}_sample)` samples from `pi(x_k | x_1=X_1_sample, ..., x_{k-1}=X_{k-1}_sample)` if `eta` has independent marginals.

### D. Sampling Conditionals of `pi`
To sample from `pi(x_{k+1:K} | x_{1:k} = x_{1:k}^*)`:
1.  Fix `x_j = x_j^*` for `j=1,...,k`.
2.  Draw `z_{k+1}, ..., z_K` from marginal reference distributions `eta_{k+1}, ..., eta_K`.
3.  Compute `x_{k+1:K}^*` using the lower block of `S^-1`:
    `x_{j}^* = S_j^-1(z_j; x_1^*, ..., x_k^*, x_{k+1}^*, ..., x_{j-1}^*)` for `j = k+1, ..., K`. (6)
    The resulting `x_{k+1:K}^*` is a sample from `pi(x_{k+1:K} | x_{1:k}^*)`.
    This is central to Bayesian inference: if `pi(x) = p(a,b)` and `x_{1:k}^* = b^*` (observed data), then `x_{k+1:K}^*` samples posterior `p(a|b*)`.
The conditional distribution factorizes as:
`pi(x_{k+1:K}|x_{1:k}^*) = pi(x_{k+1}|x_{1:k}^*) * ... * pi(x_K|x_{1:k}^*, x_{k+1:K-1}^*)` (7)

### E. Conditional Independence (Sparsification)
If, e.g., `x_3 _|_ x_1 | x_2`, then `pi(x_3|x_1,x_2) = pi(x_3|x_2)`.
The map `S_3(x_1,x_2,x_3)` simplifies to `S_3(x_2,x_3)`.
Example of full map vs. sparse map for `pi(x_1)pi(x_2|x_1)pi(x_3|x_2)pi(x_4|x_3)`:
Full `S(x_{1:4}) = [ S_1(x_1) ; S_2(x_1,x_2) ; S_3(x_1,x_2,x_3) ; S_4(x_1,x_2,x_3,x_4) ]`
Sparse `S(x_{1:4}) = [ S_1(x_1) ; S_2(x_1,x_2) ; S_3(x_2,x_3) ; S_4(x_3,x_4) ]` (8)
**Benefits:** Robustness, computational efficiency, scalability.

## III. Implementation Details

### A. General Considerations for Map `S`
*   Nonlinearity of `S` relates to "distance" between `pi` and `eta`.
*   Assume `pi`, `eta` have similar tails => `S` approx. linear far from origin.
*   **Standardization:** Pre-process data `X` for `pi` to have zero mean, unit variance if `eta = N(0,I)`. Apply `S` in standardized space, then transform back.

### B. Structuring Map Components `S_k`
Each `S_k(x_1, ..., x_k)` must be monotone in `x_k`.

**1. Monotonicity Enforcement:**

*   **Through Integration (Flexible, Costly):**
    1.  Start with `S_k^non(x_1,...,x_k) = g(x_1,...,x_{k-1}) + g_hat(x_1,...,x_k)`.
    2.  `S_k(x_1,...,x_k) = g(x_1,...,x_{k-1}) + f(x_1,...,x_k)`
        where `f(x_1,...,x_k) = integral_0^{x_k} r(g_hat(x_1,...,x_{k-1},t)) dt`. (9)
        `r: R -> R+` is a rectifier (e.g., `exp(u)`, `softplus(u)`).
    3.  Allows cross-terms (e.g., `x_1*x_k`) in `g_hat`.
    4.  Often requires numerical integration.
    5.  `g_hat` should revert to a positive constant in `x_k` tails for `f` to extrapolate linearly.

*   **Through Variable Separation (Efficient, Less Expressive):**
    1.  `S_k(x_1,...,x_k) = g(x_1,...,x_{k-1}) + f(x_k)`. (10)
    2.  `f(x_k)` is a sum of inherently monotone basis functions of `x_k` (e.g., `c_1*x_k + c_2*erf(x_k)` with `c_1,c_2 > 0`).
    3.  No numerical integration. Efficient if linear in coefficients.
    4.  Disallows cross-terms involving `x_k` and `x_1,...,x_{k-1}` in `f`.

**2. Parameterization Complexity (How `S_k` relates `x_k` to `x_1,...,x_{k-1}`):**
Affects the root-finding objective for `x_k = S_k^-1(z_k; x_1,...,x_{k-1})`, which is:
`z_k - g(x_1,...,x_{k-1}) = f(x_1,...,x_{k-1}, x_k)` (general form) (11)

*   **Marginal Maps:** `S_k(x_k) = const + f(x_k)`.
    Objective: `z_k - const = f(x_k)`.
    `f` shape and offset (relative to `z_k`) fixed. Assumes `x_k` independent of previous `x_i`.

*   **Separable Maps:** `S_k(x_1,...,x_k) = g(x_1,...,x_{k-1}) + f(x_k)`.
    Objective: `z_k - g(x_1,...,x_{k-1}) = f(x_k)`.
    `f` shape fixed; offset (from `g`) varies with `x_1,...,x_{k-1}`. Can model shifts in conditional `pi(x_k|...)`.

*   **Cross-Term Maps:** `S_k(x_1,...,x_k) = g(x_1,...,x_{k-1}) + f(x_1,...,x_{k-1},x_k)`.
    Objective: `z_k - g(x_1,...,x_{k-1}) = f(x_1,...,x_{k-1},x_k)`.
    Both `f` shape and offset `g` can vary with `x_1,...,x_{k-1}`. Most expressive. Typically uses integration for monotonicity.

**3. Basis Functions (for `g`, `g_hat`, `f`):**

*   **Polynomial Basis Functions:**
    *   Orthogonal polynomials (e.g., Hermite `He_j(x)` for Gaussian `eta`).
    *   **Edge Control:**
        *   **Hermite functions:** `H_j(x) = He_j(x) * exp(-x^2/4)`. Revert to zero in tails.
            `S_k` with only `H_j` also reverts to zero; often use unweighted linear terms + weighted higher-order terms.
        *   **Finite Support Weights:** `H_j^EC(x) = He_j(x) * SplineWeight(x, r)` (SplineWeight is zero outside `[-r,r]`).
            `S_k` becomes linear outside `(-r,r)^K` if using integrated-rectified form with these.

*   **Radial Basis Functions (RBFs):**
    *   Local influence, good for multimodality. Params `mu` (center), `sigma` (scale).
    *   Standard RBF: `(1/sqrt(tau)) * exp(-((x_k-mu)/sigma_scaled)^2)`.
    *   Monotone RBF-derived functions for `f(x_k)` or `g_hat`:
        *   Integrated RBF (iRBF): `(1/2) * (1 + erf((x_k-mu)/sigma_scaled))` (sigmoid).
        *   Left/Right Edge Terms (LET/RET): Monotone, become linear in tails.

### C. Optimization of Map Parameters (Coefficients `c`)
(This section is briefly introduced in the source text provided)

**General Approach:** Minimize an objective function measuring discrepancy between the (transformed) distributions.
1.  **Maps from Samples (`{X^i ~ pi}` available):**
    *   Objective: Make pushforward samples `Z^i = S(X^i)` resemble `eta`.
    *   E.g., Maximize likelihood of `Z^i` under `eta`, or minimize `KL(S_#pi || eta)`.

2.  **Maps from Densities (`pi(x) propto \tilde{pi}(x)`, `eta(z)` known):**
    *   Objective: Make pullback density `eta(S(x))|det(grad S(x))|` match `\tilde{pi}(x)`.
    *   E.g., Minimize `KL(pi || S_#eta)`. The map optimized here is often conceptualized as `R = S^-1` (reference-to-target).

(Specific optimization algorithms are not detailed in the provided excerpt but would typically involve gradient-based methods.)

## IV. Algorithmic Flow Summary (Conceptual)

1.  **Define Distributions:**
    *   Characterize target `pi` (samples/unnormalized density).
    *   Choose reference `eta` (e.g., `N(0,I)`).
    *   Standardize `pi` data if necessary.

2.  **Build Map `S`:**
    *   Choose variable ordering `x_1, ..., x_K`.
    *   For each `S_k`:
        *   Select parameterization (marginal, separable, cross-term).
        *   Select monotonicity enforcement (integration, separation).
        *   Select basis functions.
        *   Incorporate conditional independencies (sparsity).

3.  **Optimize Map:**
    *   Define objective function (e.g., KL divergence, likelihood).
    *   Numerically optimize map parameters `c`.
    *   *(Optional Loop):* Adapt map structure (`S_k` complexity, basis) and re-optimize if needed.

4.  **Apply Transport:**
    *   **Generative Modeling/Sampling:** `x_sample = S^-1(z_sample)` where `z_sample ~ eta`.
    *   **Density Estimation:** `pi(x_query) approx eta(S(x_query)) * |det(grad_x S(x_query))|`.
    *   **Conditional Sampling (Bayesian Inference):** Use Eq (6).

## V. Key Advantages of TTM

*   **Parsimony:** User-controlled map complexity.
*   **Sparsity:** Exploits conditional independence for efficiency and robustness in high dimensions.
*   **Numerical Convenience:** Efficient Jacobian determinant, sequential 1D inversions.
*   **Explainability:** Map components relate to statistical features (marginals/conditionals).
*   **Exact Density Evaluation:** Possible via change-of-variables.
*   **Direct Conditional Sampling:** Powerful for Bayesian inference.
