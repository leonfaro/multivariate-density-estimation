
---

## Notation (consistent with the paper, color‑free)

* Vectors/maps bold; scalars/components roman: $\mathbf{x}$, $\mathbf{z}$, $\mathbf{S}$ vs. $x_k$, $S_k$.
* $\pi$ is the target distribution; $\eta=\mathcal{N}(\mathbf{0},I_K)$ is the reference.
* $\mathbf{x}\sim\pi$, $\mathbf{z}\sim\eta$; $\mathbf{z}=\mathbf{S}(\mathbf{x})$ is the target‑to‑reference transform.
* $\mathbf{R}=\mathbf{S}^{-1}$ maps reference to target (not used for training here).
* $\partial_k S_k=\frac{\partial}{\partial x_k}S_k$ is the monotone derivative in the last argument;
  $\ell(\mathbf{x})=\sum_{k=1}^K \log \partial_k S_k(\mathbf{x})$ is the log‑Jacobian.
* Train‑only standardization: $\tilde x_k=(x_k-\mu_k)/\sigma_k$. The Jacobian includes $-\log\sigma_k$ for every $k$.
* Ranks use `ties.average`. Normal‑scores use $u\in[\frac1{n+1},\frac{n}{n+1}]$ and $z^\star=\Phi^{-1}(u)$.
* Seeds are deterministic and shared across models; Half‑Moon uses seed, seed+1, seed+2 for train, test, val‑split. The main Config‑4D pipeline runs Train/Test only.

---

## Objective (“maps‑from‑samples”, paper Eq. 36–39)

We learn a lower‑triangular, strictly monotone $\mathbf{S}$ such that $\mathbf{z}=\mathbf{S}(\mathbf{x})\sim\eta$.
The per‑sample **log‑density** under the learned model is

$$
\log \hat p(\mathbf{x})
=
\sum_{k=1}^K\big(-\tfrac12 z_k(\mathbf{x})^2 - \tfrac12\log(2\pi) + \log \partial_k S_k(\mathbf{x})\big),
\quad \mathbf{z}=\mathbf{S}(\mathbf{x}).
$$

Training minimizes the forward‑KL surrogate (up to a constant)

$$
\mathcal{J}(\mathbf{S})
=
\mathbb{E}_{\mathbf{x}\sim\pi}\Big[\tfrac12\|\mathbf{S}(\mathbf{x})\|^2 - \sum_{k=1}^K \log \partial_k S_k(\mathbf{x})\Big],
$$

with Monte‑Carlo approximation over the train set

$$
\widehat{\mathcal{J}}
=
\frac1{n_{\mathrm{tr}}}\sum_{i=1}^{n_{\mathrm{tr}}}
\Big(\tfrac12\|\mathbf{S}(\mathbf{X}^i)\|^2 - \sum_{k=1}^K \log \partial_k S_k(\mathbf{X}^i)\Big).
$$

For the separable and cross‑term parameterizations this decomposes per dimension (paper Eq. 39), enabling $K$ independent 1‑D optimizations conditional on $x_{1:k-1}$.

---

## TTM parameterizations (paper Eq. 20 / Eq. 21 / Eq. 22)

### Marginal (Eq. 20)

**Form**
$S_k(\mathbf{x}) = a_k + b_k\,\tilde x_k,\quad \tilde x_k=(x_k-\mu_k)/\sigma_k,\quad b_k>0.$

**Training (closed‑form normal‑scores)**

* Compute train‑only $\mu_k,\sigma_k$.
* Map ranks of $x_k$ to $u$, then $z^\star=\Phi^{-1}(u)$.
* Fit $z^\star \approx a_k + b_k\,\tilde x_k$ in least‑squares; clamp $b_k\ge0$.
* Store $\{\mu_k,\sigma_k,a_k,\log b_k\}$.

**Prediction and log‑density**

* $z_k=a_k+b_k\tilde x_k$.
* $\log\partial_k S_k = \log b_k - \log\sigma_k$ (constant in $x$).
* $L_k = -\tfrac12 z_k^2 - \tfrac12\log(2\pi) + \log b_k - \log\sigma_k.$

**Notes**
Very fast baseline; cannot model conditionals $x_{1:k-1}$.

---

### Separable (Eq. 21)

**Form**
$S_k(\mathbf{x}) = g_k(\mathbf{x}_{1:k-1}) + f_k(x_k)$, with $f_k'(x_k)>0$.

**Bases**

* $g_k$: polynomial features in $x_{1:k-1}$ (degrees 1–3).
* $f_k$: monotone 1‑D basis in $x_k$; derivative basis $B=\partial_{x_k}\text{basis}_f$.

**Orthogonalization (Appendix‑style)**

* Let $P_{\mathrm{non}}$ be the design for $g_k$ and $P_{\mathrm{mon}}$ for $f_k$.
* Define $M=(P_{\mathrm{non}}^\top P_{\mathrm{non}}+\lambda I)^{-1}P_{\mathrm{non}}^\top$.
* Set $A=(I-P_{\mathrm{non}}M)P_{\mathrm{mon}}$ and $D=MP_{\mathrm{mon}}$.

**Per‑k convex objective (paper Eq. 39)**

$$
J_k(c)=\tfrac12\|A c\|^2 - \sum_{i}\log(Bc)_i + \tfrac{\lambda}{2}(\|D c\|^2+\|c\|^2),
\qquad (Bc)_i>\varepsilon.
$$

**Gradient (analytic)**

$$
\nabla J_k(c) = A^\top(Ac) - \sum_i \frac{B^\top e_i}{(B c)_i} + \lambda(D^\top D c + c).
$$

**Prediction and log‑density**

* $z_k=g_k+f_k$.
* $\log\partial_k S_k=\log(Bc) - \log\sigma_k$.
* Per‑dimensional LD as in the universal formula above.

**Notes**
Adds conditional shifts via $g_k$ and preserves strict monotonicity via $f_k' > 0$. No cross‑interactions $x_k\cdot x_j$ inside the derivative beyond the monotone channel.

---

### Cross‑term (Eq. 22)

**Form**
$S_k(\mathbf{x}) = g_k(\mathbf{x}_{1:k-1}) + \int_0^{x_k} \exp(h_k(t,\mathbf{x}_{1:k-1}))\,dt.$

**Monotonicity**
$\partial_k S_k = \exp(h_k) > 0$ by construction.

**Bases and quadrature**

* $g_k$ as in separable.
* $h_k(t,\mathbf{x}_{1:k-1}) = \sum_\alpha \beta_\alpha \,\psi_\alpha(t,\mathbf{x}_{1:k-1})$ with cross terms $t^r x_j^s$.
* Training cache: build B‑spline basis $B(t)$ (df columns) and predecessor basis $\Psi(\mathbf{x}_{1:k-1})$ (m columns), then form the row‑wise Khatri–Rao design $\Pi \in \mathbb{R}^{(NQ)\times(df\cdot m)}$ with $\Pi_{(i,q),:}=B(t_{iq})\odot\Psi(\mathbf{x}_{1:k-1}^{(i)})$ (no recycling).
* Integrate via Gauss–Legendre on $[0,x_k]$; accumulate $\log\int\exp(h)$ by **log‑sum‑exp** of $\log w_q + h_q$.

**Per‑k objective (forward‑KL component)**
Mean of $0.5\,S_k^2 - h_k(x_k,\mathbf{x}_{1:k-1})$ plus regularizers.

**Prediction and log‑density**

* $z_k=g_k + \int_0^{x_k}\exp(h_k)\,dt$.
* $\log\partial_k S_k = h_k(x_k,\mathbf{x}_{1:k-1}) - \log\sigma_k$.
* Per‑dimensional LD as in the universal formula.

**Notes**
Most expressive; computationally heavier due to the 1‑D integral and stabilization of $\exp(h_k)$.

---

## TRTF as conditional $S_k$ in spirit

* TRTF fits the conditionals $p(x_k\mid x_{1:k-1})$ with transformation forests.
* The learned conditional transformation gives a monotone map in $x_k$, effectively $z_k=\Phi^{-1}(F_k(x_k\mid x_{1:k-1}))$.
* This mirrors TTM’s idea: a triangular, monotone transform from $x$ to a Gaussianized $z$.
* Evaluation is directly comparable: sum conditional log‑densities across $k$ to get joint LD for each sample.

---

## Data‑generating processes (for evaluation)

### Config‑4D (triangular, heterogeneous marginals)

* Sequential sampling:
  $X_1\sim \mathrm{Norm}(0,1)$;
  $X_2\mid X_1\sim \mathrm{Exp}(\text{rate}=\mathrm{softplus}(X_1))$;
  $X_3\mid X_{1:2}\sim \mathrm{Beta}(\alpha=\mathrm{softplus}(X_2),\ \beta=\mathrm{softplus}(X_1))$;
  $X_4\mid X_{1:3}\sim \mathrm{Gamma}(\text{shape}=\mathrm{softplus}(X_3),\ \text{scale}=\mathrm{softplus}(X_2))$.
* Optional column permutation with corresponding config re‑indexing.
* “True (Joint)” evaluates the oracle conditional log‑densities; “True (marginal)” fits independent MLEs.

### Half‑Moon‑2D (toy, controlled by $N$ and $\sigma$)

* Curves before noise: $(\cos t,\ \sin t)$ and $(1-\cos t,\ -\sin t + 0.5)$ with $t\sim\mathcal{U}[0,\pi]$.
* Add i.i.d. $N(0,\sigma^2 I_2)$; shuffle.
* Train/Test generated with seeds `seed` and `seed+1`; validation split from train with `seed+2` and `val_frac=0.2`.
* Splits saved as `results/splits_halfmoon2d_seedXXX.rds`.

---

## Predict‑API contract (shared across models)

* `predict(M, X, "logdensity_by_dim")` returns $L\in\mathbb{R}^{N\times K}$ with $L_{ik}=\log p(x_k^{(i)}\mid x_{1:k-1}^{(i)})$.
  For TTM this uses the universal per‑dim formula with $z_k$ and $\log\partial_k S_k$.
* `predict(M, X, "logdensity")` returns $L^{\mathrm{joint}}\in\mathbb{R}^N$ with $L^{\mathrm{joint}}_i=\sum_{k=1}^K L_{ik}$.
* **Invariant:** `rowSums(logdensity_by_dim) == logdensity` up to $10^{-10}$.
* **No NA/Inf:** out‑of‑support values are mapped to large negative logs; never return NA/Inf.

---

## Evaluation tables (nats)

* For each model, compute `by_dim` on $X_{\mathrm{te}}$, sum to `joint`.
* Per‑dim mean NLL is $-$colMeans(`by_dim`); SE per dim is `stderr(-by_dim[,k])`.
* SUM row reports sums of means; SE computed as `sd(rowSums(-by_dim))/sqrt(N)`.
* Columns: **True (marginal)**, **True (Joint)**, **Random Forest**, **Marginal Map**, **Separable Map**, **Cross‑term Map`.
  Note: `train_test_policy` (previously set to `train_test_only`) is omitted from the main pipeline output (`main.R`),
  but remains in auxiliary CSVs produced by Half‑Moon evaluation for bookkeeping.
* True (Joint) under permutation: evaluation uses the canonical generative order to compute per‑dimensional oracle terms and then re‑indexes them to the current display order (no NA); joint equals the sum by construction.
* Formatting: `"%.2f ± %.2f"` where the ± part is **2·SE**.

---

## Half‑Moon panels

* Grid: Achsen per Trainingssatz $\mu \pm 3\,\sigma$ je Dimension erzeugen (seed‑spezifische Splits), dann ein reguläres Gitter aufspannen.
* Auswertung: **joint** Log‑Dichte aller Modelle auf demselben Gitter berechnen (Chunking + parallel, mit Cache).
* Konturen: **Iso‑Massen‑Level** mit Massen $\{0.50,0.70,0.90\}$ pro Modell; gemeinsame Achsen in allen Panels.
* Darstellung: Datenpunkte zuletzt zeichnen (sichtbar über Konturen).
* Ausgabe: `results/moon/<timestamp>/panels_seedXXX.png` und begleitende `meta.rds` (Gitter‑ und Cache‑Metadaten).

---

## Numerical stability rules (hard constraints)

* Work in **log‑space** everywhere; never exponentiate densities during training or evaluation.
* Use **log‑sum‑exp** for all $\log\sum_i e^{u_i}$ aggregates (e.g. cross‑term integrals).
* Include $-\tfrac12\log(2\pi)$ per dimension in every TTM LD; include $-\log\sigma_k$ from train‑only standardization in each dimension’s log‑Jacobian.
* For separable, maintain $(B c)_i>\varepsilon$ (e.g. $10^{-6}$) to avoid $\log 0$.
* For marginal, ensure $\sigma_k>0$ with an $\varepsilon$ floor and clamp $b_k\ge0$.
* Support clamping for Beta/Gamma/Exp in log‑space; return large negative logs instead of NA/Inf.
* Cross‑term integrals: accumulate $\log\int \exp(h)$ via LSE of `log w + h`; never sum raw exponentials.
* Double precision throughout; avoid mixing float types.

---

## Gradients, Hessians, and fast updates

* Separable per‑k objective is strictly convex over $(B c)>0$; gradient is analytic:
  $\nabla J_k(c) = A^\top(Ac) - \sum_i \frac{B^\top e_i}{(B c)_i} + \lambda(D^\top D c + c).$
* Hessian is $A^\top A + \sum_i \frac{B^\top e_i e_i^\top B}{(B c)_i^2} + \lambda(D^\top D + I)$ and positive definite on the feasible set.
* Cross‑term gradients backprop through LSE: if $s=\mathrm{LSE}(\{a_q\})$ then $\partial a_q = \mathrm{softmax}(a)_q$; hence
  $\partial_\beta \log\int e^{h} \approx \sum_q \omega_q\,\psi_q$ with $\omega_q \propto e^{\log w_q + h_q}$.
* For “real‑time” calibration, warm‑start from previous parameters and take a few L‑BFGS steps; use diagonal preconditioning from $\mathrm{diag}(A^\top A)+\lambda$.

---

## Complexity (big‑O, high level)

* TTM‑marginal: $O(NK)$ time, $O(K)$ parameters.
* TTM‑separable: $O(N p_k)$ per dimension to build designs, plus quasi‑Newton iterations; memory $O(p_k^2)$.
* TTM‑cross: $O(N Q)$ per dimension, with $Q$ quadrature nodes; LSE dominates.
* TRTF: roughly $O(\text{ntree}\,N\log N)$ training; prediction $O(\text{ntree}\cdot\text{depth})$.

---

## Repository mapping (spec → code)

* Eq. (20) → `trainMarginalMap`, `predict.ttm_marginal`: normal‑scores, constant log‑Jacobian.
* Eq. (21) → `trainSeparableMap`, `predict.ttm_separable`: per‑k convex objective, $A,B,D$.
* Eq. (22) → `trainCrossTermMap`, `predict.ttm_cross_term`: $g_k+\int\exp(h_k)$, Gauss–Legendre + LSE.
* TRTF → `fit_TRTF`, `predict.mytrtf`: conditional transformation forests per $k$.
* TRUE (marginal) → `fit_TRUE`, `predict_TRUE`; TRUE (joint) → `true_joint_logdensity_by_dim`.
* Evaluation & formatting → `04_evaluation.R` and helpers.
* Half‑Moon splits and panels → `scripts/halfmoon_data.R`, `scripts/halfmoon_plot.R`.
* Top‑level orchestration → `main.R` and `main_moon.R`.

---

## Determinism and splits

* One global `SEED` governs all randomness; TRTF receives the same seed.
* Config‑4D shuffles once for split; optional permutation applied consistently to data and config.
* Half‑Moon: train uses `seed`, test uses `seed+1`, val‑split uses `seed+2`.
* Persist splits as `results/splits_halfmoon2d_seed%03d.rds` for exact reproducibility.
* Re‑running with the same seed reproduces tables and PNGs (modulo floating‑point roundoff).

---

## Acceptance tests (deterministic)

* `predict(..., "logdensity_by_dim")` is $N\times K$ for every model; `predict(..., "logdensity")` is length $N$.
* `rowSums(by_dim) == joint` within $10^{-10}$.
* No NA/Inf anywhere in predictions.
* TTM constants present: $-\tfrac12\log(2\pi)$ per dimension; $-\log\sigma_k$ per dimension.
* Half‑Moon: `nrow(X_tr)+nrow(X_val)==N_TRAIN` and `nrow(X_te)==N_TEST`.
* Timing table `total_sec = train_sec + test_sec` for each model.
* Permutation sanity (Config‑4D): if marginals are independent, SUM‑NLL invariant under column permutation.
* Repro test: same seed → identical CSV/PNG; different seed → different but stable metrics.

---

## Debugging cribsheet

* SUM row inflated across TTM variants → missing $-\tfrac12\log(2\pi)$.
* TTM off by a constant → missing $-\log\sigma_k$ in the log‑Jacobian.
* NA/Inf in separable → some $(B c)_i\le 0$; increase ε, add barrier, regularize.
* Exploding cross‑term gradients → not using LSE for the integral; reduce nodes; add β‑penalty.
* TRTF variance across runs → seed not propagated; fix once at entry.
* Half‑Moon panels incomparable → levels per model; switch to pooled global levels.
* RowSums mismatch → two different code paths for by‑dim vs joint; unify through the same aggregator.

---

## Pseudocode blueprints (non‑executable, repo‑aligned)

### TTM‑marginal (Eq. 20)

```
TRAIN:
  input: S with X_tr, X_val, X_te; seed
  compute μ,σ on X_tr
  for k=1..K:
    u <- rank_avg(X_tr[,k])/(n_tr+1); clamp endpoints
    z* <- qnorm(u)
    fit z* = a_k + b_k * ((X_tr[,k]-μ_k)/σ_k) with b_k >= 0
    store μ_k, σ_k, a_k, log b_k
  return model with stored params and timing

PREDICT:
  input: model, X, type
  x_std <- (X-μ)/σ
  for k: z_k = a_k + b_k * x_std[,k]
         L_k = -0.5*z_k^2 - 0.5*log(2π) + log b_k - log σ_k
  return L (by-dim) or rowSums(L)
```

### TTM‑separable (Eq. 21)

```
TRAIN:
  standardize train-only -> μ,σ
  for k=1..K:
    P_non <- features(x_{1:k-1})  # degree 1–3
    P_mon <- monotone basis in x_k
    B     <- derivative basis of P_mon
    M <- solve(P_non^T P_non + λI) P_non^T
    A <- (I - P_non M) P_mon
    D <- M P_mon
    minimize J_k(c) = 0.5||A c||^2 - sum log(B c) + 0.5 λ (||D c||^2 + ||c||^2)
         subject to (B c) > ε
    store coeffs_k
  return model

PREDICT:
  X_std <- (X-μ)/σ
  for k: g_k <- P_non(X_std_{1:(k-1)})·coef_non
         f_k <- P_mon(X_std_k)·coef_mon
         z_k <- g_k + f_k
         L_k <- -0.5*z_k^2 - 0.5*log(2π) + log(B(X_std_k)·coef_mon) - log σ_k
  return L or rowSums(L)
```

### TTM‑cross (Eq. 22)

```
TRAIN:
  standardize train-only -> μ,σ
  choose basis ψ(t, x_prev) with cross terms; choose nodes t_q and weights w_q
  for k=1..K:
    param: h_k(t, x_prev) = ψ(t, x_prev)·β_k, and g_k(x_prev) with α_k
    objective ≈ mean[ 0.5 * (g_k + ∫ exp(h_k))^2 - h_k(x_k, x_prev) ] + regs
    integral: log_int = LSE_q( log w_q(x_k) + h_k(t_q(x_k), x_prev) )
    optimize α_k, β_k with L-BFGS-B in log-space; warm-start from separable if available
    store α_k, β_k
  return model

PREDICT:
  X_std <- (X-μ)/σ
  for k: g_k <- g(X_prev; α_k)
         h*  <- h(x_k, X_prev; β_k)
         log_int <- LSE_q( log w_q(x_k) + h_k(t_q(x_k), X_prev) )
         z_k <- g_k + exp(log_int)   # computed via log-space
         L_k <- -0.5*z_k^2 - 0.5*log(2π) + h* - log σ_k
  return L or rowSums(L)
```

Implementation policy (TTM‑Cross, this repo):
- Train/Test only: no validation split is used (no CV/grid or automatic order selection).
- Hyperparameters via config: ridge weights from `getOption("mde.ctm.lambda_non")` and `getOption("mde.ctm.lambda_mon")` (fallback env `MDE_CTM_LAMBDA_NON`/`MDE_CTM_LAMBDA_MON`). Defaults: `lambda_non=3e-2`, `lambda_mon=3e-2`.
- Quadrature via config: nodes `Q` from `getOption("mde.ctm.Q")` or env `MDE_CTM_Q`; else adaptive default `min(12, 4+2·max(degree_t, degree_t_cross))`.
- Variable order: fixed by `getOption("mde.ctm.order")`/env `MDE_CTM_ORDER` ∈ {`as-is`,`x1_x2`,`x2_x1`} (2D only).
- Logging: main.R prints `[RIDGE] lambda_non=..., lambda_mon=...` after fitting TTM‑Cross.

### TRTF

```
TRAIN:
  fit marginal transformation model for k=1
  for k=2..K: fit transformation forest for x_k | x_{1:k-1}
PREDICT:
  by-dim L: log p(x_1), log p(x_2|x_1), ..., log p(x_K|x_{1:K-1})
  joint L: rowSums(by-dim)
```

---

## Expected outcomes (qualitative)

* Config‑4D: True (Joint) best; TRTF and TTM‑cross close; TTM‑separable next; TTM‑marginal worst on strong conditional structure.
* Half‑Moon: TTM‑cross best captures the arcs; TRTF shows piecewise shapes; TTM‑separable yields smooth tilted contours; marginal shows near‑circular levels.
* SUM‑row SE is computed from rowwise sums; differences beyond ±2·SE are practically meaningful.

---

## Hyperparameter ranges (practical)

* Separable: degree $g$ ∈ {1,2,3}; $\lambda\in\{10^{-6},10^{-5},10^{-4},10^{-3},10^{-2}\}$; $\varepsilon=10^{-6}$.
* Cross‑term: quadrature nodes $Q\in\{8,16,32\}$; small L2 on β.
* TRTF: `ntree ≈ n_tr`, `minsplit∈[20,100]`, `maxdepth∈[2,5]`, `minbucket∈[5,30]`.
* Half‑Moon panels: `grid_side = clamp(⌊√(100·N_tr)⌋,80,200)` unless overridden;
  contour levels from iso‑mass thresholds {0.5,0.7,0.9};
  grid evaluations cached under `results/cache/moon_<digest>.rds`.

---

## Short thesis paragraph (ready to paste)

We train lower‑triangular, strictly monotone target‑to‑reference maps $\mathbf{S}$ in the forward‑KL formulation (paper Eq. 36–39). We instantiate three parameterizations—marginal (Eq. 20), separable (Eq. 21), and cross‑term (Eq. 22)—and evaluate them against conditional transformation forests (TRTF). All computations, both training and test, run in **log‑space** with explicit Gaussian constants and train‑only standardization offsets. The per‑sample joint log‑density equals the sum over dimensions of $-\tfrac12 z_k^2 - \tfrac12\log(2\pi) + \log\partial_k S_k$. We report mean NLL (nats) with ±2·SE per dimension and in a SUM row, and we visualize Half‑Moon density shapes on a shared grid with iso-mass contour levels.

---

## Copy‑ready identities for code comments

```
# Universal per-dimension LD for TTM:
#   LD_k(x) = -0.5 * z_k(x)^2 - 0.5*log(2*pi) + log(∂_k S_k(x))
# Marginal (Eq. 20): log(∂_k S_k) = log b_k - log σ_k   # constant in x
# Separable (Eq. 21): log(∂_k S_k) = log( B(x_k) · c ) - log σ_k
# Cross-term (Eq. 22): log(∂_k S_k) = h_k(x_k, x_{1:k-1}) - log σ_k
# Joint LD: LD(x) = Σ_k LD_k(x)
# NLL(X):   -Σ_i LD(x^i)              # reported in nats
# SE (SUM row): sd(rowSums(-LD_by_dim)) / sqrt(N)
# Invariant: rowSums(LD_by_dim) ≡ LD_joint (|Δ| ≤ 1e-10)
# POLICY: training AND testing strictly in log-space; use log-sum-exp where needed.
```

---

## Final checklist

* All models produce `by_dim` of shape $N\times K$ and `joint` of length $N$.
* `rowSums(by_dim) == joint` (tolerance $10^{-10}$).
* TTM includes $-\tfrac12\log(2\pi)$ per dimension and $-\log\sigma_k$ per dimension.
* No NA/Inf anywhere in predictions.
* Seeds consistent across models; Half‑Moon uses seed, seed+1, seed+2.
* CSV and PNG written with seed in filenames.
* Panels use global levels and identical axes; points visible above contours.
* Report mean NLL in nats with ±2·SE and a timing table.
* Narrative explains why marginal < separable < cross on conditional DGPs and how TRTF compares.
* If it’s not in **log‑space**, it’s a bug.
