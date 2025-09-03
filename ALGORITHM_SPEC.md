# Algorithm Specification — Multivariate Conditional Density via Triangular Transport Maps (TTM)

This document specifies the algorithms, data flow, and evaluation protocol implemented in this repository. It is kept in lockstep with the codebase to ensure scientific reproducibility and traceability.

Last updated: 2025-09-03

## Scope and Workflow
- Problem: Estimate multivariate conditional densities p(x1,…,xK) with rigorous, reproducible experiments in R.
- Workflow: data generation → fixed train/test split → model training (TRUE, TRTF, TTM variants, optional copula) → per‑dimension and joint log‑density evaluation → CSV artifacts.
- Reference measure: standard Gaussian N(0,1). Triangular maps S push X to Z=S(X) with independent N(0,1) marginals; Jacobian diagonal enters the log‑density.

## Data Generation and Splitting
- Generator: `Generate_iid_from_config(N, cfg)` (01_data_generation.R)
  - `cfg` is a list of K entries with fields:
    - `distr`: one of `norm`, `exp`, `beta`, `gamma`.
    - `parm(prev)`: optional function mapping previously generated columns to parameter list for the current distribution.
  - Positivity is enforced for all strictly positive parameters by clamping to >= 1e-6 before sampling.
  - Returns matrix `X ∈ R^{N×K}` with columns `X1,…,XK`.
- Split: `split_data(X, seed)` (02_split.R)
  - Deterministic 80/20 train/test split using `set.seed(seed)` and a single permutation.
  - Train statistics are the only source for standardization parameters used by all models.

## Standardization
- For any fitted TTM model S:
  - `mu_k = mean(X_tr[,k])`, `sigma_k = sd(X_tr[,k]) + .Machine$double.eps`.
  - Standardized inputs: `x_std = (x - mu) / sigma` applied column‑wise.

## Models and Training

### TRUE Baselines
- TRUE (marginal): per‑dimension MLE under the univariate families in `cfg`.
  - Objective (per dim): maximize log‑likelihood; implemented via `optim` with L‑BFGS‑B and lower bounds 1e‑6 for positive params.
  - Evaluation: `logL_TRUE`, `logL_TRUE_dim` (models/true_model.R).
- TRUE (joint/oracle): evaluates the product of true conditionals induced by `cfg` (no fitting), via `true_joint_logdensity_by_dim` (models/true_joint_model.R). Handles domain clamps for exp/gamma/beta.

### TRTF (Transformation Forest) Baseline
- Implementation: `mytrtf` using `traforest` with `BoxCox` base model and `partykit::ctree_control` (models/trtf_model.R).
- Defaults: `minsplit=40`, `minbucket=5`, `maxdepth=2`, `ntree = nrow(X_tr)`, `seed` propagated.
- Predict API: `predict.mytrtf(newdata, type = 'logdensity'|'logdensity_by_dim')` returns joint or per‑dimension log‑densities. Errors if required packages (`tram`, `trtf`, `partykit`) are not installed.

### Triangular Transport Maps (TTM)
All TTM variants share: standardized inputs, triangular structure (dimension k depends only on `x1,…,xk`), and positivity of the Jacobian diagonal in `x_k` to ensure monotonicity.

- Common forward core (models/ttm/ttm_core.R):
  - For input `X`, compute standardized `Xs`, then per‑k output `Z_k` and Jacobian diagonal `J_k = ∂S_k/∂x_k_std`.
  - Per‑dimension log‑density: `LD_k = -0.5 Z_k^2 - 0.5 log(2π) + log J_k - log σ_k`.
  - Joint log‑density is the row‑sum of `LD_k`.

1) TTM‑Marginal (models/ttm/ttm_marginal.R)
   - Form: `S_k(x) = a_k + b_k x_k_std`, with `b_k>0`.
   - Closed‑form: `b_k = 1/sd(x_k_std) (≈ 1)`, `a_k = -b_k * mean(x_k_std)` which equals the MLE solution of the maps‑from‑samples objective for 1D.
   - Fit: `fit_ttm(S, algo='marginal')` returns `S` with `coeffs[[k]] = c(a,b)`.

2) TTM‑Separable (models/ttm/ttm_separable.R)
   - Form: `S_k(x) = g_k(x_<k>) + f_k(x_k)` with `f_k' > 0`.
   - Bases: `g_k` is a polynomial feature map of degree `degree_g` on predecessors; `f_k` is a monotone 1D basis with positive coefficients; derivative computed via `d_build_f`.
   - Objective per k (maps‑from‑samples, Eq. 38/39): minimize `0.5 ||g_k + f_k||^2 - sum(log f_k') + ridge`.
   - Optimization: L‑BFGS‑B on `c_mon` with positivity lower bounds; `c_non` obtained in closed form by least‑squares projection. Ridge `lambda` applies to both parts.
   - Fit: `fit_ttm(S, algo='separable', degree_g=2, lambda=0)`.

3) TTM‑Cross‑Term (models/ttm/fit_ttm_crossterm.R)
   - Form: `S_k(x) = g_k(x_<k>) + ∫_0^{x_k} exp(h_k(t, x_<k>)) dt` with `J_k = exp(h_k(x_k, x_<k>)) > 0`.
   - Bases: `h_k` uses a tensor basis built by `build_h(t, x_<k>)` = B‑splines (degree 3) in `t` with df `df_t`, optionally with symmetric interior knots and constant tail features; multiplied by polynomial features in `x_<k>` via `build_g(·,deg_g)`.
   - Quadrature: Gauss–Legendre on [0,1], `Q` nodes. The integral is implemented as `sign(x_k) * |x_k| * sum_q w_q exp(h(t_q · x_k, x_<k>))`.
   - Objective per k: `0.5 ||S_k||^2 − sum(h_k(x_k, x_<k>)) + 0.5(λ_non||α||^2 + λ_mon||β||^2)` with analytic gradient (option). `α,β` are coefficients for `g_k` and `h_k`.
   - Regularization: default shares chosen as functions of `Q` and `N` if `lambda` is `NA` and no global overrides; otherwise use provided `lambda` or options (see below).
   - Clipping: raw `h` is clamped to `[-Hmax, Hmax]` before exponentiation to stabilize tails and ensure numeric safety.
   - Optional order search: greedy adjacent‑swap search starting from either the identity or a Cholesky‑pivot heuristic (`learn_ordering` if available). Logs to `artifacts/order_search.csv`.
   - Fit: `fit_ttm(S, algo='crossterm', deg_g, df_t, Q, lambda, Hmax, maxit, order_search)`.

#### Default/Global Options (TTM Cross‑Term)
Set via `options()`:
- `cross.maxit` (int, default 200): optimizer iterations per k.
- `cross.quad_nodes` (int, default Q passed or 32): GL nodes.
- `cross.df_t` (int, default 8 but effectively `min(cross.df_t, max(4, floor(N/10)))`).
- `cross.lambda_non`, `cross.lambda_mon` (numeric or NULL): ridge weights for `α` and `β`. If `NULL` and `lambda` is `NA`, a short internal search selects a total reg share near 0.12.
- `cross.use_analytic_grad` (logical, default TRUE): use analytic gradient.
- `cross.workers` (int): parallel hints (currently not used for per‑k but respected by some utilities).
- `cross.deg_g` (int, default 3): default predecessor degree if not explicitly provided.
- `cross.sparsity_tau` (numeric in [0,1], default 0.05): predecessor feature masking by absolute correlation threshold.
- `cross.strict_monotone` (logical, default FALSE): if TRUE, fail on any monotonicity‑of‑integral violation in diagnostics.

## Prediction and Evaluation API
- Unified TTM predict: `predict_ttm(model_or_fit, X, type = 'transform'|'jac_diag'|'logdensity_by_dim'|'logdensity')`.
  - Returns `Z`, `J`, per‑dimension log‑densities, or joint log‑density respectively.
  - S3 wrappers: `predict.ttm_marginal2`, `predict.ttm_separable`, `predict.ttm_cross_term`.
- TRTF predict: `predict.mytrtf(newdata, type)`; returns per‑dimension matrix for `type='logdensity_by_dim'` and a vector for `type='logdensity'`.
- Copula NP baseline (2D labeled) or KDE product fallback (models/copula_np.R), class `copula_np` with S3 `predict`.

## Evaluation Protocol
- Primary table builder: `calc_loglik_tables(models, config, X_te, config_canonical=NULL, perm=NULL)` (04_evaluation.R):
  - Computes mean per‑dimension Negative Log‑Likelihoods (NLL; nats) for each model and a sum row; standard errors via row‑wise variance of joint NLL.
  - TRTF is wrapped in `tryCatch`; if unavailable, its columns are `NA` and a warning is emitted.
  - TRUE (joint) supports permutation‑aware evaluation by mapping test data back to canonical order.
- Invariants enforced:
  - For every TTM S and input X: `sum_k LD_k == joint_logdensity` up to 1e‑12.
  - Jacobian diagonal strictly positive for TTM (checked in cross‑term fitter; failure raises an error).
  - Cross‑term oriented integral has consistent sign with `x_k` on held‑out sets; violations are logged to CSV.
  - Monotonicity‑of‑integral on a fixed u‑grid is validated; violations logged to `artifacts/cross_integral_checks.csv` and optionally fatal if `cross.strict_monotone=TRUE`.

## Diagnostics and Artifacts
- `artifacts/order_search.csv`: order search trace (perm, NLL, acceptance, timing).
- `artifacts/cross_clip_events.csv`: per‑k counts of raw h bounds, maximum absolute raw/bounded h, and maximum exp(h_bound).
- `artifacts/cross_basis_dims.csv`: dimensions and moments of t‑basis and tensor products used at fit time.
- `artifacts/cross_sign_checks.csv`, `artifacts/cross_integral_checks.csv`: orientation and monotonicity checks on held‑out data.
- `results/*.csv`: experiment summaries (e.g., half‑moon benchmarks) written by scripts under `scripts/`.

## Reproducibility
- Seeds: set `set.seed(seed)` at every experiment start; tests fix RNGKind for determinism.
- Unit tests: `tests/testthat/test_ttm_crossterm_invariants.R` validate invariants, reproducibility, and numeric safety for cross‑term TTM.
- Logging: scripts in `scripts/` create ASCII logs and CSV artifacts; intermediate metrics are printed via `message()`.

## Numerical Safeguards
- Positivity constraints and lower bounds (1e‑6) for distribution parameters and monotone basis coefficients.
- Domain clamps for exp/gamma/beta in TRUE baselines.
- Standardization `sigma` includes `+ .Machine$double.eps` to avoid `log(0)` in evaluation.
- Cross‑term: `h` clipping to `[-Hmax, Hmax]` before exponentiation; GL integration instead of differentiating splines.
- All predict APIs check matrix shapes and finiteness; any NA/Inf triggers an error.

## Notes on Theory Alignment
- The repository references `Theory.md` for the conceptual workflow; this file is not present in the repo. The implementation follows the standard triangular transport framework (as outlined in “A Friendly Introduction to Triangular Transport”) with maps‑from‑samples objectives and monotonicity in `x_k` ensured by construction of `J_k > 0`.
- Transformation forests (TRTF) are included as a baseline model. The implemented TTM `S_k` are parametric bases (marginal/separable/cross‑term), not TRTF; this is consistent with comparing TTM vs. TRTF.

