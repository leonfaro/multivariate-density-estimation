# **Parametric Optimization – Detail-Workflow**
# *(R-Codex-ready: minimal prose, maximal formal spec)*
#
# | **S**  | **Function Name (file `03_param_baseline.R`)** | **Input (type / shape)**                  | **Output**                                                                       | **Core Logic / Formula (precise, no prose)**                                                                                                                                         |
# | ------ | ---------------------------------------------- | ----------------------------------------- | -------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
# | **-1** | **`cfg_init()`**                               | —                                         | `cfg` : list [K]                                                                | for each *k*: list(distr, n_par, parm=NULL)                                                                                                                                         |
# | **-0** | **`sample_pi_df()`**                           | `N, cfg`                                  | `X ∈ ℝ^{N×K}`                                                                    | draw `Z ~ N(0,I_K)` → `X ← S_inv(Z,cfg)`                                                                                                                                             |
# | **0**  | **`split_X()`**                                | `X`                                       | list `X_k, X_prev`                                                               | `X_prev ← X[,1:(k-1)]`, `X_k ← X[,k]`                                                                                                                                                |
# | **1**  | **`design_k()`**                               | `X_prev ∈ ℝ^{N×(k-1)}`                    | `D_k ∈ ℝ^{N×P_k}`                                                                | `D_k ← cbind(1, X_prev)` (*all* predecessors)                                                                                                                                        |
# | **2**  | **`link_k()`**                                 | `theta_k ∈ ℝ^{P_k·M_k}, D_k, distr_k`     | `pars_i`                                                                         | *Split θ_k into M_k blocks (per density):*<br>• normal: `(β_μ, β_σ)`<br>• gamma: `(β_shape, β_rate)`<br>`η_j ← D_k β_j`;<br>`param_j ← softplus(η_j)` if strictly >0 else identity |
# | **3**  | **`logpdf_k_vec()`**                           | `x_k, pars_i, distr_k`                    | `ℓ_{k,·} ∈ ℝ^N`                                                                  | element-wise `log f_{distr_k}(x_k; pars_i)`                                                                                                                                          |
# | **4**  | **`logLik_k()`**                               | `theta_k, X, k, cfg`                      | `ℓ_k ∈ ℝ`                                                                        | `D_k ← design_k(X_prev)`;<br>`pars ← link_k(theta_k,D_k,distr_k)`;<br>`ℓ_k ← sum( logpdf_k_vec(x_k,pars,distr_k) )`                                                                  |
# | **5**  | **`optim_k()`**                                | `k, start_k, X, cfg`                      | `thetâ_k`                                                                       | `optim(par=start_k, fn=function(θ) -logLik_k(θ,X,k,cfg), method="Nelder-Mead")`                                                                                                      |
# | **6**  | **`update_cfg()`**                             | `k, thetâ_k, cfg`                        | `cfg*`                                                                           | `cfg[[k]]$parm ← function(X_prev) { D←design_k(X_prev); link_k(thetâ_k,D,cfg[[k]]$distr) }`                                                                                         |
# | **7**  | **`fit_seq()`**                                | `X_train, cfg, starts`                    | list `(thetâ₁:K, cfg*)`                                                         | `for k in 1:K { thetâ_k←optim_k(...); cfg←update_cfg(...) }`                                                                                                                        |
# | **8**  | **`logLik_matrix()`**                          | `X, cfg`                                  | `L ∈ ℝ^{N×K}`                                                                    | loop k: row ℓ_{k,·} via steps 1–3                                                                                                                                                   |
# | **9a** | **`LL_train()`**                               | `X_train, cfg*`                           | `LL_train_avg`                                                                   | `mean( rowSums( logLik_matrix(X_train,cfg*) ) )`                                                                                                                                     |
# | **9b** | **`LL_test_param()`**                          | `X_test, cfg*`                            | `LL_param_avg`                                                                   | same formula on test batch                                                                                                                                                           |
# | **9c** | **`LL_test_true()`**                           | `X_test, cfg_true`                        | `LL_true_avg`                                                                    | use ground-truth cfg_true                                                                                                                                                           |
# | **10** | **`metrics()`**                                | `LL_true_avg, LL_param_avg`               | `Δ = (LL_true_avg − LL_param_avg)`                                               | expect Δ≈0 under perfect fit                                                                                                                                                         |
# | **11** | **`summary_table()`**                          | `cfg*, thetâ, LL_true_avg, LL_param_avg` | `data.frame` cols: <br>`dim, distr, ll_true_avg, ll_param_avg, delta, mle_param` | row per k: `mle_param ← list(thetâ_k)`                                                                                                                                              |
#
# ### Data-Shapes & Constants
#
# ```text
# N           : sample size per batch (train = test)
# K           : dimension of X
# P_k         : k (incl. intercept)  ⇒ P_k = k
# M_k         : #parameter-groups for density k  
#               normal→2, gamma→2, beta→2, …
# theta_k     : length = P_k · M_k
# design_k    : N × P_k
# ```
#
# ### Density Blocks (built-in)
#
# | `distr`  | param-groups *(η_j)* | pdf formula                                                       |
# | -------- | --------------------- | ----------------------------------------------------------------- |
# | `normal` | μ, σ > 0              | `log ϕ((x-μ)/σ) - log σ`                                          |
# | `gamma`  | shape > 0, rate > 0   | `shape·log(rate) + (shape-1)·log x - rate·x - log Γ(shape)`       |
# | `beta`   | α > 0, β > 0          | `log Γ(α+β) - log Γ(α) - log Γ(β) + (α-1)·log x + (β-1)·log(1-x)` |
#
# ### Overall Workflow (single line)
#
# ```
# X_train,X_test ← sample_pi_df(N,cfg_true)
# thetâ, cfg*   ← fit_seq(X_train,cfg_init(),starts)
# LL_train_avg   ← LL_train(X_train,cfg*)
# LL_param_avg   ← LL_test_param(X_test,cfg*)
# LL_true_avg    ← LL_test_true(X_test,cfg_true)
# Δ              ← metrics(LL_true_avg,LL_param_avg)
# tbl            ← summary_table(cfg*,thetâ,LL_true_avg,LL_param_avg)
# ```
#
# All steps are strictly sequential in *k*, use full predecessor set $X_{1:k-1}$, require no gradients, and deliver per-dimension MLEs plus out-of-sample log-likelihood diagnostics normalized by *N*.
