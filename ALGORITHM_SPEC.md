# Algorithm Specification for Multivariate Density Estimation

## 1. Notation
- $d$: dimensionality of the random vector $x=(x_1,\ldots,x_d)^\top$.
- $N$: number of observations.
- $\pi(x)$: target density of $x$.
- $\eta(z)$: density of the $d$-variate standard Gaussian.
- $S: \mathbb R^d \to \mathbb R^d$: lower triangular transport map, $z = S(x)$.
- $S_k$: $k$-th component of $S$ with arguments $x_{1:k}=(x_1,\ldots,x_k)$.
- $R = S^{-1}$: inverse map, $x = R(z)$.
- $\nabla S$: Jacobian matrix of $S$.
- $\ell(x)=\sum_{k=1}^{d} \log\partial_k S_k(x)$.
- $X \in \mathbb R^{N\times d}$: matrix of samples where column $j$ corresponds to $x_j$.
- $\mathcal{N}(0,I_d)$: $d$-variate standard Gaussian distribution.

All densities are evaluated in log-space. Strictly positive parameters are transformed via `softplus`.
```
########  TTM CORE (marginal)  ########
trainMarginalMap(X_or_path)
predict.ttm_marginal(S, X, type)
forwardPass(S,x)
logJacDiag(S,x)
forwardKLLoss(S,X)
inversePass(S,z)
negativeLogLikelihood(S,X)
natsPerDim(L,N,K)
stderr(v)
######################################
```

```
########  TTM CORE (separable) ########
trainSeparableMap(X_or_path)
predict.ttm_separable(S, X, type)
######################################
```

```
########  TTM CORE (cross term) ########
trainCrossTermMap(X_or_path, warmstart_from_separable=FALSE)
predict.ttm_cross_term(S, X, type)
forwardKLLoss_ct(S, X) = mean(-rowSums(predict(S,X,"logdensity_by_dim")) - 0.5*K*log(2*pi))
- uses polynomial bases with cross terms t^r x_j^s (degrees deg_t_cross, deg_x_cross)
- Gauss-Legendre quadrature on [0,1]
- analytic derivative basis `.dpsi_dt_ct` for monotonicity checks
- optimization of forward KL via L-BFGS-B
- chunked training and prediction mit gewichteter log-sum-exp-Quadratur ohne Clipping, parallelisiert über Dimensionen
- robuste Parallelisierung mit sequentiellem Retry und Null-Koeffizienten-Fallback
- optionaler Warm-Start der g_k-Koeffizienten aus der separablen Map
######################################
```

## 2. Top-Level Pipeline
The overarching routine `main()` follows the composition
$$\operatorname{mainPipeline} := f_7 \circ f_6 \circ \cdots \circ f_1,$$
where
1. $f_1 = \texttt{gen\_samples}$ – generate $X$ from configuration.
2. $f_2 = \texttt{split\_data}$ – obtain $(X_{\text{tr}}, X_{\text{val}}, X_{\text{te}})$.
3. $f_3 = \texttt{fit\_TRUE}$ – fit independent parametric marginals.
   f_3b1 = fit_TTM_marginal
   f_3b2 = fit_TTM_separable
   f_3b3 = fit_TTM_cross  # train triangular transport maps via respective trainers
   f_3b4 = fit_TRUE_JOINT    # evaluate oracle joint density
4. $f_4 = \texttt{fit\_TRTF}$ – fit transformation forests.
5. $f_5 = \texttt{logL\_\*\_dim}$ – compute dimension-wise log-likelihoods.
6. $f_6 = \texttt{add\_sum\_row}$ – append totals to tables.
7. $f_7 = \texttt{format\_loglik\_table}$ – present final evaluation table.
Optional EDA helper functions are defined in `04_evaluation.R`.

Environment variables allow selecting alternative toy datasets. Setting `DATASET=halfmoon2d`
triggers generation of a two-moons sample via `make_halfmoon_splits` with
parameters `N_TRAIN`, `N_TEST` (capped at 250), `NOISE` and `SEED`. The created
splits are stored as `results/splits_halfmoon2d_seedXXX.rds` and the RNG state is
reset afterwards to keep subsequent model training comparable across datasets.
Evaluation `eval_halfmoon` fits TRUE, TRTF and the three TTM variants on these
splits and records per-dimension and joint negative log-likelihoods (in nats)
to `results/nll_halfmoon_seedXXX.csv`.

## 3. Module Specifications
`gen_samples(G) : G \to X`
- **Description:** sequentially draws $N=G.n$ samples from distributions specified in `G.config`.
- **Pre:** each `config[[k]]$distr` matches an R sampling routine `r<distr>`.
- **Post:** matrix $X \in \mathbb R^{N\times d}$ with column names $\{X1,\ldots,Xd\}$; invalid distribution parameters are clamped to $10^{-3}$.
- **Randomness:** uses `set.seed(G.seed)` and draws via base R generators.
- **Pseudocode:**
```
function gen_samples(G)
    set seed to G.seed
    for i in 1..G.n
        for k in 1..d
            args <- if config[k].parm is null then {} else config[k].parm(previous X_i)
            replace non-finite or <=0 args by 1e-3
            X[i,k] <- r<distr_k>(1, args)
    return X
```

### split_data
`split_data(X, seed) : (\mathbb R^{N\times d}, \mathbb N) \to (X_{\text{tr}}, X_{\text{val}}, X_{\text{te}})`
- **Description:** shuffle once and split 80/10/10.
- **Randomness:** `set.seed(seed)` for permutation.
```
function split_data(X, seed)
    set seed to seed
    idx <- random permutation of 1..N
    n_tr  <- floor(0.8*N)
    n_val <- floor(0.1*N)
    return (X[idx[1:n_tr],], X[idx[n_tr+1:n_tr+n_val],], X[idx[(n_tr+n_val+1):N],])
```
### fit_TRUE
`fit_TRUE(S, config) : S \to M_{TRUE}`
- **Description:** independent maximum-likelihood estimation per dimension.
- **Pre:** distributions in `config` supported by `neg_loglik_uni`.
function fit_TRUE(S, config)
    for k in 1..d
        init <- moment_based_start(S.X_tr[,k], config[k].distr)
        theta_k <- optimize neg_loglik_uni w.r.t. init
    logL_te <- logL_TRUE((theta), S.X_te)
    return (theta, logL_te)
```

### logL_TRUE
`logL_TRUE(M_TRUE, X) : (M_{TRUE}, \mathbb R^{n\times d}) \to \mathbb R`
- **Description:** evaluate Gesamt-NLL of `X` under the independent model.
- **Pre:** fitted parameters in `M_TRUE` are valid.
- **Post:** scalar value finite.
- **Randomness:** none.
- **Pseudocode:**
```
function logL_TRUE(M, X)
    for k in 1..d
        ll_k <- log_density_vec(X[,k], distr_k, M.theta[[k]])
    return -sum(rowSums(ll_k))
```


### true_joint_logdensity_by_dim
`true_joint_logdensity_by_dim(config, X) : (config, \mathbb R^{n\times d}) \to \mathbb R^{n\times d}`
- **Description:** evaluates log $p(x_k \mid x_{1:k-1})$ row-wise using oracle parameters.
- **Pre:** each `config[[k]]$parm` accepts data frame of previous columns.
- **Post:** matrix with finite entries.
- **Pseudocode:**
```
function true_joint_logdensity_by_dim(config, X)
    for i in 1..n
        for k in 1..d
            prev <- X[i,1:(k-1)] as data.frame
            args <- if config[k].parm null then {} else config[k].parm(prev)
            sanitize positive args; map gamma shape1/shape2
            xk <- clamp(X[i,k]) to support
            ll[i,k] <- d<distr_k>(xk, args, log=TRUE)
    return ll
```

### logL_TRUE_JOINT_dim
`logL_TRUE_JOINT_dim(config, X) : (config, \mathbb R^{n\times d}) \to \mathbb R^d`
- **Description:** average negative log-likelihood per dimension under true joint density.
- **Pseudocode:** `return -colMeans(true_joint_logdensity_by_dim(config, X))`

### logL_TRUE_JOINT
`logL_TRUE_JOINT(config, X) : (config, \mathbb R^{n\times d}) \to \mathbb R`
- **Description:** Gesamt-NLL under the oracle joint density.
- **Pseudocode:** `return -mean(rowSums(true_joint_logdensity_by_dim(config, X)))`

### fit_TRUE_JOINT
`fit_TRUE_JOINT(S, config) : S \to M_{TRUE\_JOINT}`
- **Description:** no training; evaluates test set under the known joint density.
- **Pseudocode:**
```
function fit_TRUE_JOINT(S, config)
    te_dim <- logL_TRUE_JOINT_dim(config, S.X_te)
    return (config, te_dim, sum(te_dim))
```



### fit_TRTF
`fit_TRTF(S, config, seed, cores) : S \to M_{TRTF}`
- **Description:** one-shot fit of conditional transformation forest on train data, evaluation on test.
- **Pre:** training and test matrices present in `S`.
- **Post:** fitted forest with test log-likelihood.
- **Randomness:** optional `set.seed(seed)` affects forest construction.
- **Pseudocode:**
```
function fit_TRTF(S, config, seed, cores)
    if seed is not NULL:
        set seed to seed
    X_tr <- S.X_tr
    X_te <- S.X_te
    mod <- mytrtf(X_tr, ntree = nrow(X_tr), minsplit, minbucket, maxdepth, seed, cores)
    mod.logL_te <- logL_TRTF(mod, X_te, cores)
    return mod
```


### logL_TRTF
`logL_TRTF(model, X) : (M_{TRTF}, \mathbb R^{n\times d}) \to \mathbb R`
- **Description:** Gesamt-NLL for the TRTF model.
- **Pre:** forest predictions finite.
- **Post:** scalar value.
- **Pseudocode:**
```
function logL_TRTF(model, X)
    ll <- predict.mytrtf(model, X, 'logdensity')
    return -sum(ll)
```


### standardize_data
`standardize_data(X) : \mathbb R^{N\times d} \to (\tilde X, \mu, \sigma)`
- **Description:** zentriert und skaliert jede Spalte von `X`.
- **Pseudocode:**
```
function standardize_data(X)
    μ <- mean(X, axis=0)
    σ <- std(X, axis=0) + ε
    return (X-μ)/σ , μ , σ
```

### sample_reference
`sample_reference(N, d, seed) : (\mathbb N, \mathbb N, \mathbb N) \to Z`
- **Description:** erzeugt `N` Standardnormal-Samples der Dimension `d`.
- **Pseudocode:**
```
function sample_reference(N, d, seed)
    set seed to seed
    return matrix of rnorm(N*d) reshaped (N,d)
```



### fit_TTM_marginal
`fit_TTM_marginal(S) : S \to M_{TTM}`
- **Description:** Closed-form marginal map using only training data for standardization.
- **Pre:** `S` contains `X_tr`, `X_val`, `X_te`.
- **Post:** returns parameters `(\mu, \sigma, coeffA, coeffB, coeffC=0_d)` and validation metrics.
- **Pseudocode:**
```
function fit_TTM_marginal(S)
    std <- standardize_data(S.X_tr)
    X_tr <- std.X
    \mu <- std.mu; \sigma <- std.sigma
    X_val <- (S.X_val-\mu)/\sigma
    X_te  <- (S.X_te-\mu)/\sigma
    for k in 1..d
        u <- rank_avg(X_tr[,k])/(n_tr+1)
        u <- clamp(u, [1/(n_tr+1), n_tr/(n_tr+1)])
        z <- qnorm(u)
        b_k <- max(0, cov(X_tr[,k], z)/(var(X_tr[,k])+1e-12))
        a_k <- mean(z) - b_k*mean(X_tr[,k])
        coeffA_k <- log(b_k+1e-12)
        coeffB_k <- a_k
    return (\mu, \sigma, coeffA, coeffB, 0_d)
```

### MapStruct
`MapStruct(type, coeffA, coeffB, coeffC, basisF, basisG, basisH)`
- **BasisF API:** Each `basisF_k` must implement two callables: `value(x,\theta)` and `deriv(x,\theta) > 0`.
- *basisH_k muss die Signatur **h_k(t , x_{1:k-1},\,\theta)** besitzen.*
- **Description:** Container-Objekt für Koeffizienten und Basisfunktionen einer triangularen Map.
- **Pseudocode:**
```
struct MapStruct:
    type           # string
    coeffA[k]      # α-Vektor für f_k
    coeffB[k]      # β-Vektor für g_k
    coeffC[k]      # γ-Vektor für h_k
    basisF[k], basisG[k], basisH[k]   # callable handles
```
`add_sum_row(tab, label) : (data.frame, string) \to data.frame`
- **Description:** append a row with columnwise sums for numeric columns while ignoring `NA` values.
- **Pseudocode:**
```
function add_sum_row(tab, label)
    for each column in tab:
        if column == "dim": set to label
        else if numeric: sum column with na.rm = TRUE
        else: NA
    return rbind(tab, sum_row)
```

### calc_loglik_tables
`calc_loglik_tables(models, config, X_te) : (list, list, matrix) \to data.frame`
- **Description:** compute mean and standard error of negative log-likelihoods per dimension and return formatted table.
- **Pseudocode:**
```
function calc_loglik_tables(models, config, X_te)
    ll_true <- -predict_TRUE(models.true, X_te)
    ll_trtf <- -predict(models.trtf, X_te)
    ll_true_joint <- -true_joint_logdensity_by_dim(config, X_te)
    if models.ttm exists:
        ll_ttm <- -predict(models.ttm$S, X_te)
        mean_ttm <- colMeans(ll_ttm)
        se_ttm   <- apply(ll_ttm, 2, stderr)
        se_sum_ttm <- sd(rowSums(ll_ttm)) / sqrt(nrow(ll_ttm))
    else:
        mean_ttm <- rep(NA, K); se_ttm <- rep(NA, K); se_sum_ttm <- NA
    if models.ttm_sep exists:
        ll_sep <- -predict(models.ttm_sep$S, X_te)
        mean_sep <- colMeans(ll_sep)
        se_sep   <- apply(ll_sep, 2, stderr)
        se_sum_sep <- sd(rowSums(ll_sep)) / sqrt(nrow(ll_sep))
    else:
        mean_sep <- rep(NA, K); se_sep <- rep(NA, K); se_sum_sep <- NA
    if models.ttm_cross exists:
        ll_cross <- -predict(models.ttm_cross$S, X_te)
        mean_cross <- colMeans(ll_cross)
        se_cross   <- apply(ll_cross, 2, stderr)
        se_sum_cross <- sd(rowSums(ll_cross)) / sqrt(nrow(ll_cross))
    else:
        mean_cross <- rep(NA, K); se_cross <- rep(NA, K); se_sum_cross <- NA
    mean_true <- colMeans(ll_true); se_true <- apply(ll_true,2,stderr)
    se_sum_true <- sd(rowSums(ll_true)) / sqrt(nrow(ll_true))
    mean_true_joint <- colMeans(ll_true_joint)
    se_true_joint <- apply(ll_true_joint,2,stderr)
    se_sum_true_joint <- sd(rowSums(ll_true_joint)) / sqrt(nrow(ll_true_joint))
    mean_trtf <- colMeans(ll_trtf); se_trtf <- apply(ll_trtf,2,stderr)
    fmt(x,se) = sprintf("%.2f ± %.2f", round(x,2), round(2*se,2))
    tab <- data.frame(dim, distribution,
                      true = fmt(mean_true, se_true),
                      true_joint = fmt(mean_true_joint, se_true_joint),
                      trtf = fmt(mean_trtf, se_trtf),
                      ttm  = fmt(mean_ttm,  se_ttm),
                      ttm_sep = fmt(mean_sep, se_sep),
                      ttm_cross = fmt(mean_cross, se_cross))
    sum_row <- data.frame(dim="k", distribution="SUM",
                          true=fmt(sum(mean_true), se_sum_true),
                          true_joint=fmt(sum(mean_true_joint), se_sum_true_joint),
                          trtf=fmt(sum(mean_trtf), se_sum_trtf),
                          ttm =fmt(sum(mean_ttm),  se_sum_ttm),
                          ttm_sep=fmt(sum(mean_sep), se_sum_sep),
                          ttm_cross=fmt(sum(mean_cross), se_sum_cross))
    rename columns: true->"True (marginal)", true_joint->"True (Joint)",
                    trtf->"Random Forest", ttm->"Marginal Map",
                    ttm_sep->"Separable Map",
                    ttm_cross->"Cross-term Map"
    return rbind(tab, sum_row)
```


### create_EDA_report
`create_EDA_report(X, cfg, scatter_data, table_kbl, param_list)`
- **Description:** erzeugt optional Histogramme und Streuplots der Log-Dichten
  und gibt eine Liste mit `plots`, `param_plots` und `table` zurück.
- **Location:** definiert in `04_evaluation.R`.
- **Pre/Post:** keine Seiteneffekte auf Dateien.

## 4. Randomness & Reproducibility
Each module drawing random numbers sets the RNG via `set.seed` with an integer seed. Seeds are derived from the global seed (`G.seed`) with deterministic offsets. Random variables are sampled from base R distributions (`rnorm`, `rexp`, `rbeta`, `rgamma`) or via transformation forests. All optimization routines are deterministic given these seeds. To reproduce results, record `G.seed` and any hyperparameter grid.

## 5. Complexity Notes
- `gen_samples`: $\\Theta(Nd)$.
- `split_data`: $\\Theta(N)$ for shuffling.
- `fit_TRUE`: dominated by optimization per dimension; roughly $\\Theta(d n_{tr} I)$ with iteration count $I$.
- `fit_TRTF` complexity hängt von der Tiefe der Bäume ab und ist datengesteuert.

- Marginal-TTM:  \\Theta(N d)
- Separable-TTM: \\Theta(N d) + LS-Solves
- Cross-term-TTM: \\Theta(N d B) pro SGD-Batch
## 6. Appendix A – Data Schemas
- **Configuration Entry**: list with fields
  - `distr`: string naming the distribution family.
  - `parm`: optional function returning a list of parameters given previous columns.
- **Sample Matrix** `X`: numeric matrix with columns `X1` … `Xd`.
- **Model Objects**:
  - `M_TRUE`: list `theta` (per-dimension parameter vectors), `config`, `logL_te`.
  - `mytrtf`: list `ymod`, `forests`, `seed`, `varimp`, `config`, `logL_te`.
  - `ttm_model`: functions implementing triangular transport maps.
  - `logJacDiag(S,x)` → returns vector of log partial derivatives.
  - `logDetJacobian(logDiag)` → sum of log-diagonal entries.
  - `forwardKLLoss(S,X)` → mean forward-KL objective.
  - `negativeLogLikelihood(S,X)` → total NLL of dataset.
  - `natsPerDim(L,N,d)=L/(N·d)` → normalized NLL per dimension.
  - `results_table`: matrix of mean NLL per dimension for each model.

Kurzfassung:

**Wie viele vollständige TTM‑Algorithmen?**
→ **Zwei**: *Maps from samples* (Target→Reference, $S$) und *Maps from densities* (Reference→Target, $R$).&#x20;

---

### 1) *Maps from samples* (Target→Reference, $S$)

* **Idee/Objective:** Minimiere $D_{\mathrm{KL}}(\pi \,\|\, S^{\sharp}\eta)$ ≡ maximiere die Log‑Likelihood der Pullback‑Dichte auf Fix‑Samples $x\sim\pi$.
  **Gleichungen:** (31)–(39).
  **Seiten im PDF:** p. 26 (31–34), p. 27 (35–38), p. 28 (39).&#x20;
* **Kernformeln:** Zerlegung $\log\det\nabla S(x)=\sum_k \log \partial_{x_k} S_k(x)$ und Monte‑Carlo‑Kostenfunktion $\mathcal{J}(S)=\sum_{i,k}\big(\tfrac12 S_k(X_i)^2-\log\partial_{x_k}S_k(X_i)\big)$.&#x20;

### 2) *Maps from densities* (Reference→Target, $R$)

* **Idee/Objective:** Minimiere $D_{\mathrm{KL}}(\eta \,\|\, R^{\sharp}\tilde{\pi})$ (reverse KL) bei nur bis auf Konstante bekannter Zieldichte $\tilde{\pi}$; maximiere Log‑Likelihood der Pullback‑Dichte auf Referenz‑Samples $z\sim\eta$.
  **Gleichungen:** (40)–(48).
  **Seiten im PDF:** p. 28 (40–44), p. 29 (45–48).&#x20;
* **Kernformeln:** $\log\det\nabla R(z)=\sum_k \log \partial_{z_k} R_k(z)$ und MC‑Ziel $\mathcal{J}(R)=\sum_i\big(-\log \tilde{\pi}(R(Z_i))-\sum_k \log\partial_{z_k}R_k(Z_i)\big)$.&#x20;

---

## Welche davon sind im Repo umgesetzt?

* **Umgesetzt:** *Maps from samples* (Target→Reference, $S$) mit Maximierung der (negativen) Log‑Likelihood via BFGS. Das README beschreibt explizit den **lower‑triangular** Transport $z=S(x)$ und die Optimierung „on the negative log‑likelihood“. Das entspricht genau dem Objective (31)–(39) oben. ([GitHub][1])
* **Nicht umgesetzt (derzeit ersichtlich):** *Maps from densities* (Reference→Target, $R$) mit reverse‑KL‑Objective (40)–(48) – im Repo finde ich keinen Hinweis auf das Ziehen von $z\sim\eta$ und ein Training gegen $\tilde{\pi}$. ([GitHub][1])

> **Kurzantwort:** 2 vollständige TTM‑Algorithmen im Paper. *Maps from samples* (Eq. 31–39, p. 26–28) und *Maps from densities* (Eq. 40–48, p. 28–29). In eurem Repo ist der erste umgesetzt; der zweite scheint (noch) zu fehlen. ([GitHub][1])




