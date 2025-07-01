# Algorithm Specification for Multivariate Density Estimation

## 1. Notation
- $K$: dimensionality of the random vector $x=(x_1,\ldots,x_K)^\top$.
- $N$: number of observations.
- $\pi(x)$: target density of $x$.
- $\eta(z)$: density of the $K$-variate standard Gaussian.
- $S: \mathbb R^K \to \mathbb R^K$: lower triangular transport map, $z = S(x)$.
- $S_k$: $k$-th component of $S$ with arguments $x_{1:k}=(x_1,\ldots,x_k)$.
- $R = S^{-1}$: inverse map, $x = R(z)$.
- $\nabla S$: Jacobian matrix of $S$.
- $J(x)=\det(\nabla S(x))$.
- $X \in \mathbb R^{N\times K}$: matrix of samples where column $j$ corresponds to $x_j$.
- $\mathcal{N}(0,I_K)$: $K$-variate standard Gaussian distribution.

All densities are evaluated in log-space. Strictly positive parameters are transformed via `softplus`.

## 2. Top-Level Pipeline
The overarching routine `main()` follows the composition
$$\operatorname{mainPipeline} := f_9 \circ f_8 \circ \cdots \circ f_1,$$
where
1. $f_1 = \texttt{gen\_samples}$ – generate $X$ from configuration.
2. $f_2 = \texttt{train\_test\_split}$ – obtain $(X_{\text{tr}}, X_{\text{te}})$.
3. $f_3 = \texttt{fit\_TRUE}$ – fit independent parametric marginals.
4. $f_4 = \texttt{fit\_TRTF}$ – fit transformation forests.
5. $f_5 = \texttt{fit\_KS}$ – fit kernel smoother model.
6. $f_6 = \texttt{logL\_\*\_dim}$ – compute dimension-wise log-likelihoods.
7. $f_7 = \texttt{add\_sum\_row}$ – append totals to tables.
8. $f_8 = \texttt{format\_loglik\_table}$ – present final evaluation table.
9. $f_9 = \texttt{plot\_scatter\_matrix}$ – display four log-density scatter plots.
Optional EDA helper functions are defined in `04_evaluation.R`.

## 3. Module Specifications
### setup_global
`setup_global() : () \to G`
- **Description:** provide global experiment settings such as sample size $N$, configuration list and RNG seed.
- **Pre:** none.
- **Post:** returns list $G=(n,\text{config},\text{seed},\text{split\_ratio},h_{\text{grid}},\text{model\_ids},p_{\max})$.
- **Randomness:** none.
- **Pseudocode:**
```
function setup_global(cfg)
    n <- 50
    seed <- 42
    split_ratio <- 0.5
    h_grid <- 1:p_max
    model_ids <- {"TRUE"}
    return (n, cfg, seed, split_ratio, h_grid, model_ids, p_max)
```

### gen_samples
`gen_samples(G) : G \to X`
- **Description:** sequentially draws $N=G.n$ samples from distributions specified in `G.config`.
- **Pre:** each `config[[k]]$distr` matches an R sampling routine `r<distr>`.
- **Post:** matrix $X \in \mathbb R^{N\times K}$ with column names $\{X1,\ldots,XK\}$; invalid distribution parameters are clamped to $10^{-3}$.
- **Randomness:** uses `set.seed(G.seed)` and draws via base R generators.
- **Pseudocode:**
```
function gen_samples(G)
    set seed to G.seed
    for i in 1..G.n
        for k in 1..K
            args <- if config[k].parm is null then {} else config[k].parm(previous X_i)
            replace non-finite or <=0 args by 1e-3
            X[i,k] <- r<distr_k>(1, args)
    return X
```

### train_test_split
`train_test_split(X, r, seed) : (\mathbb R^{N\times K}, [0,1], \mathbb N) \to (X_{\text{tr}}, X_{\text{te}})`
- **Description:** shuffle rows of $X$ reproducibly and split into training and test matrices.
- **Pre:** $0<r<1$.
- **Post:** $X_{\text{tr}}$ has $\lfloor rN\rfloor$ rows; $X_{\text{te}}$ contains the remainder.
- **Randomness:** `set.seed(seed+1)` for permutation.
- **Pseudocode:**
```
function train_test_split(X, r, seed)
    set seed to seed+1
    idx <- random permutation of 1..N
    n_tr <- floor(r*N)
    return (X[idx[1:n_tr],], X[idx[(n_tr+1):N],])
```

### fit_TRUE
`fit_TRUE(X_tr, X_te, config) : (\mathbb R^{n_{tr}\times K}, \mathbb R^{n_{te}\times K}, config) \to M_{TRUE}`
- **Description:** independent maximum-likelihood estimation per dimension.
- **Pre:** distributions in `config` supported by `neg_loglik_uni`.
- **Post:** list with parameter estimates `theta` and test log-likelihood `logL_te`.
- **Randomness:** none once data is given.
- **Pseudocode:**
```
function fit_TRUE(X_tr, X_te, config)
    for k in 1..K
        init <- moment_based_start(X_tr[,k], config[k].distr)
        theta_k <- optimize neg_loglik_uni w.r.t. init
    logL_te <- logL_TRUE((theta), X_te)
    return (theta, logL_te)
```

### logL_TRUE
`logL_TRUE(M_TRUE, X) : (M_{TRUE}, \mathbb R^{n\times K}) \to \mathbb R`
- **Description:** evaluate mean negative log-likelihood of `X` under the independent model.
- **Pre:** fitted parameters in `M_TRUE` are valid.
- **Post:** scalar value finite.
- **Randomness:** none.
- **Pseudocode:**
```
function logL_TRUE(M, X)
    for k in 1..K
        ll_k <- log_density_vec(X[,k], distr_k, M.theta[[k]])
    return -mean(rowSums(ll_k))
```

### fit_KS
`fit_KS(X_tr, X_te, config, seed) : (...) \to M_{KS}`
- **Description:** kernel density estimate using Gaussian kernels with bandwidth `bw.nrd0`.
- **Pre:** numeric matrices; seed scalar.
- **Post:** returns training data, bandwidth vector and test log-likelihood.
- **Randomness:** `set.seed(seed)` only influences bandwidth estimation.
- **Pseudocode:**
```
function fit_KS(X_tr, X_te, config, seed)
    set seed to seed
    h <- bw.nrd0 applied columnwise on X_tr
    model <- (X_tr, h, config)
    model.logL_te <- logL_KS(model, X_te)
    return model
```

### logL_KS
`logL_KS(model, X) : (M_{KS}, \mathbb R^{n\times K}) \to \mathbb R`
- **Description:** mean negative log-likelihood via sequential KDE evaluation.
- **Pre:** bandwidths $h>0$.
- **Post:** finite scalar.
- **Pseudocode:**
```
function logL_KS(model, X)
    ll <- predict.ks_model(model, X, 'logdensity')
    return -mean(ll)
```

### fit_TRTF
`fit_TRTF(X_tr, X_te, config, grid, folds, seed) : (...) \to M_{TRTF}`
- **Description:** cross-validated conditional transformation forests.
- **Pre:** grid expands to finite hyperparameter combinations.
- **Post:** fitted forest with `best_cfg`, CV loss and test log-likelihood.
- **Randomness:** `set.seed(seed)` affects forest construction and CV splits.
- **Pseudocode:**
```
function fit_TRTF(X_tr, X_te, config, grid, folds, seed)
    set seed to seed
    best_val <- Inf
    for each cfg in expand.grid(grid)
        split data into folds
        for each fold
            m <- mytrtf(training-fold, cfg)
            val_loss <- -mean(predict(m, validation-fold, 'logdensity'))
        if mean(val_loss) < best_val
            best_cfg <- cfg
    final <- mytrtf(X_tr, best_cfg)
    final.logL_te <- logL_TRTF(final, X_te)
    return final
```

### logL_TRTF
`logL_TRTF(model, X) : (M_{TRTF}, \mathbb R^{n\times K}) \to \mathbb R`
- **Description:** mean negative log-likelihood for the TRTF model.
- **Pre:** forest predictions finite.
- **Post:** scalar value.
- **Pseudocode:**
```
function logL_TRTF(model, X)
    ll <- predict.mytrtf(model, X, 'logdensity')
    return -mean(ll)
```

### standardize_data
`standardize_data(X) : \mathbb R^{N\times K} \to (\tilde X, \mu, \sigma)`
- **Description:** zentriert und skaliert jede Spalte von `X`.
- **Pseudocode:**
```
function standardize_data(X)
    μ <- mean(X, axis=0)
    σ <- std(X, axis=0) + ε
    return (X-μ)/σ , μ , σ
```

### sample_reference
`sample_reference(N, K, seed) : (\mathbb N, \mathbb N, \mathbb N) \to Z`
- **Description:** erzeugt `N` Standardnormal-Samples der Dimension `K`.
- **Pseudocode:**
```
function sample_reference(N, K, seed)
    set seed to seed
    return matrix of rnorm(N*K) reshaped (N,K)
```

### shuffle_ordering
`shuffle_ordering(K, seed) : (\mathbb N, \mathbb N) \to \pi`
- **Description:** optionale Permutation der Spaltenindizes.
- **Pseudocode:**
```
function shuffle_ordering(K, seed)
    set seed to seed
    return random permutation of 1..K
```


### MapStruct
`MapStruct(type, coeffA, coeffB, coeffC, basisF, basisG, basisH)`
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
### evaluate_all
`evaluate_all(X_te, model_list) : (\mathbb R^{n_{te}\times K}, list) \to \text{data.frame}`
- **Description:** compute negative log-likelihoods for named models via dynamic dispatch.
- **Pseudocode:**
```
function evaluate_all(X_te, model_list)
    for each model with name id
        loss <- logL_<id>(model, X_te)
        store (id, loss)
    sort by loss
    return table
```

### add_sum_row
`add_sum_row(tab, label) : (data.frame, string) \to data.frame`
- **Description:** append a row with columnwise sums for numeric columns.
- **Pseudocode:**
```
function add_sum_row(tab, label)
    sum_row <- numeric columns summed; 'dim' <- label
    return rbind(tab, sum_row)
```

### calc_loglik_tables
`calc_loglik_tables(models, config) : (list, list) \to data.frame`
- **Description:** compute mean negative log-likelihoods per dimension for all models.
- **Pseudocode:**
```
function calc_loglik_tables(models, config)
    tab <- data.frame(dim, distribution,
                      logL_baseline, logL_trtf, logL_ks)
    tab <- add_sum_row(tab)
    return tab
```

### calc_loglik_sds
`calc_loglik_sds(models, S, config) : (list, list, list) \to data.frame`
- **Description:** standard deviations of the negative log-likelihoods across test observations.
- **Pseudocode:**
```
function calc_loglik_sds(models, S, config)
    sd_true <- sd(-log_density_true(S$X_te)) per dimension
    sd_trtf <- sd over -predict(TRTF, S$X_te)
    sd_ks   <- sd over -predict(KS, S$X_te)
    tab <- data.frame(dim, distribution, sd_true, sd_trtf, sd_ks)
    add_sum_row(tab)
```

### format_loglik_table
`format_loglik_table(tab_mean, tab_sd) : (data.frame, data.frame) \to data.frame`
- **Description:** combine mean and sd tables to strings of the form `"m ± s"`.
- **Pseudocode:**
```
function format_loglik_table(tab_mean, tab_sd)
    for each numeric column col
        tab[col] <- sprintf("%.2f ± %.2f", round(tab_mean[col],2),
                            round(2 * tab_sd[col],2))
    return tab
```


### create_EDA_report
`create_EDA_report(X, cfg, scatter_data, table_kbl, param_list)`
- **Description:** erzeugt optional Histogramme und Streuplots der Log-Dichten
  und gibt eine Liste mit `plots`, `param_plots` und `table` zurück.
- **Location:** definiert in `04_evaluation.R`.
- **Pre/Post:** keine Seiteneffekte auf Dateien.

## 4. Randomness & Reproducibility
Each module drawing random numbers sets the RNG via `set.seed` with an integer seed. Seeds are derived from the global seed (`G.seed`) with deterministic offsets. Random variables are sampled from base R distributions (`rnorm`, `rexp`, `rbeta`, `rgamma`) or via transformation forests and kernel smoothing. All optimization routines are deterministic given these seeds. To reproduce results, record `G.seed` and any hyperparameter grid.

## 5. Complexity Notes
- `gen_samples`: $\Theta(NK)$.
- `train_test_split`: $\Theta(N)$ for shuffling.
- `fit_TRUE`: dominated by optimization per dimension; roughly $\Theta(K n_{tr} I)$ with iteration count $I$.
- `fit_KS`: kernel evaluations yield $\Theta(n_{te} n_{tr} K)$ during prediction.
- `fit_TRTF` complexity hängt von der Tiefe der Bäume ab und ist datengesteuert.

## 6. Appendix A – Data Schemas
- **Configuration Entry**: list with fields
  - `distr`: string naming the distribution family.
  - `parm`: optional function returning a list of parameters given previous columns.
- **Sample Matrix** `X`: numeric matrix with columns `X1` … `XK`.
- **Model Objects**:
  - `M_TRUE`: list `theta` (per-dimension parameter vectors), `config`, `logL_te`.
  - `ks_model`: list `X_tr`, `h`, `config`, `logL_te`.
  - `mytrtf`: list `ymod`, `forests`, `seed`, `varimp`, `config`, `best_cfg`, `logL_te`.
  - `ttm_model`: functions implementing triangular transport maps.

