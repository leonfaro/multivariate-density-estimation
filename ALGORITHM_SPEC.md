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
########  TTM CORE (shared)  ########
standardizeData(X)
sampleReference(N,d)
shuffleOrdering(d)

struct MapStruct:
    type ∈ {marginal,separable,cross}
    coeffA[k], coeffB[k], coeffC[k]
    basisF[k] : {value, deriv>0}
    basisG[k], basisH[k]   # basisH_k(t, x_{1:k-1}, θ)

forwardPass(S,x)
logJacDiag(S,x)
logDetJacobian(ℓ)=Σℓ
forwardKLLoss(S,X)

monotoneIntegrator(h,t0,t)=∫exp(min(h(s),100))ds
rootFind1D(fun,target)    # interval [-10,10]
inversePass(S,z)

negativeLogLikelihood(S,X)
natsPerDim(L,N,d)=L/(N·d)
#####################################
```

## 2. Top-Level Pipeline
The overarching routine `main()` follows the composition
$$\operatorname{mainPipeline} := f_9 \circ f_8 \circ \cdots \circ f_1,$$
where
1. $f_1 = \texttt{gen\_samples}$ – generate $X$ from configuration.
2. $f_2 = \texttt{split\_data}$ – obtain $(X_{\text{tr}}, X_{\text{val}}, X_{\text{te}})$.
3. $f_3 = \texttt{fit\_TRUE}$ – fit independent parametric marginals.
   f_3b1 = fit_TTM_marginal
   f_3b2 = fit_TTM_separable
   f_3b3 = fit_TTM_cross  # train triangular transport maps via respective trainers
4. $f_4 = \texttt{fit\_TRTF}$ – fit transformation forests.
5. $f_5 = \texttt{fit\_MAF}$ – train a masked autoregressive flow.
6. $f_6 = \texttt{fit\_KS}$ – fit kernel smoother model.
7. $f_7 = \texttt{logL\_\*\_dim}$ – compute dimension-wise log-likelihoods.
8. $f_8 = \texttt{add\_sum\_row}$ – append totals to tables.
9. $f_9 = \texttt{format\_loglik\_table}$ – present final evaluation table.
Optional EDA helper functions are defined in `04_evaluation.R`.

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

### fit_KS
`fit_KS(S, config, seed) : S \to M_{KS}`
- **Description:** kernel density estimate using Gaussian kernels with bandwidth `bw.nrd0`.
- **Pre:** numeric matrices; seed scalar.
- **Post:** returns training data, bandwidth vector and test log-likelihood.
- **Randomness:** `set.seed(seed)` only influences bandwidth estimation.
- **Pseudocode:**
```
function fit_KS(S, config, seed)
    set seed to seed
    X_tr <- S.X_tr
    X_val <- S.X_val
    X_te  <- S.X_te
    h <- bw.nrd0 applied columnwise on rbind(X_tr, X_val)
    model <- (X_tr, h, config)
    model.logL_te <- logL_KS(model, X_te)
    return model
```

### logL_KS
`logL_KS(model, X) : (M_{KS}, \mathbb R^{n\times d}) \to \mathbb R`
- **Description:** Gesamt-NLL via sequential KDE evaluation.
- **Pre:** bandwidths $h>0$.
- **Post:** finite scalar.
- **Pseudocode:**
```
function logL_KS(model, X)
    ll <- predict.ks_model(model, X, 'logdensity')
    return -sum(ll)
```

### fit_TRTF
`fit_TRTF(S, config, ntree, mtry, minsplit, minbucket, maxdepth, seed, cores) : S \to M_{TRTF}`
- **Description:** one-shot fit of conditional transformation forest on train data, evaluation on test.
- **Pre:** training and test matrices present in `S`.
- **Post:** fitted forest with test log-likelihood.
- **Randomness:** `set.seed(seed)` affects forest construction.
- **Pseudocode:**
```
function fit_TRTF(S, config, ntree, mtry, minsplit, minbucket, maxdepth, seed, cores)
    set seed to seed
    X_tr <- S.X_tr
    X_te <- S.X_te
    mod <- mytrtf(X_tr, ntree, mtry, minsplit, minbucket, maxdepth, seed, cores)
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

### shuffle_ordering
`shuffle_ordering(d, seed) : (\mathbb N, \mathbb N) \to \pi`
- **Description:** optionale Permutation der Spaltenindizes.
- **Pseudocode:**
```
function shuffle_ordering(d, seed)
    set seed to seed
    return random permutation of 1..d
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
- **Description:** append a row with columnwise sums for numeric columns.
- **Pseudocode:**
```
function add_sum_row(tab, label)
    sum_row <- numeric columns summed; 'dim' <- label
    return rbind(tab, sum_row)
```

### calc_loglik_tables
`calc_loglik_tables(models, X_te) : (list, matrix) \to data.frame`
- **Description:** compute mean and standard error of negative log-likelihoods per dimension and return formatted table.
- **Pseudocode:**
```
function calc_loglik_tables(models, X_te)
    ll_true <- -predict_TRUE(models.true, X_te)
    ll_trtf <- -predict(models.trtf, X_te)
    ll_ks   <- -predict(models.ks,   X_te)
    mean_true <- colMeans(ll_true)
    se_true   <- apply(ll_true, 2, stderr)
    ... (same for trtf, ks)
    tab <- data.frame(dim, distribution,
                      true = sprintf("%.2f ± %.2f", mean_true, 2*se_true),
                      trtf = sprintf("%.2f ± %.2f", mean_trtf,2*se_trtf),
                      ks   = sprintf("%.2f ± %.2f", mean_ks,  2*se_ks),
                      ttm  = NA)
    add_sum_row(tab)
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
- `gen_samples`: $\\Theta(Nd)$.
- `split_data`: $\\Theta(N)$ for shuffling.
- `fit_TRUE`: dominated by optimization per dimension; roughly $\\Theta(d n_{tr} I)$ with iteration count $I$.
- `fit_KS`: kernel evaluations yield $\\Theta(n_{te} n_{tr} d)$ during prediction.
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
  - `ks_model`: list `X_tr`, `h`, `config`, `logL_te`.
  - `mytrtf`: list `ymod`, `forests`, `seed`, `varimp`, `config`, `logL_te`.
  - `maf_model`: fitted masked autoregressive flow from `fit_MAF`.
  - `ttm_model`: functions implementing triangular transport maps.
  - `logJacDiag(S,x)` → returns vector of log partial derivatives.
  - `logDetJacobian(logDiag)` → sum of log-diagonal entries.
  - `forwardKLLoss(S,X)` → mean forward-KL objective.
  - `negativeLogLikelihood(S,X)` → total NLL of dataset.
  - `natsPerDim(L,N,d)=L/(N·d)` → normalized NLL per dimension.


