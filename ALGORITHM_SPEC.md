
---

# Algorithm Specification for Multivariate Density Estimation 

**Date:** 2025‑08‑25
**Units:** all likelihoods in **nats** (natural log). All densities must be evaluated in log-space. Strictly positive parameters are transformed via `softplus`.
**Shapes:** `predict(..., type="logdensity_by_dim")` returns **N×K**; `predict(...,"logdensity")` returns length **N**.
**Monotonicity:** all transport components $S_k$ are strictly increasing in $x_k$.
**Train‑only standardization:** $\mu,\sigma$ from **train**; reuse on val/test; Jacobian includes $-\log\sigma_k$ per dimension.

---

## 1. Notation & Contracts

* $d$: Dimension; $x=(x_1,\dots,x_d)^\top\in\mathbb{R}^d$.
* $N$: #Samples. Matrices $X_{\text{tr}},X_{\text{val}},X_{\text{te}}\in\mathbb{R}^{n_\bullet\times d}$.
* $\pi(x)$: Zieldichte. $\eta(z)=\mathcal{N}(0,I_d)$.
* **Triangular map** $S:\mathbb{R}^d\to\mathbb{R}^d$, $z=S(x)$, Komponenten $S_k(x_{1:k})$.
* **Jacobian diag:** $\partial_k S_k(x)=\frac{\partial}{\partial x_k}S_k(x)$.
* **Log‑Jacobian:** $\ell(x)=\sum_{k=1}^d \log \partial_k S_k(x)$.
* **Predict‑API contract:**

  * `logdensity_by_dim`: $L\in\mathbb{R}^{N\times d}$ mit $L_{ik}=\log p(x^{(i)}_k\mid x^{(i)}_{1:k-1})$ bzw. für Flows $=-\tfrac12 z_{ik}^2-\tfrac12\log(2\pi)+\log\partial_k S_k(x^{(i)})$.
  * `logdensity`: Länge‑$N$ Vektor $L^{(j)}=\sum_{k=1}^d L_{jk}$.
  * **Acceptance invariant:** `rowSums(logdensity_by_dim) == logdensity` numerisch (|Δ| ≤ 1e‑10).

---

## 2. Data Generating Processes (DGPs)

### 2.1 Config‑4D (triangular, heterogeneous marginals)

For $k=1,\dots,4$, generiere sequentiell:

* $X_1\sim\mathrm{Norm}(0,1)$.
* $X_2\mid X_1\sim\mathrm{Exp}(\text{rate}= \mathrm{softplus}(X_1))$.
* $X_3\mid X_{1:2}\sim\mathrm{Beta}\big(\alpha=\mathrm{softplus}(X_2),\ \beta=\mathrm{softplus}(X_1)\big)$.
* $X_4\mid X_{1:3}\sim\mathrm{Gamma}\big(\text{shape}=\mathrm{softplus}(X_3),\ \text{scale}=\mathrm{softplus}(X_2)\big)$.

**Permutation:** Optional `perm` reorder der Spalten; `config` wird entsprechend re‑indiziert.

**Pseudocode**

```
function gen_config4d(n, seed)
    set.seed(seed)
    for i in 1..n:
        x1 <- rnorm(1,0,1)
        x2 <- rexp(1, rate = softplus(x1))
        x3 <- rbeta(1, shape1 = softplus(x2), shape2 = softplus(x1))
        x4 <- rgamma(1, shape = softplus(x3), scale = softplus(x2))
        X[i,] <- c(x1,x2,x3,x4)
    return X
```

### 2.2 Half‑Moon‑2D (toy density; N & NOISE gesteuert)

**Parameters:** `SEED`, `N_TRAIN`, `N_TEST`, `NOISE`, `val_frac=0.2`.
Train/Test werden **unabhängig** generiert (Seeds `seed`, `seed+1`), Val aus Train via (`seed+2`).

**Base curves (before noise):**
Upper arc $M_1:\ (\cos t,\ \sin t)$, lower arc $M_2:\ (1-\cos t,\ -\sin t + 0.5)$, $t\sim\mathrm{Unif}[0,\pi]$.
**Noise:** addiere iid $N(0,\sigma^2 I_2)$ pro Punkt; mische Reihenfolge (shuffle).
**Splits:** Val‑Größe $n_\text{val}=\max(10,\lfloor 0.2\cdot N_\text{train}\rfloor)$; dann `X_tr`, `X_val`, `X_te`.

**Pseudocode**

```
function make_halfmoon_splits(n_tr_total, n_te_total, noise, seed, val_frac=0.2)
    X_tr_full <- generate_two_moons(n_tr_total, noise, seed)
    X_te      <- generate_two_moons(n_te_total, noise, seed+1)
    n_val <- max(10, round(val_frac*n_tr_total))
    set.seed(seed+2)
    idx_val <- sample.int(n_tr_total, n_val)
    X_val <- X_tr_full[idx_val,]; X_tr <- X_tr_full[-idx_val,]
    return list(X_tr=X_tr, X_val=X_val, X_te=X_te, K=2, meta={seed, noise, n_train=n_tr_total, n_test=n_te_total, n_val})
```

---

## 3. Objectives (Flows: “maps from samples” – Forward‑KL)

Minimiere (bis Konstante) für triangular $S$:

$$
\mathcal{J}(S)\ =\ \mathbb{E}_{x\sim\pi}\!\Big[\tfrac12\|S(x)\|^2\ -\ \sum_{k=1}^d \log\partial_k S_k(x)\Big].
$$

Monte‑Carlo mit $X_{\text{tr}}$: $\widehat{\mathcal{J}}=\frac1{n_{\text{tr}}}\sum_i \big(\tfrac12\|S(X_i)\|^2 - \sum_k \log\partial_k S_k(X_i)\big)$.

**Log‑density under flow (per sample)**
$\log \hat p(x) = \sum_{k=1}^d \Big(-\tfrac12 z_k(x)^2 - \tfrac12\log(2\pi) + \log\partial_k S_k(x)\Big).$

---

## 4. Models

### 4.1 TRUE (marginal) – independent parametric fit

**Goal:** MLE je Dimension gemäß `config[[k]]$distr` (Norm/Exp/Beta/Gamma).
**Eval:** $\log p(x)=\sum_k \log p_k(x_k;\hat\theta_k)$.

**Pseudocode**

```
fit_TRUE(S, config):
    for k in 1..d:
        theta_k <- MLE_univariate(S.X_tr[,k], distr=config[k].distr)
    return list(theta, config)
predict_TRUE(M, X, type="logdensity_by_dim"):
    L[,k] <- d<distr_k>(X[,k], M.theta[[k]], log=TRUE)
    return L or rowSums(L)
```

### 4.2 TRUE (joint/oracle) – conditional log‑densities from config

**Goal:** nutze bekannte bedingte Dichten des DGP.
**Eval:** Matrix $L_{ik}=\log p(x^{(i)}_k\mid x^{(i)}_{1:k-1})$ durch Auswerten der `parm`‑Funktionen.

**Pseudocode**

```
true_joint_logdensity_by_dim(config, X):
    for i in 1..n:
        for k in 1..d:
            args <- (is.null(config[k].parm)) ? {} : config[k].parm(as.data.frame(X[i,1:(k-1)],...))
            L[i,k] <- d<distr_k>(X[i,k], args, log=TRUE)
    return L
```

### 4.3 TRTF – Conditional Transformation Forests

**Idea:** faktorisierte Dichte $\prod_k p(x_k\mid x_{1:k-1})$ via Transformationswälder.
**Train:**

* $k=1$: marginales Transformationsmodell.
* $k\ge2$: `traforest` (z. B. partykit/mlt) für $x_k\mid x_{1:k-1}$.
  **Predict:** `logdensity_by_dim`: $(\log p(x_1),\log p(x_2\mid x_1),\dots)$; `logdensity` = Summe.

**Pseudocode**

```
fit_TRTF(S, config, seed):
    set.seed(seed)
    model$k1 <- fit_transformation_model(S.X_tr[,1])
    for k=2..d: model$forest[[k]] <- fit_traforest(y=S.X_tr[,k], X=S.X_tr[,1:(k-1)], ctrl={minsplit,minbucket,maxdepth,ntree})
    return model
predict_TRTF(M, X, type):
    L[,1] <- logdens_k1(X[,1]); for k=2..d: L[,k] <- logdens_forest_k(M$forest[[k]], X[,1:k])
    return L or rowSums(L)
```

### 4.4 TTM – Marginal (diagonal; Eq. (20))

**Parametrization:** $S_k(x)=a_k + b_k \cdot x_{\text{std},k}$, $b_k>0$.
**Train:** *Normal‑Scores* pro Spalte (train‑only Standardisierung):

* Ränge mit `ties.average`; $u_i=\frac{\mathrm{rank}(x_i)}{n_{\text{tr}}+1}\in\big[\tfrac1{n+1},\tfrac{n}{n+1}\big]$.
* $z_i^\star=\Phi^{-1}(u_i)$.
* $b_k=\max(0,\ \mathrm{cov}(x_k,z^\star)/(\mathrm{var}(x_k)+\varepsilon))$; $a_k=\overline{z^\star}-b_k\overline{x_k}$.
* Speichere $\log b_k$ (numerisch stabil), $a_k$, $\mu_k,\sigma_k$.

**Predict:**

* Standardisiere $x_k\mapsto (x_k-\mu_k)/\sigma_k$; $z_k=a_k+b_k x_{\text{std},k}$.
* $\log\partial_k S_k = \log b_k - \log\sigma_k$ (konstant über $x$).
* **By‑dim LD:** $-\tfrac12 z_k^2 - \tfrac12\log(2\pi) + \log b_k - \log\sigma_k$.

### 4.5 TTM – Separable (Eq. (21); per‑k Objective Eq. (39))

**Parametrization:** $S_k(x_{1:k})=g_k(x_{1:k-1})+f_k(x_k)$ mit $f'_k(x_k)>0$.
**Bases (typ.):**

* $g_k$: Polynomfeatures in $x_{1:k-1}$ (Grad 1–3).
* $f_k$: monotone Basis in $x_k$ (z. B. $[x_k,\ \mathrm{erf}(x_k)]$); Ableitungsbasis $B=\frac{\partial}{\partial x_k}\text{basis}_f$.
  **Train (per k):** Orthogonalisierung und Ridge‑Regularisierung:

$$
J_k(c)=\tfrac12\|A c\|^2\ -\ \sum_{i}\log(Bc)_i\ +\ \tfrac{\lambda}{2}\big(\|Dc\|^2+\|c\|^2\big),
$$

mit $A=(I-P_\text{non}M)P_\text{mon}$, $D=MP_\text{mon}$, $M=(P_\text{non}^\top P_\text{non}+\lambda I)^{-1}P_\text{non}^\top$; Optimierung via L‑BFGS (Start $c\equiv 1$), Nebenbedingung $(Bc)_i>0$.
**Predict:** $z_k=g_k+f_k$; $\log\partial_k S_k=\log\big((B c)_i\big)-\log\sigma_k$.

### 4.6 TTM – Cross‑Term (Eq. (22); integrated exp of $h_k$)

**Parametrization:**

$$
S_k(x)=g_k(x_{1:k-1})+\int_0^{x_k}\!\exp\big(h_k(t,x_{1:k-1})\big)\,dt,
$$

Monotonie via $\exp(\cdot)$.
**Bases:**

* $g_k$: wie separabel.
* $h_k(t,x_{1:k-1})$: Polynom‑/RB‑Basen inkl. Cross‑Terms $t^r\cdot x_j^s$ (Grade steuerbar).
  **Numerik:** Gauss‑Legendre‑Quadratur (diskrete Knoten $t_q$), **stabile log‑sum‑exp** Aggregation (ohne “hard clipping”), Chunking & ggf. Parallelisierung pro $k$; Fallback‑Strategie (Null‑Koeffizienten), optionaler Warm‑Start aus separabler Map.
  **Train:** minimize Forward‑KL per $k$ über $\alpha_k$ (für $g_k$) und $\beta_k$ (für $h_k$) via L‑BFGS‑B.
  **Predict:**
* $z_k=g_k + \int_0^{x_k}\exp(h_k)\,dt$ (quadraturbasiert).
* $\partial_k S_k(x)=\exp(h_k(x_k,x_{1:k-1}))$ ⇒ $\log\partial_k S_k=h_k(x_k,x_{1:k-1})-\log\sigma_k$.
* **By‑dim LD:** $-\tfrac12 z_k^2-\tfrac12\log(2\pi)+ h_k(x)-\log\sigma_k$.

> **Gemeinsame Konstanten (alle TTM):** Pro Dimension addiert **$-\tfrac12\log(2\pi)$**. Standardisierung train‑only ⇒ **$-\log\sigma_k$** im Jacobian jeder Dim.

---

## 5. Evaluation & Tables

### 5.1 `calc_loglik_tables(models, config, X_te)`

Eingang: Liste `models` mit Einträgen (sofern vorhanden)
`true`, `true_joint`, `trtf`, `ttm` (marginal), `ttm_sep`, `ttm_cross`.
Ausgang: Formatierte Tabelle **mean NLL ± 2·SE** je Dimension + **SUM‑Zeile** (nats).

**Pseudocode**

```
compute_by_model(X_te):
    L_true_dim       <- predict_TRUE(models.true, X_te, "logdensity_by_dim")
    L_true_joint_dim <- true_joint_logdensity_by_dim(config, X_te)
    L_trtf_dim       <- predict(models.trtf, X_te, "logdensity_by_dim")
    L_ttm_dim        <- predict(models.ttm$S, X_te, "logdensity_by_dim")        (if present)
    L_sep_dim        <- predict(models.ttm_sep$S, X_te, "logdensity_by_dim")    (if present)
    L_cross_dim      <- predict(models.ttm_cross$S, X_te, "logdensity_by_dim")  (if present)
aggregate:
    mean_per_dim  <- colMeans(-L_*_dim)
    se_per_dim    <- apply(-L_*_dim, 2, stderr)
    se_sum        <- sd(rowSums(-L_*_dim))/sqrt(nrow(X_te))
format:
    "True (marginal)", "True (Joint)", "Random Forest",
    "Marginal Map", "Separable Map", "Cross-term Map".
```

### 5.2 Timing table

Zeilen: `True (marginal)`, `True (Joint)`, `Random Forest`, `Marginal Map`, `Separable Map`, `Cross-term Map`; Spalten `train_sec`, `test_sec`, `total_sec`.

### 5.3 Half‑Moon evaluation (`eval_halfmoon`)

Rückgabe‑Data‑Frame mit Spalten:
`model`, `mean_joint_nll`, `se_joint`, `per_dim_nll_1`, `per_dim_nll_2`.
Schreibt CSV: `results/nll_halfmoon_seed%03d.csv`.

**Acceptance (eval):**

* `predict(...,"logdensity_by_dim")` ist **N×K**, keine NA/Inf.
* `rowSums(LD_dim) == predict(...,"logdensity")` numerisch (|Δ|≤1e‑10).
* **NLL in nats** (inkl. $-\tfrac12\log(2\pi)$ & $-\log\sigma_k$ für TTM).

---

## 6. Standardization, Seeds & Determinism

* **Standardization (TTM):** compute $\mu,\sigma$ auf **Train**; `X_val`, `X_te` nutzen dieselben; Jacobian enthält $-\log\sigma_k$.
* **Ranks:** `ties.average`; $u\in[\frac1{n+1},\frac{n}{n+1}]$.
* **Seeds:** global `SEED` (train/eval); für Half‑Moon: Generator nutzt `seed`, `seed+1`, `seed+2` (train/test/val‑Split).
* **Determinismus:** Gleiche Seeds ⇒ identische Tabellen und PNGs (bis auf FP‑Rauschen).

---

## 7. Plotting (Half‑Moon Panels)

* **Grid:** gemeinsames Gitternetz $[x\text{-lim}]\times[y\text{-lim}]$, `grid_n ≥ 200`.
* **Global levels:** identische iso‑Log‑Dichte‑Schwellen (Quantile $\{0.90,0.70,0.50\}$) **gemeinsam** über alle Modelle.
* **Points sichtbar:** einheitliche `xlim/ylim` in allen Facetten; Punkte oberhalb der Konturen.
* **Output:** `results/halfmoon_panels_seed%03d.png`. Optional `show_plot`.

---

## 8. Numerical Stability & Edge Cases

* Additive epsilons in Varianz‑/Ableitungs‑Termen (z. B. $+1e{-12}$).
* **Cross‑term:** log‑sum‑exp bei Quadratur; kein hard clipping; Chunking großer Batches; Fallback auf stabilen Null‑Start.
* Sicherstellen: **keine** NA/Inf in `predict`; Support‑Clamping (z. B. Beta/Gamma).
* **Projection (separable):** Orthogonalisierung $A,D$ stabilisiert Optimierung und erzwingt $f'_k>0$.

---

## 9. Complexity (dominant terms)

* TRUE (marginal): $\tilde{\mathcal{O}}(d\,n_{\text{tr}}\,I)$ (per‑k Optimierung).
* TRUE (joint): $\mathcal{O}(d\,n_{\text{te}})$.
* TRTF: datenabhängig; typ. $\mathcal{O}(\text{ntree}\cdot n_{\text{tr}}\log n_{\text{tr}})$.
* TTM‑marginal: $\mathcal{O}(d\,n_{\text{tr}})$.
* TTM‑separable: $\mathcal{O}(d\,n_{\text{tr}}) +$ kleine LS/Quasi‑Newton‑Solves.
* TTM‑cross: $\mathcal{O}(d\,n_{\text{tr}}\cdot Q)$ (Quadraturknoten $Q$), plus Optim.

---

## 10. Top‑Level Pipelines

### 10.1 `main()` (default: Config‑4D; env‑switch für Half‑Moon)

```
main():
    dataset <- Sys.getenv("DATASET", "config4d")
    if dataset == "halfmoon2d":
        seed  <- as.integer(Sys.getenv("SEED"))
        ntr   <- as.integer(Sys.getenv("N_TRAIN"))
        nte   <- as.integer(Sys.getenv("N_TEST"))
        noise <- as.numeric(Sys.getenv("NOISE"))
        S <- make_halfmoon_splits(ntr, nte, noise, seed, val_frac=0.2)
        saveRDS(S, sprintf("results/splits_%s_seed%03d.rds", dataset, seed))
        df <- eval_halfmoon(S)      # writes CSV; sets results_table
        return(df)                  # early-return
    else:
        seed <- as.integer(Sys.getenv("SEED", 42)); set.seed(seed)
        X <- gen_config4d(n=50, seed)
        S <- split_data(X, seed); S <- apply_perm(S, perm); cfg <- config[perm]
        mods <- list(
            true        = fit_TRUE(S, cfg),
            true_joint  = fit_TRUE_JOINT(S, cfg),
            trtf        = fit_TRTF(S, cfg, seed=seed),
            ttm         = trainMarginalMap(S, seed=seed),
            ttm_sep     = trainSeparableMap(S, seed=seed),
            ttm_cross   = trainCrossTermMap(S, seed=seed)
        )
        tab <- calc_loglik_tables(mods, cfg, S$X_te); print(tab)
        times <- timing_table(mods); print(times)
        return(tab)
```

### 10.2 `main_moon.R` (minimal)

* Setze **nur** `SEED`, `N`, `NOISE` ganz oben.
* `Sys.setenv(DATASET="halfmoon2d", SEED=..., N_TRAIN=N, N_TEST=N, NOISE=...)`.
* `source("main.R"); tab <- main()`; prüfe CSV & PNG; drucke `tab`.

---

## 11. Tests & Acceptance

1. **Shapes:** `logdensity_by_dim` ist $N\times K$; `logdensity` Länge $N$.
2. **Consistency:** `rowSums(by_dim) == logdensity` (|Δ|≤1e‑10).
3. **Constants:** Jede TTM‑Dim enthält $-\tfrac12\log(2\pi)$ und **$-\log\sigma_k$** (train‑only standardization).
4. **Ranks:** normal scores mit `ties.average` und $u\in[\frac1{n+1},\frac{n}{n+1}]$.
5. **Repro:** identische Tabellen/PNGs bei gleichem Seed.
6. **No NA/Inf:** weder in Training noch in `predict`.
7. **Permutation sanity (Config‑4D):** Bei unabhängigen Marginals bleibt SUM‑NLL invariant bzgl. Spaltenpermutation.
8. **Timing:** `total_sec == train_sec + test_sec` je Modell.

---

## 12. Appendix: Function Stubs (Pseudocode‑Signaturen)

```
# Globals
softplus(x) = log1p(exp(x))
stderr(v)   = sd(v)/sqrt(length(v))

# Standardization
standardize_data(X): returns (X_std, mu, sigma)

# TTM-marginal core
trainMarginalMap(S, seed):
    (mu,sigma) <- from S.X_tr; compute (a_k, b_k>0) via normal scores; store coeffA=log b_k, coeffB=a_k
    return list(S=MapStruct_marginal, NLL_tr/val/te, time_train, time_pred)

predict.ttm_marginal(S, X, type):
    compute z, log|∂S/∂x|; return by-dim or sum

# TTM-separable core
trainSeparableMap(S, degree_g=2, lambda=1e-3, eps=1e-6, seed):
    for k=1..K: build P_non, P_mon, B; solve c via J_k(c); store coeffs
    return list(S=MapStruct_separable, ...)

predict.ttm_separable(S, X, type): as above, using (g_k,f_k) and Bc>0

# TTM-cross core
trainCrossTermMap(S, seed, ...):
    for k=1..K: build bases for g_k, h_k(t,x_prev); GL-quadrature; stable log-sum-exp; optimize (α,β)
    return list(S=MapStruct_cross, ...)

predict.ttm_cross(S, X, type):
    z_k = g_k + integral exp(h_k); log|∂| = h_k - log σ_k

# TRTF
fit_TRTF(S, config, seed): train marginal model for k=1 and forests for k≥2
predict.mytrtf(M, X, type): by-dim conditional log-densities and sum

# TRUE models
fit_TRUE(S, config), predict_TRUE(...)
fit_TRUE_JOINT(S, config), true_joint_logdensity_by_dim(config, X)
```

---

### One‑line Summary

**Pipeline:** Config‑4D (*oracle + forests + 3×TTM*) und Half‑Moon‑2D (N/NOISE) – **alles in nats**, **by‑dim N×K** & **joint N**, deterministisch über Seeds; TTM‑Konstanten und Standardisierung korrekt eingerechnet; Cross‑term numerisch stabil via log‑sum‑exp‑Quadratur.

