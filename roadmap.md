Nachfolgend eine **modulare Skript-Gliederung** (Dateien `.R`), wobei innerhalb jeder Datei nur **algorithmischer Pseudocode** steht.
Jeder Block beschreibt

* Zweck der Datei
* Wichtigste Funktionen (Prozeduren)
* Mathematisch präzise Schritte, gegliedert in Unterpunkte

Dadurch kann man die `config`-Liste oder die Zahl der Dimensionen $K$ ändern, ohne irgendeine andere Datei anfassen zu müssen, und jederzeit neue Modelle nachrüsten.

---

### **Script 1: 00\_globals.R**

```


FUNCTION setup_global():
    INPUT  : —
    OUTPUT : list G = (N, config, seed, split_ratio, H_grid, model_ids, P_max)

    1  N           ← 500                                 # Stichprobenumfang
    2  config      ← USER-EINGABE                        # Liste mit K Elementen
    3  seed        ← 42                                  # Reproduzierbarkeit
    4  split_ratio ← 0.70                                # Train-Anteil
    5  P_max       ← 6                                   # höchster Polynomgrad
    6  H_grid      ← {h ∈ ℕ | 1 ≤ h ≤ P_max}             # Hyperparameter für TTM
    7  model_ids   ← {"TTM", "TRUE"}                     # erweiterbar
    8  RETURN G
```

---

### **Script 2: 01\_data\_generation.R**

```
FUNCTION gen_samples(G):
    INPUT : G = (N, config, seed, …)
    OUTPUT: Matrix X ∈ ℝ^{N×K}         # gleiche Reihenfolge wie config

    1  SET rng ← seed
    2  FOR i = 1,…,N:
    3      FOR k = 1,…,K:                        # sequentiell ⇒ bedingt
    4          c_k      ← config[k]
    5          params_k ←
                 IF is.null(c_k.parm)
                 THEN {}
                 ELSE c_k.parm( X[i, 1:(k−1)] )  # greift auf frühere Spalten
    6          X[i,k]  ← draw_from( c_k.distr , params_k )
    7  RETURN X
```

*`draw_from`* ist eine generische Hilfsroutine (z. B. Norm, Gamma usw.).

---

### **Script 3: 02\_split.R**

```
FUNCTION train_test_split(X, split_ratio, seed):
    INPUT : Datenmatrix X, split-Ratio r, seed
    OUTPUT: list S = (X_tr, X_te)

    1  SET rng ← seed + 1              # anderer Stream
    2  idx   ← random_permutation(1,…,N)
    3  N_tr  ← ⌊r·N⌋
    4  X_tr  ← X[ idx[1:N_tr],  ]
    5  X_te  ← X[ idx[N_tr+1:N], ]
    6  RETURN S
```

---

### **Script 4: models/ttm\_model.R**  

# 4.1 Helfer für stabilen Log-Raum-Integrations- und Exponential-Umgang

```
FUNCTION log_phi_K(z):
    # log-Dichte der K-dimensionalen Standardnormalen
    RETURN −0.5 * (K * log(2π) + ∥z∥²)

FUNCTION logsumexp(v):
    # stabiler log( Σ exp(v_j) )
    m ← max(v)
    RETURN m + log( Σ exp(v − m) )

FUNCTION log_integrate_exp(f, a, b, n=32):
    # log ∫_a^b exp( f(t) ) dt  via Gauss-Legendre-Quadratur
    # 1. transformiere (a,b) → (−1,1), 2. wandle Summation per logsumexp
    (w, s) ← gauss_legendre_nodes_weights(n)         # Gewichte w_j, Stütz­stellen s_j ∈ (−1,1)
    t      ← 0.5*(b−a)*s + 0.5*(b+a)                 # Rücktransformation
    v      ← log(w) + log(0.5*(b−a)) + f(t)          # additiv, kein expl. exp()
    RETURN logsumexp(v)
```

# 4.2 Transport-Map `S`, Jacobi‐Logdet und Neg-Loglikelihood

```
FUNCTION S(x ; θ, h):
    INPUT  : x=(x₁,…,x_K),     θ = {β_k, α_k}_k,     Grad h
    OUTPUT : z=(z₁,…,z_K)

    FOR k = 1,…,K:
        1  g_k ← Polynomial_deg(h)( x₁,…,x_{k−1} ; β_k )        # rein additiv
        2  P_k(t, x₁:_{k−1}) ← Polynomial_deg(h)( t, x₁,…,x_{k−1} ; α_k )
        3  logI_k ← log_integrate_exp( λ t: P_k(t, x₁:_{k−1}) , 0 , x_k )
        4  z_k ← g_k + exp( logI_k )                             # nur hier exp()
    RETURN z
```

```
FUNCTION logJ(x ; θ, h):
    INPUT  : x, θ, h
    OUTPUT : Σ_{k=1}^{K}  P_k( x_k , x₁:_{k−1} )
```

```
FUNCTION ℓ(θ | x, h):
    # log-Dichte‐Pullback
    z     ← S(x ; θ, h)
    logJx ← logJ(x ; θ, h)
    RETURN log_phi_K(z) + logJx
```

# 4.3 Training
```
FUNCTION 𝔏_train(θ | h):
    RETURN − |X_tr|^{-1} Σ_{x∈X_tr} ℓ(θ | x, h)       # Neg-Loglikelihood‐Mittel
```

```
FUNCTION fit_TTM(X_tr, X_te, H_grid):
    INPUT  : X_tr, X_te, H_grid
    OUTPUT : M_TTM = (θ*, h*, logL_te)

    FOR each h ∈ H_grid:
        1  θ^(0) ← 0
        2  θ̂(h) ← argmin_θ  𝔏_train(θ | h)      # L-BFGS-B
               stopping:  ‖∇𝔏_train‖_∞ < 10^{−6}
        3  logL_te(h) ← − |X_te|^{-1} Σ_{x∈X_te} ℓ( θ̂(h) | x , h )
        4  MESSAGE "h={h}, logL_te={logL_te(h)}"

    5  h* ← argmin_h logL_te(h)
    6  θ* ← θ̂(h*)
    7  RETURN (θ*, h*, logL_te(h*))
```

```
FUNCTION logL_TTM(M_TTM, X):
    INPUT  : (θ*, h*),  X
    OUTPUT : − |X|^{-1} Σ_{x∈X} ℓ( θ*, x , h* )
```
```
FUNCTION sample_TTM(M_TTM, Z):
    INPUT  : Modell M_TTM, Z ~ N(0, I_K)
    OUTPUT : X ~ π̂  via sequentielle Inversion der S_k
```

---

### **Script 5: models/true\_model.R**  („wahrer“ Baseline-MLE)

```
FUNCTION fit_TRUE(X_tr, X_te, config):
    INPUT : X_tr, X_te, config
    OUTPUT: M_TRUE = (Θ̂, logL_te)

    FOR k = 1,…,K:
        distr_k ← config[k].distr
        Θ̂_k    ← argmax_θ_k Σ_{x∈X_tr[,k]} log f_{distr_k}(x | θ_k)   # wie bisher
                 # in der Implementierung ist log f() ohnehin im Log-Raum
    logL_te ← − |X_te|^{-1} Σ_i Σ_k log f_{distr_k}( X_te[i,k] | Θ̂_k )
    RETURN (Θ̂, logL_te)

FUNCTION logL_TRUE(M_TRUE, X):
    RETURN − |X|^{-1} Σ_i Σ_k log f_{distr_k}( X[i,k] | Θ̂_k )


```

*Anmerkung:* Jede Dimension wird unabhängig behandelt, denn „wahrer“ Mechanismus kennt die bedingten Dichten.

---

### **Script 6: 04\_evaluation.R**

```
FUNCTION evaluate_all(X_te, model_list):
    INPUT : Testdaten X_te, Liste model_list = {M_TTM, M_TRUE, …}
    OUTPUT: Tabelle P = (model_id, −logL_te)

    1  INIT empty table P
    2  FOR i = 1,…,|model_list|:
    3      id   ← names(model_list)[i]
    4      M    ← model_list[[i]]
    5      loss ← logL_<id>(M, X_te)            # dynamischer Funktionsaufruf
    6      append_row(P, (id, loss))
    7  sort_by loss ascending
    8  RETURN P
```

`model_specific_logL` ruft intern `logL_TTM`, `logL_TRUE`, oder eine Analogie für künftige Modelle.

---

### **Script 7: 05\_main.R**  (Orchestrator)

```
config <- list(
  list(distr = "norm", parm = NULL),
  list(distr = "exp",  parm = function(d) list(rate = d$X1)),
  list(distr = "beta", parm = function(d) list(shape1 = softplus(d$X2), shape2 = 1)),
  list(distr = "gamma", parm = function(d) list(shape = softplus(d$X3), scale = 1))
)

FUNCTION main():
    1  G        ← setup_global()                    # Script 1
    2  X        ← gen_samples(G)                    # Script 2
    3  S        ← train_test_split(X, G.split_ratio, G.seed)   # Script 3
    4  M_TTM    ← fit_TTM(S.X_tr, S.X_te, G.H_grid) # Script 4
    5  M_TRUE   ← fit_TRUE(S.X_tr, S.X_te, G.config) # Script 5
    6  results  ← evaluate_all(S.X_te, {M_TTM, M_TRUE})   # Script 6
    7  print(results)
    8  # einfache Möglichkeit, weitere Modelle:
       # Quellcode in models/<neues>.R mit
       #   fit_<NAME>()  und  logL_<NAME>()
       # anschließend hier einfach laden & anhängen.
```

---
