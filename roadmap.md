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
    6  H_grid      ← {h ∈ ℕ | 1 ≤ h ≤ P_max}             # Gradbeschränkung
    7  model_ids   ← {"TRUE"}                             # erweiterbar
    8  RETURN G
```

---

### **Script 2: 01\_data\_generation.R**

```
FUNCTION gen_samples(G):
    INPUT : G = (N, config, seed, …)
    OUTPUT: Matrix X ∈ ℝ^{N×K}         # gleiche Reihenfolge wie config

    1  SET rng ← seed
    2  K ← |config|
    3  INIT X ← matrix(NA, N, K)
    4  FOR i = 1,…,N:
    5      FOR k = 1,…,K:                        # sequentiell ⇒ bedingt
    6          c_k      ← config[k]
    7          params_k ←
                 IF is.null(c_k.parm)
                 THEN {}
                 ELSE c_k.parm( X[i, 1:(k−1)] )
    8          IF c_k.distr = "gamma" ∧ {shape1, shape2} ⊂ names(params_k):
                 params_k ← {shape = params_k.shape1,
                              scale = params_k.shape2}
    9          FOR p ∈ params_k:
                 IF ¬finite(p) ∨ p ≤ 0: p ← 1e−3
   10          fun ← get("r" ++ c_k.distr)
   11          X[i,k]  ← fun(1, params_k)
   12  SET colnames(X) ← {"X1",…,"XK"}
   13  RETURN X
```

Die Generierung nutzt die jeweiligen Basisfunktionen `r<distr>` aus R und
ersetzt ungültige Parameter (≤0 oder nicht endlich) durch `1e−3`.

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


### **Script 4 (neu): models/**`triangular_transport_map.R`

> **Ziel** Trianguläre **Transport Map**
>
> $$
> S_\theta(x)=\bigl(S_1(x_1),\,S_2(x_{1:2}),\dots,S_K(x_{1:K})\bigr)^\top,\qquad  
> S_k = g_k(x_{1:k-1};\beta_k) \;+\; f_k(x_k;\alpha_k)
> $$
>
> *Alle* numerisch kritischen Schritte laufen **im Log-Raum**.

---

#### 4.1 Hilfsfunktionen (Log-Raum)

```
FUNCTION log_phi_K(z):
    # log-Dichte der K-dimensionalen Standardnormalen
    RETURN −0.5 * (K * log(2π) + ||z||²)

FUNCTION logsumexp(v):
    m ← max(v)
    RETURN m + log( Σ exp(v − m) )

FUNCTION poly_deg(h)(t ; γ):
    # schlichtes Polynom   Σ_{j=0}^h γ_j t^j
    RETURN Σ_{j=0}^h γ_j * t^j

FUNCTION f_k (x_k ; α_k, h):
    # monotone 1-D Komponente   f_k = ∫₀^{x_k} exp( poly_deg(h)(t;α_k) ) dt
    # (=> immer streng steigend, kein Overflow dank Log-Raum)
    logI ← log_integrate_exp( λ t: poly_deg(h)(t;α_k), 0, x_k )
    RETURN exp( logI )

FUNCTION log_fprime_k (x_k ; α_k, h):
    # log f'_k(x_k) = poly_deg(h)( x_k ; α_k )
    RETURN poly_deg(h)( x_k ; α_k )
```

*(`log_integrate_exp` wie im alten Skript; wird nur für $f_k$ gebraucht.)*

---

#### 4.2 Map $S_\theta$, Jacobian-Logdet und Log-Likelihood

```
FUNCTION S(x ; θ=(β,α), h):
    INPUT  : x ∈ ℝ^K
    OUTPUT : z ∈ ℝ^K

    FOR k = 1,…,K:
        g_k ← poly_deg(h)( x[1:(k−1)] ; β_k )        # nur frühere Coordinates
        f_k_val ← f_k( x_k ; α_k , h )                # 1-D, streng steigend
        z_k ← g_k + f_k_val
    RETURN z
```

```
FUNCTION logJ(x ; θ, h):
    # log |det ∇S| = Σ_k log f'_k(x_k)
    logdet ← 0
    FOR k = 1,…,K:
        logdet ← logdet + log_fprime_k( x_k ; α_k , h )
    RETURN logdet
```

```
FUNCTION ℓ(θ | x, h):
    z      ← S(x ; θ , h)
    RETURN log_phi_K(z) + logJ(x ; θ , h)      # sample-LogLikelihood
```

---

#### 4.3 Training (empirische KL / –log L)

```
FUNCTION 𝔏_train(θ | h, X_tr):
    RETURN − (1/|X_tr|) * Σ_{x ∈ X_tr} ℓ(θ | x , h)
```

```
FUNCTION fit_SEPAR(X_tr, X_te, H_grid):
    INPUT  : Trainings-/Testdaten, Polynomgrade
    OUTPUT : M_SEP = (θ*, h*, logL_te*)

    FOR h ∈ H_grid:
        θ⁰      ← 0
        θ̂(h)   ← argmin_θ  𝔏_train(θ | h, X_tr)     # L-BFGS-B
                  stop wenn ||∇𝔏_train||_∞ < 1e−6
        logL_te(h) ← − (1/|X_te|) Σ_{x ∈ X_te} ℓ( θ̂(h) | x , h )
        MESSAGE "h={h}, logL_te={logL_te(h)}"

    h*  ← argmin_h logL_te(h)
    θ*  ← θ̂(h*)
    RETURN (θ*, h*, logL_te(h*))
```

```
FUNCTION logL_SEPAR(M_SEP, X):
    (θ*, h*) ← M_SEP
    RETURN − (1/|X|) Σ_{x ∈ X} ℓ( θ*, x , h* )
```

---

#### 4.4 Sampling und Dichteschätzung

```
FUNCTION sample_SEPAR(M_SEP, Z):
    INPUT  : Z ~ N(0,I_K)
    OUTPUT : X  via sequentielle Inversion

    (θ*, h*) ← M_SEP
    FOR k = 1,…,K:
        g_k ← poly_deg(h*)( X[1:(k−1)] ; β*_k )
        # löse   f_k(x_k) = Z_k − g_k   nach x_k  (1-D root-finding, z.B. Newton)
        x_k ← invert_f_k( Z_k − g_k ; α*_k , h* )
        X_k ← x_k
    RETURN X
```

```
FUNCTION density_SEPAR(M_SEP, x):
    z      ← S(x ; θ*, h*)
    logdet ← logJ(x ; θ*, h*)
    RETURN exp( log_phi_K(z) + logdet )        # π̂(x)
```

`invert_f_k` benutzt z.B. Newton–Raphson mit Startwert $x_k^{(0)} = Z_k$; nur $f_k$ und $f'_k$ nötig.



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
    INPUT : Testdaten X_te, Liste model_list = {M_TRUE, …}
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

`model_specific_logL` ruft intern `logL_TRUE` oder eine Analogie für künftige Modelle.

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
    4  M_TRUE   ← fit_TRUE(S.X_tr, S.X_te, G.config) # Script 5
    5  results  ← evaluate_all(S.X_te, {M_TRUE})      # Script 6
    7  print(results)
    8  # einfache Möglichkeit, weitere Modelle:
       # Quellcode in models/<neues>.R mit
       #   fit_<NAME>()  und  logL_<NAME>()
       # anschließend hier einfach laden & anhängen.
```

---
