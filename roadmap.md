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

### **Script 4: models/ttm\_model.R**  (Referenz-Modell aus dem Paper)


#  Hilfsroutinen 

FUNCTION log_phi_K(z):
    # Log-Dichte der K-dimensionalen Standardnormalverteilung
    RETURN −0.5 * ( K * log(2π) + ||z||² )

FUNCTION logsumexp(v):
    # stabil: log( Σ_j exp(v_j) )
    m ← max(v)
    RETURN m + log( Σ exp(v − m) )

FUNCTION log_integrate_exp(f, a, b, n = 32):
    # log ∫_a^b exp(f(t)) dt  mittels Gauss-Legendre-Quadratur (n Knoten)
    (w, s) ← gauss_legendre_nodes_weights(n)          # Gewichte w_j, Stützstellen s_j∈(−1,1)
    t      ← 0.5*(b−a) * s + 0.5*(b+a)               # Rücktransformation
    v      ← log(w) + log(0.5*(b−a)) + f(t)          # alles additiv
    RETURN logsumexp(v)

#  Transport-Map und Logdet 

FUNCTION S(x ; θ, h):
    INPUT  : x = (x₁,…,x_K),  θ = {β_k, α_k}_k,  Grad h
    OUTPUT : z = (z₁,…,z_K)

    FOR k = 1,…,K:
        1  g_k ← Polynomial_deg(h)( x₁,…,x_{k−1} ; β_k )
        2  P_k(t, x₁:_{k−1}) ← Polynomial_deg(h)( t, x₁,…,x_{k−1} ; α_k )
        3  logI_k ← log_integrate_exp( λ t: P_k(t, x₁:_{k−1}) , 0 , x_k )
        4  z_k ← g_k + exp( logI_k )                  # einziger exp-Schritt
    RETURN z

FUNCTION logJ(x ; θ, h):
    # Log-Jacobi-Determinante
    RETURN Σ_{k=1}^{K} P_k( x_k , x₁:_{k−1} )

#  (Negativ-)Log-Likelihood 

FUNCTION ℓ(θ | x, h):
    z     ← S(x ; θ, h)
    logJx ← logJ(x ; θ, h)
    RETURN log_phi_K(z) + logJx                      # log π̂(x)

FUNCTION 𝔏_train(θ | h):
    RETURN − |X_tr|^{-1} Σ_{x∈X_tr} ℓ(θ | x, h)      # zu minimieren

#  Training / Hyperparameter-Sweep 

FUNCTION fit_TTM(X_tr, X_te, H_grid):
    INPUT  : Trainings-/Test­daten, H_grid
    OUTPUT : M_TTM = (θ*, h*, logL_te)

    FOR each h ∈ H_grid:
        1  θ⁰ ← 0                                   # alle Koeffizienten = 0
        2  θ̂(h) ← argmin_θ  𝔏_train(θ | h)          # L-BFGS-B
               stopping: ‖∇𝔏_train‖_∞ < 1e−6
        3  logL_te(h) ← −|X_te|^{-1} Σ_{x∈X_te} ℓ(θ̂(h) | x, h)
        4  MESSAGE "h={h}, logL_te={logL_te(h)}"

    5  h* ← argmin_h logL_te(h)
    6  θ* ← θ̂(h*)
    7  RETURN (θ*, h*, logL_te(h*))

#  Auswertung & optionales Sampling 

FUNCTION logL_TTM(M_TTM, X):
    INPUT  : (θ*, h*),  X
    OUTPUT : −|X|^{-1} Σ_{x∈X} ℓ(θ*, x, h*)

FUNCTION sample_TTM(M_TTM, Z):
    # optional für Diagnosen
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

## **Wie man die Flexibilität nutzt**

1. **Dimensionalität ändern**
   *Passe nur* `config` in `00_globals.R` an.
   Das Sampling in `01_data_generation.R` sowie alle Schleifen über $k=1,\dots,K$ adaptieren automatisch.

2. **Andere Bedingungsstrukturen**
   Ersetze in `config[[k]]$parm` die Funktion, die aus den bereits generierten Spalten Parameter ableitet.
   Kein weiterer Code muss geändert werden.

3. **Neues Modell hinzufügen**

   * Datei `models/<neues_modell>.R` anlegen
   * Zwei Funktionen implementieren: `fit_<NAME>()`, `logL_<NAME>()`
   * In `00_globals.R` den Namen ins Set `model_ids` und in `05_main.R` nach dem Fitting in `evaluate_all` übergeben.

4. **Hyperparameter-Sweep verändern**
   Ausschließlich `H_grid` in `00_globals.R` anpassen.

Damit bleibt das Gesamtsystem **skript-modular** und **konfigurationsgetrieben**: sämtliche Pipeline-Änderungen erfolgen, ohne andere Dateien „anzufassen“ oder internen Code umschreiben zu müssen.

