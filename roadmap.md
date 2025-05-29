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

```
FUNCTION fit_TTM(X_tr, X_te, H_grid):
    INPUT  : Trainings­daten X_tr, Test­daten X_te, Hyperparameter-Menge H_grid
    OUTPUT : model M_TTM = (θ*, h*, logL_te)

    # -------- Hilfsdefinitionen --------
    DEFINE  S(x ; θ, h)         # Triangular Transport-Map
            logJ(x ; θ, h)      # Log-Jacobi-Determinante
            ℓ(θ | x, h)         = log φ_K ( S(x; θ,h) ) + logJ(x; θ,h)
            𝔏_train(θ | h)      = − |X_tr|^{-1} Σ_{x∈X_tr} ℓ(θ | x, h)

    FOR each h ∈ H_grid:
        1  INIT θ^(0) ← 0                      # alle Koeffizienten = 0
        2  θ̂(h) ← argmin_θ  𝔏_train(θ | h)    # L-BFGS-B  ohne Box-Constraints
                    stopping rule:  ‖∇𝔏_train‖_∞ < 10^{−6}
        3  logL_te(h) ← − |X_te|^{-1} Σ_{x∈X_te} ℓ( θ̂(h) | x , h )
        4  MESSAGE "h={h}, logL_te={logL_te(h)}"

    5  h* ← argmin_h  logL_te(h)
    6  θ* ← θ̂(h*)
    7  RETURN (θ*, h*, logL_te(h*))

FUNCTION S(x ; θ, h):
    INPUT  : Beobachtung x=(x₁,…,x_K), Parameter θ=(θ₁,…,θ_K), Polynomgrad h
    OUTPUT : z = (z₁,…,z_K)

    FOR k = 1,…,K:
        1  g_k ← Polynom_{Grad=h}( x₁,…,x_{k−1} ; β_k )
        2  P_k(t, x₁:_{k−1}) ← Polynom_{Grad=h}( t, x₁,…,x_{k−1} ; α_{k} )
             # P_k linear in den Koeffizienten, keine t-Konstante
        3  z_k ← g_k  +  ∫_{0}^{x_k}  exp( P_k( t, x₁:_{k−1} ) ) dt
    RETURN z

FUNCTION logJ(x ; θ, h):
    INPUT  : x, θ, h
    OUTPUT : Σ_{k=1}^{K}  P_k( x_k , x₁:_{k−1} )

FUNCTION logL_TTM(M_TTM, X):
    INPUT  : Modell M_TTM = (θ*, h*, …), Datenmatrix X
    OUTPUT : − |X|^{-1} Σ_{x∈X} ℓ( θ*, x , h* )

FUNCTION sample_TTM(M_TTM, Z):
    # optional für Diagnosen
    INPUT  : Modell M_TTM, Z ~ N(0, I_K)
    OUTPUT : X ~ π̂  via sequentielle Inversion der S_k
```


---

### **Script 5: models/true\_model.R**  („wahrer“ Baseline-MLE)

```
FUNCTION fit_TRUE(X_tr, X_te, config):
    INPUT : Trainingsdaten X_tr, Testdaten X_te, vollständige Verteilungskonfiguration
    OUTPUT: model M_TRUE = (Θ̂, logL_te)

    1  FOR k = 1,…,K:
    2      distr_k ← config[k].distr
    3      # univariate MLE für jeden k getrennt
    4      Θ̂_k   ← argmax_θ_k   Σ_{x∈X_tr[,k]}  log f_{distr_k}(x | θ_k)
                       via optim()
    5  logL_te ← − |X_te|^{-1} Σ_{i} Σ_{k} log f_{distr_k}(X_te[i,k] | Θ̂_k)
    6  RETURN (Θ̂, logL_te)

FUNCTION logL_TRUE(M_TRUE, X):
    INPUT : (Θ̂,…), Datenmatrix X
    OUTPUT: −log-Likelihood pro Beobachtung
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

