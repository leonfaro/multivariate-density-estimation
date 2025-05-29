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
    INPUT: —
    OUTPUT: list G = (N, config, seed, split_ratio, H_grid, model_ids)

    1  N            ← 500                                # Stichprobenumfang
    2  config       ← USER-EINGABE                       # Liste mit K Elementen
    3  seed         ← 42                                 # Reproduzierbarkeit
    4  split_ratio  ← 0.70                               # Train-Anteil
    5  H_grid       ← {p ∈ ℕ | 1 ≤ p ≤ P_max}            # Hyperparameter für TTM
    6  model_ids    ← {"TTM", "TRUE"}                    # erweiterbar
    7  RETURN G
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
FUNCTION fit_TTM(X_tr, H_grid):
    INPUT : Trainingsdaten X_tr, Hyperparameter-Grid H_grid
    OUTPUT: model M_TTM = (θ*, h*, logL_te)

    1  FOR each h ∈ H_grid:                          # z. B. Polynomiellgrad
    2      θ̂(h) ← argmin_θ  𝔏_train(θ | h)           # s. Block 4 oben
                       subject to θ_{k,j}^f ≥ 0
    3      logL_te(h) ← −  |X_te|^{-1}   Σ_{x ∈ X_te}  ℓ_{θ̂(h)}(x)
    4  h*    ← argmin_h  logL_te(h)
    5  θ*    ← θ̂(h*)
    6  RETURN (θ*, h*, logL_te(h*))

FUNCTION logL_TTM(M_TTM, X):
    INPUT : Modell (θ*, …), Datenmatrix X
    OUTPUT: −log-Likelihood pro Beobachtung
```

*Subschritte*:

* $ℓ_{θ}(x)$ und $𝔏_\text{train}$ exakt wie im vorherigen Post.
* Optimierung via LBFGS-B (weil lineare Ungleichungen).

---

### **Script 5: models/true\_model.R**  („wahrer“ Baseline-MLE)

```
FUNCTION fit_TRUE(X_tr, config):
    INPUT : Trainingsdaten X_tr, vollständige Verteilungskonfiguration
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
    2  FOR each M in model_list:
    3      id ← model_name(M)
    4      loss ← model_specific_logL(M, X_te)
    5      append_row(P, (id, loss))
    6  sort_by loss ascending
    7  RETURN P
```

`model_specific_logL` ruft intern `logL_TTM`, `logL_TRUE`, oder eine Analogie für künftige Modelle.

---

### **Script 7: 05\_main.R**  (Orchestrator)

```
FUNCTION main():
    1  G        ← setup_global()                    # Script 1
    2  X        ← gen_samples(G)                    # Script 2
    3  S        ← train_test_split(X, G.split_ratio, G.seed)   # Script 3
    4  M_TTM    ← fit_TTM (S.X_tr, G.H_grid)        # Script 4
    5  M_TRUE   ← fit_TRUE(S.X_tr, G.config)        # Script 5
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

