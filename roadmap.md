#Roadmap

• **Block A (K = 3 Festlegung)**
 ◦ `config3` unverändert (norm, exp, gamma)
  ◦ Notation strikt wie *notation.txt* (Suffix `_pi`, `_eta`, etc.) 

• **Block B (Daten-Generating-Prozess)**
  ◦ Funktionen `eta_sample()`, `S_inv()` beibehalten
  ◦ Generiere zwei unabhängige Samples
   ▪ Train $N_{\text{train}} = 100$
   ▪ Test $N_{\text{test}} = 100$
  ◦ Speichere `X_pi`, `U_eta`, `Z_eta`, `logd` separat als CSV für Debug

• **Block C (Parametrische Baseline ≙ “ohne Forest”)**
  ◦ MLE-Schätzung jeder Konditional­verteilung $f_{d_k}\bigl(x_k\mid x_{<k}\bigr)$ via `optim` schon umgesetzt 
  ◦ Berechne **Out-of-Sample**‐Loglikelihood
   $\ell_{\text{param}}^{\text{test}} = \sum_{i,k}\log f_{d_k}(x_{ik}\mid x_{i,<k};\hat\theta_k)$
  ◦ Vergleiche mit wahrer Loglikelihood
   $\Delta\ell_k = \ell_{\text{true},k} - \ell_{\text{param},k}$ (table already printed)

• **Block D (Transformation Forests)**
  ◦ Für jedes $k=2,\dots,3$:
   ▪ Formular `y_k ~ x_1 + … + x_{k-1}`
   ▪ Basismodell `mlt()` linearer shift-scale; danach `traforest()` 
   ▪ Parameter: `ntree = 200`, `mtry = ⌈(k-1)/3⌉`, `minbucket ≥ 20`
   ▪ Keine Informationen aus `config` (Forest sieht **nur Daten**)
  ◦ $k=1$: einfaches marginales `mlt()` ohne Forest
  ◦ Vorhersage auf Testdaten:
   `ld_k <- predict(obj_k, newdata = X_test, type = "logdensity")`
◦ Baue Matrix `LD_hat` ($N_{\text{test}}\times3$)

• **Block E (Joint Likelihood-Evaluation)**
  ◦ Addiere Beiträge:
   $\hat\ell_i = -\tfrac12‖Z_{i,\eta}‖^2 -2\log(2\pi) + \sum_{k} ld_{ik}$
  ◦ Vergleiche gegen wahres $\ell_i$ (aus DGP)
   ▪ `all.equal(sum(LD_hat), sum(ll_true), tol = 1e-1)`
  ◦ Streudiagramm pro Dimension:
    x-Achse $\ell_{ik}^{\text{true}}$, y-Achse $\ell_{ik}^{\text{forest}}$
    Ideal ≈ 45°-Linie

• **Diagnostik-Plots**
  ◦ Histogramme `Δℓ_k`
  ◦ QQ-Plot $Z_{\eta}$ gegen $\mathcal N(0,1)$
  ◦ Heatmap empirische Copula der Residuen

• **Dateistruktur**
  ◦ `run_all.R` ruft `00_setup.R` bis `05_joint_evaluation.R` auf
  ◦ Output-Verzeichnis `results/` für CSV und PDF-Plots

• **To-Do-Checkliste**
  ▢ Code refaktorieren: Funktionen `fit_param()`, `fit_forest()`
  ▢ Seed setzen `set.seed(…)` für Reproduzierbarkeit
  ▢ Variable-Order $x_1,x_2,x_3$ **nicht ändern** (Professor­vorgabe) 
  ▢ Sample-Größen im Skript als Parameter
  ▢ Laufzeit-Kommentar: Forests \~ minuten­lang, akzeptabel

• **Nächster Schritt nach Proof-of-Concept**
  ◦ Skalieren auf $K=25$ (nur Rechenzeit)
  ◦ Evtl. Sparsity-Flags $f_3$ aktivieren für bedingte Unabhängigkeit 

• **Lieferobjekte fürs Meeting**
  ◦ Scripts, CSV-Daten, Plots
  ◦ Tabelle “true vs forest” inkl. $\Delta\ell$
  ◦ Kurze Slide mit Bullet-Summary der Ergebnisse

✔️ Erfüllt alle vom Prof geforderten Punkte; fehlende Aufgaben oben abhaken.
