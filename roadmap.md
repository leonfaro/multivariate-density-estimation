---

### **Block 0 (Setup & Reproduzierbarkeit) – neu**

* **00\_setup.R**

  * Setze `set.seed(SEED)` mit SEED = 24
  * Definiere `EPS`, `SAFE_CLIP`, `logsumexp()`, usw. wie bereits vorhanden („numerically safe utility layer“).
---

### Block A (K = 3 Festlegung)

* **config3** bleibt (norm, exp, pois), doch:
  * **Statt `lambda = d$X2`** verwende `lambda = clip(exp(d$X2), EPS, 1e6)`; so ist Positivität garantiert und identisch zum Clip-Schema anderer Verteilungen.
  * Hinterlege `K <- length(config)` global; alle subsequent scripts lesen `K`, nicht umgekehrt.
---

### Block B (Daten-Generating-Prozess)

* **02\_generate\_data.R**
   Erzeuge `samp_train`, `samp_test` exakt wie jetzt, aber
    * **speichere den Seed in der CSV-Datei** (`attr(..., "seed") <- SEED`), damit später rekonstruierbar.
    * **speichere auch `det_J` und die wahren Log-Likelihood‐Spalten** in den CSVs; das erleichtert unit tests.
  * Die Hilfsfunktionen `eta_sample()`, `S_inv()` usw. bleiben unverändert .

---

### Block C (Parametrische Baseline)

* **03\_param\_baseline.R**

  * **Generischer Ansatz:** baue `nll_fun_from_cfg(k, cfg)` und `eval_ll_from_cfg()` – so verschwindet der manuelle Switch-Block und die Gefahr des Gamma/Poisson-Vertauschers (Fehler aus dem ursprünglichen Code).
  * Führe nach jedem Fit `stopifnot(abs(Δℓ_k) < 1e2)` aus und addiere `testthat`-Files in `tests/testthat/`.
  * **Out-of-sample Δℓ** wird in `ll_delta_df_test` berechnet, aber die Spalte heißt nun strikt `delta_ll_param` wie in run3.R.
---

### Block D (Evaluierung & Reporting) – neu

* **run3.R**
* Für ein strikt monoton-trianguläres S müssen alle Diagonal-Ableitungen ∂<sub>x k</sub>S<sub>k</sub>>0 gelten; damit ist auch det (J)>0 .
    ---

### Weitere Rigorositäts-Details

1. **Notation**: Jede Variable erhält `_pi` (target) oder `_eta` (reference) Suffix gemäß Tabelle 1 des Tutorials .
2. **Log-Determinant**: Stelle sicher, dass `det_J(logd)` summiert **vor** dem Bias-Term im Log-Likelihood-Call auftaucht .
