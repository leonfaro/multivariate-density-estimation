# Ablaufdiagramm der Skripte und Modelle

Das folgende Mermaid-Diagramm veranschaulicht die Interaktion der Skripte
`01_data_generation.R`, `02_split.R`, `04_evaluation.R` und `05_main.R` sowie der
vier Modellimplementierungen **TRTF**, **TrueModel**, **KS-Modell** und
**TTM-Modell**. Alle Berechnungen der Log-Likelihood finden konsequent im
Log-Raum statt, wie in `README.md` beschrieben. Parameter, die strikt positiv
sein m\u00fcssen, werden \u00fcber `softplus()` transformiert.

```mermaid
graph TD
    subgraph Script01["01_data_generation.R"]
        A1["gen_samples(G)"] --> A2["Generate_iid_from_config"]
        A2 -->|"nutzt config"| A3["r<distr>-Aufrufe"]
        A3 --> A4["X (N\u00d7K)"]
    end

    subgraph Script02["02_split.R"]
        B1["train_test_split(X, split_ratio, seed)"] -->|"Zufallspermutation"| B2["X_tr"]
        B1 --> B3["X_te"]
    end

    subgraph TrueModel["models/true_model.R"]
        C1["fit_TRUE(X_tr, X_te, config)"] --> C2["neg_loglik_uni (optim)"]
        C2 --> C3["theta_list"]
        C1 --> C4["logL_TRUE(X_te)"]
    end

    subgraph TRTFModel["models/trtf_model.R"]
        D1["fit_TRTF(X_tr, X_te, config)"] --> D2["mytrtf(data)"]
        D2 --> D3["traforest" ]
        D1 --> D4["logL_TRTF(X_te)"]
    end

    subgraph KSModel["models/ks_model.R"]
        E1["fit_KS(X_tr, X_te, config)"] --> E2["bw.nrd0 pro Spalte"]
        E1 --> E3["logL_KS(X_te)"]
    end

    subgraph TTMModel["models/ttm_model.R"]
        F1["fit_TTM(X_tr, X_te)"] --> F2["S_forward / R_inverse"]
        F2 --> F3["LogDet_Jacobian"]
        F1 --> F4["logL_TTM(X_te)"]
    end

    subgraph Script04["04_evaluation.R"]
        G1["evaluate_all(X_te, model_list)"] -->|"ruft logL_<id> auf"| G2["Tabelle"]
    end

    subgraph Script05["05_main.R"]
        H1["main()"] -->|"erstellt G"| Gsetup["G = list(N, config, seed, split_ratio)"]
        H1 -->|"ruft"| A1
        H1 -->|"ruft"| B1
        H1 -->|"fit"| C1
        H1 -->|"fit"| D1
        H1 -->|"fit"| E1
        H1 -->|"fit"| F1
        H1 -->|"evaluiert"| G1
        G2 --> H2["Ausgabe / Plots"]
    end

    A4 --> B1
    B2 --> C1
    B2 --> D1
    B2 --> E1
    B2 --> F1
    B3 --> G1
```
```
