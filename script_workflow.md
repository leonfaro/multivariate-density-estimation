# Script Workflow

The following diagram details how the analysis scripts `00_setup.R` through `07_dvine_copula.R` interact when executed via `run_all.R`.  Boxes refer to scripts and arrow labels indicate which objects are passed.  All names follow [Notation.md](Notation.md).

```mermaid
flowchart TD
    Setup["00_setup.R\nclip(), safe_*(), config"] -->|"source"| Utils["01_transport_utils.R\nS_inv(), pi_sample(), pdf_k" ]
    Utils -->|"S_inv"| Data["02_generate_data.R\nX_{\pi,train}, X_{\pi,test}, Z_{\eta,train}, Z_{\eta,test}, logd"]
    Data -->|"X_{\pi,train}"| Baseline["03_param_baseline.R\nparam_est, ll_delta_df_test"]
    Data -->|"X_{\pi,train}"| Forest["04_forest_models.R\nmodel, LD_hat"]
    Data -->|"X_{\pi,train}"| Kernel["06_kernel_smoothing.R\nks_model, KS_hat"]
    Forest -->|"U_{\hat{\pi}}"| Copula["07_dvine_copula.R\nloglik_dvine"]
    Data -->|"ll_test"| Eval["05_joint_evaluation.R"]
    Baseline --> Eval
    Forest --> Eval
    Kernel --> Eval
    Copula --> Eval
```

Scripts are sourced in this order. `00_setup.R` defines numerical safeguards and the distribution list `config`. `01_transport_utils.R` constructs the triangular map via `S_inv()` and helper densities. `02_generate_data.R` draws samples from the target distribution, computes the log-determinant and stores CSV files. In `03_param_baseline.R` parametric log-likelihoods are fitted and compared to the truth, producing `ll_delta_df_test`. `04_forest_models.R` trains transformation forests and yields component log-density estimates `LD_hat`. `06_kernel_smoothing.R` implements a sequential kernel estimator `KS_hat`. Using the forest CDFs, `07_dvine_copula.R` fits a D-vine copula and outputs the joint log-likelihood `loglik_dvine`. Finally `05_joint_evaluation.R` gathers all log-likelihood contributions, writes diagnostics and aggregates them in `results/evaluation_summary.csv`.
