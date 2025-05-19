# Script Workflow

The following Mermaid diagram summarises how the main scripts interact. Each box denotes an R script. Arrows indicate which objects are passed forward.

```mermaid
flowchart TD
    A["00_setup.R"] -->|"config, K"| B["01_transport_utils.R"]
    B -->|"pi_sample(), S_inv()"| C["02_generate_data.R"]
    C -->|"X_pi_train, X_pi_test\nU_eta_train, U_eta_test\nZ_eta_train, Z_eta_test"| D["03_param_baseline.R"]
    C -->|"X_pi_train, X_pi_test"| E["04_forest_models.R"]
    D -->|"ll_delta_df_test"| F["05_joint_evaluation.R"]
    E -->|"LD_hat"| F
    C -->|"ll_test, Z_eta_test, true_ll_mat_test"| F
```

This workflow mirrors `run_all.R`.  The setup script defines global constants, followed by utilities for the triangular transport map. Data generation draws from the reference distribution and transforms it into the target space. The parametric baseline and forest models fit conditional densities using the generated training data. Finally, the evaluation script compares the estimated densities against the true likelihood.
