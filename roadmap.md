
### **Block 0 (Setup & Reproducibility)**

* **00_setup.R**
  * Call `set.seed(SEED)` with `SEED <- 24`.
  * Define `EPS`, `SAFE_CLIP`, `logsumexp()`, etc. as in the existing
    numerically safe utility layer.



### Block A (K = 3 setting)

* **config3** stays as (norm, exp, pois), but:
  * **Instead of `lambda = d$X2`** use  
    `lambda = clip(exp(d$X2), EPS, 1e6)`; this guarantees positivity and
    matches the clipping scheme of the other distributions.
  * Store `K <- length(config)` as a global so that all subsequent
    scripts *read* `K` rather than redefine it.



### Block B (Data generating process)

* **02_generate_data.R**  
  Create `samp_train` and `samp_test` exactly as now, but
  * **write the seed to the CSV file** with  
    `attr(csv_object, "seed") <- SEED` so runs are reproducible later.
  * **also store `det_J` and the true log-likelihood columns** in the CSVs
    to simplify unit tests.
* Helper functions `eta_sample()`, `S_inv()`, etc. remain unchanged.



### Block C (Parametric baseline)

* **03_param_baseline.R**
  * **Generic approach:** implement `nll_fun_from_cfg(k, cfg)` and
    `eval_ll_from_cfg()`.  
    This removes the manual switch block and prevents the
    Gamma/Poisson mix-up that was in the original code.
  * After every fit run  
    `stopifnot(abs(delta_l_k) < 1e2)`  
    and add matching `testthat` files in `tests/testthat/`.
  * **Out-of-sample delta_l** is computed in `ll_delta_df_test`, but the
    column is now strictly called `delta_ll_param` as in `run3.R`.



### Block D (Evaluation & reporting) - new

* **run3.R**
* For a strictly monotone triangular map `S`, every diagonal derivative  
  `partial_{x_k} S_k` must be `> 0`; this implies `det(J) > 0`.



### Further rigor details

1. **Notation:** Every variable gets the suffix `_pi` (target domain) or
   `_eta` (reference domain) exactly as in Table 1 of the tutorial.
2. **Log determinant:** Make sure `det_J(logd)` is added **before** the
   bias term inside the log-likelihood call.

