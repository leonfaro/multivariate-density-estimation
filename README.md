# Multivariate Density Estimation via Triangular Transport Maps
All densities are evaluated in log-space and strictly positive parameters use softplus transforms, as detailed in [ALGORITHM_SPEC.md](ALGORITHM_SPEC.md).

* * *

## Folder Structure
- [data/](data/) — Raw and generated datasets for experiments.
- [docs/](docs/) — Thesis material and supplementary explanations.
- [models/](models/) — Estimation models and supporting code.
- [results/](results/) — Output tables and trained objects.
- [tests/](tests/) — Test suite using testthat.

* * *

## Top-level Scripts
- [00_globals.R](00_globals.R) — Global settings and helper utilities.
- [01_data_generation.R](01_data_generation.R) — Conditional sampling via triangular transport mapping.
- [02_split.R](02_split.R) — Reproducible 80/20 Train/Test split (Train/Test only).
- [04_evaluation.R](04_evaluation.R) — Model evaluation helpers and table summaries.
- [demo_script_do_not_change.R](demo_script_do_not_change.R) — Demonstration of conditional transformation forests.
- [main.R](main.R) — Orchestrates full pipeline with global seed.
- [main_moon.R](main_moon.R) — Run Half‑Moon (2D) experiments and plotting pipeline.
- [miniboone.R](miniboone.R) — Preprocesses MiniBooNE data into splits.
- [replicate_code.R](replicate_code.R) — Concatenates sourced scripts for reproducibility.

* * *

## Notes
- Developed within a master's thesis project.
- Reproducibility: global RNG seed configured in [main.R](main.R).
- Cross-term (Eq. 22) policy: trainer uses train/test only (no validation/CV); ridge weights `lambda_non`/`lambda_mon` and quadrature nodes `Q` can be configured via R options (`mde.ctm.*`) or env vars (`MDE_CTM_*`).
- Half‑Moon usage: use `main_moon.R` for Half‑Moon runs; `main.R` handles the Config‑4D pipeline only.
