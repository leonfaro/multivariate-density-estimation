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
- [02_split.R](02_split.R) — Reproducible 80/10/10 data partitioning.
- [04_evaluation.R](04_evaluation.R) — Model evaluation helpers and table summaries.
- [demo_script_do_not_change.R](demo_script_do_not_change.R) — Demonstration of conditional transformation forests.
- [main.R](main.R) — Orchestrates full pipeline with global seed.
- [miniboone.R](miniboone.R) — Preprocesses MiniBooNE data into splits.
- [replicate_code.R](replicate_code.R) — Concatenates sourced scripts for reproducibility.

* * *

## Notes
- Developed within a master's thesis project.
- Reproducibility: global RNG seed configured in [main.R](main.R).
