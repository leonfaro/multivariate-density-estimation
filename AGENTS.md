# Repository Guidelines

## Project Structure & Module Organization
Core R modules in `code/R`: `loader.R` initialises `globals.R`, data generation, splitting, evaluation, and model backends (TRUE, TRTF, TTM, optional `copula_np`). Experiments in `code/experiments` (e.g. `4d`, `NF/main`, `halfmoon`) call those helpers and expect artefacts in neighbouring `results/` folders. Reporting lives in `report` with knitr-driven `.Rnw`, LaTeX resources, caches, and the `Makefile`; keep generated PDFs and diagnostics inside `report/` and `results/`.

## Build, Test, and Development Commands
- `Rscript code/experiments/4d/main.R` — runs the synthetic pipeline; override sample size via `N_OVERRIDE=200`.
- `TRTF_DATASET=power TRTF_N=1000 Rscript code/experiments/NF/main/main_NF.R` — reproduces normalizing-flow benchmarks and writes logs under `code/experiments/NF/main/results/`.
- `make -C report all` — knits `MSc_Report.Rnw` and compiles the PDF; `make clean` clears stale LaTeX artefacts.
- `Rscript -e "source('code/R/loader.R'); initialize_repo()"` — loads modules and options for interactive R sessions.

## Coding Style & Naming Conventions
Use two-space indents, explicit guards (`stopifnot`), and Roxygen headers for exported helpers. Follow existing names: externally called helpers often use `Camel_Snake` (`Generate_iid_from_config`), while internal utilities stay `snake_case`. Log informative messages instead of silently ignoring missing packages. In the report, reuse macros from `header.sty` and avoid hard-coded formatting.

## Environment & Dependencies
Install dependencies before running experiments: `install.packages(c('tram','trtf','partykit','mlt','dplyr','parallel'))`. Scripts default to single-threaded BLAS; raise `options(trtf.train_cores = 4)` only when coordinating with other users. Dataset CSVs are tracked via Git LFS—run `git lfs pull` after cloning. Optional packages such as `cli` enhance progress output but should stay behind availability checks.

## Testing Guidelines
There is no dedicated test harness; rerun the relevant experiment scripts after model changes and compare `results/*.txt` outputs or console NLL summaries. For report edits, rebuild with `make -C report all` and inspect the LaTeX log for warnings before committing regenerated PDFs. Capture manual verification steps in the commit body when outputs remain unchanged.

## Commit & Pull Request Guidelines
Recent history favours short, imperative subjects in English or German (e.g. `Make experiment scripts locate loader dynamically`) with occasional Conventional-Commit prefixes (`chore: track power CSVs`). Keep the first line under 72 characters, reference affected datasets or seeds, and note whether artefacts were regenerated. Pull requests should summarise experiment outcomes, link issues, attach result snippets or tables, and confirm the report still builds cleanly after LaTeX changes.
