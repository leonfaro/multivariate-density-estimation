# Thesis in Multivariate Conditional Density Estimation: Comparison with Transformation Forest Models, Copulas and Normalizing Flows

This repository collects the R notebooks and scripts used in the associated thesis.  All experiments are run directly in R to study multivariate conditional densities with mathematical rigour.

**Important:** Read `Theory.md` and 'roadmap.md' before starting any task. Always analyzing outputs for numerical numerical validation. You must always stay within the algorithmic workflow of `Theory.md`. If that is not possible, tell where there might be conflicts.

## Scientific approach
This repository is a research notebook rather than a software project. All scripts are executed directly in R to study problems in mathematical statistics and probability theory. The focus is on precise and reproducible experiments with **Mathematical Rigour**, not on standard software engineering workflows or packaging. The choice of R reflects the statistical setting; Python is intentionally avoided. The triangular transport methodology is described in `Theory.md`. Our \(S_k(x_1,\ldots,x_k)\) functions are simply transformation forests for the regression \(x_k \sim x_1,\ldots,x_{k-1}\). Monotonicity in \(x_k\) is therefore automatically ensured.

# Mathematical Writing and Statistical Algorithms Guidelines

- Read `Theory.md` before running any scripts.
- Strictly adhear to notations outlined in `Theory.md`.
- Adjust `ALGORITHM_SPEC.md` after every task if necessary, so that `ALGORITHM_SPEC.md` always reflects the current state of the codebase.
- Always analyzing outputs for numerical numerical validation.
- You must always stay within the algorithmic workflow of `Theory.md`. If that is not possible, tell where there might be conflicts.
- Keep Code as simple as possible. Goal is a Proof-of-Concept.
- Check numeric results for unexpected NaN or NA values and inspect the full output for numerical plausibility.
- Before sending PR-Message, check if change changes lead to conflicts with algorithmic workflow of `Theory.md`.
- Pure documentation or comment changes do not require running tests or lint checks.
- Jeder Commit-Kommentar muss auf Deutsch verfasst sein.
- Bei Bugfixes eine kurze Zusammenfassung im Titel, z.B. "Clamp probabilities in qtf_k to avoid beta support violations." schreiben.
- Anschließend einen Body mit Erläuterung der Änderung und eventuellen Nebenwirkungen anfügen.
  
### Testing, Debugging, Reproducibility

* Write unit tests with **testthat**; place them in `tests/testthat/`.
* Use `set.seed()` at every experiment start and record the seed in your report.
* Log intermediate metrics with `futile.logger`, `logger`, or simple `message()`.
* Track experiments in a plain CSV or YAML file; include hyper-parameters.
