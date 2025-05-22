# PhD Thesis in Multivariate Conditional Density Estimation: Comparison with Transformation Forest Models, Copulas and Normalizing Flows

This repository collects the R notebooks and scripts used in the associated thesis.  All experiments are run directly in R to study multivariate conditional densities with mathematical rigour.

**Important:** Read `AGENTS.md` and `Theory.md` before running any scripts. Always analyzing outputs for numerical numerical validation. You must always stay within the algorithmic workflow of `Theory.md`. If that is not possible, tell where there might be conflicts.

## Scientific approach
This repository is a research notebook rather than a software project. All scripts are executed directly in R to study
problems in mathematical statistics and probability theory. The focus is on
precise and reproducible experiments with **Mathematical Rigour**, not on
standard software engineering workflows or packaging. The choice of R reflects
the statistical setting; Python is intentionally avoided.

The triangular transport methodology is described in `Theory.md`. Our \(S_k(x_1,\ldots,x_k)\) functions are simply transformation forests for
the regression \(x_k \sim x_1,\ldots,x_{k-1}\). Monotonicity in \(x_k\) is therefore automatically ensured.
