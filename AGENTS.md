# Mathematical Writing and Statistical Algorithms Guidelines

- Read `Theory.md` before running any scripts.
- Always analyzing outputs for numerical numerical validation.
- You must always stay within the algorithmic workflow of `Theory.md`. If that is not possible, tell where there might be conflicts.
- Keep Code as simple as possible. Goal is a Proof-of-Concept.
- Check numeric results for unexpected NaN or NA values and inspect the full output for numerical plausibility.
- Before sending PR-Message, check if change changes lead to conflicts with algorithmic workflow of `Theory.md`.
- Pure documentation or comment changes do not require running tests or lint checks.
  
### Testing, Debugging, Reproducibility

* Write unit tests with **testthat**; place them in `tests/testthat/`.
* Use `set.seed()` at every experiment start and record the seed in your report.
* Log intermediate metrics with `futile.logger`, `logger`, or simple `message()`.
* Track experiments in a plain CSV or YAML file; include hyper-parameters.
