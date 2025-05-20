# Mathematical Writing and Statistical Algorithms Guidelines

- Keep Code as simple as possible. Goal is a Proof-of-Concept.
- Check numeric results for unexpected NaN or NA values and inspect the full output for numerical plausibility.
- Pure documentation or comment changes do not require running tests or lint checks.
  
### Testing, Debugging, Reproducibility

* Write unit tests with **testthat**; place them in `tests/testthat/`.
* Use `set.seed()` at every experiment start and record the seed in your report.
* Log intermediate metrics with `futile.logger`, `logger`, or simple `message()`.
* Track experiments in a plain CSV or YAML file; include hyper-parameters, commit hash, date, hostname, CPU or GPU type.
* Keep the code in a Git repository; cite the commit id that produced each table or figure.
