# Mathematical Writing and Statistical Algorithms Guidelines

- Always adhere to Notation.md.
- Check numeric results for unexpected NaN or NA values and inspect the full output for numerical plausibility.
- Pure documentation or comment changes do not require running tests or lint checks.

## 1. Mathematical Rigor and Precise Notation

- Define every symbol exactly once before use, and never reuse it with a different meaning.  
- Make every statement answer both questions: *what does it mean?* and *why is it true?*  
- Replace logic symbols in running text with words.
- Keep notation simple. Prefer set-builder notation over complicated indices.  


## 2. Implementation of Statistical Algorithms (R plus Ubuntu 24)

### Numerical Stability and Conditioning

* Inspect the condition number of every design or covariance matrix.
  If it is large, the problem is ill-conditioned and results may vary wildly with tiny data changes.
* Always center and scale features before numeric work. Scaling reduces the condition number and avoids loss of precision when very small and very large numbers meet.
* Prefer numerically stable algorithms: For solving linear systems use `qr()` or `chol()` instead of naive `solve()`. Work in log space for likelihoods; apply the log-sum-exp trick via `matrixStats::logSumExp()`.
* Guard against underflow and overflow; for example compute softmax as

  ```r
  softmax <- function(x) {
    exp_x <- exp(x - max(x))
    exp_x / sum(exp_x)
  }
  ```

### Accuracy versus Compute Time

* Seek a pragmatic optimum: sufficient accuracy at acceptable runtime.
* Tune Monte Carlo sample sizes, integration grids, or step sizes until incremental precision gains flatten.
* Choose model complexity that balances interpretability, overfitting risk, and computation.
* Document sensitivity tests that compare coarse and fine approximations.

### Testing, Debugging, Reproducibility

* Write unit tests with **testthat**; place them in `tests/testthat/`.
* Use `set.seed()` at every experiment start and record the seed in your report.
* Log intermediate metrics with `futile.logger`, `logger`, or simple `message()`.
* Track experiments in a plain CSV or YAML file; include hyper-parameters, commit hash, date, hostname, CPU or GPU type.
* Keep the code in a Git repository; cite the commit id that produced each table or figure.

## 3. Clarity, Conciseness, Logical Structure

* Every sentence must advance the argument; delete filler such as "it is obvious that".
* Follow the rule **as detailed as necessary, as brief as possible**.
* Present ideas in complete grammatical sentences; avoid telegram style.
* Do **not** start a sentence with a formula.
* Embed every displayed equation grammatically and end the surrounding sentence with a period.

## 4. Proof Construction and Documentation

* Start each result in a numbered `Lemma`, `Theorem`, or `Proposition` environment, then begin the proof with **Proof.**
* The statement must be self-contained; no variables may appear that are defined only later.
* As soon as a variable appears in the proof, specify its type and domain.
* Make the proof strategy explicit up front.
* Inside the proof, explain every nontrivial inference or cite a known result.
* Break long calculations with guiding sentences that announce what comes next and why.
* Never introduce new definitions inside the proof that are needed later; place them before the proof.
* Conclude with a q.e.d. symbol or simply end the proof; avoid redundant closing phrases.

## 5. Revision and Proofreading

* Read the entire text aloud; fix awkward formulations and punctuation.
* Check that the prose is understandable even if formulas are skipped.
* Iterate: clarity often emerges only after several rewrites.
