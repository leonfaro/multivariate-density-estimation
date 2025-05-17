- Keep the entire README context in mind when analyzing outputs or generating code.
- Provide precise answers with mathematical rigour.
- Run R scripts using `Rscript` (e.g. `Rscript Code/00_setup.R`).
- Execute any tests or lint checks inside this environment if they exist.
  Use `./run_checks.sh` to run the repository's basic tests and linting.
- Follow the variable names and mathematical symbols in [Notation.md](Notation.md)
  when writing new code.
- Check new numeric results for unexpected `NaN` or `NA` values and inspect the
  full output for numerical plausibility from a mathematical statistics
  perspective. Apply *Mathematical Rigour* in all analyses.
- This project does not follow a classic software engineering workflow and will
  never be released as an R package. The code base is a vehicle for scientific
  experimentation in probability theory using R, so focus on reproducibility and
  precision rather than packaging.
- Pure documentation or comment changes do not require running tests or lint
  checks.
