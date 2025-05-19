#!/bin/sh

# Run basic tests and lint. Exit with non-zero if anything fails.

Rscript basic_tests.R || exit 1
Rscript -e "testthat::test_dir('tests/testthat')" || exit 1
bash lint.sh || exit 1

echo "Checks completed successfully."
