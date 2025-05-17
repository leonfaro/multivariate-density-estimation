#!/bin/bash
# Wrapper to execute the R workflow if R is available.

set -e

if ! command -v Rscript > /dev/null 2>&1; then
    echo "R not available; skipping R-based scripts." >&2
    exit 0
fi

Rscript Code/run_all.R "$@"
