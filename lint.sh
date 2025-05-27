#!/bin/sh
# Simple lint: check for trailing whitespace in R files

fail=0
for f in *.R; do
  case "$f" in
    demo* ) continue ;;
  esac
  if grep -nE "[[:space:]]$" "$f" >/dev/null; then
    echo "Trailing whitespace found in $f" >&2
    fail=1
  fi
  if grep -n $'\t' "$f" >/dev/null; then
    echo "Tab character found in $f" >&2
    fail=1
  fi
done
if [ $fail -ne 0 ]; then
  exit $fail
fi
Rscript extra_lint.R || exit 1
exit 0
