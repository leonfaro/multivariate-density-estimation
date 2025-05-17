#!/bin/sh
# Simple lint: check for trailing whitespace in R files

fail=0
for f in Code/*.R tests/*.R; do
  case "$f" in
    Code/demo* ) continue ;;
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
exit $fail
