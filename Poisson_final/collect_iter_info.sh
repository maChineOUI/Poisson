#!/usr/bin/env bash
set -euo pipefail

echo "== (2.3) Extract tol/maxit from src/tp_poisson1D_iter.c =="
grep -nE 'double[[:space:]]+tol[[:space:]]*=|int[[:space:]]+maxit[[:space:]]*=' -n src/tp_poisson1D_iter.c || true
echo

echo "== Rebuild (clean + make) =="
# 让 clean 不因为文件不存在而失败
make clean || true
mkdir -p bin
make
echo

run_one () {
  local m="$1"
  echo "=============================="
  echo "== (2.1) Run tpPoisson1D_iter method=$m =="
  ./bin/tpPoisson1D_iter "$m"
  echo
  echo "== (2.2) RESVEC stats (method=$m) =="
  wc -l results/RESVEC.dat || true
  echo "-- head -n 3 --"
  head -n 3 results/RESVEC.dat || true
  echo "-- tail -n 5 --"
  tail -n 5 results/RESVEC.dat || true
  echo
}

run_one 0
run_one 1
run_one 2

echo "Done."
