#!/usr/bin/env bash
set -euo pipefail

echo "== Clean =="
make clean

echo "== Build =="
make

echo "== Show binary timestamp =="
ls -l bin/tpPoisson1D_iter

echo
echo "== Run ALPHA (0) =="
./bin/tpPoisson1D_iter 0
echo "RESVEC lines:"; wc -l /app/results/RESVEC.dat || true
echo "RESVEC tail:"; tail -n 5 /app/results/RESVEC.dat || true
echo

echo "== Run JAC (1) =="
./bin/tpPoisson1D_iter 1
echo "RESVEC lines:"; wc -l /app/results/RESVEC.dat || true
echo "RESVEC tail:"; tail -n 5 /app/results/RESVEC.dat || true
echo

echo "== Run GS (2) =="
./bin/tpPoisson1D_iter 2
echo "RESVEC lines:"; wc -l /app/results/RESVEC.dat || true
echo "RESVEC tail:"; tail -n 5 /app/results/RESVEC.dat || true
echo

echo "== Sanity: ensure the new printf string exists inside binary =="
strings bin/tpPoisson1D_iter | grep "\[ITER VERSION\]" || (echo "ERROR: new version string not found in binary" && exit 1)

echo
echo "All done."
