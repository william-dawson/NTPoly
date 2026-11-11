#!/bin/bash
# v4 GPU-only benchmark: cuSPARSE hybrid with persistent buffers + warm calls.
set -euo pipefail
cd "$HOME/agent/work/ntpoly-gemm"
module load nvhpc-hpcx
export UCX_TLS=sm,sysv,posix
export UCX_NET_DEVICES=""
export LD_LIBRARY_PATH=/shared/software/hpc_sdk/Linux_aarch64/26.3/math_libs/lib64:$LD_LIBRARY_PATH

mkdir -p results
BIN=NTPoly/Build/bin/yukawa_gemm

echo "=== v4 cuSPARSE hybrid GPU benchmark ==="
date

echo "--- 64k --cusparse ---"
OMP_NUM_THREADS=1 mpirun -np 1 $BIN \
  --input matrices/yukawa_sc_64000.bin --mode fixed --rounds 10 --warmup 3 --cusparse \
  --out results/v4_cusparse_64k.json 2>&1 | grep -E '(round|gemm|trace|nnz:|asym)'
echo "exit=$?"

echo "--- 32k rcut=6 --cusparse ---"
OMP_NUM_THREADS=1 mpirun -np 1 $BIN \
  --input matrices/fill32k_rcut6.bin --mode fixed --rounds 10 --warmup 3 --cusparse \
  --out results/v4_cusparse_32k_rcut6.json 2>&1 | grep -E '(round|gemm|trace|nnz:|asym)'
echo "exit=$?"

echo "--- 32k rcut=10 --cusparse ---"
OMP_NUM_THREADS=1 mpirun -np 1 $BIN \
  --input matrices/fill32k_rcut10.bin --mode fixed --rounds 10 --warmup 3 --cusparse \
  --out results/v4_cusparse_32k_rcut10.json 2>&1 | grep -E '(round|gemm|trace|nnz:|asym)'
echo "exit=$?"

echo "--- 32k rcut=15 --cusparse --direct (fallback) ---"
OMP_NUM_THREADS=1 mpirun -np 1 $BIN \
  --input matrices/fill32k_rcut15.bin --mode fixed --rounds 5 --warmup 2 --cusparse --direct \
  --out results/v4_cusparse_32k_rcut15.json 2>&1 | grep -E '(round|gemm|trace|nnz:|asym)'
echo "exit=$?"

echo "--- 32k rcut=20 --cusparse --direct (fallback) ---"
OMP_NUM_THREADS=1 mpirun -np 1 $BIN \
  --input matrices/fill32k_rcut20.bin --mode fixed --rounds 5 --warmup 2 --cusparse --direct \
  --out results/v4_cusparse_32k_rcut20.json 2>&1 | grep -E '(round|gemm|trace|nnz:|asym)'
echo "exit=$?"

echo "=== done ==="
date
