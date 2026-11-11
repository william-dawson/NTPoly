#!/bin/bash
# Final cuSPARSE hybrid benchmark: low fill uses cuSPARSE, high fill uses
# --cusparse --direct for OpenACC fallback when cuSPARSE fails.
set -euo pipefail
cd "$HOME/agent/work/ntpoly-gemm"
module load nvhpc-hpcx
export UCX_TLS=sm,sysv,posix
export UCX_NET_DEVICES=""
export LD_LIBRARY_PATH=/shared/software/hpc_sdk/Linux_aarch64/26.3/math_libs/lib64:$LD_LIBRARY_PATH

mkdir -p results

echo "=== cuSPARSE hybrid final benchmark ==="
date

# --- 64k: cuSPARSE (low fill, cuSPARSE works) ---
echo "--- 64k --cusparse ---"
OMP_NUM_THREADS=1 mpirun -np 1 NTPoly/Build/bin/yukawa_gemm \
  --input matrices/yukawa_sc_64000.bin --mode fixed --rounds 10 --warmup 3 --cusparse \
  --out results/cusparse_final_64k.json 2>&1
echo "exit=$?"

# --- 32k rcut=6: cuSPARSE (low fill) ---
echo "--- 32k rcut=6 --cusparse ---"
OMP_NUM_THREADS=1 mpirun -np 1 NTPoly/Build/bin/yukawa_gemm \
  --input matrices/fill32k_rcut6.bin --mode fixed --rounds 10 --warmup 3 --cusparse \
  --out results/cusparse_final_32k_rcut6.json 2>&1
echo "exit=$?"

# --- 32k rcut=10: cuSPARSE (medium fill) ---
echo "--- 32k rcut=10 --cusparse ---"
OMP_NUM_THREADS=1 mpirun -np 1 NTPoly/Build/bin/yukawa_gemm \
  --input matrices/fill32k_rcut10.bin --mode fixed --rounds 10 --warmup 3 --cusparse \
  --out results/cusparse_final_32k_rcut10.json 2>&1
echo "exit=$?"

# --- 32k rcut=15: cusparse+direct (cuSPARSE fails, OpenACC fallback) ---
echo "--- 32k rcut=15 --cusparse --direct (OpenACC fallback) ---"
OMP_NUM_THREADS=1 mpirun -np 1 NTPoly/Build/bin/yukawa_gemm \
  --input matrices/fill32k_rcut15.bin --mode fixed --rounds 10 --warmup 3 --cusparse --direct \
  --out results/cusparse_final_32k_rcut15.json 2>&1
echo "exit=$?"

# --- 32k rcut=20: cusparse+direct (cuSPARSE fails, OpenACC fallback) ---
echo "--- 32k rcut=20 --cusparse --direct (OpenACC fallback) ---"
OMP_NUM_THREADS=1 mpirun -np 1 NTPoly/Build/bin/yukawa_gemm \
  --input matrices/fill32k_rcut20.bin --mode fixed --rounds 10 --warmup 3 --cusparse --direct \
  --out results/cusparse_final_32k_rcut20.json 2>&1
echo "exit=$?"

echo "=== done ==="
date
