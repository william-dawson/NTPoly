#!/bin/bash
# Test cuSPARSE fallback at high fill (rcut=15, 20) -- both --cusparse alone
# (distributed fallback) and --cusparse --direct (OpenACC fallback).
set -euo pipefail
cd "$HOME/agent/work/ntpoly-gemm"
module load nvhpc-hpcx
export UCX_TLS=sm,sysv,posix
export UCX_NET_DEVICES=""
export LD_LIBRARY_PATH=/shared/software/hpc_sdk/Linux_aarch64/26.3/math_libs/lib64:$LD_LIBRARY_PATH

mkdir -p results

echo "=== cuSPARSE fallback test ==="
date

# --- 32k rcut=15: cuSPARSE alone (distributed fallback) ---
echo "--- 32k rcut=15 --cusparse (distributed fallback) ---"
OMP_NUM_THREADS=1 mpirun -np 1 NTPoly/Build/bin/yukawa_gemm \
  --input matrices/fill32k_rcut15.bin --mode fixed --rounds 3 --warmup 1 --cusparse \
  --out results/cusparse_32k_rcut15_fallback.json 2>&1
echo "exit=$?"

# --- 32k rcut=20: cuSPARSE alone (distributed fallback) ---
echo "--- 32k rcut=20 --cusparse (distributed fallback) ---"
OMP_NUM_THREADS=1 mpirun -np 1 NTPoly/Build/bin/yukawa_gemm \
  --input matrices/fill32k_rcut20.bin --mode fixed --rounds 3 --warmup 1 --cusparse \
  --out results/cusparse_32k_rcut20_fallback.json 2>&1
echo "exit=$?"

# --- 32k rcut=15: cuSPARSE + direct (OpenACC fallback) ---
echo "--- 32k rcut=15 --cusparse --direct (OpenACC fallback) ---"
OMP_NUM_THREADS=1 mpirun -np 1 NTPoly/Build/bin/yukawa_gemm \
  --input matrices/fill32k_rcut15.bin --mode fixed --rounds 3 --warmup 1 --cusparse --direct \
  --out results/cusparse_32k_rcut15_direct_fb.json 2>&1
echo "exit=$?"

# --- 32k rcut=20: cuSPARSE + direct (OpenACC fallback) ---
echo "--- 32k rcut=20 --cusparse --direct (OpenACC fallback) ---"
OMP_NUM_THREADS=1 mpirun -np 1 NTPoly/Build/bin/yukawa_gemm \
  --input matrices/fill32k_rcut20.bin --mode fixed --rounds 3 --warmup 1 --cusparse --direct \
  --out results/cusparse_32k_rcut20_direct_fb.json 2>&1
echo "exit=$?"

echo "=== done ==="
date
