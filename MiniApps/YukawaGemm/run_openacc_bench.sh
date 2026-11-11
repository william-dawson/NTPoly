#!/bin/bash
# OpenACC GPU benchmark — test the ported NTPoly on GPU.
# Compare against the CPU baseline and the cuSPARSE results.
set -euo pipefail

WORK="$HOME/agent/work/ntpoly-gemm"
cd "$WORK"
module load nvhpc-hpcx

export UCX_TLS=sm,sysv,posix
export UCX_NET_DEVICES=""

mkdir -p results

echo "=== OpenACC GPU benchmark ==="
date

# --- 64k matrix at np=1 (single rank, pure GPU) ---
M64="matrices/yukawa_sc_64000.bin"
echo "--- 64k np=1 (GPU) ---"
OMP_NUM_THREADS=1 mpirun -np 1 NTPoly/Build/bin/yukawa_gemm \
  --input "$M64" --mode fixed --rounds 10 --warmup 3 \
  --out results/openacc_64k_np1.json 2>&1
echo "exit=$?"

# --- 64k at np=4 (distributed, each rank on GPU) ---
echo "--- 64k np=4 omp=1 (GPU) ---"
OMP_NUM_THREADS=1 mpirun --oversubscribe -np 4 NTPoly/Build/bin/yukawa_gemm \
  --input "$M64" --mode fixed --rounds 10 --warmup 3 \
  --out results/openacc_64k_np4.json 2>&1
echo "exit=$?"

# --- Fill sweep: 32k matrices at np=1 ---
# This tests whether OpenACC handles high fill (where cuSPARSE crashed)
for rcut in 6 10 15 20; do
  M="matrices/fill32k_rcut${rcut}.bin"
  if [ ! -f "$M" ]; then
    echo "Skipping rcut=$rcut (matrix not found)"
    continue
  fi
  echo "--- 32k rcut=$rcut np=1 (GPU) ---"
  OMP_NUM_THREADS=1 mpirun -np 1 NTPoly/Build/bin/yukawa_gemm \
    --input "$M" --mode fixed --rounds 5 --warmup 2 \
    --out "results/openacc_32k_rcut${rcut}_np1.json" 2>&1
  echo "exit=$?"
done

# --- 175k matrix at np=1 (was OOM on CPU) ---
M175="matrices/yukawa_sc_175616.bin"
if [ -f "$M175" ]; then
  echo "--- 175k np=1 (GPU) ---"
  OMP_NUM_THREADS=1 mpirun -np 1 NTPoly/Build/bin/yukawa_gemm \
    --input "$M175" --mode fixed --rounds 5 --warmup 2 \
    --out results/openacc_175k_np1.json 2>&1
  echo "exit=$?"
fi

echo "=== done ==="
ls -la results/openacc_*.json
date
