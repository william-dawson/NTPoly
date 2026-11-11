#!/bin/bash
# Benchmark the --direct path (bypass distributed wrapper at np=1).
# Compare against the original distributed path.
set -euo pipefail

WORK="$HOME/agent/work/ntpoly-gemm"
cd "$WORK"
module load nvhpc-hpcx

export UCX_TLS=sm,sysv,posix
export UCX_NET_DEVICES=""

mkdir -p results

echo "=== direct-path benchmark ==="
date

# --- 32k rcut=6 (fastest, 0.12% fill) ---
M32="matrices/fill32k_rcut6.bin"

echo "--- 32k rcut=6 np=1 DIRECT ---"
OMP_NUM_THREADS=1 mpirun -np 1 NTPoly/Build/bin/yukawa_gemm \
  --input "$M32" --mode fixed --rounds 10 --warmup 3 --direct \
  --out results/direct_32k_rcut6_np1.json 2>&1
echo "exit=$?"

echo "--- 32k rcut=6 np=1 DISTRIBUTED (for comparison) ---"
OMP_NUM_THREADS=1 mpirun -np 1 NTPoly/Build/bin/yukawa_gemm \
  --input "$M32" --mode fixed --rounds 10 --warmup 3 \
  --out results/dist_32k_rcut6_np1.json 2>&1
echo "exit=$?"

# --- 64k (the main benchmark matrix) ---
M64="matrices/yukawa_sc_64000.bin"

echo "--- 64k np=1 DIRECT ---"
OMP_NUM_THREADS=1 mpirun -np 1 NTPoly/Build/bin/yukawa_gemm \
  --input "$M64" --mode fixed --rounds 10 --warmup 3 --direct \
  --out results/direct_64k_np1.json 2>&1
echo "exit=$?"

echo "--- 64k np=1 DISTRIBUTED (for comparison) ---"
OMP_NUM_THREADS=1 mpirun -np 1 NTPoly/Build/bin/yukawa_gemm \
  --input "$M64" --mode fixed --rounds 10 --warmup 3 \
  --out results/dist_64k_np1.json 2>&1
echo "exit=$?"

# --- 32k fill sweep with DIRECT ---
for rcut in 10 15 20; do
  M="matrices/fill32k_rcut${rcut}.bin"
  if [ ! -f "$M" ]; then
    echo "Skipping rcut=$rcut (matrix not found)"
    continue
  fi
  echo "--- 32k rcut=$rcut np=1 DIRECT ---"
  OMP_NUM_THREADS=1 mpirun -np 1 NTPoly/Build/bin/yukawa_gemm \
    --input "$M" --mode fixed --rounds 5 --warmup 2 --direct \
    --out "results/direct_32k_rcut${rcut}_np1.json" 2>&1
  echo "exit=$?"
done

echo "=== done ==="
date
