#!/bin/bash
# Benchmark with insertion-sort fix + direct path.
set -euo pipefail
cd "$HOME/agent/work/ntpoly-gemm"
module load nvhpc-hpcx
export UCX_TLS=sm,sysv,posix
export UCX_NET_DEVICES=""
mkdir -p results

echo "=== sort-fix benchmark ==="
date

# --- 32k fill sweep (direct path) ---
for rcut in 6 10 15 20; do
  M="matrices/fill32k_rcut${rcut}.bin"
  [ -f "$M" ] || continue
  echo "--- 32k rcut=$rcut np=1 DIRECT ---"
  OMP_NUM_THREADS=1 mpirun -np 1 NTPoly/Build/bin/yukawa_gemm \
    --input "$M" --mode fixed --rounds 10 --warmup 3 --direct \
    --out "results/sortfix_32k_rcut${rcut}_np1.json" 2>&1
  echo "exit=$?"
done

# --- 64k (direct path) ---
M64="matrices/yukawa_sc_64000.bin"
echo "--- 64k np=1 DIRECT ---"
OMP_NUM_THREADS=1 mpirun -np 1 NTPoly/Build/bin/yukawa_gemm \
  --input "$M64" --mode fixed --rounds 10 --warmup 3 --direct \
  --out results/sortfix_64k_np1.json 2>&1
echo "exit=$?"

# --- 64k distributed (for comparison) ---
echo "--- 64k np=1 DISTRIBUTED ---"
OMP_NUM_THREADS=1 mpirun -np 1 NTPoly/Build/bin/yukawa_gemm \
  --input "$M64" --mode fixed --rounds 5 --warmup 2 \
  --out results/sortfix_64k_np1_dist.json 2>&1
echo "exit=$?"

echo "=== done ==="
date
