#!/bin/bash
# Benchmark the cuSPARSE hybrid backend (--cusparse flag).
set -euo pipefail
cd "$HOME/agent/work/ntpoly-gemm"
module load nvhpc-hpcx
export UCX_TLS=sm,sysv,posix
export UCX_NET_DEVICES=""
export LD_LIBRARY_PATH=/shared/software/hpc_sdk/Linux_aarch64/26.3/math_libs/lib64:$LD_LIBRARY_PATH

mkdir -p results

echo "=== cuSPARSE hybrid benchmark ==="
date

# --- 64k (the primary benchmark) ---
M64="matrices/yukawa_sc_64000.bin"
echo "--- 64k np=1 cuSPARSE ---"
OMP_NUM_THREADS=1 mpirun -np 1 NTPoly/Build/bin/yukawa_gemm \
  --input "$M64" --mode fixed --rounds 10 --warmup 3 --cusparse \
  --out results/cusparse_64k_np1.json 2>&1
echo "exit=$?"

# --- 32k fill sweep ---
for rcut in 6 10 15 20; do
  M="matrices/fill32k_rcut${rcut}.bin"
  [ -f "$M" ] || continue
  echo "--- 32k rcut=$rcut np=1 cuSPARSE ---"
  OMP_NUM_THREADS=1 mpirun -np 1 NTPoly/Build/bin/yukawa_gemm \
    --input "$M" --mode fixed --rounds 10 --warmup 3 --cusparse \
    --out "results/cusparse_32k_rcut${rcut}_np1.json" 2>&1
  echo "exit=$?"
done

echo "=== done ==="
date
