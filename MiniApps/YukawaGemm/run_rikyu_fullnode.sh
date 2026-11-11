#!/bin/bash
# NTPoly YukawaGemm full-node CPU baseline — v10
# Two matrices:
#   sc 40^3 = 64,000 atoms: full sweep (all configs fit, no INT32 overflow)
#   sc 56^3 = 175,616 atoms: np=144 only (pool=5.1 GB, fits 11 GB/rank)
set -euo pipefail

WORK="$HOME/agent/work/ntpoly-gemm"
cd "$WORK"
module load nvhpc-hpcx

export UCX_TLS=sm,sysv,posix
export UCX_NET_DEVICES=""

mkdir -p matrices results

# --- Generate sc 40^3 (64,000 atoms) ---
M64="matrices/yukawa_sc_64000.bin"
if [ ! -f "$M64" ]; then
  echo "=== generating sc 40^3 (64,000 atoms) ==="
  NTPoly/Build/bin/generate_yukawa --side 40 --lattice 2.7 --screening 1.0 \
    --rcut 6.0 --jitter 0.05 --seed 42 --diag 1.0 \
    --out matrices/yukawa_sc_64000.mtx
  mpirun -np 1 NTPoly/Build/bin/yukawa_gemm \
    --input matrices/yukawa_sc_64000.mtx --convert --out "$M64"
  rm -f matrices/yukawa_sc_64000.mtx
fi

# --- Generate sc 56^3 (175,616 atoms) if not already ---
M175="matrices/yukawa_sc_175616.bin"
if [ ! -f "$M175" ]; then
  echo "=== generating sc 56^3 (175,616 atoms) ==="
  NTPoly/Build/bin/generate_yukawa --side 56 --lattice 2.7 --screening 1.0 \
    --rcut 6.0 --jitter 0.05 --seed 42 --diag 1.0 \
    --out matrices/yukawa_sc_175616.mtx
  mpirun -np 1 NTPoly/Build/bin/yukawa_gemm \
    --input matrices/yukawa_sc_175616.mtx --convert --out "$M175"
  rm -f matrices/yukawa_sc_175616.mtx
fi

echo "=== sc 40^3 (64k) full sweep ==="
# Pure-MPI: 4, 8, 16, 36, 72, 144 ranks
for np in 4 8 16 36 72 144; do
  echo "--- np=$np omp=1 ---"
  OMP_NUM_THREADS=1 mpirun --oversubscribe -np "$np" \
    NTPoly/Build/bin/yukawa_gemm \
    --input "$M64" --mode fixed --rounds 10 --warmup 3 \
    --out "results/sc64k_np${np}_omp1.json" 2>&1
  echo "exit=$?"
done

# Hybrid MPI×OpenMP (total=144 cores)
for np in 4 8 16 36; do
  omp=$((144 / np))
  echo "--- np=$np omp=$omp ---"
  OMP_NUM_THREADS="$omp" mpirun --oversubscribe -np "$np" \
    NTPoly/Build/bin/yukawa_gemm \
    --input "$M64" --mode fixed --rounds 10 --warmup 3 \
    --out "results/sc64k_np${np}_omp${omp}.json" 2>&1
  echo "exit=$?"
done

echo "=== sc 56^3 (175k) np=144 only ==="
echo "--- np=144 omp=1 ---"
OMP_NUM_THREADS=1 mpirun --oversubscribe -np 144 \
  NTPoly/Build/bin/yukawa_gemm \
  --input "$M175" --mode fixed --rounds 10 --warmup 3 \
  --out results/sc175k_np144_omp1.json 2>&1
echo "exit=$?"

echo "=== done ==="
ls -la results/sc64k_*.json results/sc175k_*.json 2>&1
date
