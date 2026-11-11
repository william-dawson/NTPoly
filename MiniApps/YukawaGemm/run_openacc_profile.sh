#!/bin/bash
# Profile the OpenACC GPU code to see where time is spent.
# PGI_ACC_TIME=1 makes the OpenACC runtime print per-region timing to stderr.
set -euo pipefail
cd "$HOME/agent/work/ntpoly-gemm"
module load nvhpc-hpcx
export UCX_TLS=sm,sysv,posix
export UCX_NET_DEVICES=""

# Enable OpenACC per-region profiling (stderr)
export PGI_ACC_TIME=1
export NVACC_TIME=1

mkdir -p results

# --- 32k rcut=6 np=1 (fastest case, ~0.7s/round) ---
M32="matrices/fill32k_rcut6.bin"
echo "=== profiling 32k rcut=6 np=1 ==="
OMP_NUM_THREADS=1 mpirun -np 1 NTPoly/Build/bin/yukawa_gemm \
  --input "$M32" --mode fixed --rounds 3 --warmup 1 \
  --out results/profile_32k_rcut6_np1.json 2>&1

# --- 64k np=1 (slowest case, ~5.6s/round) ---
M64="matrices/yukawa_sc_64000.bin"
echo "=== profiling 64k np=1 ==="
OMP_NUM_THREADS=1 mpirun -np 1 NTPoly/Build/bin/yukawa_gemm \
  --input "$M64" --mode fixed --rounds 2 --warmup 1 \
  --out results/profile_64k_np1.json 2>&1

echo "=== done ==="
