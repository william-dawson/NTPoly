#!/bin/bash
# Re-profile the optimized code (sort fix + direct path) with PGI_ACC_TIME.
set -euo pipefail
cd "$HOME/agent/work/ntpoly-gemm"
module load nvhpc-hpcx
export UCX_TLS=sm,sysv,posix
export UCX_NET_DEVICES=""

# Enable OpenACC per-region profiling (stderr)
export PGI_ACC_TIME=1

mkdir -p results

echo "=== profile 64k np=1 direct (optimized) ==="
OMP_NUM_THREADS=1 mpirun -np 1 NTPoly/Build/bin/yukawa_gemm \
  --input matrices/yukawa_sc_64000.bin --mode fixed --rounds 3 --warmup 1 --direct \
  --out results/profile_opt_64k_np1.json 2>&1

echo "=== profile 32k rcut=6 np=1 direct (optimized) ==="
OMP_NUM_THREADS=1 mpirun -np 1 NTPoly/Build/bin/yukawa_gemm \
  --input matrices/fill32k_rcut6.bin --mode fixed --rounds 3 --warmup 1 --direct \
  --out results/profile_opt_32k_rcut6_np1.json 2>&1

echo "=== done ==="
