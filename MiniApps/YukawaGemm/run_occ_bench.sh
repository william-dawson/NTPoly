#!/bin/bash
# Benchmark occupancy-increased MultiplyBlock with PGI_ACC_TIME profiling.
set -euo pipefail
cd "$HOME/agent/work/ntpoly-gemm"
module load nvhpc-hpcx
export UCX_TLS=sm,sysv,posix
export UCX_NET_DEVICES=""
export PGI_ACC_TIME=1
mkdir -p results

echo "=== occupancy test: num_gangs(4096) vector_length(256) ==="
date

echo "--- 64k np=1 direct ---"
OMP_NUM_THREADS=1 mpirun -np 1 NTPoly/Build/bin/yukawa_gemm \
  --input matrices/yukawa_sc_64000.bin --mode fixed --rounds 5 --warmup 2 --direct \
  --out results/occ4096_64k_np1.json 2>&1

echo "--- 32k rcut=6 np=1 direct ---"
OMP_NUM_THREADS=1 mpirun -np 1 NTPoly/Build/bin/yukawa_gemm \
  --input matrices/fill32k_rcut6.bin --mode fixed --rounds 5 --warmup 2 --direct \
  --out results/occ4096_32k_rcut6_np1.json 2>&1

echo "=== done ==="
date
