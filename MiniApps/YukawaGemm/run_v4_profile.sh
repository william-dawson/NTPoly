#!/bin/bash
# Profile cuSPARSE hybrid v4: per-phase timing with warm-call optimization.
set -euo pipefail
cd "$HOME/agent/work/ntpoly-gemm"
module load nvhpc-hpcx
export UCX_TLS=sm,sysv,posix
export UCX_NET_DEVICES=""
export LD_LIBRARY_PATH=/shared/software/hpc_sdk/Linux_aarch64/26.3/math_libs/lib64:$LD_LIBRARY_PATH

BIN=NTPoly/Build/bin/yukawa_gemm

echo "=== cuSPARSE v4 per-phase profiling ==="
echo "--- 64k, 4 rounds, 1 warmup (default alg = ALG1) ---"
OMP_NUM_THREADS=1 mpirun -np 1 $BIN \
  --input matrices/yukawa_sc_64000.bin --mode fixed --rounds 4 --warmup 1 --cusparse \
  --out results/v4_phase_64k.json 2>&1
echo "exit=$?"

echo ""
echo "--- 32k rcut=6, 4 rounds, 1 warmup ---"
OMP_NUM_THREADS=1 mpirun -np 1 $BIN \
  --input matrices/fill32k_rcut6.bin --mode fixed --rounds 4 --warmup 1 --cusparse \
  --out results/v4_phase_32k_rcut6.json 2>&1
echo "exit=$?"

echo "=== done ==="
