#!/bin/bash
# Profile the cuSPARSE hybrid path: per-phase timing inside the C wrapper.
set -euo pipefail
cd "$HOME/agent/work/ntpoly-gemm"
module load nvhpc-hpcx
export UCX_TLS=sm,sysv,posix
export UCX_NET_DEVICES=""
export LD_LIBRARY_PATH=/shared/software/hpc_sdk/Linux_aarch64/26.3/math_libs/lib64:$LD_LIBRARY_PATH

echo "=== cuSPARSE hybrid profiling ==="
echo "--- 64k, 3 rounds, 1 warmup ---"
OMP_NUM_THREADS=1 mpirun -np 1 NTPoly/Build/bin/yukawa_gemm \
  --input matrices/yukawa_sc_64000.bin --mode fixed --rounds 3 --warmup 1 --cusparse \
  --out results/profile_cusparse_64k.json 2>&1
echo "exit=$?"

echo ""
echo "--- 32k rcut=6, 3 rounds, 1 warmup ---"
OMP_NUM_THREADS=1 mpirun -np 1 NTPoly/Build/bin/yukawa_gemm \
  --input matrices/fill32k_rcut6.bin --mode fixed --rounds 3 --warmup 1 --cusparse \
  --out results/profile_cusparse_32k_rcut6.json 2>&1
echo "exit=$?"

echo "=== done ==="
