#!/bin/bash
set -euo pipefail
cd "$HOME/agent/work/ntpoly-gemm"
module load nvhpc-hpcx
export UCX_TLS=sm,sysv,posix
export UCX_NET_DEVICES=""
export LD_LIBRARY_PATH=/shared/software/hpc_sdk/Linux_aarch64/26.3/math_libs/lib64:$LD_LIBRARY_PATH

echo "=== debug: 32k rcut=15 --cusparse ==="
OMP_NUM_THREADS=1 mpirun -np 1 NTPoly/Build/bin/yukawa_gemm \
  --input matrices/fill32k_rcut15.bin --mode fixed --rounds 1 --warmup 0 --cusparse \
  --out results/debug_rcut15.json 2>&1
echo "exit=$?"
