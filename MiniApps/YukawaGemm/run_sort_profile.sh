#!/bin/bash
# Quick profile: measure sort/CSR time in PruneList.
set -euo pipefail
cd "$HOME/agent/work/ntpoly-gemm"
module load nvhpc-hpcx
export UCX_TLS=sm,sysv,posix
export UCX_NET_DEVICES=""

echo "=== 32k rcut=6 direct profile ==="
OMP_NUM_THREADS=1 mpirun -np 1 NTPoly/Build/bin/yukawa_gemm \
  --input matrices/fill32k_rcut6.bin --mode fixed --rounds 3 --warmup 1 --direct \
  --out /dev/null 2>&1

echo "=== 64k direct profile ==="
OMP_NUM_THREADS=1 mpirun -np 1 NTPoly/Build/bin/yukawa_gemm \
  --input matrices/yukawa_sc_64000.bin --mode fixed --rounds 2 --warmup 1 --direct \
  --out /dev/null 2>&1

echo "=== done ==="
