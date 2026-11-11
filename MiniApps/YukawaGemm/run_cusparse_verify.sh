#!/bin/bash
# Run cuSPARSE SpGEMM benchmark on 64k matrix with thresholding.
set -euo pipefail
cd "$HOME/agent/work/ntpoly-gemm"
module load nvhpc-hpcx

export CUDA_HOME=/shared/software/hpc_sdk/Linux_aarch64/26.3/cuda
export LD_LIBRARY_PATH=/shared/software/hpc_sdk/Linux_aarch64/26.3/math_libs/lib64:$LD_LIBRARY_PATH

mkdir -p results

echo "=== cuSPARSE 64k ==="
./bench_cusparse matrices/yukawa_sc_64000.bin --rounds 10 --warmup 3 --threshold 1e-4 2>&1

echo "=== cuSPARSE 32k rcut=6 ==="
./bench_cusparse matrices/fill32k_rcut6.bin --rounds 10 --warmup 3 --threshold 1e-4 2>&1

echo "=== cuSPARSE 32k rcut=10 ==="
./bench_cusparse matrices/fill32k_rcut10.bin --rounds 10 --warmup 3 --threshold 1e-4 2>&1

echo "=== done ==="
