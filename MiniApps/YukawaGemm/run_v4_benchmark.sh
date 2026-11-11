#!/bin/bash
# Profile cuSPARSE hybrid with all persistent buffers + algorithm selection.
set -euo pipefail
cd "$HOME/agent/work/ntpoly-gemm"
module load nvhpc-hpcx
export UCX_TLS=sm,sysv,posix
export UCX_NET_DEVICES=""
export LD_LIBRARY_PATH=/shared/software/hpc_sdk/Linux_aarch64/26.3/math_libs/lib64:$LD_LIBRARY_PATH

BIN=NTPoly/Build/bin/yukawa_gemm

echo "=== 1. Verify correctness (default algorithm chain) ==="
echo "--- 64k, 3 rounds, 1 warmup ---"
OMP_NUM_THREADS=1 mpirun -np 1 $BIN \
  --input matrices/yukawa_sc_64000.bin --mode fixed --rounds 3 --warmup 1 --cusparse \
  --out results/v4_default_64k.json 2>&1
echo "exit=$?"

echo ""
echo "=== 2. Benchmark ALL algorithms (CUSPARSE_ALG=-2) ==="
echo "--- 64k, 1 round, 0 warmup (alg-bench runs all 4 internally) ---"
CUSPARSE_ALG=-2 OMP_NUM_THREADS=1 mpirun -np 1 $BIN \
  --input matrices/yukawa_sc_64000.bin --mode fixed --rounds 1 --warmup 0 --cusparse \
  --out results/v4_algbench_64k.json 2>&1
echo "exit=$?"

echo ""
echo "=== 3. Per-algorithm timing (64k, 3 rounds each) ==="
for alg in 0 1 2 3; do
  echo "--- ALG=$alg ---"
  CUSPARSE_ALG=$alg OMP_NUM_THREADS=1 mpirun -np 1 $BIN \
    --input matrices/yukawa_sc_64000.bin --mode fixed --rounds 3 --warmup 1 --cusparse \
    --out results/v4_alg${alg}_64k.json 2>&1 | grep -E '(round|prof.*SpGEMM|prof.*mark|prof.*compact|prof.*D2H|trace|nnz:|asym)'
  echo "exit=$?"
done

echo ""
echo "=== 4. 32k rcut=6 (default + best alg) ==="
echo "--- default ---"
OMP_NUM_THREADS=1 mpirun -np 1 $BIN \
  --input matrices/fill32k_rcut6.bin --mode fixed --rounds 3 --warmup 1 --cusparse \
  --out results/v4_default_32k_rcut6.json 2>&1 | grep -E '(round|trace|nnz:|asym)'
echo "exit=$?"

echo ""
echo "=== done ==="
