#!/bin/bash
# Fill sweep: benchmark CPU (NTPoly) vs GPU (cuSPARSE) at different sparsity.
# Matrix: sc 32^3 = 32,768 atoms. Sweep rcut to vary fill from ~0.1% to ~5%.
#
# rcut=6:  ~46 nnz/row   → 0.14% fill  (our current regime)
# rcut=10: ~200 nnz/row  → 0.65% fill
# rcut=15: ~700 nnz/row  → 2.2%  fill
# rcut=20: ~1600 nnz/row → 5.2%  fill  (realistic DFT density matrix)
set -euo pipefail

WORK="$HOME/agent/work/ntpoly-gemm"
cd "$WORK"
module load nvhpc-hpcx

export UCX_TLS=sm,sysv,posix
export UCX_NET_DEVICES=""

mkdir -p matrices results

for rcut in 6 10 15 20; do
  MATRIX="matrices/fill32k_rcut${rcut}.bin"

  # Generate if not present
  if [ ! -f "$MATRIX" ]; then
    echo "=== generating rcut=$rcut ==="
    NTPoly/Build/bin/generate_yukawa \
      --side 32 --lattice 2.7 --screening 1.0 --rcut $rcut \
      --jitter 0.05 --seed 42 --diag 1.0 \
      --out matrices/fill32k_rcut${rcut}.mtx
    mpirun -np 1 NTPoly/Build/bin/yukawa_gemm \
      --input matrices/fill32k_rcut${rcut}.mtx --convert --out "$MATRIX"
    rm -f matrices/fill32k_rcut${rcut}.mtx
  fi
  echo "=== rcut=$rcut: $(du -h "$MATRIX" | cut -f1) ==="

  # CPU benchmark (best config: np=36 omp=4)
  echo "--- CPU np=36 omp=4 ---"
  OMP_NUM_THREADS=4 mpirun --oversubscribe -np 36 \
    NTPoly/Build/bin/yukawa_gemm \
    --input "$MATRIX" --mode fixed --rounds 5 --warmup 2 \
    --out "results/fill_cpu_rcut${rcut}.json" 2>&1
  echo "cpu_exit=$?"

  # GPU benchmark (cuSPARSE)
  echo "--- GPU cuSPARSE ---"
  ./bench_cusparse "$MATRIX" --rounds 5 --warmup 2 --threshold 1e-4 \
    > "results/fill_gpu_rcut${rcut}.json" 2> "results/fill_gpu_rcut${rcut}.log"
  cat results/fill_gpu_rcut${rcut}.log
  echo "gpu_exit=$?"
done

echo "=== done ==="
ls -la results/fill_*.json
date
