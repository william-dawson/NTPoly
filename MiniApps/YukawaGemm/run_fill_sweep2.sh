#!/bin/bash
# GPU fill sweep: rcut=15 and rcut=20 (the ones that failed/crashed before)
set -euo pipefail

WORK="$HOME/agent/work/ntpoly-gemm"
cd "$WORK"
module load nvhpc-hpcx

export UCX_TLS=sm,sysv,posix
export UCX_NET_DEVICES=""

mkdir -p matrices results

# --- rcut=15: GPU only (CPU data already collected) ---
M15="matrices/fill32k_rcut15.bin"
echo "=== rcut=15 GPU ==="
./bench_cusparse "$M15" --rounds 5 --warmup 2 --threshold 1e-4 \
  > results/fill_gpu_rcut15.json 2> results/fill_gpu_rcut15.log
cat results/fill_gpu_rcut15.log
echo "gpu15_exit=$?"

# --- rcut=20: generate + CPU + GPU ---
M20="matrices/fill32k_rcut20.bin"
if [ ! -f "$M20" ]; then
  echo "=== generating rcut=20 ==="
  NTPoly/Build/bin/generate_yukawa \
    --side 32 --lattice 2.7 --screening 1.0 --rcut 20 \
    --jitter 0.05 --seed 42 --diag 1.0 \
    --out matrices/fill32k_rcut20.mtx
  mpirun -np 1 NTPoly/Build/bin/yukawa_gemm \
    --input matrices/fill32k_rcut20.mtx --convert --out "$M20"
  rm -f matrices/fill32k_rcut20.mtx
fi
echo "=== rcut=20: $(du -h "$M20" | cut -f1) ==="

echo "--- CPU np=36 omp=4 ---"
OMP_NUM_THREADS=4 mpirun --oversubscribe -np 36 \
  NTPoly/Build/bin/yukawa_gemm \
  --input "$M20" --mode fixed --rounds 3 --warmup 1 \
  --out "results/fill_cpu_rcut20.json" 2>&1
echo "cpu20_exit=$?"

echo "--- GPU cuSPARSE ---"
./bench_cusparse "$M20" --rounds 5 --warmup 2 --threshold 1e-4 \
  > results/fill_gpu_rcut20.json 2> results/fill_gpu_rcut20.log
cat results/fill_gpu_rcut20.log
echo "gpu20_exit=$?"

echo "=== done ==="
date
