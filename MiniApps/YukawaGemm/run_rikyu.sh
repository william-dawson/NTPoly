#!/bin/bash
# NTPoly YukawaGemm CPU baseline on RIKYU (RIKEN AI4S / GB200).
# Builds NTPoly (Fortran-only, nvfortran + HPC-X MPI) and runs the
# fixed-mode sparse GEMM benchmark sweep. Results -> results/*.json
set -euo pipefail

WORK="$HOME/agent/work/ntpoly-gemm"
cd "$WORK"

module load nvhpc-hpcx

echo "=== environment ==="
date
which mpif90 mpicc nvfortran cmake make mpirun
nvfortran --version 2>&1 | head -1 || true

echo "=== unpack source ==="
if [ ! -d NTPoly ]; then
  tar xzf ntpoly-bundle.tar.gz
fi

echo "=== BLAS/LAPACK probe ==="
# The toolchain needs a working -llapack -lblas; probe where it resolves.
cat > /tmp/blas_probe_$$.f90 <<'EOF'
program p
  implicit none
  real(8) :: a(2,2), b(2,2), c(2,2)
  a = 1.0d0; b = 2.0d0; c = 0.0d0
  call dgemm('N', 'N', 2, 2, 2, 1.0d0, a, 2, b, 2, 0.0d0, c, 2)
  print *, 'dgemm ok', sum(c)
end program p
EOF
BLASFLAGS=""
if nvfortran -o /tmp/blas_probe_$$.x /tmp/blas_probe_$$.f90 -llapack -lblas 2>/tmp/blas_probe_$$.err; then
  BLASFLAGS="-llapack -lblas"
  echo "plain -llapack -lblas works"
else
  echo "plain link failed:"
  head -5 /tmp/blas_probe_$$.err
  SDKLIB="$(dirname "$(dirname "$(which nvfortran)")")/lib"
  echo "trying SDK lib dir: $SDKLIB"
  ls "$SDKLIB" 2>/dev/null | grep -iE 'blas|lapack' | head || true
  if ls "$SDKLIB"/liblapack.* >/dev/null 2>&1; then
    BLASFLAGS="-L$SDKLIB -llapack -lblas"
  else
    echo "FATAL: no BLAS/LAPACK found"
    exit 2
  fi
fi
rm -f /tmp/blas_probe_$$.f90 /tmp/blas_probe_$$.x /tmp/blas_probe_$$.err
echo "BLASFLAGS=$BLASFLAGS"

echo "=== generate toolchain ==="
cat > NTPoly/Targets/rikyu.cmake <<EOF
set(CMAKE_SYSTEM_NAME Linux)
set(CMAKE_C_COMPILER mpicc)
set(CMAKE_Fortran_COMPILER mpif90)
set(CMAKE_CXX_COMPILER mpicxx)
set(TOOLCHAIN_LIBS "$BLASFLAGS")
set(F_TOOLCHAINFLAGS_RELEASE "-O3 -mp")
set(CXX_TOOLCHAINFLAGS_RELEASE "-O3")
set(F_TOOLCHAINFLAGS_DEBUG "-O0")
set(CXX_TOOLCHAINFLAGS_DEBUG "-O0")
EOF

echo "=== build ==="
mkdir -p NTPoly/Build
cd NTPoly/Build
cmake -DCMAKE_TOOLCHAIN_FILE=../Targets/rikyu.cmake \
      -DCMAKE_BUILD_TYPE=Release -DFORTRAN_ONLY=Yes ..
make -j8
cd "$WORK"

mkdir -p results

echo "=== smoke: tiny matrix (MatrixMarket path), np=2 ==="
export OMP_NUM_THREADS=1
mpirun -np 2 NTPoly/Build/bin/yukawa_gemm \
  --input NTPoly/MiniApps/YukawaGemm/yukawa_sc_4096.mtx \
  --mode fixed --rounds 4 --warmup 2 \
  --out results/tiny_fixed_np2.json

echo "=== sweep: small (32768, binary input), fixed mode ==="
for np in 1 2 4 8 16; do
  echo "--- np=$np ---"
  mpirun -np "$np" NTPoly/Build/bin/yukawa_gemm \
    --input NTPoly/MiniApps/YukawaGemm/yukawa_sc_32768.bin \
    --mode fixed --rounds 8 --warmup 3 \
    --out "results/small_fixed_np${np}.json"
done

echo "=== done ==="
ls -la results/
date
