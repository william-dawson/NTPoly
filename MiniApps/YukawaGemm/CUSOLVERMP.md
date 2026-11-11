# cuSOLVERMp density-matrix benchmark

`cusolvermp_bench.cu` is a standalone MPI/CUDA benchmark of FP64 distributed
`cusolverMpSyevd` followed by `cublasMpGemm` for
`D = Q(:,1:nelec) Q(:,1:nelec)^T`. It computes all eigenvectors.
The occupation convention is trace(D)=nelec, without a spin factor.
This compares with ordinary TRS2 using identity overlap. Generalized
eigenproblems and nonidentity overlap are outside this benchmark.

Build on RIKYU after `module load nvhpc-hpcx`:

```sh
SDK=/shared/software/hpc_sdk/Linux_aarch64/26.3
MATH="$SDK/math_libs"
NCCL="$SDK/comm_libs/13.1/nccl"
MPI_FLAGS=""
for inc in $(mpicxx --showme:incdirs); do MPI_FLAGS="$MPI_FLAGS -I$inc"; done
nvcc -c -O3 -std=c++17 -arch=sm_100 $MPI_FLAGS \
  -I"$MATH/include" -I"$NCCL/include" cusolvermp_bench.cu -o bench.o
mpicxx bench.o -o cusolvermp_bench -L"$MATH/lib64" -L"$NCCL/lib" \
  -L"$SDK/cuda/lib64" -lcusolverMp -lcublasmp -lnccl -lcudart
export LD_LIBRARY_PATH="$MATH/lib64:$NCCL/lib:$SDK/cuda/lib64:${LD_LIBRARY_PATH:-}"
OMP_NUM_THREADS=1 mpirun -np 4 ./cusolvermp_bench input.bin 2 2 32000 3 512
```

Arguments are input, process rows, process columns, occupied states,
repetitions, and optional block size (default 256). Each node uses one MPI
rank per visible GPU, in local-rank order. Matrices use row-major process
mapping and column-major block-cyclic local storage. All ranks read the sparse
input and import only locally owned entries; no full host dense matrix is
needed for large inputs. The input must contain unique coordinates and both
triangles of a real symmetric matrix in NTPoly binary format (20-byte header,
then little-endian column:int32,row:int32,value:float64 records).

The program reuses overwritten A as D and shares the maximum workspace between
the sequential library calls. It queries device and host workspace sizes,
checks free device memory before major allocations, and aborts collectively
on allocation/API/numerical errors. The printed explicit allocation estimate
excludes internal CUDA/NCCL/library allocations and temporary sparse output;
it is not measured peak HBM. No recovery from arbitrary library OOM is claimed.

Each phase is synchronized and reports the maximum elapsed time over ranks:

- `setup`: CUDA/NCCL/library initialization, sparse file reading, descriptors,
  workspace queries and allocations. This is outside repeated solves.
- `import`: restoring the original H from cached sparse device coordinates.
- `evd`: the full eigensolve including eigenvectors.
- `density`: occupied-state GEMM. `solve` is `evd+density`.
- `export`: count/compact entries with abs(D)>1e-4 and transfer distributed
  coordinate records to host, including allocation. No disk output, sorting,
  NTPoly CSR construction or redistribution is included.
- `complete`: import+evd+density+export; validation is separate. This is an
  optimistic comparison with native TRS2 returning a normal NTPoly matrix.

Repetition zero captures first-use costs; report it separately from later
repetitions. A fresh executable invocation initializes a fresh context.
Native TRS2 currently creates contexts per call, so compare complete native
API times as well as resident loop times and state the different reuse scope.
File reading, external input permutation, and disk writes are outside the
native TRS2 API timing. The eigensolver uses the original ordering; TRS2 uses
its fixed load-balancing permutation. Energy and trace are invariant to that
permutation; no full entrywise parity between these two large outputs is claimed.

Checks include Syevd info=0, trace(D), Frobenius norm squared of D, and
Tr(H D) versus the occupied eigenvalue sum, with relative tolerance 1e-7.
Inputs N<=1024 also receive exhaustive CPU checks of H symmetry, HQ-Q Lambda,
Q^T Q-I and the GEMM projector, with residual tolerances 1e-10.
The large checks are sanity checks, not a complete residual bound.
Diagonalization is unpruned internally; TRS2 truncates intermediates at 1e-4.
Expect accuracy differences and report them alongside timing.

Tested installation: HPC SDK 26.3, cuSOLVERMp 0.7.2, cuBLASMp 0.8.0.
The installed cuSOLVERMp header provides neither `cusolverMpMalloc` nor
`cusolverMpBufferRegister`; workspace uses `cudaMalloc`. Newer releases
describe improved Blackwell SYTRD performance with registered workspace.
Results from this installation do not establish the performance of those
newer versions. See NVIDIA's [API documentation](https://docs.nvidia.com/cuda/cusolvermp/usage/functions.html)
and [release notes](https://docs.nvidia.com/cuda/cusolvermp/release_notes/index.html).

The outer experiment repository preserves the exact scripts, job IDs and raw
output in `validation/rikyu-cusolvermp-20261111/`.
