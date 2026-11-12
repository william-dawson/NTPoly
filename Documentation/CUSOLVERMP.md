# Optional cuSOLVERMp dense eigendecomposition

Configure with `-DNTPOLY_ENABLE_CUSOLVERMP=ON` and set `NTPOLY_GPU_EIGEN=1`
at runtime. The ordinary distributed `EigenDecomposition` API then tries the
real/complex FP64 GPU backend before the configured CPU LAPACK/EigenExa implementation.
It is independent of `NTPOLY_ENABLE_GPU_TRS2`; both can be enabled together.
Custom Cholesky is unchanged.

The shared entry point also serves dense matrix functions, dense density/Fermi
construction and SVD's Hermitian eigenproblem. Those callers keep their other
operations and their normal Matrix_ps return types. The low-level local
DMatrix LAPACK routines are unchanged.

## Interface and numerical conventions

- Eigenvalues remain ascending and are returned in a distributed diagonal
  matrix. Eigenvectors remain columns of a square distributed matrix.
- `nvals` keeps the lowest requested columns/values, leaving the rest zero.
  An omitted vector output still works. The installed0.7.2 implementation
  computes vectors internally even for that case, then skips vector export.
- The threshold applies separately to returned eigenvalues/eigenvector entries,
  as in the CPU path. Thresholding can reduce orthogonality; use zero threshold
  for residual/orthogonality checks.
- Inputs use the CPU LAPACK upper-triangle convention. The adapter forms
  the symmetric/Hermitian distributed matrix from that triangle, conjugating
  the mirror for complex input and ignoring the diagonal imaginary part.
  Eigenvalues remain real; complex eigenvectors retain both components.
  Vector pruning uses complex magnitude, rather than pruning components.
- Inputs are redistributed from Matrix_ps sparse triplets into block-cyclic
  dense GPU arrays. Eigenvectors are compacted on GPU before host return and
  redistribution into NTPoly's normal layout. These are explicit boundary
  operations, not a claim of resident dense storage or zero-copy conversion.

## Placement and fallback

`NTPOLY_GPU_EIGEN_RANK_POLICY=auto|root|all` controls active ranks. Automatic
mode coalesces dimensions up to16384 to one GPU, reflecting the measured10k
baseline; larger problems use all admitted application ranks. The solver uses
its own balanced2D grid and redistributes back to the caller's grid. A1x4
application grid can therefore use a2x2 dense solver grid.

Real and complex FP64 matrices on one NTPoly process slice are supported.
Multiple process slices retain CPU fallback. Multiple GPUs and
multiple slices are different axes: a2x2x1 application uses four ranks on one
slice;2x2x2 uses eight ranks across two slices.

The default all-rank path permits at most four GPUs on one node with unique
GPU identities. Cross-node NCCL requires the existing explicit
`YUKAWA_NCCL_ALLOW_MULTINODE=1`; that is not a claim that multi-node performance
or reliability has been established for this new entry point.

Capacity admission includes A, Q, eigenvalues, packing storage, solver device
workspace and worst-case compacted vector output, leaving headroom in live
free device memory. `NTPOLY_GPU_EIGEN_MEMORY_MB` optionally lowers the budget;
zero forces CPU fallback. Invalid layouts/capacity decline collectively.
Runtime or numerical failures are errors, not a silent CPU rerun after a
failed GPU collective. Empty-owner layouts are declined: installed0.7.2
returned an incorrect spectrum with info0 for N129/block64/1x4, while the
nonempty and balanced controls passed. No cluster-wide fault is diagnosed.

The backend owns its cuSOLVERMp/NCCL resources and releases them after each
call. Handle/workspace caching is future performance work, rather than an
extra persistent ownership system in this first integration.

## Build and tests

The library owns this backend under Source/GPU. Set CUDA_HOME,
CUSOLVERMP_ROOT and NCCL_HOME when installed outside the configured SDK search
paths. `NTPOLY_CUDA_ARCH` defaults to `sm_100` for RIKYU and is configurable.
The CPU-only build keeps all CUDA dependencies optional.

`eigen_gpu_test` compares CPU/GPU eigenvalues and spectral projectors, checks
residuals and orthogonality, partial/zero/oversized nvals, omitted vectors,
thresholding, changed inputs, upper-triangle semantics, complex Hermitian matrices and sliced
fallback, and dense callers. Tests use the ordinary public API.
`eigen_gpu_bench` measures complete public calls; spectrum trace/Frobenius
fingerprints and repeated-output checks are outside its timer. A values-only
benchmark excludes host eigenvector export and is not a complete density build.

The RIKYU toolchain exports BLAS/LAPACK as transitive library dependencies,
so a separately configured CMake client can link the installed target and
retain the CPU fallbacks. Configure its install prefix explicitly.
