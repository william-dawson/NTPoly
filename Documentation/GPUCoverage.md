# GPU solver coverage and implementation boundaries

The optional backend is built with `NTPOLY_ENABLE_GPU_TRS2=ON`. It currently
provides real and complex FP64 shared products on one process slice. The
specialized resident TRS2 iteration remains real FP64. GPU dispatch is opt-in:

- `NTPOLY_GPU_TRS2=1` selects the resident TRS2 implementation, including
  working-H/basis products, bounds, initial guess and its iterative algebra.
- `NTPOLY_GPU_ALGEBRA=1` selects GPU products through the ordinary distributed
  `MatrixMultiply` API. This reaches the other iterative solver families without
  adding a separate GPU implementation to every solver. Other operations still
  use the ordinary CPU representation between products. This is an operation
  coverage path; whole-solver residency and speedups must be measured separately.
- Both flags can be used together. The resident TRS2 loop remains the specialized
  fast path. Existing rank-policy, memory/index admission and safe CPU fallback
  apply to the shared products. Multiple process slices use the CPU path.

The caller can build the matrix on either processor. Real imports borrow accessible
CSC arrays synchronously and create the backend CSR representation. Multiple
local blocks and root coalescing still involve CPU composition/gathering. There
is no pointer-address cache: changed values at the same address are re-imported.

Complex shared products use native complex accumulation, magnitude-based
pruning of `alpha*A*B`, then addition of `beta*C`. Mixed real/complex inputs
are supported. Coordinate import/export remains a host boundary; row/column
panel payloads use the leased NCCL GPU context. Temporary buffers and panel
sorting are per product. This provides solver coverage, not the residency or
performance of the specialized real TRS2 loop. Complex TRS2 uses shared GPU
products when `NTPOLY_GPU_ALGEBRA=1`; `NTPOLY_GPU_TRS2=1` alone retains its CPU
path for complex input. Complex admission uses larger payload and wire-count
bounds based on actual input nnz and the buffers this path owns; the ranking
model has not been calibrated with complex workloads. All-rank output is already
partitioned and skips redistribution. The64k sparse benchmark now admits four GPUs.

## Coverage map

| Family | Shared GPU product coverage | Remaining CPU work |
|---|---|---|
| PM, TRS2, TRS4, HPCP, ScaleAndFold | Iteration and basis products | Control, nonresident matrix operations |
| Hotelling inverse, pseudoinverse | Iterative products | Norms, scaling, transpose, increments |
| Square/inverse square roots, general roots | Iterative products, including higher-order polynomials | Control and linear operations |
| Sign, polar | Iterative products | Transpose/conjugation and linear operations |
| Polynomial, Chebyshev, Hermite | Recurrences and factorized products | Coefficient assembly and linear operations |
| Exponential/logarithm, sine/cosine | Powers, recurrences and internal iterative solver products | Scaling, convergence, linear operations |
| CG | A*P and initial residual product | Direct Frobenius dot products and updates |
| Canonical/grand-canonical Fermi expansion | Expansion and density/basis products | Scalar occupation/control and updates |
| Energy density, McWeeny, gap estimate, geometry extrapolation | Their shared products | Control and linear operations |
| Dense functions, dense density/FOE | Products; optional real/complex cuSOLVERMp eigendecomposition | Control and unsupported eigenproblem fallback |
| SVD | Polar/products; optional real/complex cuSOLVERMp Hermitian eigenproblem | Unsupported eigenproblem fallback |
| Cholesky | No product-based acceleration | Existing column/pivot implementation |

Real resident GPU matrices already provide owned copy, scale, thresholded increment,
trace, dot, product and bounds operations. Basic public operations such as scale,
copy, diagonal scaling, Hadamard, transpose, permutation and symmetrization still
run on CPU. Keep these simple until whole-call measurements identify a material
cost or a residency requirement; do not add import/export solely to move a cheap
operation to a GPU.

The optional `NTPOLY_ENABLE_CUSOLVERMP` backend connects real/complex FP64 dense
EigenDecomposition to cuSOLVERMp with `NTPOLY_GPU_EIGEN=1`; dense wrappers and
SVD inherit that eigenproblem path. See [CUSOLVERMP.md](CUSOLVERMP.md) for its
placement, conversion, fallback and validation boundaries. Custom Cholesky
remains on CPU pending discussion.

## Validation scope

The API suite compares49 real or48 complex entry points and18 matrix-product cases on1/2/4
ranks, root/all GPU policies, zero-memory CPU fallback and a two-slice CPU
fallback. It checks full matrix differences, finite payloads, slot lifetime,
context teardown and complex products. Dense/SVD/Cholesky comparisons include
CPU stages and do not establish GPU implementation of those stages. Small
coverage tests are correctness tests, not performance benchmarks.

## Numerical semantics

The product is `prune(alpha*A*B, threshold)` followed by addition of `beta*C`
with the ordinary zero-threshold cancellation rule. Root-coalesced shared
algebra preserves the entire logical dimension, including permutation entries
in padding; its memory/index admission uses that same dimension. The final sum is not
thresholded at the product threshold. The CPU dense local GEMM branch formerly
pruned before scaling, disagreeing with its sparse branch; both now prune the
scaled product.

CG uses direct Frobenius inner products rather than building thresholded Gram
matrices for their traces. This preserves small nonzero residuals, avoids
spurious zero denominators caused by pruning, and removes three products per
iteration. An exactly solved initial residual exits before dividing by zero. The residual
is checked before taking another step; an erased search product is retried
without pruning rather than used as a false zero denominator.

The optional row/column selector of FillMatrixPermutation is tested safely
only when present, retaining the documented row-permutation default.

## Measured general-solver performance

On the64k SPD benchmark I+0.2*H/||H||1, with1e-8 pruning and1e-5 convergence,
shared GPU inverse takes1.420s warm versus1.842s on128 CPU threads; inverse
square root takes1.174s versus1.616s (1.30x/1.38x). Output differences are below
8.91e-16; validated residual norms are about1e-6. This includes per-product
boundaries. It is not the physical-H purification workload. Tiny problems are
correctness tests and generally favor CPU. The generic flag remains opt-in.

## Measured complex multiplication

Matched complete MatrixMultiply calls, including import/export and allocation,
use a diagonal-unitary phase transform to preserve real/complex product patterns.
After removing redundant output redistribution, warm4GPU medians are15.710/37.616ms
real/complex at10k (2.39x),35.020/105.858ms at20k (3.02x), and126.794/551.916ms at64k
(4.35x). The64k complex call is1.10x faster than CPU128 (609.117ms); the smaller
CPU controls are faster than the shared GPU path. One-GPU10k is32.296/84.610ms
(2.62x). All full-matrix errors are below4.83e-16. These are products, not
complete purification calls. This does not establish complex solver speedups.

A warm rank-zero profile measures34.1us real and52.4us complex multiply kernels
(1.54x), but allocation/free API time totals71.11ms across3 real and6 complex
calls, versus12.31ms of all GPU kernels. Those aggregated scopes are distinct.
Owned complex storage, reusable scratch and avoiding redundant conversion are
more useful next steps than PTX tuning. See Finding36 and the benchmark evidence.

## Cleanup direction

1. Keep GPU admission, context leasing, import/export and ownership in the shared
   backend. Solver code should express mathematics and convergence conditions.
2. Move GPU build ownership and implementation from the mini-app into a library
   source directory. The mini-app should link the same library implementation.
   Preserve the CPU-only build and explicit optional dependencies.
3. Separate owned storage/basic operations, resident products and benchmark-only
   legacy kernels in the large CUDA source. Retain ABI and regression coverage
   while moving code; avoid changing algorithms in that cleanup commit.
4. Prefer a small number of explicit ownership boundaries to long-lived pointer
   caches. Retain managed-memory compatibility and fresh-input semantics.
5. Add ordinary build/test integration, bindings and precision/grid coverage.
   Do not call an untested sliced path GPU-supported because its CPU fallback
   succeeds. Complex whole-solver performance and residency remain separate work.
