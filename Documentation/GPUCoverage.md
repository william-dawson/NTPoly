# GPU solver coverage and implementation boundaries

The optional backend is built with `NTPOLY_ENABLE_GPU_TRS2=ON`. It currently
provides real FP64 products on one process slice. GPU dispatch is opt-in:

- `NTPOLY_GPU_TRS2=1` selects the resident TRS2 implementation, including
  working-H/basis products, bounds, initial guess and its iterative algebra.
- `NTPOLY_GPU_ALGEBRA=1` selects GPU products through the ordinary distributed
  `MatrixMultiply` API. This reaches the other iterative solver families without
  adding a separate GPU implementation to every solver. Other operations still
  use the ordinary CPU representation between products. This is an operation
  coverage path; whole-solver residency and speedups must be measured separately.
- Both flags can be used together. The resident TRS2 loop remains the specialized
  fast path. Existing rank-policy, memory/index admission and safe CPU fallback
  apply to the shared products. Complex and sliced products use the CPU path.

The caller can build the matrix on either processor. Imports borrow accessible
CSC arrays synchronously and create the backend CSR representation. Multiple
local blocks and root coalescing still involve CPU composition/gathering. There
is no pointer-address cache: changed values at the same address are re-imported.

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
| Dense functions, dense density/FOE | Products; optional real cuSOLVERMp eigendecomposition | Control and unsupported eigenproblem fallback |
| SVD | Polar/products; optional real cuSOLVERMp Hermitian eigenproblem | Unsupported eigenproblem fallback |
| Cholesky | No product-based acceleration | Existing column/pivot implementation |

Resident GPU matrices already provide owned copy, scale, thresholded increment,
trace, dot, product and bounds operations. Basic public operations such as scale,
copy, diagonal scaling, Hadamard, transpose, permutation and symmetrization still
run on CPU. Keep these simple until whole-call measurements identify a material
cost or a residency requirement; do not add import/export solely to move a cheap
operation to a GPU.

The optional `NTPOLY_ENABLE_CUSOLVERMP` backend now connects real FP64 dense
EigenDecomposition to cuSOLVERMp with `NTPOLY_GPU_EIGEN=1`; dense wrappers and
SVD inherit that eigenproblem path. See [CUSOLVERMP.md](CUSOLVERMP.md) for its
placement, conversion, fallback and validation boundaries. Custom Cholesky
remains on CPU pending discussion.

## Validation scope

The API suite compares49 entry points and18 matrix-product cases on1/2/4
ranks, root/all GPU policies, zero-memory CPU fallback and a two-slice CPU
fallback. It checks full matrix differences, finite payloads, slot lifetime,
context teardown and complex fallback. Dense/SVD/Cholesky comparisons include
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
   Do not call an untested complex or sliced path GPU-supported because its
   CPU fallback succeeds.
