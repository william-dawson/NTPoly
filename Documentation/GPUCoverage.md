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
panel payloads use the leased NCCL GPU context. Device buffers and Thrust
sorting/scan scratch are retained by that process-grid context and reused.
Input values and patterns are imported afresh, including changes at the same
address. Memory admission can reclaim idle scratch from any local grid before
retrying; cap changes and grid teardown also release it. Scratch reuse preserves
communicators and does not retain matrix contents between operations.
This provides solver coverage, not the residency or
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

## Measured resident multiplication (primary GPU performance scope)

The requested GPU performance scope starts after input placement and ends with
a GPU sparse result. Import/conversion and final host output construction are
excluded. `resident_product_bench` loads once, repeatedly calls the production
GPU product and exports/validates the final result once. Timers include NCCL,
GPU panel assembly/sorting, clear, arithmetic, pruning/counting, compaction,
scalar metadata and completion synchronization. This is a distributed resident
operation, not just an isolated multiply kernel. Warm medians on4 GPUs:

| Size | CPU128 real ms | GPU real ms | CPU128 complex ms | GPU complex ms |
|---|---|---|---|---|
| 10k | 26.253 | 0.536 | 30.637 | 1.271 |
| 20k | 60.388 | 0.803 | 73.577 | 3.020 |
| 64k | 457.110 | 2.995 | 608.403 | 21.486 |

64k multiplication speedups versus matched CPU128 are152.6x real and28.3x
complex. Local clear/arithmetic/prune/compaction after panels exist is2.488ms
real and19.672ms complex; that subset excludes communication/assembly. These
measurements do not establish complete purification speedups or physical HBM
utilization. CPU and GPU algorithms move different amounts of data.

Complex64k profiling identifies11.05ms dense pruning,2.29ms count and3.80ms
packing versus0.39ms arithmetic multiply. PTX/SASS shows vector loads/stores,
28 registers, no spills, and an inlined robust magnitude calculation even on
zero cells. The useful next GPU optimization is avoiding full dense complex
passes through bitmap/touched-entry processing and pruning fusion. Public
input/output conversion is outside the user's current performance target.
See Finding38 and validation/rikyu-resident-products-20261112.

## Measured public-boundary complex multiplication

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

Reusable grid-owned scratch subsequently lowers matched warm complex calls
from34.979/102.223/535.151ms to27.757/84.445/360.017ms at10k/20k/64k (Finding37).
The new complex/real ratios are1.87x/2.46x/2.88x.64k complex is1.71x faster than
matched CPU128 (614.200ms);20k CPU128 remains faster (74.842ms).10k is roughly
tied. Separate datatype processes prevent retained caches from distorting the
primary comparison. Interleaved64k checks also keep all real/complex calls on
GPU by reclaiming idle scratch before admission retry. Teardown is outside warm
product timers. These remain full product measurements, not resident solvers.

## Full real public-call stage diagnostics

`NTPOLY_GPU_ALGEBRA_TIMING=1` records optional real shared-product stages;
`GetGPUProductTimings` returns42 local values for the benchmark to gather after
its public timer. No extra CUDA synchronization or inter-stage MPI barriers
are inserted. See Finding39 and validation/rikyu-product-stages-20261112.

At64k/4GPUs the full real call averages126.664ms over14 warm calls: checks/
context0.501ms, A/B imports6.102ms, product8.195ms, host result108.222ms,
cleanup3.596ms and dispatch/reporting0.047ms. The host result comprises CPU
CSR-to-CSC35.396ms and local-block population72.091ms plus construction/freeing.
The general block splitter performs two additional sparse transposes even for
one local block. These costs remain outside the user's resident GPU target.
The same-allocation retained product is2.992ms. Fresh slots/panel buffers make
the product stage inside the public call8.195ms; subtracting separate scopes
does not isolate conversion alone. Baseline/timer-disabled controls have
126–127ms medians. No algorithm or one-block shortcut is changed here.

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
