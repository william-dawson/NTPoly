# NTPoly GPU Port — Findings Log

Tracking key discoveries and decisions for the NTPoly GPU port to RIKYU (GB200).

---

## Finding 1: NTPoly's local GEMM uses a dense accumulator (memory bottleneck)

**Date:** 2026-11-03
**File:** `Source/Fortran/dense_includes/ConstructMatrixMemoryPool.f90`

The local sparse GEMM kernel (`SMatrixAlgebraModule` → `GemmMatrix_lsr` →
`SparseBranch_lsr` → `MultiplyBlock_lsr` + `PruneList_lsr`) uses a memory pool
that allocates **five dense arrays** sized `columns × rows` per rank:

```fortran
ALLOCATE(this%pruned_list(columns * rows))        ! INTEGER   (4 bytes)
ALLOCATE(this%value_array(columns, rows))         ! REAL(NTREAL) (8 bytes)
ALLOCATE(this%dirty_array(columns, rows))        ! LOGICAL   (4 bytes)
ALLOCATE(this%hash_index(columns, rows))          ! INTEGER   (4 bytes)
ALLOCATE(this%inserted_per_bucket(columns, rows)) ! INTEGER   (4 bytes)
```

Total: ~24 bytes × (N/√np)² per rank. This is the dominant memory cost and
limits matrix size on a single node:

| Matrix   | N        | np=4    | np=36  | np=72  | np=144 |
|----------|----------|--------:|-------:|-------:|-------:|
| 64k      | 64,000   | 2.7 GB  | 0.3 GB | 0.1 GB | 0.03 GB |
| 175k     | 175,616  | 162 GB  | 20 GB  | 10 GB  | 5 GB   |
| 885k     | 884,736  | 4.1 TB  | 510 GB | 130 GB | 65 GB  |

Also: `columns * rows` is computed as 32-bit INTEGER, which overflows at
local block > ~46k×46k (np=4 on 175k).

**Implication:** The dense accumulator is the primary porting bottleneck. Any
GPU port must either keep it (limited by HBM) or replace it with a sparse
accumulator.

---

## Finding 2: CPU full-node baseline (64k matrix)

**Date:** 2026-11-04
**Job:** 277784
**File:** `results/rikyu-fullnode/SUMMARY.md`

Best CPU config: **np=36 omp=4 at 128 ms median (105 ms min)** on sc 40³ =
64,000 atoms. Trace = 64287.048350011 (bit-identical across all 10 configs).

- Best pure-MPI: np=72 at 153 ms median
- np=144 slower (200 ms) — communication overhead at this size
- Hybrid wins because OpenMP threads parallelize the local block multiply
  without MPI communication
- Scaling is sublinear (4→72 ranks = 18× ranks, only 3.6× speedup)

175k and 885k matrices OOM on CPU at all single-node configs.

---

## Finding 3: cuSPARSE is 9× faster at very low fill (0.06%)

**Date:** 2026-11-04
**Jobs:** 280288, 280289, 280291
**File:** `results/rikyu-fullnode/GPU_VS_CPU.md`

Standalone cuSPARSE SpGEMM benchmark (`bench_cusparse.cu`) on the same
Yukawa matrices:

| Matrix   | GPU (cuSPARSE) | CPU best   | Speedup |
|----------|---------------:|-----------:|--------:|
| 64k      | 13.9 ms        | 128 ms     | 9.2×    |
| 175k     | 36.6 ms        | OOM        | —       |
| 885k     | 185 ms         | OOM        | —       |

GPU eliminates the memory bottleneck entirely (198 GB HBM). Trace matches.
GPU achieves ~0.85 ms per million output nonzeros — nearly constant throughput.

---

## Finding 4: cuSPARSE speedup collapses at realistic DFT fill

**Date:** 2026-11-04
**Jobs:** 281261, 281266
**File:** `results/rikyu-fullnode/FILL_SWEEP.md`

**This is the critical finding.** Real DFT density matrices have ~1-5% fill.
A fill sweep (varying rcut on sc 32³ = 32,768 atoms) shows:

| Input fill | Input nnz/row | CPU median | GPU (cuSPARSE) | Speedup |
|-----------:|--------------:|-----------:|---------------:|--------:|
| 0.12%      | 40            | 64.6 ms    | 7.0 ms         | 9.2×    |
| 0.58%      | 189           | 234 ms     | 118.6 ms       | 2.0×    |
| 1.83%      | 599           | 553 ms     | **CRASH**      | —       |
| 4.0%       | 1314          | 1747 ms    | **CRASH**      | —       |

**Root cause:** cuSPARSE computes the **full product** (no thresholding
during multiply). At high fill, the full product approaches dense — cuSPARSE's
internal data structures can't handle it. All four algorithms (DEFAULT,
ALG1, ALG2, ALG3) return `CUSPARSE_STATUS_INSUFFICIENT_RESOURCES`.

NTPoly's CPU code keeps working because it thresholds **during** the multiply
— the output stays at ~336 nnz/row regardless of input fill. This
thresholding-during-multiply is the key algorithmic feature that cuSPARSE
lacks.

**Decision: cuSPARSE is NOT viable for the GPU port.** Must port NTPoly's
hash-accumulator-with-thresholding algorithm to the GPU via OpenACC.

---

## Finding 5: NTPoly's GEMM algorithm structure

**Date:** 2026-11-04
**Files:** `Source/Fortran/sparse_includes/{MultiplyBlock,PruneList,SparseBranch,GemmMatrix}.f90`

The local sparse GEMM has three phases:

### Phase 1: MultiplyBlock (the inner loop)
```fortran
DO II = 1, matAT%columns           ! outer: columns of A^T (= rows of A)
  DO AA = 1, elements_per_inner_a  ! middle: nnz in column II of A^T
    val_a = matAT%values(...)
    idx_a = matAT%inner_index(...)  ! row index in B^T
    DO BB = 1, elements_per_inner_b ! inner: nnz in row idx_a of B^T
      idx_b = matBT%inner_index(...)
      val_b = matBT%values(...)
      ! Hash-based scatter-add into dense accumulator:
      IF (.NOT. dirty_array(idx_b, II)) THEN
        dirty_array(idx_b, II) = .TRUE.
        ! Insert idx_b into hash bucket
      END IF
      value_array(idx_b, II) += val_a * val_b
    END DO
  END DO
END DO
```

The outer loop over `II` (columns) is **independent** — each column's
accumulation is independent. This is the parallelization target.

### Phase 2: PruneList (threshold + extract)
```fortran
DO RII = 1, mat_c_rows              ! rows
  DO CII = 1, num_hash_buckets      ! hash buckets
    DO HII = 1, inserted_per_bucket(CII, RII)
      working_col = hash_index(HII + ..., RII)
      working_value = value_array(working_col, RII)
      value_array(working_col, RII) = 0  ! clear
      dirty_array(working_col, RII) = .FALSE.
      IF (ABS(alpha * working_value) > threshold) THEN
        pruned_list(PII) = ...  ! keep
        PII = PII + 1
      END IF
    END DO
  END DO
END DO
```

Also parallelizable over rows.

### Phase 3: Sort + construct CSR
```fortran
CALL SortTripletList(unsorted_pruned_list, ..., sorted_pruned_list, bubble_in=.TRUE.)
CALL ConstructMatrixFromTripletList(matAB, sorted_pruned_list, ...)
```

The GemmMatrix wrapper also has a sparsity-based dispatch: if both inputs
are >10% dense, it calls DenseBranch (dense BLAS GEMM) instead.

### Parallelization strategy for OpenACC
- **MultiplyBlock:** `!$acc parallel loop independent` over `II` (columns).
  Each column accumulates into its own slice of the dense arrays — but the
  arrays are `value_array(idx_b, II)`, so column `II` writes to column `II`
  only. **No conflicts between columns.** This is embarrassingly parallel.
- **PruneList:** `!$acc parallel loop independent` over `RII` (rows).
  Each row extracts its own entries. No conflicts.
- **Memory:** Use `!$acc data copyin(matA, matB) create(memorypool) copyout(matC)`
  or `-gpu=managed` for zero-copy. The dense accumulator lives on GPU.

### The data flow
```
matA (CSR) → TransposeMatrix → matAT (CSR, columns = rows of A)
matB (CSR) → TransposeMatrix → matBT (CSR, columns = columns of B)
                           ↓
              MultiplyBlock(matAT, matBT, memorypool)
                           ↓
              memorypool (dense hash accumulator on GPU)
                           ↓
              PruneList(memorypool, threshold, matAB)
                           ↓
                       matAB (CSR)
```

---

## Finding 6: OpenACC port — correct but slow (dense accumulator bottleneck)

**Date:** 2026-11-04
**Jobs:** 281314, 281318
**Files:** `Source/Fortran/sparse_includes/MultiplyBlock.f90`, `PruneList.f90`,
`Source/Fortran/SMatrixAlgebraModule.F90`, `dense_includes/ConstructMatrixMemoryPool.f90`

### What was done

Ported NTPoly's local sparse GEMM kernel to GPU via OpenACC directives:
- **MultiplyBlock**: `!$acc parallel loop independent` over II (columns) —
  each II writes to an independent column of the 2D accumulator arrays
- **PruneList**: `!$acc parallel loop independent` over RII (rows) —
  atomic capture for the shared PII counter
- Removed `PURE` from 6 subroutines (OpenACC compute constructs not allowed in PURE)
- Fixed INT32 overflow in `pruned_list(columns * rows)` allocation (use NTLONG)
- Built with `-acc -gpu=mem:managed -gpu=cc100 -Minfo=accel`
- Compiler confirms: II loop → `gang, vector(128)` on GPU, inner loops → `seq`

### Results — correct at ALL fill levels (unlike cuSPARSE)

**64k matrix:**

| Config          | Median   | Trace          |
|-----------------|----------|----------------|
| GPU OpenACC np=1 | 5600 ms  | 64287.048350011 ✓ |
| GPU OpenACC np=4 | 3581 ms  | 64287.048350011 ✓ |
| CPU np=4 omp=1  | 550 ms   | 64287.048350011 ✓ |
| CPU np=36 omp=4 | 128 ms   | 64287.048350011 ✓ |
| cuSPARSE        | 13.9 ms  | 64287.048350011 ✓ |

**32k fill sweep (GPU np=1 vs CPU np=36 omp=4):**

| rcut | Fill  | GPU OpenACC | CPU (np=36) | cuSPARSE  | GPU/CPU |
|------|-------|------------|------------|-----------|---------|
| 6    | 0.12% | 724 ms     | 64.6 ms    | 7.0 ms    | 11.2× slower |
| 10   | 0.58% | 2175 ms    | 234 ms     | 118.6 ms  | 9.3× slower |
| 15   | 1.83% | 3134 ms    | 553 ms     | **CRASH** | 5.7× slower |
| 20   | 4.0%  | 6348 ms    | 1747 ms    | **CRASH** | 3.6× slower |

### Why the GPU is slow

The OpenACC port is correct and handles all fill levels (where cuSPARSE crashes
at 1.8%+). But it's 3.6–11× slower than the CPU. Three root causes:

1. **Dense accumulator = random memory access.** The memory pool allocates
   `columns × rows` dense arrays (~24 bytes/element). For 32k at np=1, that's
   24 GB of random-access memory. GPU HBM has ~100 MB L2 cache — the working
   set is 240× larger. Each scatter-write to `value_array(idx_b, II)` is a
   random access into a 24 GB array → guaranteed L2 miss → ~150 ns latency.

2. **Low GPU occupancy.** The II loop has 32k parallel iterations (one per
   column). The compiler uses `gang, vector(128)` → 32k/128 = 256 gangs. On
   152 SMs, that's ~1.7 gangs/SM × 128 threads = ~219 threads/SM (max is ~2048).
   The inner AA/BB loops are `seq` (each thread processes ~1600 product terms
   sequentially). With only 219 threads/SM, the GPU can't hide the 150 ns
   memory latency.

3. **Managed memory page migration.** With `-gpu=mem:managed`, pages migrate
   between CPU and GPU. PruneList's sort/CSR construction runs on CPU, pulling
   the accumulator data back from GPU, then the next MultiplyBlock round pushes
   it back to GPU. This thrashing adds overhead.

### The gap narrows at higher fill

The GPU/CPU ratio improves from 11.2× at rcut=6 to 3.6× at rcut=20. At higher
fill, the CPU's hash accumulator has more collisions and more sequential
overhead, while the GPU's dense accumulator access pattern doesn't change
(the array size is the same regardless of fill). But the GPU is still slower
at all tested fill levels.

### Next optimization steps

1. **Flatten the computation** — pre-enumerate all (II, idx_b, val) product
   triplets, then process them in parallel with atomic additions. This gives
   millions of independent operations (vs 32k columns), dramatically increasing
   GPU occupancy.

2. **Replace the dense accumulator** — use a per-column sparse hash table
   instead of dense `columns × rows` arrays. This reduces memory usage and
   improves cache locality. Each column's hash table only stores the nonzeros
   it actually generates.

3. **Avoid managed memory** — use explicit `!$acc data` regions to keep the
   accumulator on GPU and avoid page migration.

4. **Move sort/CSR to GPU** — keep the entire pipeline on GPU to avoid
   CPU↔GPU data transfers between MultiplyBlock and PruneList.

---

## Finding 7: CPU overhead eliminated — sort fix + INT32 overflow fix + direct path

**Date:** 2026-11-04
**Jobs:** 281331 (profile), 281497 (sort profile), 281505 (final)
**Files:** `Source/Fortran/triplet_includes/SortTripletList.f90`,
`Source/Fortran/sparse_includes/PruneList.f90`, `MiniApps/YukawaGemm/yukawa_gemm.f90`

### Root cause: 93% of GPU time was CPU overhead, not GPU compute

OpenACC profiling (`PGI_ACC_TIME=1`, job 281331) on 64k np=1:

| Phase                  | Per-round time |
|------------------------|---------------:|
| MultiplyBlock (GPU)    | 336 ms         |
| PruneList (GPU)         | 31 ms          |
| **GPU compute total**  | **367 ms**     |
| **Measured gemm total**| **5,447 ms**   |
| **CPU overhead**        | **5,080 ms (93%)** |

Instrumented PruneList with `system_clock` (job 281497) to find where the
5,080 ms CPU overhead went:

| Phase                  | 32k (2.4M nnz) | 64k (4.8M nnz) |
|------------------------|---------------:|---------------:|
| GPU prune kernels      | 16 ms          | 35 ms          |
| copy (GPU→CPU)         | 6 ms           | 13 ms          |
| **sort (bubble sort)** | **640 ms**     | **5,020 ms**   |
| CSR construction        | 5 ms           | 12 ms          |

The **bubble sort** in `SortTripletList` was 97% of the CPU overhead. After
a bucket sort by column, it scanned the **entire list** (2.4M–4.8M entries)
up to 76 times (max bucket size), each pass touching every element —
guaranteed L3 cache thrashing.

### Fix 1: Insertion sort replaces bubble sort

Replaced the O(max_bucket × nz) bubble sort with a **per-bucket insertion
sort**. Each column bucket is ~75 entries (nz/column), which fits entirely in
L1 cache. The insertion sort does O(bucket_size²) work per bucket, but the
constant factor is tiny since the data is already in L1.

Result on 32k rcut=6: **716 ms → 105 ms (6.8×)**. But 64k didn't improve —
the sort dispatch condition had an INT32 overflow.

### Fix 2: INT32 overflow in sort dispatch

The dispatch condition `list_length .GT. matrix_rows*matrix_columns * 0.1`
computes `matrix_rows * matrix_columns` as 32-bit INTEGER. For 64k:
`64000 * 64000 = 4.1B > INT32_MAX (2.1B)` → **overflow** → negative value →
condition true → `SortDenseTripletList` called instead of the bucket sort.

`SortDenseTripletList` allocates and zeroes **48 GB of dense buffers**
(32 GB `value_buffer` + 16 GB `dirty_buffer`, both `matrix_rows × matrix_columns`).
The zeroing alone takes ~5 seconds — that was the entire 64k bottleneck.

Fixed by comparing in REAL: `REAL(list_length) .GT. REAL(matrix_rows)*REAL(matrix_columns)*0.1`.

### Fix 3: Direct path bypasses distributed wrapper

Added `--direct` flag to the mini-app. At np=1, the distributed
`MatrixMultiply_ps` wrapper does gather/compose/transpose on CPU every
round — pure overhead. The `--direct` path calls `GemmMatrix_lsr` directly
(pre-transposing the input once), eliminating per-round wrapper overhead.

### Final results (job 281505)

**64k matrix (the primary benchmark):**

| Config                      | Median   | Min    | Trace              | Speedup |
|-----------------------------|----------|--------|--------------------|--------:|
| GPU np=1 (before fixes)     | 5,513 ms | —      | 64287.048350011 ✓  | 1.0×    |
| GPU np=1 direct (after)     | **536 ms** | 536 ms | 64287.048350011 ✓  | **10.3×** |
| GPU np=1 distributed (after)| 611 ms   | 610 ms | 64287.048350011 ✓  | 9.0×    |
| CPU np=36 omp=4 (baseline)  | 128 ms   | 105 ms | 64287.048350011 ✓  | —       |

**32k fill sweep (direct path, after fixes):**

| rcut | Fill  | Before    | After    | Speedup | Trace                |
|------|-------|-----------|----------|---------|----------------------|
| 6    | 0.12% | 724 ms    | 105 ms   | 6.9×    | 32913.93065162921 ✓  |
| 10   | 0.58% | 2,175 ms  | 307 ms   | 7.1×    | 32914.01527298696 ✓  |
| 15   | 1.83% | 3,134 ms  | 848 ms   | 3.7×    | 32914.01529141428 ✓  |
| 20   | 4.0%  | 6,348 ms  | 3,198 ms | 2.0×    | 32914.01529141477 ✓  |

All traces bit-identical to the CPU baseline. All asymmetry = 0.

### Current state

GPU compute is now ~371 ms (MultiplyBlock ~336 ms + PruneList ~35 ms on 64k).
Total measured is 536 ms — the remaining ~165 ms is CPU↔GPU managed-memory
page migration and CSR construction. The GPU is now **4.2× slower** than the
CPU best (536 ms vs 128 ms), down from 43× slower before the fixes.

The speedup is smaller at higher fill (2× at rcut=20) because the GPU
MultiplyBlock kernel time grows with fill while the CPU sort was already a
smaller fraction at high fill.

### Files changed

- `SortTripletList.f90`: bubble sort → per-bucket insertion sort; INT32
  overflow in dispatch condition fixed (compare in REAL)
- `PruneList.f90`: profiling instrumentation removed (clean version restored)
- `yukawa_gemm.f90`: `--direct` flag added (calls `GemmMatrix_lsr` directly
  at np=1, pre-transposing input once instead of every round)

---

## Finding 8: cuSPARSE hybrid backend — 2.1× faster than CPU

**Date:** 2026-11-04
**Jobs:** 281525–281551 (segfault debug → calling convention fix → final benchmark)
**Files:** `MiniApps/YukawaGemm/cusparse_gemm.cu`, `MiniApps/YukawaGemm/CMakeLists.txt`,
`MiniApps/YukawaGemm/yukawa_gemm.f90`

### Approach

Compute the full sparse product C = A*A with cuSPARSE SpGEMM on GPU, then
threshold on CPU. NTPoly's Matrix_lsr stores CSC (column pointers, row
indices). Passing these arrays to cuSPARSE as CSR computes C^T in CSR = C
in CSC for symmetric C = A*A — **no format conversion needed**, only
index-base adjustment (NTPoly inner_index is 1-based, cuSPARSE is 0-based).

The C wrapper (`cusparse_gemm.cu`) is callable from Fortran via ISO_C_BINDING.
CMake compiles the `.cu` with nvcc at configure time and links the `.o` into
the Fortran executable.

### Bugs found and fixed

1. **Fortran-C calling convention mismatch** (segfault at address 0xfa00):
   The Fortran interface declared `n`, `nnz_a`, `threshold` with `VALUE`
   (pass by value), but the C function expected pointers and dereferenced
   them — receiving the literal value 64000 as a pointer address. Fixed by
   changing the C function to accept these by value.

2. **NTPoly outer_index is 0-based, not 1-based**: The code subtracted 1
   from `outer_index` (column pointers), producing -1 for the first entry.
   Fixed by copying `outer_index` directly (no subtraction). Only
   `inner_index` (row indices) is 1-based and needs -1.

3. **GPU threshold kernel produced all zeros**: The custom `threshold_mark_kernel`
   silently failed (all `keep` entries were 0 despite values > threshold).
   Root cause unclear (possibly managed-memory context conflict). Replaced
   with CPU-side thresholding — simpler, reliable, and the data is already
   downloaded to CPU for the compaction step.

4. **cuSPARSE ALG2 silently returns C_nnz=0 at high fill**: Instead of
   returning `CUSPARSE_STATUS_INSUFFICIENT_RESOURCES` (like ALG0/ALG1),
   ALG2 and ALG3 return `CUSPARSE_STATUS_SUCCESS` with 0 nonzeros. Added a
   check: if `C_nnz == 0`, treat as failure and try the next algorithm.

5. **Fallback to distributed path was broken**: When cuSPARSE failed and
   `use_direct` was false, the result went to `matC` (distributed) but the
   code read from the empty `matC_local`. Fixed by copying the distributed
   result to `matC_local` for consistent handling.

### Results (job 281551)

**64k matrix (the primary benchmark):**

| Backend                      | Median   | Min    | vs CPU | Trace              |
|-----------------------------|----------|--------|--------|--------------------|
| CPU np=36 omp=4 (baseline)  | 128 ms   | 105 ms | 1.0×   | 64287.048350011    |
| OpenACC np=1 direct          | 536 ms   | 536 ms | 0.24×  | 64287.048350011 ✓  |
| **cuSPARSE hybrid np=1**     | **61 ms** | **61 ms** | **2.1×** | **64287.048350011 ✓** |

cuSPARSE hybrid is **2.1× faster than the full-node CPU baseline** and
**8.8× faster than the OpenACC port**.

**32k fill sweep:**

| rcut | Fill  | Backend used     | Median   | nnz       | Trace                | Asymmetry |
|------|-------|------------------|----------|-----------|----------------------|-----------|
| 6    | 0.12% | cuSPARSE         | 31 ms    | 2,444,530 | 32913.930651629  ✓   | 1.4e-16   |
| 10   | 0.58% | cuSPARSE         | 235 ms   | 3,747,210 | 32914.015272987  ✓   | 1.8e-16   |
| 15   | 1.83% | OpenACC fallback | 883 ms   | 3,751,584 | 32914.015291414  ✓   | 0         |
| 20   | 4.0%  | OpenACC fallback | 3258 ms  | 3,751,584 | 32914.015291414  ✓   | 0         |

All traces match the CPU/OpenACC baselines. The `--cusparse --direct`
combination gives the best of both worlds: cuSPARSE at low fill (2–9× faster
than OpenACC), automatic fallback to OpenACC at high fill where cuSPARSE
crashes.

### Performance breakdown (64k, estimated from standalone benchmark)

| Phase                          | Time  |
|--------------------------------|------:|
| cuSPARSE SpGEMM (GPU)           | ~13 ms |
| Download full product to CPU    | ~20 ms |
| CPU threshold + compact         | ~10 ms |
| Copy result to Fortran arrays   | ~18 ms |
| **Total**                       | **~61 ms** |

The GPU SpGEMM is only 13 ms — the remaining ~48 ms is data transfer and
CPU-side thresholding. Optimizing the thresholding (GPU kernel or avoiding
the full download) could bring the total closer to 13–20 ms.

### Limitations

- cuSPARSE SpGEMM crashes (`INSUFFICIENT_RESOURCES` or silent C_nnz=0) at
  ≥1.8% fill (rcut≥15 on 32k). The hybrid backend automatically falls back
  to OpenACC when this happens.
- The `cusparse_gemm: all algorithms failed` messages on stderr are noisy
  but informative during fallback.
- Only works at np=1 (single GPU). Multi-GPU would need a different approach.

---

## Finding 9: cuSPARSE hybrid v4 — 16.9× speedup via persistent buffers + warm calls

**Date:** 2026-11-04
**Jobs:** 281609 (alg sweep), 281610 (per-phase profiling), 281613 (final benchmark)
**Files:** `MiniApps/YukawaGemm/cusparse_gemm.cu`, `MiniApps/YukawaGemm/yukawa_gemm.f90`

### What changed (v1 → v4 optimization progression)

| Version | Key change | 64k median | vs CPU |
|---------|-----------|-----------:|------:|
| v1 (Finding 8) | Per-round cudaMalloc/cudaFree everywhere | 61 ms | 2.1× |
| v3 | Persistent A on GPU, skip per-round C→Fortran copy | 18.2 ms | 7.0× |
| **v4** | **Persistent C/work buffers + warm-call (skip size queries)** | **7.58 ms** | **16.9×** |

### Root cause of v1's 61 ms: GPU idle time during cudaMalloc

Profiling with CUDA events revealed that the cuSPARSE SpGEMM GPU event time
was 16.3 ms in v3, but the **actual GPU kernel time was only 6.3 ms**. The
difference was GPU idle time caused by `cudaMalloc` calls inside the event
window:

- `cudaMalloc(d_C_col_ind, 60 MB)` — allocates the full product column indices
- `cudaMalloc(d_C_val, 120 MB)` — allocates the full product values
- `cudaMalloc(d_work1/buf2, variable)` — cuSPARSE internal work buffers

Each `cudaMalloc` of a large buffer can cause a GPU synchronization point,
leaving the GPU idle while the host allocator runs. With 180 MB of per-round
allocations, this added ~10 ms of GPU idle time per round.

### v4 optimizations

1. **Persistent C (full product) buffers** — `d_C_row_ptr`, `d_C_col_ind`,
   `d_C_val` are allocated once and reused. Realloc only if `C_nnz` exceeds
   capacity. Eliminates 180 MB of per-round `cudaMalloc`/`cudaFree`.

2. **Persistent cuSPARSE work buffers** — `d_work1`, `d_work2` (for
   `workEstimation` and `compute` phases) are reused across calls. Realloc
   only if the required size grows.

3. **Warm-call optimization** — After the first successful SpGEMM, the size
   query calls (`cusparseSpGEMM_workEstimation` and `cusparseSpGEMM_compute`
   with `nullptr` buffer) are skipped. The cached buffer sizes are reused
   directly. This saves ~1.5 ms per round (7.8 ms → 6.3 ms compute).

4. **Persistent threshold-pipeline buffers** — `d_keep`, `d_row_count`,
   `d_out_row_ptr`, `d_out_col_ind`, `d_out_val` are all persistent. Only
   realloc if capacity is exceeded.

5. **Persistent host output buffers** — `h_out_row_ptr`, `h_out_col_ind`,
   `h_out_val` avoid per-round `malloc`/`free`.

6. **Algorithm selection** — `CUSPARSE_ALG` env var (0=DEFAULT, 1=ALG1,
   2=ALG2, 3=ALG3, -2=benchmark all). Benchmark mode runs all 4 algorithms
   and prints per-algorithm GPU timing. Default: try DEFAULT first, then
   fallback to ALG1/ALG2/ALG3.

### Algorithm benchmark results (64k, job 281609)

| Algorithm | First call (cold) | Steady-state | Status |
|-----------|----------------:|-------------:|--------|
| DEFAULT   | 30.8 ms         | 6.5 ms       | ✓ (maps to ALG1) |
| ALG1      | 30.7 ms         | **6.5 ms**   | ✓ Best |
| ALG2      | 20.2 ms         | 0.13 ms      | ✗ Silent C_nnz=0 failure |
| ALG3      | 17.6 ms         | 0.13 ms      | ✗ Silent C_nnz=0 failure |

ALG2/ALG3 silently return `C_nnz=0` at this fill level. DEFAULT maps to ALG1
internally. ALG1 is the best algorithm for this problem.

### Per-phase GPU timing (64k steady-state, job 281610)

| Phase                        | GPU time | % of total |
|------------------------------|---------:|-----------:|
| cuSPARSE compute (warm)      | 6.3 ms   | 84%        |
| cuSPARSE copy                | 0.11 ms  | 1.5%       |
| cuSPARSE workEstimation      | ~0.1 ms  | 1.3%       |
| GPU threshold mark + count   | 0.15 ms  | 2%         |
| CPU prefix scan (n elements) | 0.04 ms  | 0.5%       |
| GPU compact                  | 0.50 ms  | 6.7%       |
| D2H download (thresholded)   | 0.32 ms  | 4.3%       |
| **Total**                    | **7.58 ms** | 100%    |

The cuSPARSE compute kernel is 84% of the total — it is the dominant cost and
the theoretical floor for this approach. The threshold pipeline (mark + count +
compact) is only 1.0 ms total.

### Final benchmark results (job 281613)

**64k matrix (primary benchmark, 0.06% fill):**

| Backend                      | Median   | Min    | vs CPU | Trace              | Asymmetry |
|-----------------------------|----------|--------|--------|--------------------|-----------|
| CPU np=36 omp=4 (baseline)  | 128 ms   | 105 ms | 1.0×   | 64287.048350011    | —         |
| OpenACC np=1 direct          | 536 ms   | 536 ms | 0.24×  | 64287.048350011 ✓  | 0         |
| cuSPARSE hybrid v1           | 61 ms    | 61 ms  | 2.1×   | 64287.048350011 ✓  | 1.6e-16   |
| **cuSPARSE hybrid v4**       | **7.58 ms** | **7.56 ms** | **16.9×** | **64287.048350011 ✓** | **1.6e-16** |

**32k fill sweep:**

| rcut | Fill  | Backend      | GPU median | CPU (np=36) | Speedup | nnz       | Trace                | Asymmetry |
|------|-------|-------------|-----------:|------------:|--------:|-----------|----------------------|-----------|
| 6    | 0.12% | cuSPARSE    | 4.05 ms    | 64.6 ms     | **16.0×** | 2,444,530 | 32913.930651629  ✓ | 1.4e-16   |
| 10   | 0.58% | cuSPARSE    | 50.4 ms    | 234 ms      | **4.6×**  | 3,747,210 | 32914.015272987  ✓ | 1.8e-16   |
| 15   | 1.83% | OpenACC     | 850 ms     | 553 ms      | 0.65×  | 3,751,584 | 32914.015291414  ✓ | 0         |
| 20   | 4.0%  | OpenACC     | 3203 ms    | 1747 ms     | 0.55×  | 3,751,584 | 32914.015291415  ✓ | 0         |

All traces match the CPU baseline. cuSPARSE achieves 4.6–16.9× speedup at
low fill (≤0.58%). At high fill (≥1.83%), cuSPARSE fails and the OpenACC
fallback is slower than CPU.

The 32k rcut=10 case (50.4 ms) is slower than rcut=6 (4.05 ms) because the
input has ~6.1M nnz (vs 1.3M), producing a much larger full product that
cuSPARSE must compute and store.

### Bandwidth analysis: why SpGEMM can't hit peak HBM bandwidth

The GB200's B200 GPU has 198 GB HBM3e at ~8 TB/s bandwidth. The user's goal of
"hitting peak memory bandwidth" is **not achievable for SpGEMM** — and here's
the quantitative model showing why:

**64k matrix, cuSPARSE compute = 6.3 ms:**

| Data              | Size (CSR)   |
|-------------------|-------------|
| A row_ptr         | 256 KB      |
| A col_ind         | 10 MB       |
| A val             | 20 MB       |
| A read twice      | 60 MB       |
| C row_ptr         | 256 KB      |
| C col_ind (full)  | 60 MB       |
| C val (full)      | 120 MB      |
| **Total traffic** | **~241 MB** |

| Metric                          | Value       |
|--------------------------------|------------:|
| Streaming BW floor (241 MB / 8 TB/s) | 0.030 ms |
| Compute floor (1.2 GFLOP / 90 TFLOP/s) | 0.013 ms |
| **Actual cuSPARSE compute**     | **6.3 ms**  |
| Gap from streaming floor        | 210×        |
| Achieved bandwidth              | 38 GB/s     |
| **% of peak HBM bandwidth**     | **0.48%**   |

SpGEMM is **not a streaming operation**. It requires:
- **Hash table operations** — each output element is inserted into a hash
  table, with probing and atomic operations. This is latency-bound, not
  bandwidth-bound.
- **Random memory access** — merging row patterns requires indirect indexing
  into A's CSR arrays, causing L2 cache misses.
- **Sequential dependencies** — the symbolic phase (computing C's sparsity
  pattern) must complete before the numeric phase (computing values).

Dense GEMM (cuBLAS) can hit >90% of peak bandwidth because it streams
predictable data. SpGEMM cannot — the irregular access pattern means most of
the 6.3 ms is spent on data structure management, not memory transfers.

**Comparison to dense GEMM:** If the 64k matrix were dense, cuBLAS would
achieve ~7 TB/s (90% of peak) for a 64k × 64k × 64k FP64 GEMM. But at 0.06%
fill, dense GEMM would do 99.94% wasted work — the sparse format is essential
for correctness, at the cost of bandwidth efficiency.

### What would it take to go faster?

The cuSPARSE compute kernel (6.3 ms) is 84% of the total time. To go faster:

1. **Custom SpGEMM kernel** — exploit the specific structure of this problem
   (symmetric A, moderate fill, known output pattern in steady-state). A
   custom merge-based kernel with per-row hash tables could potentially be
   faster, but this is a major engineering effort and unlikely to beat
   cuSPARSE's highly tuned implementation.

2. **Exploit symmetry** — C = A*A is symmetric. Computing only the upper
   triangle would halve the output. But cuSPARSE doesn't support triangular
   output, and a custom approach would need a different algorithm entirely.

3. **Symbolic-numeric reuse** — in steady-state purification, the sparsity
   pattern of C is similar to A. If cuSPARSE exposed a "reuse symbolic phase"
   API, we could skip the symbolic computation on subsequent rounds. The
   CUDA 13.1 API does not expose this for ALG1 (ALG2/ALG3 have separate phases
   but fail at our fill level).

4. **Different library** — cuSPARSKIT or GraphBLAS might have different
   trade-offs, but cuSPARSE is generally the fastest on NVIDIA GPUs.

The 7.58 ms total (16.9× over CPU) is likely close to the practical limit for
cuSPARSE-based SpGEMM on this problem. The remaining 1.3 ms of overhead
(threshold pipeline + D2H) could be reduced with pinned memory and better
kernel parallelism, but the 6.3 ms cuSPARSE kernel is the floor.

---

## Finding 10: Multi-GPU scaling — 92× speedup over CPU with 8 GPUs

**Date:** 2026-11-05
**Jobs:** 281962 (np=8), plus login-node runs (np=1/2/4)
**Files:** `cusparse_gemm.cu` (generalized A×B path), `yukawa_gemm.f90` (`--cusparse-mpi` mode)

### Approach: 1D row decomposition

Instead of modifying NTPoly's 3D/2.5D SUMMA (which has a complex task scheduler
with non-blocking MPI and OpenMP), we use a simpler **1D row decomposition**:

1. Each rank reads the full matrix A directly from the NTPoly binary file
   (replicated — 41 MB for 64k, fits easily in 173 GB HBM).
2. Each rank extracts its row-block [start_row : end_row) in CSR format.
3. Each rank calls cuSPARSE SpGEMM for C_local = A_local × A on its own GPU.
4. Each rank thresholds the result on-GPU (same pipeline as v4).
5. Fingerprint computed via MPI_AllReduce (trace, grand_sum, nnz).

**No communication during the multiply** — the computation is embarrassingly
parallel. The only MPI calls are `MPI_Barrier` (synchronization before each
round, ~0.01 ms) and `MPI_Allreduce` for the fingerprint (3 reductions at the
end, outside the timed loop).

### Implementation

**C wrapper (`cusparse_gemm.cu`):**
- `init_gpu()` selects GPU via `OMPI_COMM_WORLD_LOCAL_RANK` / `SLURM_LOCALID`
  environment variable, so each MPI rank gets a different GPU.
- New persistent state for A_local (`d_Aloc_*`) and B (`d_B_*`) matrices.
- `ensure_aloc_uploaded()` and `ensure_b_uploaded()` — upload once, reuse
  across rounds (same pattern as `ensure_a_uploaded()` in v4).
- `run_spgemm_general()` — refactored from `run_spgemm_one()` to accept two
  arbitrary descriptors (A and B) instead of using `g_desc_A` for both.
- `do_spgemm()` — extracted algorithm selection (env var `CUSPARSE_ALG`).
- `do_threshold_pipeline()` — extracted threshold pipeline (mark + count +
  scan + compact + download), shared by both A×A and A_local×B paths.
- New entry point `cusparse_gemm_ab_threshold_()` — accepts A_local (m×n) and
  B (n×n) separately, uploads both persistently, calls `do_spgemm` with
  `g_desc_Aloc` and `g_desc_B`, then `do_threshold_pipeline` with `c_rows=m`.

**Fortran mini-app (`yukawa_gemm.f90`):**
- New `--cusparse-mpi` flag.
- `ReadBinaryCSR()` subroutine — reads NTPoly binary format directly (header:
  3×int32 + int64, then (col, row, val) triplets), builds CSR arrays
  (row_ptr 0-based, col_ind 1-based — matching NTPoly's format).
- 1D row decomposition: `start_row = rank * N/np`, `end_row = (rank+1) * N/np`
  (last rank gets the remainder).
- Extracts A_local row-block by slicing `row_ptr`, `col_ind`, `val` arrays.
- Fingerprint: local trace (sum of diagonal elements in row-block), local
  grand_sum (sum of all values), then `MPI_Allreduce(..., MPI_SUM, ...)`.

### Scaling results (64k matrix, 0.06% fill)

| np | GPUs | Nodes | median (ms) | min (ms) | vs np=1 | vs CPU (128ms) | efficiency |
|----|------|-------|------------:|---------:|--------:|---------------:|-----------:|
| 1  | 1    | 1     | 7.76        | 7.75     | 1.00×   | 16.5×          | 100%       |
| 2  | 2    | 1     | 4.33        | 4.31     | 1.79×   | 29.6×          | 89.7%      |
| 4  | 4    | 1     | 2.35        | 2.33     | 3.30×   | 54.5×          | 82.6%      |
| 8  | 8    | 2     | 1.39        | 1.37     | 5.58×   | **92.1×**      | 69.8%      |

All traces match the CPU baseline (64287.048350011) to 13 significant figures
(differences at the 13th digit are from floating-point summation order).

### Scaling analysis

The scaling efficiency drops from 90% (np=2) to 70% (np=8). The main causes:

1. **Fixed overhead per rank** — cuSPARSE handle creation, descriptor setup,
   and the threshold pipeline (kernel launches, CPU prefix scan) have fixed
   costs that don't shrink with fewer rows per rank. At np=8, each rank
   processes only 8000 rows (vs 64000 at np=1), so the fixed overhead is a
   larger fraction of the total.

2. **cuSPARSE sublinear scaling** — the SpGEMM compute time drops from 6.3 ms
   (np=1) to ~1.0 ms (np=8, estimated), a 6.3× speedup with 8× fewer rows.
   cuSPARSE has per-call overhead (work estimation, hash table initialization)
   that doesn't scale with matrix size.

3. **Load imbalance** — the 1D row decomposition is not perfectly balanced.
   The nnz/row varies from 70.35 (np=8) to 76.02 (np=1), because outer rows
   of the Yukawa matrix have fewer neighbors than inner rows. The last rank
   (with the remainder) may have slightly different work.

4. **MPI_Barrier synchronization** — the barrier waits for the slowest rank.
   At np=8, the barrier cost is ~0.01 ms (negligible vs 1.39 ms compute).

### Communication vs compute breakdown

The multiply is **embarrassingly parallel** — no data exchange between ranks
during the SpGEMM. The only communication is:

| Communication               | Count/round | Cost       | % of 1.39ms |
|----------------------------|------------:|-----------:|------------:|
| MPI_Barrier (sync)         | 1           | ~0.01 ms   | <1%         |
| MPI_Allreduce (timing)     | 1           | ~0.01 ms   | <1%         |
| MPI_Allreduce (fingerprint) | 3 (end only)| ~0.03 ms   | N/A (outside loop) |

Communication is **not the bottleneck** for this problem size. The scaling is
compute-bound, limited by cuSPARSE's per-call overhead and the threshold
pipeline's fixed costs.

### Comparison to single-GPU `--cusparse` path

| Path              | np=1 median | Difference |
|-------------------|------------:|-----------:|
| `--cusparse` (v4)  | 7.58 ms     | —          |
| `--cusparse-mpi`   | 7.76 ms     | +2.4%      |

The multi-GPU path at np=1 is 2.4% slower than the single-GPU v4 path, because
it uploads A_local and B separately (two uploads vs one) and uses separate
warm flags. This small overhead is acceptable given the scaling benefit.

### What would it take to scale further?

1. **Larger matrices** — at 64k, each GPU processes only 8000 rows at np=8.
   For larger matrices (175k, 885k), each GPU would have more work, improving
   the compute-to-overhead ratio. The 175k/885k matrices OOM on a single GPU,
   but multi-GPU scaling provides aggregate HBM (4×173 GB = 692 GB per node).

2. **Better load balancing** — a nnz-weighted partition (instead of uniform
   row partition) would distribute work more evenly. The current partition
   gives the last rank slightly different nnz/row.

3. **Overlap threshold with next round's compute** — pipeline the threshold
   pipeline of round N with the SpGEMM of round N+1. This would hide the
   ~1.0 ms threshold overhead behind the ~1.0 ms compute.

4. **Reduce D2H transfer** — the fingerprint (trace, grand_sum) could be
   computed on-GPU and reduced via MPI_Allreduce, avoiding the D2H transfer
   of the full thresholded matrix. This would save ~0.5 ms per round at np=8.

5. **cuSPARSE multi-GPU API** — CUDA 13.1's cuSPARSE does not have a native
   multi-GPU SpGEMM. If available in future CUDA versions, it could handle
   the decomposition internally with better load balancing and overlap.

---

## Finding 11: 2D SUMMA + symmetry exploitation — 47× over CPU at np=16

**Date:** 2026-11-05
**Jobs:** 284212 (per-rank diag), 284226 (alg bench), 284253 (np=4 sym),
284260 (np=16 deadlocked → fixed), 284261 (32k np=4), 284329 (np=4 regress),
284330 (np=16), 284336 (np=9), 284337 (175k np=16), 284342 (175k np=1),
284343 (175k np=4)
**Files:** `yukawa_gemm.f90` (2D SUMMA + `--symmetric`), `cusparse_gemm.cu`
(`*nnz_c = 0` on failure)

### Background: why 1D row decomposition was rejected

The 1D row decomposition (Finding 10) achieved 92× speedup at np=8, but was
rejected for production use because:
1. **Redistribution each round kills purification** — NTPoly's purification
   iterates C = f(C) many times. With 1D row decomposition, each round
   redistributes the full matrix, adding O(N²/np) communication.
2. **Replicating the full matrix OOMs at scale** — each rank stores a full
   copy of A (41 MB for 64k, but 116 MB for 175k, 593 MB for 885k). At np=16,
   that's 9.5 GB of replicated A for 885k — wasteful.

The 2D SUMMA approach stores A in a 2D block layout (no full replication),
and the output stays in 2D block layout (no redistribution for purification).

### 2D SUMMA design

Process grid: pr × pc × 1 (rows × cols × 1 slice). Requires square grid
(pr = pc). Each rank stores block (my_row, my_col) of A in CSC format
(NTPoly's native storage).

Per-round computation (for fixed mode, gather/transpose done once):
1. **AllGather** row-block along `row_comm` → A_rows (m_loc × N) in CSC
2. **Bcast** row-block from rank my_col along `col_comm` → B source
3. **CSC→CSR transpose** of A_rows → A_local (m_loc × N) in CSR
4. **B aliasing** — B is A[my_col,:] in CSC (m_loc×N). For symmetric A,
   CSC of m_loc×N = CSR of N×m_loc (the transpose), so B can be passed
   directly as CSR of N×m_loc. No transpose needed.
5. **cuSPARSE SpGEMM**: C = A_local × B (m_loc × kk), then threshold

Key constraints:
- `OMP_NUM_THREADS=1` required — with multiple OpenMP threads, NTPoly creates
  multiple local blocks (`number_of_blocks > 1`), and AllGather only gathers
  block (1,1), missing data.
- NTPoly's `outer_index` (col_ptr) is 0-based, `inner_index` (row_ind) is
  1-based. The GPU kernel subtracts 1 from col_ind at upload time.

### Discovery: cuSPARSE SpGEMM pathology (120× throughput gap)

Per-rank diagnostics on np=4 (2×2 grid) revealed a structural cuSPARSE issue:

| Rank | Grid pos | flop_count | max_row_flop | C_nnz | SpGEMM time |
|------|----------|------------|---------------|-------|-------------|
| 0 | (0,0) diag | 51,894,159 | 2197 | 7,235,332 | 3.4 ms |
| 1 | (0,1) off-diag | 1,493,894 | 881 | 284,710 | **0.49 ms** |
| 2 | (1,0) off-diag | 1,493,894 | 860 | 284,710 | **58 ms** |
| 3 | (1,1) diag | 51,889,465 | 2224 | 7,236,080 | 3.4 ms |

The two off-diagonal blocks have **identical** flop_count (1,493,894) and
similar max_row_flop, yet one takes 0.49 ms and the other 58 ms — a **120×
throughput gap**. All cuSPARSE algorithms (DEFAULT/ALG1/ALG2/ALG3) give the
same 57 ms for the slow block. ALG2/ALG3 fail on ALL blocks at this fill.

The diagonal blocks have 35× more flops but are 17× faster. cuSPARSE
throughput: 15.3M flops/ms (diagonal) vs 26K flops/ms (slow off-diagonal) —
a **600× throughput gap**. This is a structural cuSPARSE issue, not a
workload issue.

### Fix: symmetry exploitation (`--symmetric` flag)

For symmetric A, C = A² is symmetric: C[my_row,my_col] = C[my_col,my_row]^T.
This means:
- **Upper triangle** (my_row < my_col): compute SpGEMM, transpose result,
  send C^T to lower-triangle partner via MPI_Send
- **Lower triangle** (my_row > my_col): skip SpGEMM entirely, receive C^T
  from upper-triangle partner via MPI_Recv
- **Diagonal** (my_row == my_col): compute normally

Partner rank = my_col * pc + my_row (transpose position in grid).

This eliminates the slow 58 ms off-diagonal block entirely. Only 3 SpGEMM
calls per round at np=4 (2 diagonal + 1 upper-triangle), down from 4.

**Result: np=4 from 58 ms → 4.0 ms (14.4× improvement, 32× over CPU).**

### Zero-nnz handling (distant blocks)

At np=16 (4×4 grid), blocks with |my_row - my_col| > 1 have zero product
(flop_count = 0) because the Yukawa potential decays with distance. cuSPARSE
returns C_nnz=0 for these blocks, which our wrapper treats as failure.

This caused a **deadlock**: upper triangle didn't send (treated as failure),
lower triangle waited for recv. Fix:
- Upper triangle **always** sends ct_nnz (even if 0); data arrays only sent
  if ct_nnz > 0
- Lower triangle only receives data arrays if recv_nnz > 0
- Fingerprint checks `cs_ab_nnz_c > 0` before using pointer data
- C wrapper sets `*nnz_c = 0` on failure (prevents stale values)

### Scaling results — 64k matrix

| np | Grid | GPUs | Nodes | Median (ms) | Min (ms) | vs CPU (128ms) | vs np=1 (7.6ms) |
|----|------|------|-------|------------:|---------:|---------------:|----------------:|
| 1  | 1×1  | 1    | 1     | 7.57       | 7.56     | 16.9×          | 1.0×            |
| 4  | 2×2  | 4    | 1     | 4.13       | 4.10     | 31.0×          | 1.83×           |
| 9  | 3×3  | 12   | 3     | 4.18       | 4.17     | 30.6×          | 1.81×           |
| 16 | 4×4  | 16   | 4     | **2.73**   | 2.71     | **46.9×**      | **2.77×**       |

All traces match to 13 significant figures (64287.0483500112–113) —
differences are floating-point round-off from different summation order.

**Note: np=9 is NOT faster than np=4** (4.18 vs 4.13 ms). The bottleneck is
the diagonal SpGEMM, which doesn't scale linearly:

| Grid | Diagonal block | flop_count | SpGEMM time | flops/ms |
|------|---------------|------------|-------------|----------|
| 2×2  | 32k × 32k     | 52M        | 3.4 ms      | 15.3M    |
| 3×3  | 21k × 21k     | 33M        | 3.8 ms      | 8.7M     |
| 4×4  | 16k × 16k     | 24M        | 2.4 ms      | 10.0M    |

cuSPARSE SpGEMM has highly variable throughput (8.7M–15.3M flops/ms) that
doesn't correlate with problem size — the 21k blocks are actually slower than
32k blocks. This is the same pathology seen in the off-diagonal blocks.

### Scaling results — 175k matrix (OOM on CPU)

The 175k matrix (175,616 × 175,616, 7.2M nnz) OOM'd on CPU at all single-node
configs (Finding 2). With 2D SUMMA, it runs on 16 GPUs across 4 nodes:

| np | Grid | GPUs | Nodes | Median (ms) | vs np=1 | Trace |
|----|------|------|-------|------------:|--------:|------:|
| 1  | 1×1  | 1    | 1     | 20.60      | 1.0×    | 176409.89522371584 |
| 4  | 2×2  | 4    | 1     | 10.65      | 1.93×   | 176409.89522371738 |
| 16 | 4×4  | 16   | 4     | **6.18**   | **3.33×** | 176409.89522371633 |

Traces match to 13 significant figures. The 175k scaling is better than 64k
(more work per GPU → lower overhead fraction). At np=16, the 175k matrix
(2.75× larger than 64k) takes only 2.26× longer — superlinear from better
cuSPARSE efficiency on larger blocks.

### Scaling results — 32k matrix

| np | Grid | Median (ms) | Trace |
|----|------|------------:|-------|
| 4  | 2×2  | 2.2         | 32913.93 ✓ |

### Scaling results — 885k matrix (OOM on CPU, OOM on single GPU in standalone)

The 885k matrix (884,736 × 884,736, 37M nnz) is the largest test matrix.
The standalone cuSPARSE benchmark (Finding 3) ran it in 185 ms on one GPU.
The mini-app with v4 persistent buffers runs it in 104.88 ms on one GPU (1.76×
faster than standalone). CPU OOM'd at all single-node configs.

| np | Grid | GPUs | Nodes | Median (ms) | vs np=1 | Trace |
|----|------|------|-------|------------:|--------:|------:|
| 1  | 1×1  | 1    | 1     | 104.88      | 1.0×    | 888767.955692882211 |
| 4  | 2×2  | 4    | 1     | 53.70       | 1.95×   | 888767.955692885611 |
| 16 | 4×4  | 16   | 4     | **28.00**   | **3.75×** | 888767.95569288160 |

Traces match to ~13 significant figures. The scaling improves with matrix
size (np=16: 2.77× for 64k, 3.33× for 175k, 3.75× for 885k) because larger
matrices have more work per GPU, reducing the per-call overhead fraction.

### Communication vs compute

Per-round communication (symmetry exchange):
- Upper→lower: 1 MPI_Send(ct_nnz) + 3 MPI_Send(row_ptr, col_ind, vals)
- Lower←upper: 1 MPI_Recv(ct_nnz) + 3 MPI_Recv(data arrays)
- Data size: ~2-3 MB per pair (285K nnz × 16 bytes)
- At IB XDR (800 Gbps = 100 GB/s): ~0.03 ms per pair — negligible

The MPI_Send/MPI_Recv are blocking but overlap across independent pairs.
The AllGather and Bcast (steps 1-2) happen once for fixed mode, not per round.

### What limits scaling?

1. **cuSPARSE SpGEMM variable throughput** — the diagonal blocks have 8-15M
   flops/ms throughput (600× below peak). Different block sizes hit different
   cuSPARSE performance regimes. The slowest block determines wall time.

2. **Per-call overhead** — cuSPARSE has fixed costs (handle setup, work
   estimation, hash table init) that don't shrink with block size. At
   np=16, each diagonal block is 16k × 16k, and the 2.4 ms SpGEMM includes
   ~1 ms of overhead.

3. **RIKYU GPU limit** — max 16 GPUs per job (4 nodes × 4 GPUs). The 5×5
   grid (np=25) would need 25 GPUs, but RIKYU only allows [1,2,3,4,8,12,16].

4. **Blocking MPI** — the Send/Recv for symmetry exchange is blocking. With
   many pairs, the slowest pair can delay all ranks. Non-blocking MPI
   (MPI_Isend/MPI_Irecv) would overlap communication with the threshold
   pipeline.

### Comparison: 1D row decomposition vs 2D SUMMA + symmetry

| Approach | np=4 median | np=8 median | np=16 median | OOM risk | Redistribution |
|----------|------------|------------|-------------|----------|---------------|
| 1D row (Finding 10) | 2.35 ms | 1.39 ms | — | High (full A replicated) | Per-round |
| 2D SUMMA + sym | 4.13 ms | — | 2.73 ms | Low (block layout) | None (fixed mode) |

The 1D approach is faster at np=4 (2.35 vs 4.13 ms) because it has no
communication and simpler cuSPARSE calls. But it doesn't scale beyond 8 GPUs
on RIKYU (GPU limit) and can't handle large matrices in purification
(redistribution cost). The 2D SUMMA + symmetry approach is the production
choice.

---

## Finding 12: Load balancing — 6.0× scaling on 16 GPUs (101× over CPU)

**Date:** 2026-11-05
**Jobs:** 286603 (np1 LB), 286604 (np16 LB sym), 286605 (np9 LB sym),
286877 (np16 LB sym instrumented), 286942/286943/286944 (LB no-sym),
286982/287008 (np4 LB determinism), 287009 (np4 LB over IB), 287055 (np4 LB
no-sym over IB)
**Files:** `yukawa_gemm.f90` (`--load-balance` flag via
`PermutationModule` + `LoadBalancerModule`)

### Motivation: the diagonal-block lottery

The timer breakdown (Finding 11 analysis) showed the 2D SUMMA wall time is
set by the **slowest diagonal block's cuSPARSE SpGEMM**, because the Yukawa
matrix's spatial structure puts all near-diagonal work on diagonal blocks
(24.5M flops vs 0–1.5M on off-diagonals at np=16). Worse, one diagonal block
at np=9/np=16 fell into a cuSPARSE slow regime (~10.6 GF/s vs 14.8 GF/s for
identical flops) — the "pathological block lottery" — and the MAX-over-ranks
wall always pays for the worst block.

NTPoly has native load balancing (`LoadBalancerModule`/`PermuteMatrix` +
`PermutationModule`/`ConstructRandomPermutation`), which replaces A with
P·A·Pᵀ (fixed-seed random permutation, Knuth-shuffled once, broadcast).
Properties that make it safe for a benchmark: **trace/grand-sum/nnz of
A² are invariant** under symmetric permutation, and symmetry is preserved.

### Result 1: load balancing eliminates the lottery — but breaks the exchange

| np | Config | Diagonal blocks | Off-diag blocks | Wall |
|----|--------|----------------|-----------------|------|
| 16 | no-LB + sym | 1.81, 1.81, 1.82, **2.42** ms (lottery) | 0.47 ms (+ 0.3 zero-product) | 2.73 ms |
| 16 | LB + sym | 1.24, 1.25, 1.25, 1.44 ms (**uniform ✓**) | 0.66–0.74 ms (all nonzero) | 5.36 ms ✗ |

Load balancing worked exactly as designed — the diagonals fell into a tight
1.24–1.44 ms spread. But the wall time DOUBLED, because:

1. **The exchange data exploded.** Without LB, Yukawa's spatial decay makes
   off-diagonal products nearly all below threshold → the C^T exchange
   messages are tiny. With random shuffling every block produces ~265K nnz
   output (m_loc=16000 blocks at np=16) → ~3.2 MB per Cᵀ message.
2. **Lower-triangle ranks idle 3.8–5.4 ms** (measured via per-rank
   instrumentation) waiting on the exchange; the critical path moved from
   diagonal SpGEMM to the exchange.
3. **cuSPARSE is ~5% slower on shuffled matrices** even at np=1 (7.57 →
   7.92 ms) — spatial locality helps the kernel too.

### Result 2: the big win — load balancing with NO exchange at all

With uniform blocks from LB, `--symmetric` is no longer needed to dodge the
pathology, so we can drop the exchange entirely: every rank computes its own
C block. This wastes the symmetry (upper and lower both compute) but removes
ALL communication from the per-round critical path:

| np | no-LB + sym | LB + sym | **LB + no-sym (best)** | LB no-sym vs no-LB sym |
|----|------------|----------|------------------------|------------------------|
| 1  | 7.57      | —        | 7.92 (5% worse)        | —                      |
| 4  | 4.13      | 13.6     | **3.08**               | 1.34×                  |
| 9  | 4.18      | 8.69     | **1.78**               | 2.35×                  |
| 16 | 2.73      | 5.36     | **1.27**               | 2.15×                  |

**64k scaling (LB + no-sym):** 7.57(no-LB np=1) → 3.08 → 1.78 → 1.27 ms =
**6.0× on 16 GPUs = 101× over the 128 ms CPU baseline.**

All traces match (64287.048350011…) and grand sums match (98232.41…) for
every correct run.

### RIKYU transport landmine: UCX_TLS=sm corrupts NTPoly's distributed path

The np=4 single-node runs (which set `UCX_TLS=sm,sysv,posix` /
`UCX_NET_DEVICES=""` to dodge mlx5 UAR exhaustion) produced **wrong,
deterministic** grand sums with `--load-balance` (111441, 124650 — vs
reference 98232), BUT correct results with IB/default transport, and correct
results at np=1/9/16 anywhere. Root cause: NTPoly's distributed multiply
(used by `PermuteMatrix`) uses **non-blocking collectives**
(MPI_IAllGather/MPI_IAllGatherv) via the UCX shared-memory transport, which
deterministically corrupts some data at np≥2 on a single node. Our cuSPARSE
2D SUMMA path uses only **blocking** collectives (MPI_AllGatherv, MPI_Bcast,
MPI_Send/Recv) and is unaffected — all no-LB np=4 results matched exactly.

**Rule going forward:** single-node runs that touch NTPoly's distributed
solver/PermuteMatrix path must NOT set `UCX_TLS=sm,sysv,posix`; use default
(IB loopback) transport. The blocking-collective cuSPARSE 2D SUMMA path can
use either. At np=4 single-node, omitting `UCX_TLS` gives correct results
(np=4 LB no-sym: 3.08 ms, grand_sum 98232.41 ✓, job 287055).

### Updated production recommendation

For multi-GPU cuSPARSE 2D SUMMA: **`--load-balance` WITHOUT `--symmetric`**
under default (IB) transport. The permutation is computed once (fixed seed
424242, broadcast), stays applied for the entire purification chain (results
map back with `UndoPermuteMatrix` at the end, exactly as InverseSolvers does),
and trades duplicated symmetric compute for zero per-round exchange — a clear
win at today's GPU:IB bandwidth ratio.

### Remaining scaling gap (why 6.0× and not 16×)

Even with uniform blocks and no exchange, np=16 is 6.0× not 16×: cuSPARSE
per-block throughput still degrades as blocks shrink (16.9 → 14.8 GF/s at
16k blocks) from fixed per-call overhead, and the LB + no-sym config does
2× the flops of the symmetric config (both triangles computed). The next
lever is the 2.5D k-split (Finding 13): splitting the shared dimension
across slices reduces per-call work further, at the cost of a partial-product
reduction.

## Finding 13: 2.5D k-split — correct but 47–2400× slower; closed for the
## ultra-low-fill regime

**Date:** 2026-11-05
**Jobs:** 287309/287316 (ksplit pre-fix), 287330/287331 (ksplit final)
**Files:** `yukawa_gemm.f90` (`--ksplit` flag, `MergeLocalBlocksToCSC`
helper, within-slice fingerprint), `cusparse_gemm.cu`
(`cusparse_gemm_ab_partial_`, `cusparse_ksplit_merge_`)

### The experiment (agreed item 2)

Split the shared dimension of C = A·A across process slices: grid
pr×pc×ps, slice s computes the partial C_ij^(s) = A_i·(:,K_s)·B_·j(K_s,:)
over the s-th slab K_s of the inner dimension (kfull = N/ps), then partials
are summed across `between_slice_comm` and thresholded. On RIKYU's 16-GPU
cap the only feasible 16-rank k-split is (2,2,4); (2,2,2) gives the np=8
point. This mirrors NTPoly's own 3D multiply, which uses exactly this
structure (slice s owns sub-block column/row s, reduced threshold
thr/(ps·1000), `ReduceAndSumMatrix` across `blocked_between_slice_comm`).

Implementation: `--ksplit` (with `--cusparse-2d --load-balance`). Per rank:
AllGather the full row panel (as in 2D), slice the k-slab out of A_local
and B on the CPU once, then per round: `cusparse_gemm_ab_partial_`
(SpGEMM, no threshold, partial downloaded to host), `MPI_Allgatherv` over
`between_slice_comm`, and `cusparse_ksplit_merge_` (H2D, scatter-add into
a dense m_loc×kk GPU buffer via atomicAdd, threshold, compact back to
CSR, D2H). The dense buffer for (2,2,·) blocks at 64k is 32000² × 8 B =
8.2 GB — fits B200 HBM easily; the Fortran side rejects configs whose
dense buffer exceeds 64 GiB (which rules out k-split at 175k/885k — the
scalable merge is future work that this finding makes moot).

### Three porting bugs found on the way (all fixed)

1. **Sub-block gather.** With ps>1, NTPoly splits `local_data_r` into a
   (ps × ps) array of CSC sub-blocks (uniform tiles of the local block,
   remainder in the last; the same mechanism as the OMP-threads landmine,
   driven by `column/row_block_multiplier = (pr/pc)·ps`). Reading only
   `local_data_r(1,1)` gathers 1/ps² of the block and segfaults in the
   CSC→CSR transpose. Fix: `MergeLocalBlocksToCSC` merges all sub-blocks
   into one full-block CSC before the AllGather (also makes the gather
   OMP-thread-count independent).
2. **Index base of the partial.** The unthresholded partial product
   downloaded from cuSPARSE has **0-based** col indices (the 1-based
   conversion normally happens in the threshold pipeline's compact
   kernel, which the partial path skips). The merge scatter assumed
   1-based and computed col−1 → illegal memory access. Fix: scatter
   consumes 0-based directly; compact still emits 1-based.
3. **Fingerprint double-count.** After the k-split merge, every (i,j)
   block exists once **per slice**; reducing the fingerprint over
   MPI_COMM_WORLD counts each block ps times (trace came out exactly
   2×/4× the reference). Fix: reduce over `within_slice_comm`, which
   holds each block exactly once (identical to the full grid at ps=1).

### Results (64k Yukawa, fixed mode, threshold 1e-4, LB, 5 rounds)

| Config | partial (SpGEMM+D2H) | exchange | merge (GPU) | round wall | vs 2D |
|--------|----------------------|----------|-------------|------------|-------|
| np=8 (2,2,2) k-split | 1.4–5.0 ms | 4.3–47 ms | 32–39 ms | **84.6 ms** | 47× slower than 1.78 ms |
| np=16 (2,2,4) k-split | 0.9–4.2 ms | **2.71–2.99 s** | 37–54 ms | **3040 ms** | 2400× slower than 1.27 ms |
| np=8 (2,2) 2D SUMMA+LB | — | — | — | 1.78 ms | — |
| np=16 (4,4) 2D SUMMA+LB | — | — | — | 1.27 ms | — |

Fingerprints (fresh runs, within-slice reduce) match the reference
(trace 6.42870483500114e+04, grand_sum 9.82324134997852e+04) to 13
significant figures for both ps=2 and ps=4 — the k-split is numerically
correct, including cancellation and threshold behavior after summation.

### Why it loses (the physics)

At 0.06% input fill the unthresholded partial products are **~40–90×
denser than A itself** (measured: 2.1–2.7M nnz per partial vs 648k per
input panel). Each round must move and merge every partial triplet:

- exchange: 3 × 32 MB per rank per round at np=16 (pageable host buffers,
  cross-node 4-way fan-in) collapsed to ~36 MB/s per rank → seconds.
  Even the best measured exchange (np=8, ~5 ms) plus merge (33–39 ms:
  8.2 GB memset + scan + compact + unpinned H2D/D2H) is ~25× the entire
  2D wall time.
- the SpGEMM itself did get faster (0.5–4 ms vs ~0.9–1.4 ms per block at
  (4,4)) — the compute lever works, but it is buried under 10–1000× of
  data-movement tax.

Threshold cannot be applied before the sum (values cancel and
accumulate), so the partial inflation is unavoidable for exact results.
A GPU-direct/NVLink peer merge would only cut the exchange; the merge
tax alone keeps k-split uncompetitive. **2.5D k-split for thresholded
SpGEMM at ultra-low fill is closed**: it pays to communicate partial
products that are ~100× the input data, to save at most ~2× on
latency-bound SpGEMM compute. The corollary: 2D + load balancing is the
right algorithmic point at this fill/GPU-count, and further gains must
come from the SpGEMM kernel itself (F8/F9) or from exploiting the
matrix's spatial decay without global reshuffles.

### Side-discovery: per-block density differs between diagonal and
### off-diagonal blocks (correction of an earlier claim)

An earlier version of this finding claimed cuSPARSE ALG1 emits ~28%
unmerged duplicate column entries. **That was wrong** — it rested on
assuming all 16 blocks hold 422,300 entries because rank 0's round print
said so. In fact rank 0 is a *diagonal* block: measured unthresholded
per-block counts are ~1.25–1.30M for the four diagonal blocks vs
~0.79–0.86M for the twelve off-diagonal blocks (job 287415, per-rank
`C_nnz` lines), and the thresholded counts follow the same pattern
(~422K diag, ~272K off-diag). The totals are exactly consistent:
4·422K + 12·272K = 4,865,120 = the np=1 hybrid's `final_nnz`, and the
k-split (2,2) blocks (2 diagonal-type + 2 off-diagonal-type quarters)
come out at 1.386M ≈ 2·422K + 2·272K. cuSPARSE ALG1's output is clean.

The real effect worth noting: **even under load balancing, diagonal
blocks are ~1.5× denser in C than off-diagonal blocks** (flops are
uniform per F12, but C-density is not). The guaranteed-nonzero diagonal
of C (64k entries, |C_rr| ≥ Σ_k A_rk² > threshold) plus the on-site
energy d on A's diagonal (every C entry gets a 2d·A'_rc term) bias
surviving-entry counts toward diagonal blocks, and this — not A-nnz
imbalance — is the residual per-rank SpGEMM spread under LB (compute
0.53–0.89 ms at np=16, diagonal ranks slowest, Finding 15).

## Finding 14: Load balancing at 175k/885k — wins at scale; runtime LB OOMs
## at 885k, pre-permuted file fixes it

**Date:** 2026-11-05
**Jobs:** 287266 (175k LB np16 ✓), 287267 (885k runtime-LB np16 OOM),
287296 (permute 885k), 287311 (885k pre-permuted LB np16 ✓)
**Files:** `permute_binary.py`

### Result

- **175k, np=16, (4,4), runtime `--load-balance`, no-sym: 2.849 ms**
  (trace 1.7640989522371601e+05, grand_sum 2.7062020900536800e+05; vs
  no-LB sym 6.18 ms → **2.17×**).
- **885k runtime LB OOMs**: `PermuteMatrix` builds the permutation via
  NTPoly's distributed multiply, whose dense local accumulator (F1,
  ~24 B/elem of the 221184² local block ≈ 1.1 TB) cannot fit node RAM at
  (4,4) — ranks die in the LB phase itself (job 287267, oom_kill). Note
  the 2D multiply itself runs fine at 885k (cuSPARSE path, 28 ms no-LB) —
  it is specifically the LB step's NTPoly multiply that OOMs.
- **Workaround: pre-permute the file.** `permute_binary.py` (uv/numpy,
  seed 424242) applies a random symmetric permutation to the NTPoly
  binary directly (COO index remap, sorted output); runs then use the
  permuted file WITHOUT `--load-balance`. Same distribution effect, no
  dense accumulator. **885k np=16: 10.457 ms** (trace
  8.8876795569288661e+05, grand_sum 1.3690607691073134e+06; vs no-LB
  sym 28.0 ms → **2.68×**).

Production table after this finding:

| Matrix | CPU full node | 2D SUMMA+LB no-sym np=16 | vs CPU |
|--------|---------------|--------------------------|--------|
| 64k    | 128 ms        | **1.267 ms**             | 100.8× |
| 175k   | OOM           | **2.849 ms**             | —      |
| 885k   | OOM           | **10.457 ms**            | —      |

(175k uses runtime LB; 885k needs the pre-permuted file until/unless
PermuteMatrix's accumulator is ever bypassed — e.g. by teaching the
mini-app to permute COO triplets directly.)

## Finding 15: Detailed timings of the production config

**Date:** 2026-11-05
**Jobs:** 287415 (64k np=16 LB), 287417 (64k np=4 LB), 287416 (885k np=16,
pre-permuted)
**Files:** `yukawa_gemm.f90` (SPROF setup timers), `cusparse_gemm.cu`
(rank-tagged `[prof r<rank>]` lines, CUSPARSE_PROFILE)

Per-phase timings from the rank-tagged CUSPARSE_PROFILE instrumentation
(GPU phases via cudaEvent, host phases via wall clock). "med/max" are over
all ranks; the round wall is set by the **max** per-rank total (MAX-reduce).

### Steady-state per-round phases (median of 5 timed rounds)

64k, np=16 (4×4, 4 nodes) — round wall 1.256 ms (max 1.297 ms):

| Phase | med (ms) | max (ms) | share of max path |
|-------|---------:|---------:|------------------:|
| SpGEMM compute (GPU, warm) | 0.567 | 0.893 | ~62% |
| SpGEMM copy | 0.050 | 0.053 | |
| full SpGEMM call (incl. work-est + descriptors) | 0.699 | 1.034 | ~67% |
| GPU threshold mark + count | 0.034 | 0.036 | |
| CPU prefix scan + row_ptr | 0.023 | 0.025 | |
| GPU compact | 0.063 | 0.080 | |
| D2H download (thresholded block) | 0.047 | 0.064 | |
| **round total (max path)** | | **1.26** | 100% |
| H2D upload A_local + B (first call only) | ~12 | ~41 | amortized |

64k, np=4 (2×2, 1 node) — round wall 3.087 ms:

| Phase | med (ms) | max (ms) |
|-------|---------:|---------:|
| SpGEMM compute (GPU, warm) | 2.454 | 2.481 |
| full SpGEMM call | 2.612 | 2.645 |
| GPU mark + count | 0.061 | 0.063 |
| CPU prefix scan | 0.029 | 0.032 |
| GPU compact | 0.165 | 0.166 |
| D2H download | 0.121 | 0.157 |
| **round total (max path)** | | **3.09** |

885k, np=16 (4×4, 4 nodes, pre-permuted) — round wall 10.441 ms:

| Phase | med (ms) | max (ms) |
|-------|---------:|---------:|
| SpGEMM compute (GPU, warm) | 7.774 | 8.537 |
| full SpGEMM call | 8.106 | 8.888 |
| GPU mark + count | 0.348 | 0.351 |
| CPU prefix scan | 0.114 | 0.126 |
| GPU compact | 0.542 | 0.566 |
| D2H download | 0.452 | 0.855 |
| **round total (max path)** | | **10.44** |
| H2D upload A_local + B (first call only) | ~40 | ~67 | amortized |

Readings:

1. **The SpGEMM is 80–84% of the round** at every scale; the GPU threshold
   pipeline is small (0.14–0.20 ms at 64k np=16; 1.5 ms at 885k) and the
   CPU prefix scan is negligible (0.02–0.12 ms).
2. **Fixed per-call cuSPARSE overhead** (full-call span − compute − copy)
   is ~0.08–0.15 ms per round regardless of size — 10% of the 64k np=16
   round, 1% at 885k.
3. The compute spread across ranks under LB is real and tracks block
   *type*, not A-nnz: diagonal blocks compute ~1.5× denser C products
   (0.53–0.89 ms compute at np=16, diagonal ranks slowest — see the
   corrected side-discovery in Finding 13). The MAX-reduce wall pays for
   the diagonal blocks.
4. np=4 → np=16 compute scaling is 4.3× (2.454 → 0.567 ms) on 4× the GPUs.

### One-time setup phases (as measured, first call)

| Stage | 64k np=4 (1 node) | 64k np=16 (4 nodes) | 885k np=16 (4 nodes) |
|-------|------------------:|--------------------:|---------------------:|
| read + LB permutation | 37.8 s | 8.9–10.2 s | 19.3 s (read only) |
| gather (AllGather row panel) | 13.7 ms | 183–247 ms | **2.85–3.04 s** |
| transpose (CSC→CSR) | 14–25 ms | 6.6–11.5 ms | 99–186 ms |
| bcast B along column | 3–13 ms | 2–67 ms | 26–292 ms |

Notes:

- read+LB at np=4 is 3.7× slower than at np=16 (PermuteMatrix's NTPoly
  multiply shrinks with rank count); at 885k there is no LB multiply
  (pre-permuted file), so it is pure MPI-IO read of the 593 MB binary.
- The np=16 gather (0.18–0.25 s at 64k) includes MPI first-touch connection
  setup across nodes and absorption of the read skew (read_lb spread is
  1.3 s; early ranks wait inside the AllGather). Single-node np=4 measures
  13.7 ms for the identical code path.
- **Purification-chain implication**: in fixed mode the gather/transpose/
  bcast run once, but a real purification chain recomputes them every
  round (the input matrix changes). At 885k the gather alone moves
  ~111 MB/rank over IB and measured 3.0 s — that would add ~30% per round
  on top of the 10.4 ms multiply unless it is overlapped, fused, or the
  per-round matrix update is expressed in the 2D block layout directly
  (no re-gather). At 64k np=16 the equivalent first-gather cost (~0.2 s)
  is dominated by first-touch effects and needs a steady-state measurement
  before quoting a per-round number. This is the main open question for
  end-to-end purification performance (handoff next-steps item 2).

## Finding 16: P1 — managed-memory residency + combine ops (scale /
## increment / copy) on GPU

**Date:** 2026-11-05
**Jobs:** 287633 (scale), 287634 (increment), 287632 (copy), plus debug
runs 287592/287605/287612/287614
**Files:** `yukawa_gemm.f90` (`--op scale|increment|copy`, `--reps`,
`--alpha`, `MergeLocalBlocksToCSC` reuse), `cusparse_gemm.cu`
(`gpu_ops_init_/alloc_/scale_/increment_/copy_/grandsum_`)

### Design (user-directed): CUDA MANAGED memory, not explicit buffers

Per GPU_ROADMAP.md P1, each rank's local block lives in
`cudaMallocManaged` CSR arrays (slot model, 2 slots): Fortran fills them
directly via `C_F_POINTER` (no upload step), GPU kernels operate in
place, and any op can fall back to CPU code touching the same arrays.
`cudaDeviceSynchronize()` before returning keeps CPU/GPU handoff simple.
`gpu_ops_init_` (cudaSetDevice per local rank) must run **after** any
NTPoly OpenACC multiply — the OpenACC runtime binds its device lazily
and an early cudaSetDevice races it (`CUDA_ERROR_INVALID_HANDLE` in
`MultiplyBlock.f90`).

Increment is a per-row two-pointer merge of the two sorted CSR rows
(duplicates summed on the fly, C + α·B), with the row_ptr prefix scan
done by the CPU **directly on managed memory** (no explicit D2H), then a
merge-write kernel, then a pointer ping-pong into slot C. Thresholding
is not applied in v1 (NTPoly's solver calls use threshold 0 here).

### Results (64k, np=16, (4,4) grid, LB, 16k×16k blocks, 232K nnz/block)

| Op | GPU median | CPU median (NTPoly) | Speedup | Parity |
|----|-----------:|--------------------:|--------:|-------:|
| scale (C ← αC) | **15.2 µs** | 50.6 µs | 3.3× | PASS 1.9e-14 |
| increment (C += αB) | **330.9 µs** | 920.4 µs | 2.8× | PASS 3.1e-14 |
| copy | **45.8 µs** | 1,613.9 µs | 35.2× | PASS 9.2e-16 |

Grand-sum parity between the GPU resident path and NTPoly's CPU path
after 23 identical applications (3 warmup + 20 timed): 13+ significant
figures on all three ops. The P1 exit criteria are met (combine ops
≪ 0.5 ms with resident data; fingerprint parity; git baseline).

### Bugs found on the way (all fixed, all educational)

1. **Parity double-reduce (the big one):** `MatrixGrandSum` already
   returns the GLOBAL sum on every rank (it allreduces over
   `within_slice_comm` internally); allreducing it again over
   MPI_COMM_WORLD multiplied the CPU reference by the rank count
   (16.02× at np=16). The GPU slots were correct from the first managed
   run — two full debug cycles were spent chasing a GPU bug that was a
   double-reduce in the *checker*. Rule for all future reductions: NTPoly
   reduction results are global; never reduce them again.
2. **Warmup is stateful:** warmup rounds apply the op (needed to warm
   the buffers), so the final state reflects warmup+reps applications
   (grand scaled by α²³, not α²⁰). The expected-final formula and the
   CPU reference must run the identical history.
3. **Early cudaSetDevice vs OpenACC:** see above — init order matters.
4. **Copy mode used an uninitialized `matB_ps`** (only built for
   increment) → hung runs. Both non-scale modes now build it.

### Observations

- Scale at 15 µs is dominated by kernel-launch + sync latency (~10 µs
  floor); the arithmetic itself is ~1 µs of memory traffic (232K × 8 B).
- Increment at 331 µs is dominated by the counts D2H/prefix/H2D round
  trip through managed memory; a CUB device-wide scan or persistent
  row_ptr capacity would cut it further (future optimization, not
  needed for P1 exit).
- The 35× copy win reflects how heavy NTPoly's triplet-level CopyMatrix
  is; in solver chains, copy should become free (pointer swap).
- Per-block data after LB is *not* uniform (232K nnz on diagonal blocks
  vs 135–142K off-diagonal; grand sums 17.3K vs 0.83K — the on-site
  diagonal concentrates in diagonal blocks). This does not affect
  correctness (grand sums are permutation-invariant and fingerprints
  match) but matters for per-block load expectations in P3+.

### Next (P2)

Reductions (trace / norm / dot) as per-block GPU reduces + MPI allreduce
of scalars — the grandsum machinery here is already the template.

## Finding 17: P2 — GPU reductions (trace / dot / norm / sigma /
## grandsum), exact parity

**Date:** 2026-11-05
**Jobs:** 288266 (trace), 288264 (dot), 288265 (norm), 288268 (sigma),
288267 (grandsum), 288270 (trace re-run after the diag-offset fix)
**Files:** `yukawa_gemm.f90` (`--op trace|dot|norm|sigma|grandsum`),
`cusparse_gemm.cu` (`gpu_ops_trace_`, `gpu_ops_dot_`,
`gpu_ops_colsum_`)

### Semantics (verified against NTPoly source)

- **trace**: Σ of entries with global row == global col. Block-local
  condition: `col_ind == r + 1 + diag_offset` with
  `diag_offset = (my_row - my_col) * m_loc`; off-diagonal blocks get an
  offset outside the local column range and contribute exactly 0.
- **dot**: Σ over common (r,c) of A·B — per-row two-pointer walk over the
  two sorted rows (NTPoly computes it as GrandSum(PairwiseMultiply)).
- **norm**: NTPoly's `MatrixNorm` is the **max column 1-norm**:
  per-column Σ|v| summed over row-blocks (`column_comm`), max over
  column-blocks (`row_comm`).
- **sigma**: `1/(max column abs-sum)²` — same two-stage reduce, used for
  eigenvalue bounds (Chebyshev/PM).
- MPI combines stay on the Fortran side (it owns the grid communicators);
  only scalars / a kk-length vector cross to the host.

### Results (64k, np=16, (4,4) LB, 20 reps, median)

| Reduction | GPU | CPU (NTPoly) | Speedup | Parity |
|-----------|----:|-------------:|--------:|-------:|
| trace | 960 µs | 3,328 µs | 3.5× | 0.00 (exact) |
| dot | 661 µs | 1,140 µs | 1.7× | 0.00 (exact) |
| norm | 1,243 µs | 3,457 µs | 2.8× | 0.00 (exact) |
| sigma | 1,075 µs | 3,695 µs | 3.4× | 0.00 (exact) |
| grandsum | 602 µs | 77 µs | 0.1× | 9.2e-16 |

trace = exactly 64,000 (the on-site diagonal: 64k entries of 1.0); dot =
64,287.04835001135 = trace(A²)'s famous companion value; norm/sigma
bit-identical to NTPoly.

### Notes

1. **GrandSum is slower on GPU (0.1×)**: NTPoly's CPU grand sum is a
   tight local loop (77 µs) while the GPU version pays kernel + managed-
   memory + allreduce latency. In a resident GPU chain this is moot — the
   alternative (CPU GrandSum) requires downloading the whole matrix
   (444 MB at 885k). Use the GPU reduction in the resident regime.
2. **The "<5% of a round" convergence-check target is NOT yet met**: a
   trace+norm check costs ~2.2 ms vs the 1.267 ms multiply round. The
   per-op costs are dominated by fixed overheads (kernel launch, the
   2–3 MPI collectives at 4-node latency, D2H). The fix is to fuse
   trace+norm+grandsum into ONE kernel pass with ONE packed allreduce
   (and ultimately ride the multiply's existing allreduce) — deferred to
   P3, where the convergence check is actually consumed.
3. **Bug found in the bench harness (fixed)**: the ops branch never
   assigned `my_row`/`my_col` before computing the trace's diag_offset —
   garbage offsets let off-diagonal blocks match a few accidental
   entries (trace = 64,000.0613 instead of 64,000). After assigning from
   the process grid: exact.
4. `MeasureAsymmetry` (= Norm(A − Aᵀ)) needs the distributed transpose →
   deferred to P4 as planned.

## Finding 18: P3 correctness review and local repairs (2026-11-10)

The remote Fortran/CUDA files were byte-identical to the local files.
Cluster availability is not a sufficient explanation for the failed P3
runs: the source had several independent defects:

- Chunk packing used row numbers/row-pointer-array offsets as value-array
  offsets. Even a single slab could overwrite entries and read uninitialized
  data. Slab assignment and bounds also disagreed for ragged dimensions.
- Global symmetry was incorrectly used to alias each off-diagonal block's
  CSR as CSC. Slots now always contain true CSR, with explicit conversion
  to/from NTPoly CSC; the solver gathers the actual B column panel.
- Upload caches reused old values/patterns whenever shape and nnz matched.
  A's shared dimension was missing from its cache key. Allocations are now
  reused while the supplied arrays are refreshed; solver workspace sizes
  are queried afresh for changing products.
- Increment swapped buffers without swapping allocation capacities, causing
  later writes beyond the actual scratch allocation. Capacities now follow
  their pointers. Its merge also applies the reference threshold after
  summation, removing cancellations and small entries.
- Failed SpGEMM partials were silently skipped, and failed merges became
  empty slots with stale/uninitialized pointers. Errors now stop the MPI
  run; valid empty products/slots have explicit empty CSR representations.
- GPU convergence now uses NTPoly's monitor. Final trace is recomputed
  after the last update; the driver checks full-matrix and scalar relative
  errors against the reference, with a 1e-12 acceptance tolerance.
- Duplicate Fortran declarations prevented a fresh local syntax check.
  The existing remote executable's provenance was therefore uncertain.
  CUDA object compilation now has a source dependency instead of occurring
  only during CMake configuration.

`solver_csr.f90` provides bounds-tested transpose, slab planning, and packing.
Slabs are planned using actual intermediate-product counts, bounded by
`--chunk-work` (default 4 million). A single shared index exceeding that
budget causes an explicit error; this is not a vendor-resource guarantee.
The merge still concatenates all partials and uses a dense local accumulator,
so these fixes do not establish arbitrary high-fill scalability.

Local regression command (macOS, existing CPU `BuildCheck` library required):

```sh
python3 NTPoly/MiniApps/YukawaGemm/tests/run_host_checks.py
```

The suite checks actual operation and upload-helper bodies with serial
CUDA emulation under ASan/UBSan, bounds-checked CSR helpers, and the actual
Fortran driver with deterministic CPU substitutes for cuSPARSE. TRS2 cases
exercise 1/4 MPI ranks, load balancing, padded odd dimensions, forced
multiple slabs, and empty off-diagonal blocks. Matrix/scalar parity and
reference iteration counts are checked; an injected multiply error must
abort. These are host tests, **not CUDA validation or performance results**.

The complete suite also passes all eight P1/P2 operation benchmarks on four
MPI ranks and verifies CUDA rebuild dependencies with a mock compiler.
Those checks exposed an additional pre-existing P1 harness bug: stateful
scale/increment/copy never assigned the CPU comparison scalar. The harness
now computes the final NTPoly grand sum after timing. Bounds checks also
exposed an insertion-sort guard that assumed Fortran short-circuits `.AND.`;
the median helper and triplet sorter now check bounds before indexing.
Historical P1 parity claims should be re-established with the corrected
harness on RIKYU.

### Real CUDA validation (2026-11-10/11)

Build 311070 succeeded in 66 s in the isolated
`~/agent/work/ntpoly-p3-tiny-20261110` directory. Tiny suite 311076 passed
all seven 1/4-rank solver cases and the actual CUDA operation regression;
Compute Sanitizer reported **0 errors**. Worst tiny scalar/matrix relative
error was 2.34e-16. A tiny `-Ktrap=inv` diagnostic (311087) completed without
trapping, although nonfatal IEEE-invalid exit flags remain unexplained.

The first 64k attempt (311090) stalled after the first partial products and
timed out; its gated full run never started. Instrumented 311105 reproduced
the stall on different nodes and localized it to Fortran partial-buffer
growth, before the merge. A standalone initialized-buffer growth test passed,
so the precise runtime trigger is not established. Replacing repeated
growth/copying of the entire allocated capacity with `ReservePartials` (one
allocation sized for the required entries; copy only the initialized prefix)
resolved it. The host suite now includes growth beyond the 65,536-entry seed.

Build 311150 succeeded; job **311151** completed both 16-GPU checks in 58 s:

| Check | Iterations | Energy rel. error | Trace rel. error | Matrix rel. error |
|---|---:|---:|---:|---:|
| Short, capped | 10 | 2.28e-16 | 1.59e-15 | 8.49e-16 |
| Full, converged | 24 | 5.70e-16 | 1.14e-16 | 1.61e-14 |

Reference and GPU iteration counts match (NTPoly reports 11 after exhausting
a cap of 10; this represents 10 executed iterations). Full energy is
15947.315148941969; trace is 32000.000143710073. This validates correctness
for the gapped 64k chain; it does not establish arbitrary high-fill coverage.

**Performance remains open:** summing root iteration timings gives 14.760 s
for the GPU loop versus 3.7 s for the NTPoly reference. The reference in this
build uses OpenACC, so it is not a pure CPU baseline. GPU multiply accounts
for 13.512 s, reductions 1.138 s, updates 0.035 s, and host synchronization
0.076 s. These are root timings, not independently reduced per-phase maxima,
and omit setup. No solver speedup is demonstrated.

The original remote `ntpoly-gemm` workspace was preserved. Tested source
and executable are under `ntpoly-p3-tiny-20261110/source-v2` and
`build-v2/bin/yukawa_gemm`. Local evidence is archived under
`validation/rikyu-p3-tiny-20261110/`, including source hashes, build logs,
raw tiny/64k output, `verified-summary.json`, and `big-summary.json`.

## Finding 19: Solver profile and genuine CPU comparison (2026-11-11)

All comparable cases use the gapped 64k matrix, LB seed 424242, grid
4×4×1, 16 MPI ranks on four nodes, threshold 1e-4 and convergence 1e-6.
Medians of three full 24-iteration runs:

| Implementation | Threads/rank | Median seconds | Range seconds |
|---|---:|---:|---:|
| Pure CPU NTPoly | 1 | 4.897 | 4.832–4.898 |
| Pure CPU NTPoly | 32 | **1.577** | 1.552–1.704 |
| OpenACC dense-accumulator NTPoly reference | 1 | 3.746 | 3.679–3.765 |
| CUDA/cuSPARSE custom solver loop | 1 | 14.735 | 14.450–15.315 |

The CPU library uses `-O3 -mp`, with neither `-acc` nor `-gpu`.
Linkage was checked for absence of CUDA/cuSPARSE/OpenACC libraries. The
new `trs2_reference` target runs the existing NTPoly TRS2, without calling
the CUDA wrappers. CPU MPI ranks were bound to 1 or 32 cores, respectively,
using four ranks per node; the 32-thread case uses 512 CPU cores overall.
The allocation reserves 16 GPUs because RIKYU bundles its CPU resources
with GPU allocations; these GPUs are unused by the CPU benchmark.

CPU 1/32-thread cases have the same unpadded logical dimension (64000),
Gershgorin bounds (0.1942489437, 1.8087116049), 24 iterations and final
fingerprints. Job 311220 additionally wrote both density matrices and
compared them entrywise via relative matrix norm: **0.0**.

Timing scope: both reference timings include TRS2's internal initialization
and exclude file reading/permutation. The custom CUDA loop starts with X0
already prepared, so it excludes that initialization. Its global MPI-MAX
elapsed time includes the new profile reductions/output. Verbose CUDA
event logging was disabled (`YUKAWA_CUDA_PROFILE=OFF`) for these runs.
The custom path is still about 9.34× slower than the valid threaded CPU
baseline, despite its narrower timed scope.

### Where the custom path spends its time

Per-iteration timings now separately record trace, row gathering, row
transpose, column gathering, planning/packing, chunk preparation,
SpGEMM wrapper, partial storage, merge, cleanup, product-slot upload,
update, host CSC reconstruction, energy and wall time. Each CSV reports
MPI-MAX per phase. Independent maxima can come from different ranks, so
their sums must **not** be interpreted as a wall-time decomposition.
Medians across three runs of the sums over 24 iterations:

| Phase | Seconds, sum of per-iteration maxima |
|---|---:|
| Row gather | 6.871 |
| Column gather/assembly | 8.260 |
| SpGEMM wrapper (uploads, allocation, compute and download) | 1.235 |
| Dense merge | 0.282 |
| Row transpose | 0.157 |
| Planning/packing | 0.110 |
| Partial storage | 0.097 |

A separate Nsight Systems capture covers **only the custom loop on rank 0**,
using cudaProfilerStart/Stop markers. It took 16.528 s (instrumented),
versus a 14.735 s median without Nsight. Its trace establishes:

- `MPI_Allgatherv`: 96 calls, 10.861 s; `MPI_Allgather`: 96 calls, 1.654 s.
- All **23,586 `cuMemcpyAsync` API calls** are contained inside MPI
  collective intervals on the same thread (verified from SQLite).
  23,106 calls / **9.627 s** are inside Allgatherv; the remaining
  480 calls / 0.073 s are inside Allgather. These times overlap the MPI
  times and must not be added to them.
- GPU kernels total only **0.213 s** on this rank. The largest is dense
  compact: 0.144 s across 24 launches. GPU UM-fault intervals total 0.045 s.

**Conclusion:** MPI's CUDA copy/staging path during managed-memory panel
gathering is the main profiling target. First investigate ordinary host
communication buffers or explicit host staging while preserving managed
operation slots. Tuning cuSPARSE kernels or the dense compact kernel alone
cannot recover the observed gap. This is a profiling diagnosis, not a
validated staging fix.

### Excluded configurations and limitations

The 36-thread layout pads the logical dimension to 64224. It changes bounds
and fingerprints; one diagnostic run diverged to NaN/-Inf while NTPoly's
convergence monitor still reported success. It is **not a valid performance
comparison**. The standalone reference driver now rejects nonfinite outputs.
This threaded/padded NTPoly issue remains open; it was not fixed by changing
the benchmark to the valid 32-thread layout.

The first standalone OpenACC reference harness stalled before the solver
and hit its 120 s timeout. OpenACC timings above therefore come from the
validated mini-app's embedded reference, not that failed harness.

Jobs: CPU build 311191; initial CPU 311194; audited CPU 311207;
GPU builds 311195/311201; failed standalone harness 311204;
successful CUDA/profile 311208; full CPU matrix check 311220.
All successful CUDA runs preserved 24 iterations and full-matrix relative
error 1.61e-14 against their embedded NTPoly reference.

Evidence: `validation/rikyu-profile-20261111/summary.json`, `analyze.py`,
raw phase CSVs, CPU linkage/build logs, and the Nsight `.nsys-rep`/SQLite
under `results/nsys/`. The remote working directory is
`~/agent/work/ntpoly-profile-20261111`.

## Finding 20: NCCL on one node, 1/2/4 GPUs (2026-11-11)

Added optional `--nccl` solver communication, linked against RIKYU's NCCL
2.29.3+cuda13.1. Device communication uses persistent `cudaMalloc` buffers:
one bulk host-to-device transfer, NCCL exchange, one bulk device-to-host
transfer per payload. Fixed-count data uses AllGather; variable-count CSR
data uses grouped Send/Recv. Solver scalar sums also use NCCL. MPI handles
bootstrap IDs and profiling reductions; the embedded reference remains
the existing MPI/OpenACC implementation. CPU panel construction remains.
The grouped-P2P design follows NVIDIA's documented matching-count/group
requirements: https://docs.nvidia.com/deeplearning/nccl/user-guide/docs/usage/p2p.html

This backend explicitly rejects more than four ranks or multiple nodes.
Two-GPU grids are now supported by explicit merge/split of NTPoly's local
sub-blocks, true CSR/CSC conversion, and correct rectangular-block diagonal
offsets. The square-only restriction is no longer needed for this solver.
The fixed-GEMM paths retain their own pre-existing grid requirements.

### Validation and measurements

Build 311255 succeeded. Job 311260 passed fixed/variable/empty/scalar
transport checks on 1/2/4 GPUs, plus tiny MPI/NCCL solver comparisons,
including both two-rank orientations, load balancing, and empty blocks.
NCCL INFO logs show P2P/CUMEM peer communication. NVLS multicast was not
available; it was not forced or assumed.

Job **311265** ran on a single node, **c117**, with four GPUs allocated.
Every 64k gapped-matrix case converged in **24 iterations** and passed
energy, trace and full-matrix parity against its reference. Worst relative
error across all cases: **1.616e-14**. Threshold 1e-4, convergence 1e-6,
LB seed 424242; grids 1×1, 1×2 and 2×2:

| GPUs | MPI solver loop (s) | NCCL solver loop (s) | Samples/backend |
|---|---:|---:|---:|
| 1 | 14.427 | 12.483 | 1 |
| 2 | 7.921 | 6.834 | 1 |
| 4 | **5.217** | **5.521** | 3 (median) |

One-GPU NCCL mode is a local-copy path with no inter-GPU collective; do not
attribute that difference to NCCL communication speed. The 1/2-GPU figures
are single samples. Four-GPU ranges were 5.114–5.307 s (MPI) and
5.501–5.560 s (NCCL). All values are for the same new executable/build and
include the first timed communication calls. Reading, permutation, X0
construction and communicator bootstrap are outside the timed loop.

### First-use cost versus later iterations

The initial NCCL communication calls are expensive despite communicator
bootstrap occurring before timing. Four-GPU medians:

| Timed portion | MPI (s) | NCCL (s) |
|---|---:|---:|
| Iteration 1 | 0.286 | 1.758 |
| Iterations 2–24 | 4.932 | 3.761 |

Thus NCCL is about **24% faster for iterations 2–24**, but roughly **6%
slower for the complete 24-iteration loop**. The later-iteration values
sum per-iteration MPI-MAX wall timings; they are not a measured warm-start
full solve. First-use costs must be reported separately rather than
silently subtracting them from a claimed full-solver time.

For the first NCCL four-GPU sample, iteration-1 trace, row gathering and
column gathering cost 0.754 / 0.385 / 0.372 s; subsequent trace reductions
are around 0.0002 s. The earlier impression that scalar reductions were
steadily expensive was first-use overhead, not steady-state behavior.
For iterations 2–24, the median sum of row/column phase maxima was 1.054 s
for MPI and 0.411 s for NCCL. As before, phase maxima need not belong to
the same rank and must not be added into a wall-time decomposition.

**Result:** the single-node NCCL prototype is correct and improves later
iterations. The complete short solve has not gained at four GPUs because
first-use cost dominates the savings. Next experiment: explicitly warm
and/or reuse communicators before measuring repeated solver chains, while
reporting that setup cost. CPU panel conversion/staging and the dense merge
are still present; this is not a fully device-resident solver.

Local artifacts: `validation/rikyu-nccl-20261111/summary.json`, `analyze.py`,
raw outputs/phase CSVs, NCCL transport logs, build hashes and benchmark
scripts. Remote source/executable: `~/agent/work/ntpoly-nccl-20261111/source`
and `build/bin/yukawa_gemm`. Local host regressions also passed, including
the newly added two-rank rectangular-grid test. Existing workspaces and
executables were preserved.

## Finding 21: GPU accumulation, code-generation audit and timing model (2026-11-11)

All solver runs below use one node, four GPUs/ranks, grid 2×2×1, gapped
N=64000, LB seed 424242, threshold 1e-4, convergence 1e-6, and 24
iterations. Every full case passed energy/trace/full-matrix parity;
worst relative error across these experiments is 1.818e-14. Three-sample
medians unless explicitly labeled otherwise.

### Implemented paths and validation

- `--device-partials`: scatter each cuSPARSE partial directly from device
  output into the persistent dense GPU accumulator. Eliminate raw-product
  downloads, CPU concatenation and re-upload; threshold once after all
  chunks. Legacy callers keep their previous interface and behavior.
- Replace the dense compactor's one active thread per row with one complete
  warp. Ballot/popcount preserves sorted CSR output, including tail tiles.
- Sequential device chunks have unique columns within each partial row;
  their accumulation needs no atomicAdd. The direct sparse multiply below
  has concurrent contributions and correctly retains atomicAdd.
- `--gpu-dense`: direct sparse-panel multiplication into the dense GPU
  accumulator, removing chunk planning/packing, symbolic SpGEMM, workspace
  queries, intermediate partial arrays and their kernel launches.
- Compute energy before the CPU CSC reconstruction; binary-search the
  diagonal instead of linearly scanning the sparse row.
- `--nccl-warmup` explicitly warms collectives/P2P before timing. Report
  communicator/warmup setup and setup-plus-loop separately. Panel buffers
  and CPU CSR/CSC conversions are still present.
- `YUKAWA_NCCL_SHARED_COMM=1` uses one world communicator and bootstrapped
  peer maps for row/column P2P, reducing duplicated communicator setup.
  Default separate communicators remain available for comparison.
- Buffer coherent per-rank phase vectors and gather/write `.ranks.csv`
  after the timed loop. This supports models without subtracting unrelated
  phase maxima from different ranks.

Host MPI tests pass for all three multiply paths on 1/2/4 ranks, including
ragged and empty blocks, load balancing and injected multiply failure.
Actual CUDA tests cover cancellation, sum-before-threshold, empty products,
reset, sorted compaction across 32-column boundaries and a 67-column tail.
Compute Sanitizer reports zero errors. NCCL fixed/variable/empty/scalar
tests pass on 1/2/4 GPUs with both separate and shared communicators. Tiny
solver grids 1×1, 1×2, 2×1 and 2×2 pass real CUDA parity.

### Measurements and timing scope

| Experiment | GPU loop (s) | NCCL setup/warmup (s) | Setup + loop (s) |
|---|---:|---:|---:|
| Tuned chunked cuSPARSE, 64M work budget, c242 | 3.096 | 5.311 | 8.408 |
| Paired cuSPARSE on c116 | 2.978 | 5.438 | 8.417 |
| Direct GPU dense accumulator, c116 | **2.468** | 5.035 | 7.685 |
| Direct dense, separate communicators, c117 | 2.666 | 5.293 | 7.959 |
| Direct dense, shared communicator, c117 | 2.660 | **3.584** | **6.243** |

Medians of sums need not equal sums of medians. The earlier cold-loop
paired test on c117 improved 5.547 → 5.158 s with device partials alone.
64M is a tested tuning choice for this matrix; default remains 4M.

A genuine pure-CPU one-node baseline, four ranks ×32 threads (128 cores),
is **3.238 s** (3.167–3.333 s), same 24 iterations, unpadded dimension and
matching fingerprints. Four ranks ×1 thread takes approximately 15.6 s.
The previous **1.577 s CPU** result uses four nodes and 512 cores.

The GPU loop excludes spectral bounds, X0 construction/uploads and
communicator initialization. CPU TRS2 timing includes its internal
initialization. Thus the warm loop is faster, but this is **not a claim
of cold end-to-end GPU speedup**: setup alone is more than the one-node
CPU solve. Cold GPU runtime/device binding is also outside these figures.
Repeated solver chains should reuse communicators and device state.

### PTX/SASS and communication audit

Actual sm_100 PTX/SASS and ptxas resource reports were generated on RIKYU.
All emitted kernels have **zero spill stores/loads**. Original compaction
had a serial, compiler-unrolled scalar column scan without a lane index.
The replacement emits warp ballot/popcount and coalesced lane-indexed
loads. Streaming partial accumulation has no atomic instruction; direct
multiply retains the required FP64 atomic addition. No local-memory
declarations were found in these kernels. These checks do not imply that
all instruction-level optimization opportunities have been exhausted.

Rank-0 instrumented direct-dense MPI capture: 24 multiply kernels total
**0.006884 s**; compaction **0.078566 s**, versus **0.479397 s** in the
original compactor capture. Count scan takes 0.044778 s; 24 dense clears
write 196.608 GB in 0.070656 s. These are separate instrumented runs,
not a same-run wall-time decomposition. Dense buffer zero/count/compact
move approximately 589.824 GB/GPU over 24 iterations; the captured device
reports nominal memory bandwidth 7.928 TB/s, giving an idealized bandwidth
floor of 0.0744 s for those passes alone.

MPI loop counts match implementation exactly: **96 Allgathers, 96
Allgathervs, 48 scalar Allreduces, 24 profiling Reduces**. MPI wall times
include waiting and cannot be labeled pure network transfer time.

The successful shared-NCCL loop capture matches independent predictions:
24×8 = **192 panel SendRecv kernels**, 24×2 = **48 scalar AllReduce
kernels**. Exactly **240 uploads +240 downloads +192 local device copies
=672 staging copies** appear on the NCCL stream. Total copied data is
1.926 GB on rank 0. The stage-copy GPU durations total about 0.132 s;
NCCL kernel durations include inter-rank waits. Peer application payload
is 514.183 MB sent and 382.976 MB received, excluding scalar allreduce
protocol traffic and wire framing. MPI performs only **24 profiling
Reduces** inside this NCCL loop; bulk MPI collectives are absent.

A separate uninstrumented upload probe **rejects** a generic pinned-bounce
fix: at 64 MiB, host-built managed input uploads at 22.49 GB/s, versus
20.24 GB/s including memcpy to pinned staging. Do not infer raw bandwidth
from fault-instrumented copy intervals or assume an extra bounce helps.

Early rank-0-only MPI/NCCL profiling wrappers omitted the documented
`NSYS_MPI_STORE_TEAMS_PER_RANK=1` setting. NCCL communicator initialization
failed only under those profiling configurations; unprofiled jobs passed.
See NVIDIA's [subset-rank profiling instructions](https://docs.nvidia.com/nsight-systems/UserGuide/).
The corrected wrapper and follow-up capture are in
`validation/rikyu-profile-fixed-20261111`: job 312018 completed, including
the transport toy and full MPI/NCCL captures with separate communicators.
The corrected NCCL trace has the predicted 96 AllGather +96 SendRecv +48
AllReduce kernels and 576 staging copies (240 uploads, 240 downloads,
96 P2P local copies). Shared mode instead copies its own segment for all
192 panel calls, accounting for its 672 copies. Do not call these failures a
RIKYU outage or an established NCCL library defect.

### Counterfactual model and next target

The model removes phases **within each rank's vector**, then takes the
maximum residual for each iteration. It preserves other measured waits
and assumes removed phases vanish; it is an optimistic static sensitivity
model, not a real critical-path simulation or measured resident solve.

| Direct-dense c116 scenario | Predicted loop (s) |
|---|---:|
| Reconstructed per-iteration wall maxima | 2.497 |
| Zero row/column panel exchange time | 1.908 |
| Also remove CPU row transpose and host CSC reconstruction | **0.773** |

Shared-communicator c117 gives approximately **0.877 s** for the last
scenario. New GPU panel assembly/transpose and real communication will
cost nonzero time. A sub-1.5 s reused-device loop is a plausible target,
which could beat both the one-node CPU baseline and the earlier four-node
CPU result. This has not been demonstrated. A cold single solve still
cannot beat CPU while NCCL setup alone exceeds CPU runtime.

**Next:** keep updated X, CSC/CSR conversion, row/column panel assembly
and NCCL payloads on device; retain only counts/scalars on host initially.
Measure real repeated solves in one MPI process and include X0/upload
setup before claiming full solver speedup. Kernel arithmetic is already
a small fraction of loop time; further isolated SpGEMM tuning has limited
upside. The dense local buffer still scales as O(N²/GPUs); no arbitrary
high-fill or large-N scalability claim follows from these results.

Artifacts and replay: `validation/rikyu-gpu-audit-20261111/analyze.py` and
`summary.json`, plus per-experiment archives in `rikyu-device-partials-*`,
`rikyu-gpu-opt-*`, `rikyu-gpu-dense-*`, `rikyu-nccl-shared-*` and
`rikyu-staging-*`. Build/source hashes and scripts are preserved.
Jobs: 311533 (device partials); 311544 (all optimization benchmarks passed,
later profiling failed); 311536/311542 (CPU and initial audit, profiling
wrapper issues); 312009 (direct dense, memory check and MPI trace);
312010 (shared communicator and full NCCL trace); 312014 (staging probe);
312018 (corrected profiler wrapper). All runs use at most four GPUs.

## Finding 22: Resident panels, touched bitmaps and measured reuse (2026-11-11)

Same gapped N=64000, seed 424242, threshold 1e-4, convergence 1e-6.
Every full case below takes 24 iterations and passes full matrix/scalar
parity; observed errors are approximately 5–8e-15. At most four GPUs are
active per solve, on one node. Default comparison paths remain available.

### Payload residency

`--resident-panels` requires `--gpu-dense --nccl`. The local GPU CSR slot
is copied to explicit device communication storage, row/column components
are exchanged on GPU, and CUDA kernels/scans construct horizontal A and
vertical B CSR panels. Multiplication, threshold counts, prefix scan and
compaction write the next operation slot on GPU. Per-iteration host CSR↔CSC
conversion, packing and product downloads are bypassed. A host matrix is
materialized after timing for validation. Counts/scalars and the existing
update's small row-count scan can still touch CPU; this is payload residency,
not a claim of zero host activity. Events order the default compute stream
against the nonblocking NCCL stream. No symmetry alias is used.

Independent CUDA tests compare every local CSR entry for nonsymmetric
squares, cancellations, empty peers/all-empty subgroups, reset and sorted
output. Grids 1×1, 1×2, 2×1 and 2×2 pass; sizes 6/66 exercise partial/full
32-column bitmap words. Tiny padded TRS2 and repeated-solve tests pass,
as do local host MPI/bounds/ASan regressions.

The initial sanitizer job stopped on handled driver API probe errors
inside OpenMPI/UCX (context queries and allocation-handle queries on host
buffers). Backtraces are retained in `slurm-312028.out`; no invalid kernel
access was reported. Repeating with `--tool memcheck --report-api-errors no`
passes with zero memory-access errors on all four ranks, including bitmap
tails. Backend CUDA statuses remain checked. This disables API-error
reporting, not memory instrumentation; see [NVIDIA's reporting modes](https://docs.nvidia.com/compute-sanitizer/ComputeSanitizer/index.html).

Three paired fresh-process samples on c116:

| Four-GPU path | First loop median (s) |
|---|---:|
| Direct multiply, CPU panels | 2.532 |
| Resident panels and result slot | **0.279** |

This is 9.07× faster. Finding 21's static model preserved measured waits
and upload costs inside the remaining phases; those also disappear with
resident payloads. That explains why the real result beats its ~0.8 s
forecast. This is a loop comparison, not a cold executable comparison.

The capture matches **192 panel +48 scalar NCCL kernels** and only 24
profiling MPI Reduces. NCCL-stream host staging drops from ~1.411 GB to
**576 bytes uploaded +768 bytes downloaded** for the whole solve. Its 192
local GPU copies still transfer 514.184 MB device-to-device. Other library
count/scalar accesses are not included in that communication-stream total.

NCCL panel kernels total 6.434 ms instead of 144.295 ms in the earlier
shared-communicator capture; scalar kernels total 0.292 ms instead of
155.804 ms. These include inter-rank waits: CPU preparation/arrival skew
previously made communication look slow. Captures are instrumented rank-0
runs; their overlapping durations are not additive wall decompositions.

### Touched bitmap

Optional `YUKAWA_BITMAP_ACCUM=1` marks every contributed position with an
atomic OR. Count/compaction inspect touched values; cancellations are
thresholded after all contributions. Later iterations clear previously
touched entries, including ones removed by threshold. Shape/mode changes
invalidate the bitmap and force full initialization. It adds one bit per
dense element: 128 MB/GPU beside 8.192 GB dense storage at 64k/2×2. Dense
O(N²/GPUs) storage remains; other matrix/fill regimes may have a different
crossover. The option is off by default.

Three paired five-solve processes on c116, same binary:

| Mode | First loop (s) | Later loop (s) | Five loops + resets (s) |
|---|---:|---:|---:|
| Full scans | 0.274 | 0.1825 | 1.009 |
| Touched bitmap | **0.190** | **0.0963** | **0.582** |

First-loop medians use three processes; later-loop medians use solves 2–5.
The bitmap improves later throughput 1.90×, not cold application latency
1.90×. Full-buffer memsets drop from 24 /196.608 GB to two /8.320 GB
(dense plus bitmap once), with 23 additional selective-clear shaders.
Those shader writes are not included in memset-byte totals. Selective
clearing itself takes 27.3 ms, approximately the old clear time; main gains
are count scan 40.8→4.13 ms and compaction 80.0→25.2 ms. Product time rises
6.88→10.77 ms from the extra atomic OR, a favorable net tradeoff here.

Actual panel/bitmap PTX and SASS have zero register spill loads/stores.
Bitmap multiply emits FP64 atomic add and integer atomic OR; compaction
uses vote/popcount and sorted writes. No local-memory loads/stores were
found in these kernels.

### Real repetition and setup scope

`--solver-reps N` restores cached GPU X0, resets convergence state, and
reuses communicators/allocations across independently converged solves.
Every repetition receives full reference parity and its own `.repN`
profiles. This is identical-H reuse; changing H needs new bounds/X0.
The CPU library recomputes its internal initialization each call.

On c117, five pure-CPU solves (4 ranks ×32 threads, 128 cores) total
**15.815 s**. Three fresh GPU processes, each with five resident solves,
give **5.796 s** median including measured binding, numerical setup and
resets: 2.73× for the numerical chain before the bitmap change. Bitmap
batches on c116 give 5.645 s median with setup despite much faster loops;
setup noise dominates short batches.

`gpu_chain_with_setup_seconds` includes explicit runtime binding,
communicator/warmup, spectral bounds, X0 construction/upload/snapshot,
state resets and loop time. MPI startup, I/O/permutation, input validation,
embedded oracle, parity checks and CPU output materialization are excluded.
The result is ready in GPU slots. This is not total executable wall time.

### One/two/four active GPUs

The final bitmap binary was swept on c116. The allocation reserved four
GPUs; each solve actively used 1/2/4. One/two GPU results use three fresh
five-solve processes; four GPUs is one control, with earlier three-process
bitmap comparisons independently confirming its ~0.096 s later-loop time.

| Active GPUs | First loop (s) | Later loop (s) | First solve + measured setup (s) | Five solves + measured setup (s) |
|---|---:|---:|---:|---:|
| 1 | 0.383 | 0.2844 | **1.949** | **3.098** |
| 2 | 0.244 | 0.1550 | 3.283 | 3.911 |
| 4 | 0.189 | **0.0963** | 4.973 | 5.365 |

First-solve prefixes include the measured reset within the five-solve
process; they are not separate standalone one-shot measurements. One GPU
uses local copies and no inter-GPU collectives. It wins the measured short
batch, while four GPUs gives best reused throughput. The CPU five-solve
baseline uses 128 cores on one node; the older 1.577 s per-solve CPU result
uses 512 cores across four nodes.

A constant-startup plus constant-later-solve fit predicts one GPU best
through roughly 11 identical solves, two for roughly 12–29, four thereafter.
Those crossover counts are unmeasured beyond five solves and sensitive to
setup noise, especially the single four-GPU control. Coefficients and
assumptions are preserved in the replay script, not universal GPU-count rules.

Artifacts: `validation/rikyu-resident-20261111/{analyze.py,summary.json}`
and archives under `rikyu-resident-repeat-*`, `rikyu-bitmap-*`, and
`rikyu-resident-scale-*`. Jobs: 312028 (build/tests/API probe reports),
312030 (memory/tiny/full/profile), 312034 (repetition/CPU), 312039 (bitmap),
312053 (sweep). Latest executable:
`~/agent/work/ntpoly-bitmap-20261111/build/bin/yukawa_gemm`.
Its tested-source manifest matches the current implementation.

**Next:** profile/move bounds and X0 setup to GPU; optimize bitmap traversal
if throughput matters; replace dense storage before larger-N claims; test
changing Hamiltonians and integrate this backend into NTPoly proper.


## Finding 23: Two-node NCCL and paired CPU scaling (2026-11-11)

The previous NCCL solver had a deliberate Fortran guard rejecting more than
four ranks or more than one shared-memory node. This was an implementation
restriction, not evidence of a RIKYU network failure. It remains the default;
`YUKAWA_NCCL_ALLOW_MULTINODE=1` explicitly enables experiments. The C++ NCCL
transport already bootstraps communicators through MPI and uses node-local
GPU binding in the solver. The transport correctness test now uses node-local
ranks too, and both CUDA tests accept eight ranks (resident test N=68).

### Matched CPU/GPU placement

Same gapped 64k input, seed 424242, threshold 1e-4, convergence 1e-6;
all solves converge in 24 iterations. Job 312261 runs CPU and GPU serially
in the same two-node allocation, c392/c393. Two fresh processes/configuration,
three solves/process; table pools solves 2 and 3, taking their median (four
observations). CPU uses the pure-CPU reference executable, 32 OpenMP threads
per MPI rank, `ppr:<ranks/node>:node:PE=32 --bind-to core`, OMP_PLACES=cores,
OMP_PROC_BIND=close. GPU uses `ppr:<ranks/node>:node --bind-to none`, one
thread/rank, shared NCCL communicator, resident panels and touched bitmap.
The one-node controls run only on the first node of the allocation.

| Nodes used | CPU cores | GPUs | Process grid | CPU TRS2 (s) | Initialized GPU loop (s) |
|---|---:|---:|---|---:|---:|
| 1 | 128 | 4 | 2×2 | 3.1444 | 0.09763 |
| 2 | 128 | 4 | 2×2 | 3.0436 | 0.09777 |
| 2 | 256 | 8 | 2×4 | 2.4528 | 0.06931 |
| 2 | 256 | 8 | 4×2 | 2.2207 | 0.06968 |

Four ranks split across two nodes retain 128 CPU cores/four GPUs total.
Neither backend suffers an abrupt penalty: CPU 3.144→3.044 s, GPU
0.09763→0.09777 s. Small differences should not be interpreted beyond this
sample size; placement also changes memory/cache resources, not only MPI.
Doubling resources improves the best CPU layout 1.416× and the best GPU
layout 1.409×. CPU favors 4×2; GPU 2×4 is slightly faster than 4×2 here.
This is healthy but sublinear strong scaling on this particular matrix.

CPU TRS2 includes its internal initialization. The GPU loop excludes binding,
communicator/warmup, bounds and cached-X0 setup; resets/output validation are
outside each loop. These are throughput comparisons, not executable wall
speedups. In the paired job, three GPU solves including measured binding,
numerical setup and resets have median totals 5.504 s local4, 3.020 s split4,
5.139 s eight2×4, 5.383 s eight4×2. They exclude MPI startup, file input,
permutation, oracle/parity and CPU output materialization. Identical-H reuse
is required. Setup variability and the CPU's repeated initialization remain
important when comparing short chains.

### Transport evidence and phase timing

Independent job 312254, c388/c391, first builds the opt-in binary and passes
both shared and separate communicator tests: fixed/variable/all-empty/scalar
NCCL and nonsymmetric/empty/reset distributed squares at N=68. Then two
fresh processes/configuration, three solves each, give later medians
0.09755 s local4, 0.09988 s split4, 0.06781 s eight2×4, 0.06964 s eight4×2.
All 48 GPU solves across the two jobs pass full matrix/scalar parity, with
errors below 8e-15. All 24 pure-CPU solves have matching energy/trace/grand-sum
fingerprints across layouts; their separate executable does not export a
full matrix in this benchmark. Host MPI/bounds/ASan regressions also pass.

NCCL 2.29.3 logs select `IBext_v11`, report GPU Direct RDMA via DMA-BUF,
and report four 800000 Mb/s virtual HCA devices. Across the first job's 80
NCCL logs there are no NCCL WARN lines or `Using network Socket` selections.
Bootstrap/OOB network traffic is separate from data transport. Log evidence
shows RDMA selection, not a direct measurement of wire bandwidth or a claim
that every access bypasses host memory. Default facility/module networking
is used; no forced NCCL HCA, protocol or network settings were needed.

Rank-0 phase totals, medians over the four later solves of job 312254 (ms):

| Grid / placement | Row exchange | Column exchange + assembly | Local multiply/count/compact | Rank-0 loop wall |
|---|---:|---:|---:|---:|
| 2×2 local4 | 4.54 | 6.21 | 70.09 | 96.26 |
| 2×2 split4 | 4.83 | 9.63 | 68.77 | 98.60 |
| 2×4 eight | 7.75 | 7.49 | 36.65 | 66.44 |
| 4×2 eight | 4.10 | 10.55 | 38.51 | 68.21 |

These include counts, synchronization and arrival skew; they are not pure
NCCL transfer durations or maximum-rank wall times. Independent medians
need not add. Local multiply/count/compact nearly halves on eight GPUs;
remaining exchanges, reductions, updates and synchronization limit total
scaling. Row panels and column panels are still serialized before multiply.
For rank0, a rough Amdahl model using 70.09 ms scalable local work plus
26.17 ms remaining time predicts ~61.2 ms on eight GPUs before extra panel
costs; observed rank0 time is ~66.4 ms. This is consistent with a modest
network/overhead penalty, not transport collapse. It does not identify the
CPU's bottleneck; CPU communication phases were not separately profiled.

Artifacts and replay: `validation/rikyu-multinode-20261111/analyze.py`,
`analysis.json`, both raw results directories, verified archives, scripts,
Slurm logs, source hashes and run metadata. Job 312254 completed in 6:51,
312261 in 7:01, each allocating eight GPUs/two nodes. The second allocation
runs CPU and GPU controls in sequence. Latest experimental executable:
`~/agent/work/ntpoly-multinode-20261111/build/bin/yukawa_gemm`.
The prior four-GPU executable/workspaces are preserved. No active jobs remain.

**Limits / next:** these are two node pairs and one 64k input. No 3/4-node
or larger/fill-heavier claim follows. Profile CPU communication separately,
reduce scalar/metadata synchronization or overlap panels to improve strong
scaling; retain default network selection until measurements justify tuning.

## Finding 24: Native TRS2 integration and protection against overdistribution (2026-11-11)

The fast resident backend now runs inside the ordinary library TRS2 routine,
not a separate copied solver loop. `NTPOLY_ENABLE_GPU_TRS2=ON` (default OFF)
adds the Fortran backend and shared CUDA/NCCL objects to `libNTPoly`.
`NTPOLY_GPU_TRS2=1` selects it at runtime. The `Matrix_ps` layout and TRS2
signature remain unchanged; the CPU branch stays available in the same
build. CMake >=3.20 policies/dependencies support CUDA objects attached to
the library across directories. The mini-app no longer defines duplicate
objects/modules when the library feature is enabled. Both native reference
and original mini-app link successfully. CUDA objects are reused, not copied.

### Ownership and coherence boundary

`PSMatrixGPUBackendModule.F90` provides GPU matrices with private dynamic
handles, type-bound deep-copy assignment and elemental finalization. The
CUDA slot table now uses a deque; new library owners are allocated beyond
the four reserved mini-app slots. Releasing an owner frees its CSR buffers
and recycles its handle. Existing references survive table growth. Scalar
and array assignment are checked with 41 simultaneous GPU owners, mutation
of the source after copying, finalization and a live-owner count. Host
ASan/UBSan tests additionally exercise 32 simultaneous owners/recycling.
The multiply and merge scratch remain process-wide caches; this is not a
fully concurrent GPU matrix backend yet.

The existing TRS2 loop, convergence monitor, sigma history and chemical
potential calculation are retained. Its trace, square, scale/increment/copy
and energy dot use GPU objects within an explicit region. WH/X0 setup and
initial/final basis transforms remain CPU. The result is materialized in
X_k before any host printing, undo-permutation or final transform. Thus the
normal return value is a usable distributed host matrix, with no stale
host/device aliases. Other library solvers and ordinary algebra entry points
are not yet dispatched to this backend. General resident AB multiplication,
complex matrices and pslices>1 remain outside this integration.

### Small matrices on many ranks

`NTPOLY_GPU_RANK_POLICY=auto` initially coalesces N<=10000 onto one active
GPU (the application process grid's rank zero), while larger matrices use
all ranks. `all` and `root` provide explicit controls. Input WH/X0 are gathered
once; nonparticipants wait outside the iteration, and output is redistributed
once using the original grid. Iteration count, energy and sigma history are
broadcast so every rank can compute consistent chemical potentials. Padding
stays zero and the normal matrix distribution is preserved. The one-GPU
region bootstraps NCCL on MPI_COMM_SELF; the full-grid region uses the
matrix's communicator. `InitSolverNCCL` now accepts an optional bootstrap
communicator, preserving MPI_COMM_WORLD for the old mini-app. Two independent
MPI subgroups are tested concurrently, each with its own active root.

This addresses the stated BigDFT regime (roughly 10k matrices with many
allocated nodes). It is an initial dimension-based policy, not an optimal
work/communication/memory model. Forced root is still limited by O(N²)
accumulator memory; all-rank mode requires one rank per GPU and experimental
multi-node opt-in above four ranks. No universal performance guarantee follows
from one matrix family or four nodes.

Job 312452: four nodes c338–c341, sixteen GPUs reserved, application rank
counts 1/4/8/16. A principal N=10000 submatrix is extracted from the existing
gapped 64k binary. Seed 424242, threshold 1e-4, convergence 1e-6, nelec=5000.
All cases take 18 iterations. Each configuration uses two fresh processes,
three solves/process; later medians pool solves 2/3 (four observations).
CPU uses one OpenMP thread per MPI rank, not the earlier 128-core baseline.
GPU automatic mode activates one GPU even when the application spans nodes.

| Application ranks / nodes | CPU TRS2 (s) | GPU auto (s) | GPU forced all (s) |
|---|---:|---:|---:|
| 1 / 1 | 3.39792 | 0.13791 | 0.10885 |
| 4 / 1 | 1.15398 | 0.08314 | 2.16636 |
| 8 / 2 | 0.65587 | 0.07539 | 2.12158 |
| 16 / 4 | 0.35520 | 0.07757 | 2.13115 |

These are normal TRS2 API times: internal CPU preparation, CUDA/NCCL setup,
CSR import, iteration, export, communicator teardown and final transform
are included. MPI launch, reading/permuting the input outside TRS2, scalar
validation and result file writing are excluded. Unlike the preceding
initialized mini-app numbers, the returned result is on CPU. Each native
call rebuilds bounds/X0 and communicators; it does not reuse cached X0.
At one rank, automatic coalescing still does redundant gather/redistribution,
explaining some overhead versus the direct all-rank path; avoiding these
identity transfers is a straightforward remaining optimization.

Median totals for three native calls, including the first call:

| Application ranks | CPU (s) | GPU auto (s) | GPU all (s) |
|---|---:|---:|---:|
| 1 | 10.16252 | 0.60271 | 0.55675 |
| 4 | 3.46554 | 0.43452 | 9.67776 |
| 8 | 1.96997 | 0.41679 | 9.33968 |
| 16 | 1.07185 | 0.41640 | 9.69984 |

At sixteen ranks automatic coalescing is 27.5× faster than forced-all for
later API calls and 4.58× faster than this 16-core CPU configuration. The
three-call batch is 2.57× faster than CPU. It stays nearly flat between
8/16 ranks; it does not promise latency independent of arbitrary rank count.
Much of the forced-all penalty is communicator setup, not network payload:
its initialized GPU loop remains ~0.03–0.04 s. Retaining communicators will
change the policy crossover and must be measured before tuning the cutoff.

Twelve complete saved density matrices (first process of every rank/policy)
are independently compared entrywise against the one-rank CPU reference;
worst relative column-1-norm error is **6.398e-16**. The replay validates
finite values, index ranges and duplicate-free entries. All 72 large-case
solves agree in iteration counts and scalar fingerprints. Four additional
16-rank tiny runs exercise all/auto mode on the world grid and two concurrent
eight-rank subgroups, including nonidentity diagonal inverse overlap,
load balancing, changed H, repeated calls, full matrix/energy/chemical
potential parity and owned-array cleanup. The subgroup logger can interleave;
numerical assertions run independently in each subgroup.

### Native 64k result and remaining setup cost

Job 312449, c333: 1/2/4-rank tiny native all/auto tests pass, 36 CPU/GPU
pairs in total, with full matrix/energy/chemical-potential and matching
iteration counts. The ordinary `trs2_reference` binary then calls native
TRS2 three times on gapped 64k. CPU uses four ranks ×32 threads; GPU four
ranks ×one thread, all-rank mode. Every solve takes 24 iterations and
matches the established energy/trace/grand-sum fingerprints.

| Native path | First API call (s) | Later API calls (s) |
|---|---:|---:|
| CPU, 128 cores | 3.3317 | 3.4364 / 3.1708 |
| GPU, four ranks | 6.7370 | 3.9922 / 3.8906 |

GPU phase counters: setup 5.781/3.106/2.952 s, loop
0.1955/0.1526/0.1368 s, export 0.01886/0.01849/0.01835 s. Setup counters
cover runtime/NCCL/import only, excluding preceding CPU bounds/X0 work.
They are per-phase maximum-rank times; coalesced counters may include
idle-rank arrival skew, so do not sum them into API wall time. The native
large-matrix path is not yet faster end-to-end than the threaded CPU:
rebuilding communicators each call consumes most of the advantage. Earlier
mini-app reuse results remain valid but have a different setup/output scope.

The initial build attempts 312434/312435 stopped on CMake target-policy/
custom-command dependency issues. Attempt 312436 then found a bad test
permutation: it covered actual N=5 rather than padded N=6. The fixture now
allocates the logical size, permutes actual indices, and fixes padding.
The corrected suite passes. Host regression checks (actual kernel bodies
under ASan/UBSan, MPI 1/2/4 legacy parity, bounds and repeated solves) pass.

Artifacts: `validation/rikyu-native-trs2-20261111` holds the native build,
first source snapshot, scripts, Slurm log and results. The final test adds
subcommunicator cases; its source manifest is also in
`validation/rikyu-native-small-20261111`, alongside verified archives, twelve
binary matrices, raw logs, `analyze.py` and `analysis.json`. Job 312449
completed in 1:48 after incremental build; 312452 completed in 5:17.
Executable: `~/agent/work/ntpoly-native-trs2-20261111/build/bin/trs2_reference`.
No active jobs remain.

**Next:** retain GPU communicator contexts with clear lifetime ownership,
remove redundant one-rank redistribution, calibrate active-rank selection
using workload/fill as well as N, and move setup onto GPU. Then generalize
resident AB and wire more native solver families into the backend.

## Finding 25: Work-aware selection and conservative memory admission (2026-11-11)

The fixed N<=10000 cutoff is replaced for native TRS2. Selection now runs
after CPU WH and X0 construction, so it can inspect the actual working
pattern without merging/copying payloads: sum local sub-block nnz counts
with 64-bit MPI reductions. `GPUPlanningModule.F90` isolates the arithmetic
from MPI/CUDA so ranking, memory budgets and overflow checks are testable
on the host. The native solver still uses its existing loop and API.

### Candidate ranking

Candidates remain one active GPU or the application's complete 2D grid.
For each, the ranking model uses dimension/layout, initial average degree
nnz(X0)/N, an N*degree² product-work estimate, bitmap/CSR traversal, launch
latency, panel overhead and native communicator startup. Initial coefficients
are fitted to the existing 10k/64k regimes, with up to 24 expected iterations
(capped by the requested maximum) and ~2 s full-grid communication setup.
They are a selection heuristic, not measured timings or an exact roofline.
It does not yet score CPU execution, select intermediate subgrids, learn
online, or anticipate all irregular/fill-growth behavior. Communicator
reuse will require recalibrating these costs.

Sparse N=12000 now selects one GPU, demonstrating removal of the old size
boundary. Pure unit tests also change the ranking by raising fill while
keeping N=10000. A one-rank application in automatic mode uses the direct
layout, eliminating Finding 24's identity gather/redistribution.

### Memory and index admission

Use actual N for the root candidate and padded logical N for the full grid.
For m=N/pr, k=N/pc, E=m*k, the estimate contains:

- 8E bytes for the numerical dense accumulator;
- 4m*ceil(k/32) bytes for the bitmap (also reserved when the option is off);
- 36E bytes for X, X² and merge scratch, allowing full local CSR fill;
- 12E*(1+2pc+pr) bytes for local CSR, row receive, column receive and A panels;
- 12*global_nnz(WH), a conservative bound on fixed WH storage on every rank;
- 512*(N+m+k+1) auxiliary bytes and a 1 GiB runtime/protocol allowance.

These are worst-case payload bounds plus a runtime allowance, not a formal
upper bound on all opaque vendor allocations. They avoid trusting the
initial sparsity to remain unchanged during purification. E, pc*E and
pr*E must fit the current 32-bit CSR counts. Invalid grids and overflowing
estimates are rejected rather than wrapped. Input WH storage does not
assume perfect load balance. Scratch already cached in the process reduces
reported free memory, making admission more conservative.

Budget defaults to **0.8 * currently free GPU bytes**, limited further by
`NTPOLY_GPU_MEMORY_MB` (nonnegative integer MiB). A supplied cap cannot
increase availability; the lowest supplied rank cap is used for selection.
Overflowing/negative configuration values are rejected. Automatic mode
queries only root when it is safe and cheaper, avoiding CUDA contexts on
inactive ranks. All-rank mode queries active GPUs and uses the smallest
budget. UUIDs, not device ordinals alone, detect duplicate GPU mappings
within a node; this accommodates per-rank visibility of distinct devices.
All-rank candidates also honor the existing multi-node NCCL opt-in.

A rejected candidate does not enter its GPU iteration. `root`/`all` are
layout requests, not overrides of safety. If no permitted layout fits,
TRS2 executes its existing CPU iteration using the already-built WH/X0.
Admission is rechecked collectively after CUDA/NCCL warmup, before importing
matrix payloads; a failed recheck destroys the GPU communication region
before continuing on CPU. Rejection flags are reduced across the complete
application communicator so ranks make the same decision. The post-warmup
recheck's allowed path is exercised; the GPU tests below force pre-bootstrap
rejection, not an injected post-bootstrap allocation failure.

There is no live VRAM reservation, hard allocator quota, mid-iteration
repartition or arbitrary-OOM recovery. Another process can allocate after
admission; runtime internals can exceed the allowance. Existing CUDA/NCCL
allocation/error checks remain active. CPU setup/fallback memory and raw
mini-app helper calls are outside this admission policy. The model may
reject sparse layouts that would fit in practice; bounded/chunked workspaces
are the next way to relax this conservatism safely.

### GPU and CPU-fallback validation

Validation uses eight-GPU/two-node allocations with one, four or eight
application ranks.
The first job, 312579 (c080/c081), passed twelve two-call cases before its
eight-GPU capped case timed out during startup, before iteration payloads.
Job 312609 (c392/c393) then tested the old native binary and the new binary
in the same allocation. Both eight-GPU runs passed; the original timeout
was not reproduced. No cause is established and this is not evidence of a
RIKYU outage. Retry validation completed before packaging failed separately
with `Disk quota exceeded` in the remote home directory.

Observed choices for gapped principal submatrices / the established 64k input:

| Case | Selected execution |
|---|---|
| 10k or 12k, normal budget, four application ranks | One GPU |
| 10k, one application rank | Direct one-GPU path, no coalescing copies |
| 10k, auto, 6 GiB cap, four application ranks | Four GPUs |
| 10k, auto, 6 GiB cap, eight application ranks | Eight GPUs |
| 10k, zero budget | CPU |
| 10k, forced root, 6 GiB cap | CPU |
| 10k, forced all, 1 MiB cap | CPU |
| Eight ranks mapped to four physical GPUs on one node | CPU |
| CUDA_VISIBLE_DEVICES empty | CPU |
| Eight ranks across nodes, 6 GiB cap, no multi-node opt-in | CPU |
| 64k, normal budget, four application ranks | Four GPUs |

For 10k, estimated root peak is **10,305,974,640 B**; full-grid peaks are
**4,291,474,640 B (four ranks)** and **3,288,634,640 B (eight ranks)**.
Thus a 6 GiB cap excludes root but admits distribution. These are admission
estimates, not observed HBM traffic or measured usage. For 64k/2×2, the
estimate is **132,370,380,176 B/rank**, against a minimum live budget of
~157.45 GB in the retry. The root 64k plan is rejected by CSR index capacity.

All 40 completed native calls across the twenty large cases agree with CPU
fingerprints: 10k takes 18 iterations, 12k takes 20, 64k takes 24. Nineteen
complete saved 10k/12k output matrices are independently compared against
the CPU outputs, with worst relative column-1-norm error **2.372e-15**.
The 64k output is checked by scalar fingerprints, not a new full-matrix
comparison in this experiment. Zero/unsafe-cap fallback logs show no NCCL
initialization and no GPU phase counters. Two corrected tiny native suites
also pass full matrix/energy/chemical-potential/iteration parity, nonidentity
overlap, load balance, changed H, repeated calls and owner cleanup.

Illustrative second-call native API times (one process, two calls each;
not a statistical speedup claim): 12k root **0.0992 s**, CPU reference
**1.4245 s** from the initial allocation (one thread/rank, four ranks).
New eight-GPU capped control **1.9714 s**, old control **2.0712 s** in the
retry allocation. Setup still dominates full-grid native calls. The 64k
retry's second API call is **2.6247 s**, but no matched threaded-CPU run
was made in that allocation. Selection/query overhead is included in API
time; queries can initialize CUDA before the narrower GPU-setup counter.
Do not compare these as initialized-kernel timings or claim a new paired
128-core CPU speedup.

Host tests pass with Fortran bounds/overflow checks and the existing actual
operation-body ASan/UBSan/MPI regressions. The pure planning tests cover sparse
versus dense ranking at fixed N, constrained budgets, caps that cannot raise
live availability, invalid dimensions and large integer/index limits. NVHPC
builds/runs the same planning test on RIKYU.

### Artifacts and status

`validation/rikyu-memory-policy-20261111` contains both stage summaries,
first/final source manifests, the first module snapshot, scripts, Slurm logs,
43 raw transferred files and a transfer manifest with verified SHA-256.
`analyze.py` independently rechecks source/transfer hashes, fingerprints,
plan admission inequalities and full matrices, producing `analysis.json`.
An incomplete remote archive remains: the filesystem deletion tool requires
explicit exact-path confirmation, so it was not removed. All raw files are
preserved locally and remotely; archive failure did not invalidate solver
checks. Neither scheduler job is marked successfully completed: the first
timed out in NCCL; the retry's solver validation passed but packaging failed.
No active jobs remain. Latest executable:
`~/agent/work/ntpoly-memory-policy-20261111/build/bin/trs2_reference`.

**Next:** communicator lifetime/reuse and cost recalibration, intermediate
active grids, broader N/fill/imbalance measurements, and a bounded sparse or
chunked workspace to reduce the pessimistic full-fill admission requirement.

## Finding 26 — cuSOLVERMp diagonalization comparison (2026-11-11)

The standalone `cusolvermp_bench.cu` benchmarks installed cuSOLVERMp 0.7.2
Syevd with eigenvectors, cuBLASMp 0.8.0 occupied-state FP64 GEMM, and sparse
host coordinate export at 1e-4. The input is the same gapped 64k Hamiltonian
and its principal 10k submatrix; half filling and identity overlap match TRS2.
TRS2 uses threshold 1e-4, convergence 1e-6, and the fixed random permutation.
Diagonalization uses the original ordering; energy/trace comparisons are
permutation invariant. Dense setup/file reading is outside repeated timings.

Later-call means (two calls after repetition zero, seconds):

| Matrix | cuSOLVERMp active GPUs | Dense EVD + density + host export | Native GPU TRS2 API, same allocation | Ratio |
|---|---:|---:|---:|---:|
| 10k | 1 | 0.7873 | 0.08124 (1 active GPU, 4 app ranks) | 9.69 |
| 10k | 4 | 4.6891 | 0.08124 | 57.7 |
| 64k | 4 | 58.4877 | 2.6062 | 22.4 |

These are workload-specific, different-accuracy comparisons. Dense output is
host COO; native output is a normal distributed NTPoly matrix. Dense COO
sorting/CSR construction/redistribution is excluded, making its comparison
optimistic. Dense contexts/workspaces persist across repetitions; native
TRS2 still initializes contexts per call. The native timings include normal
setup, loop and export; **do not use its ~0.134 s loop as the headline API
speedup**. Small-matrix coalescing is important for both algorithms.

Job 312738 on c396 used block 512. At 64k, later EVD, density GEMM and export
average 56.5214, 1.9451 and 0.02008 s. The density GEMM performs
2*N*N*nelec operations, about 134.8 TFLOP/s aggregate; the eigensolve is the
main expense. Explicit A/Q/workspace/coordinate allocations are about
50.22 GB per GPU, excluding internal allocations and temporary output. Dense
setup is 5.1822 s; first complete solve 62.5899 s versus first native API
5.6022 s. At 10k/1 GPU, setup is 1.0949 s and first complete solve 0.83066 s;
the paired native first call is 0.50505 s. No setup is silently amortized
into the reported later means. Initial job 312731 on c388 used block 256:
64k later EVD+density 59.0543 s, native API 2.5234 s. It used accepted jobz Z
and omitted sparse export; the final source uses documented jobz V and full
export timing. Preserve the earlier source/hash separately.

Bound CPU128 control (4 MPI x 32 threads, core binding, job 312746/c393)
passes the 64k fingerprint and averages 3.1516 s, with paired native GPU
2.6065 s. Unbound CPU controls are slower and are not the selected baseline.
The 10k CPU32 controls are **numerically invalid**, despite exit status zero:
trace 4978.6163 without explicit binding, 4932.4405 with binding, instead of
5000.00002794. Exclude their timings. Valid one-thread CPU10k controls
(job 312752/c191) give 3.4339 s on one rank and 1.1463 s on four ranks.
This extends the padded threaded-CPU correctness concern beyond OMP36;
root cause is unresolved and requires a separate minimal regression/fix.

Every dense solve has info=0, trace/nelec and norm-squared/nelec within
roundoff, and Tr(H D) agrees with the occupied eigenvalue sum to roundoff.
An independent CPU check at N512 on 1/4 GPUs gives relative eigen-residual
6.74e-17, orthogonality error 5.03e-15, and zero projector entry error.
These exhaustive checks are limited to the small input. At 64k, exact dense
energy is 15947.3141381810 versus TRS2 15947.3151489420: absolute difference
0.0010107610, relative 6.34e-8. At 10k the energy difference is 0.0001195884
(4.80e-8 relative). TRS2 prunes intermediates and is not exact EVD; no full
large-matrix norm parity or equal-accuracy speedup is claimed.

`CUSOLVERMP.md` documents building, memory and timing scopes. Raw text,
scripts, input/source/binary hashes, statuses and the replayable analyzer
are in outer `validation/rikyu-cusolvermp-20261111`. Jobs 312731, 312738,
312746 and 312752 all completed; 312730 failed immediately on nvcc's rejection
of mpicxx's -pthread flag, fixed by separating compilation and MPI linking.
One intermediate diagnostic line in 312746 contains malformed binary bytes;
the exact raw log is preserved locally/remotely with SHA-256, while the
tracked decoded text escapes control bytes. Final fingerprints are intact.
Do not interpret the filesystem tool's UTF-8 decoding error as SSH failure.

Installed 0.7.2 lacks registered-workspace APIs. NVIDIA's newer
[release notes](https://docs.nvidia.com/cuda/cusolvermp/release_notes/index.html)
describe improved SYTRD/STEDC and Blackwell registered-workspace performance;
this experiment does not measure those releases. Testing a newer cuSOLVERMp,
larger GPU grids, representative BigDFT matrices, and matched accuracy remain
follow-ups. Native communicator reuse and the padded CPU bug remain open.

## Finding 27 — grid-owned NCCL reuse in native TRS2 (2026-11-11)

Native solves now borrow communication contexts from their owning process
grid. Full-grid and coalesced-root contexts have separate lazy slots; the
first use initializes/warms NCCL, later uses reactivate it. Each grid owns
its communicators, stream, ordering events and small staging buffers;
`DestructProcessGrid` releases them before MPI handles, and `CopyProcessGrid`
starts an independent lifetime. This avoids matching recycled MPI handles or
silently borrowing a different application's communicator. Concurrent/nested
GPU calls remain unsupported; multiply scratch is still process-wide.
Matrix owners and panel payloads are freed at each solve boundary. Cached
communication has no H, bounds, X0, numerical matrix or dimension payload.

### Measured native API costs

At gapped 64k/half filling, threshold 1e-4 and convergence 1e-6, five solves
on one node/four GPUs (job 312787, c183, OMP1) give:

| Scope | Previous code | Grid-owned reuse |
|---|---:|---:|
| First call | 5.6599 s | 5.3479 s |
| Later four calls, mean | 2.7842 s | 0.67184 s |
| Later range | 2.6275–3.1021 s | 0.5654–0.9193 s |
| Sum of five API calls | 16.7966 s | 8.0353 s |
| Later narrow GPU-setup counter | 1.8153 s | 0.02909 s |
| Later resident loop | 0.13262 s | 0.13809 s |

The measured gain is **4.14x for later API calls**, with essentially the same
GPU loop; this removes repeated bootstrap/destruction, not numerical work.
The same allocation's CPU128 control (4 MPI x 32 bound threads) averages
3.2369 s on later calls, 4.82x the new GPU API cost. These are full TRS2 API
times, not loop-only timings; external input read/permutation, MPI launch,
matrix validation and disk writes are excluded. First-call cost is still
large and must not be silently amortized. Final communication teardown now
occurs at grid destruction, outside the API sums. Recorded five-solve executable
wall times, including launch, input, density write and teardown, are 26.60 s
old versus 18.09 s new; these include filesystem costs.
At 10k/coalesced one GPU, later
means are 0.08048 s old and 0.08221 s new; no useful speed improvement there.

A separate bound-affinity control (312802, also c183) gives later API means
0.73371 s for GPU/OMP1, 0.91378 s for GPU/OMP32, and 3.13395 s for CPU128.
OMP32 is numerically valid at this unpadded 64k size, but slower for GPU;
retain OMP1 for these GPU runs. The measured five-call API sums are 9.8774 s
GPU/OMP1 versus 15.9290 s CPU128, a 1.61x batch gain after including first use.
This does not resolve the separate padded 10k threaded-CPU bug from Finding 26.

### Correctness, lifetime and transport checks

One-node tests cover auto/root coalescing, full 2x2 and 1x4 grids, interleaved
independent application communicators, nonidentity overlap, load balancing,
changed H, repeated CPU/GPU alternation, deep GPU matrix ownership,
independent grid copies, alternating grids, reshaped grids and teardown.
Each process asserts zero owned contexts after final grid destruction.
Full 10k/64k old/new density comparisons give maximum entry differences
6.67e-16/8.11e-15. Separate CPU128 versus GPU/OMP1 and GPU/OMP32 64k matrices
agree within 6.67e-15. The large matrices are retained on group storage;
the exact comparison scripts and raw error maxima are in the text records.

Two-node/eight-GPU validation succeeds in job 312812 on c387/c390, using
`NCCL_IB_MERGE_NICS=0`: old/new later native API means 2.3825/0.35887 s,
resident loops 0.10636/0.10597 s. Both use the same setting/allocation.
First new API call is 5.0276 s. Full-grid and split-application tiny suites
also pass full CPU/GPU matrix parity, grid switching and zero-context teardown.
Large two-node solves pass energy/trace/grandsum fingerprints; they did not
save large matrices for a separate full entrywise comparison in this run.

Initial multi-node job 312792 timed out during first bootstrap on c186/c188,
before reuse. A paired control on the same nodes (312805) times out with
**both unchanged old and new executables** (code 110, no completed solves).
Its scheduler job exits zero because the driver records these failures and
continues; do not confuse that with numerical validation. Logs show asymmetric
merged NIC descriptions, but that is not a proven cause. The successful run
changes both node pair and merging setting, so it **does not establish a
transport fix or a causal effect of disabling NIC merging**. No cluster outage
or regression from reuse is diagnosed. Library network defaults are unchanged;
multi-node remains opt-in. NVIDIA documents the test setting in its
[NCCL 2.29 environment reference](https://docs.nvidia.com/deeplearning/nccl/archives/nccl_2297/user-guide/docs/env.html).

Legacy compatibility job 312829 passes actual NCCL fixed/variable/empty/scalar
transport checks at four ranks, and resident nonsymmetric/empty/reset tests
on both 1x2 and 2x1 grids, with shared communication disabled and enabled.
CPU-only NTPoly/reference targets build, and the full host regression suite
passes. An all-target CPU build still encounters the existing unlinked GPU
mini-app symbols without CUDA; this is not claimed as a successful build.

`validation/rikyu-context-20261111` preserves source/transfer hashes, specs,
scheduler statuses, raw outputs, both failed/successful transport logs, and
an independent analyzer checking fingerprints, context reuse and timing scopes.
All six jobs finished. Generated matrices/builds stay on group storage;
home quota is not consumed by new builds or density output. Next: bound the
quadratic workspace, recalibrate placement using cache state, and move CPU
bounds/X0/setup to GPU before generalizing resident operations/solver families.

## Finding 28: OpenACC cold multiplies dominate; PTX has no spills (2026-11-11)

Nsight Systems job312960 on c119 captures the validated OpenACC mini-app's
embedded NTPoly reference and the separate native four-GPU library backend.
The OpenACC path really launches GPU kernels: 26 real MultiplyBlock launches
total47.953 s, but three cold launches alone take47.686 s (99.44%).
Their durations are11.312/17.977/18.398 s; median across26 is10.849 ms.
PruneList totals0.4332 s. The reference API reports47.135 s, with correct
energy/trace. The capture includes input-permutation products and a later
legacy CUDA solve; these must not all be assigned to the reference API.
OpenACC Compute and Wait events overlap and must not be summed.

Cold launches coincide with thousands of managed-memory GPU fault groups.
Whole-capture CPU faults number1,287,988 versus453 for native4.
This motivates explicit workspace residency/first-touch investigation; fault
counts do not measure unique pages or prove a migration-time breakdown.
A separate OpenACC fixed H-squared test, one warmup plus three timed rounds,
has0.58353 s median. Its runtime reports a13.1519 s cold multiply, versus
0.36774 s minimum, but warns CUPTI is unavailable: device-time subfields are
unreliable. Fixed H-squared differs from late TRS2 squares.

Recompilation with NVHPC26.3 keep/lineinfo/ptxinfo and inspection of extracted
PTX/SASS show112 physical registers for real multiply (118 complex), and
zero stack/spill loads/stores. Runtime confirms112regs,608x128 launch.
The outer column loop is parallel; both inner sparse loops remain sequential.
Threads varying the second index of Fortran value_array(idx_b,II) access
columns512 KiB apart at64k. Source/PTX therefore suggest poor locality and
coalescing; achieved bandwidth/occupancy were not measured. Static PTX has
seven signed-divide sites for hash insertion, FP64 FMA, descriptor/indirect
loads, no local-memory accesses and no calls. Static unrolled sites are not
dynamic counts per update. There is no evidence supporting spill tuning.

Native4 has no OpenACC launches. Baseline rank-zero kernel time over three
solves is62.94% bitmap clear+compact,12.3% product and6.7% NCCL send/receive.
Eliminating clear+compact entirely could improve kernel time at most2.70x;
API/startup costs remain. MPI initialization/input I/O dominate whole-process
MPI totals; resident payload communication uses NCCL. Detailed reports,
per-launch SQL audit, source/resource interpretation and replayable analyzer
are in outer validation/rikyu-hotspots-20261111.

## Finding 29: resident AB, GPU bounds/X0 and a useful CPU setup fix (2026-11-11)

The resident panel product now accepts independent A/B, retaining the square
ABI. Both inputs must have matching square global layouts, with rectangular
local blocks supported. Either output/input alias is rejected. Nonidentity
working-H and final density transforms use two resident products when
collective memory/index admission succeeds; unsafe cases retain CPU fallback.
Transform boundaries still import/export normal host matrices.

After working-H construction/balancing, admission uses WH plus an upper
bound for the identity-shifted X0 pattern. GPU computes Gershgorin column
contributions; only O(local columns) vectors go to host MPI reductions.
Root coalescing explicitly retains the CPU padded-zero-column convention.
X0 is copied/scaled/shifted on GPU, avoiding CPU construction and a second
full payload upload. A sparse identity is still imported. CPU balancing,
host public matrices and other solver families remain.

New optional NTPOLY_TRS2_PROFILE=1 stage counters exposed a CPU identity fill
that scanned every local row/column pair looking for diagonal entries.
FillMatrixIdentity now visits only the owned diagonal intersection, clipped
to actual dimension. It preserves real/complex, padded and empty-block
semantics. In the paired 64k GPU control, the allocation stage fell from
0.328/0.593 s (identity/banded, third call) to 1.90/2.04 ms.
This was materially more useful than tuning a small GPU kernel.

Job313129 completed on c117 with four GPUs. CPU is four MPI ranks x32 bound
threads; GPU uses one thread/rank. Each fresh executable makes three calls.
Threshold1e-4, convergence1e-6, actualN64000:

| Full TRS2 API | Identity overlap | Banded inverse-overlap |
|---|---:|---:|
| New CPU, later mean | 2.8100 s | 5.3639 s |
| Previous GPU, later mean | 0.5881 s | 1.0403 s |
| New GPU, later mean | 0.19795 s | 0.37275 s |
| GPU first call | 6.2316 s | 6.5839 s |
| New CPU/GPU later speedup | 14.20x | 14.39x |
| New CPU/GPU three-call speedup | 1.30x | 2.22x |
| Previous/new GPU later speedup | 2.97x | 2.79x |

The previous GPU binary in this table already has the new AB transforms,
but retains CPU bounds/X0 and quadratic identity fill; it is the separately
hashed initial revision, not an unchanged 7fb3231b executable. Both binaries
run on the same allocation. API max-rank times include internal setup,
transfers, iteration, export and final transform; exclude launch/input,
external permutation, full-matrix checking and output writes. Final grid
communication teardown is outside the API. Cold startup remains expensive.

The banded input factor is SPD tridiagonal I+0.01T in the externally permuted
input ordering, constructed outside timing. It represents a supplied
inverse-overlap factor, not a separately timed inverse square root.
Physical trace(K) need not equal electron count in this basis.
Full CPU/GPU density matrices agree within6.22e-15 identity and3.34e-16 banded;
energy/trace/grandsum agree. Tiny tests also compare chemical potentials and
full matrices with nonidentity overlap, reverse balancing, changed H,
auto/all placement, split application communicators, CPU alternation,
zero-cap fallback, ownership and zero-context teardown.

Resident A/B tests pass with shared NCCL off/on at1/2/4 ranks, both1x2/2x1,
independent nonsymmetric matrices, zero A/B, pruning, reuse after empty and
both alias rejections. Identity tests pass real/complex, N1/5/8,1/2/4 ranks
and both grid orientations. GPU positive padded bounds match CPU conventions
in full/root modes. CPU library/reference/identity targets and the full host
regressions pass. The existing padded threaded CPU10k failure is not fixed;
its timings remain excluded. Latest code has single-node validation.

Updated native Nsight job313171 on c365 confirms bounds takes1.5–2.2 ms per
call, under1% of warm API time. It uses26 registers with no spills.
Clear+compact still dominate kernel time, but eliminating them entirely gives
only an estimated1.34x warm identity API ceiling and1.22–1.24x banded,
assuming all that work is on the critical path. The cold-call ceiling is
about1%. Optimize useful end-to-end headroom: prioritize startup and remaining
transform/host handoffs, defer bounds microtuning. No further speculative
kernel changes were made. Full scripts, hashes, raw text, per-call audit and
replayable analyses are in outer validation/rikyu-basis-20261111 and
validation/rikyu-hotspots-20261111.

## Finding 30: borrow accessible CSC input and skip redundant working-H copies (2026-11-11)

The caller may build H on CPU or GPU. Native import now checks CUDA pointer
accessibility and borrows existing CSC arrays during a synchronous GPU
conversion to the backend CSR layout. It supports managed allocations,
current-device allocations, mapped host memory, and system allocations when
the device reports pageableMemoryAccess. It does not touch source payloads
on CPU or copy them into input staging buffers. Source pointers do not escape
the call; conversion output is independently owned. Address-based numerical
caching is not introduced. Unsupported allocations retain host conversion.

For one local block, import skips CPU block composition as well. Multiple
local blocks still use CPU merge; root coalescing still uses host gather.
With identity overlap and no balancing, GPU TRS2 borrows H as its working-H
input rather than copying it onto CPU. If admission selects CPU fallback,
WH is materialized before the CPU algorithm. Identity with balancing retains
a separate working matrix and has a new regression case.

NTPOLY_GPU_IMPORT=auto is default; host is the layout-conversion control.
NTPOLY_GPU_IMPORT_REPORT=1 emits diagnostic borrowed-import markers; final
performance runs disable it. The existing stage profiler also reports
identity-input borrowing. Runtime pointer-query errors remain fatal rather
than being mistaken for unsupported memory.

### Matched one-node measurements

Job313383 completed on c118, four GPUs, GPU4xOMP1 and bound CPU4xOMP32.
Each GPU executable makes five API calls; warm means calls2–5, retaining
outliers. Old is unchanged b5fe68ce. Host control includes merge/identity-copy
removal but retains CPU CSC-to-CSR conversion.

| Full API | Old warm mean | New host-conversion control | New auto warm mean |
|---|---:|---:|---:|
| 10k identity | 60.66 ms | 59.83 ms | 58.42 ms |
| 64k identity | 207.25 ms | 183.54 ms | 183.31 ms |
| 64k banded | 379.21 ms | 330.89 ms | 322.71 ms |

The 64k gains are 1.13x identity and 1.18x banded for the full API.
Warm medians are 198.95→179.29 ms and 379.10→311.69 ms; both means and
medians are retained because isolated slower calls affect these short runs.
The host/auto comparison shows most mean improvement comes from removing
redundant CPU composition/copies, rather than the GPU conversion alone.
The new path additionally provides required device-accessible input coverage.
Bound CPU banded mean is 5.39385 s, 16.71x the new GPU warm mean.

Cold 64k calls remain expensive: old/new 6.021/6.034 s identity,
6.163/6.258 s banded. Five-call sums are 6.850/6.767 s and 7.680/7.549 s:
only 1–2% batch gain after including startup. At10k the five-call sum
increases from0.5045 to0.5162 s despite slightly faster warm calls.
Do not claim a robust small-matrix or cold-start improvement.

Full old/new density differences are at most6.11e-15 identity64k,
4.45e-16 banded64k and5.56e-16 identity10k. CPU/GPU banded scalar
fingerprints agree. Tiny full matrices/chemical potentials pass auto/all,
1x4, independent communicators, identity/nonidentity balancing, changed H,
CPU alternation, zero-cap fallback and teardown.

### Memory checks and coverage limits

The real CUDA test exercises24 allocation/dimension cases: managed, device,
system and mapped host memory; rectangular shapes, zero nnz/empty rows and
columns, sorted output, slot reuse, and repeated device-only value changes
on a separate nonblocking stream. All pass; source values are not copied
back to CPU before conversion. Compute Sanitizer reports zero errors when
NV_COMPUTE_SANITIZER_SHARED_ADDRESSING_SUPPORT=force is used.

Initial job313375 stopped on system-memory accesses under sanitizer auto
detection; job313380 used obsolete --hmm-support and exited before tests.
The force setting resolved the instrumented run without implementation
changes. Both failed logs are preserved; these are not cluster outages.
CPU-only library/reference targets and the full host regressions pass.

Two-node/eight-GPU job313385 passes small/all/split suites and full10k/64k
matrix comparisons on c116/c117, with NIC merging disabled for both binaries.
Later old/new means are55.19/54.75 ms at10k,1.1032/1.1012 s identity64k
and0.9563/0.9523 s banded64k: no meaningful large-case improvement.
The identity iteration stage alone takes1.059 s, so import optimization
does not explain or remove the dominant cost in that allocation/layout.
Keep the earlier unresolved default-NCCL transport caveat.

The direct input path is not a whole-library no-copy guarantee. It still
performs required layout conversion, leaves multi-block merging/root gather
on CPU and exports normal public matrices. Further work should follow
measured full-API headroom and required coverage. Evidence, specs, source
hashes, all failures/successes and replayable analysis are in outer
validation/rikyu-import-20261111.

## Finding 31: matched eight-GPU slowdown is NCCL/wait dominated (2026-11-11)

Follow-up job313441 completed on c117/c118 in one two-node allocation.
The same new executable runs four GPUs on one node, then eight GPUs in
2x4 and4x2 process grids, with explicit placement maps saved. All use
OMP1, NIC merging disabled and the same64k identity-overlap input.

| Layout | Warm full API mean |
|---|---:|
| Four GPUs, 2x2 | 0.17623 s |
| Eight GPUs, 2x4 | 1.08754 s |
| Eight GPUs, 4x2 | 1.02818 s |

All energy/trace/grandsum fingerprints pass. The two-node numerical tests
in Finding30 also save full matrices; this follow-up checks fingerprints
and timing rather than another full matrix dump. Reversing grid orientation
does not remove the slowdown. Prior old/new eight-GPU controls are equally
slow, so this is not evidence of a regression from borrowed imports.

The rank-zero eight-GPU Nsight capture over three solves records2.50445 s
in584 NCCL send/receive kernels and0.46803 s in145 scalar allreduce kernels:
about96% of GPU kernel time is communication/wait. Product totals16.1 ms,
clear42.3 ms and compact38.9 ms. Kernel lifetime includes waiting for peers;
it is not a direct measurement of physical network transfer bandwidth.
No transport debug log was captured, so IB/RDMA versus socket selection
and rank-arrival imbalance are not established by this trace.

The next performance investigation must inspect transport selection and
per-rank arrival skew on a matched allocation. Local kernel/PTX tuning
cannot address most of this runtime. The placement model still chooses only
root or all ranks; at64k the single-GPU index bound rejects root, so it
cannot retain four active GPUs when the application uses eight ranks.
Intermediate active grids are a concrete remaining scalability feature.
This result does not invalidate the earlier successful faster two-node
runs on other nodes/layouts, nor establish a cluster-wide transport failure.
Exact maps, raw profiler reports and analysis are in the scaling subdirectory
of outer validation/rikyu-import-20261111.

## Finding 32: cuSOLVERMp strong scaling and native profiling (2026-11-11)

Job313506 completed collection on c102/c117 with eight allocated GPUs.
All runs use the same binaries, input matrices, OMP1, explicit placement and
NIC merging disabled. Dense uses installed cuSOLVERMp0.7.2/cuBLASMp0.8.0,
block512, FP64 eigenvectors, occupied-density GEMM and sparse host export.
Warm means exclude the first of three solves.

| N | Application ranks | Dense complete | Native TRS2 API | Native active GPUs |
|---|---:|---:|---:|---:|
| 10k | 1 | 0.7871 s | 0.05126 s | 1 |
| 10k | 2 | 3.0773 s | 0.06190 s | 1 |
| 10k | 4 | 5.6551 s | 0.05663 s | 1 |
| 10k | 8 | Timeout, no result | 0.05392 s | 1 |
| 64k | 1 | Workspace budget rejection | Not admitted by current index bounds | — |
| 64k | 2 | 94.9443 s | Not admitted by current index bounds | — |
| 64k | 4 | 59.7541 s | 0.17346 s | 4 |
| 64k | 8 | Timeout, no result | Bootstrap timeout, no result | — |

Dense64k two→four gives1.589x speedup,79.4% efficiency relative to ideal
doubling. EVD91.1579→57.7893 s, density3.7460→1.9439 s and export
0.0380→0.0197 s. The dense four-GPU result is97% EVD; host export is
negligible. Dense10k worsens with more GPUs; four are7.18x slower than one.
Native10k is a coalescing/CPU/redistribution curve, not GPU strong scaling.

The eight-GPU512 dense correctness test passes, but no large dense result
completes within585 s, and no10k dense result within175 s. Native64k
times out during NCCL bootstrap after135 s. These cases are preserved as
failures and excluded from speedup calculations; collection exit zero does
not mean every solver case passed. The dense harness prints phase times
only after a complete pass, so its timeout does not identify the failing
phase or establish the same cause as native bootstrap.

The prior completed matched native4→8 result from Finding31 remains a
separate allocation:0.1762→1.03–1.09 s, about8% scaling efficiency with
NCCL/wait dominance. Thus four GPUs/one node is the useful observed64k
configuration in this deployment; no universal16-GPU or larger-N limit
is inferred. Root/all selection cannot choose four active GPUs from an
eight-rank64k application; intermediate placement remains important.

Dense complete sums import/EVD/density/export phase maxima; native complete
is the normal maximum-rank API wall. Persistent dense setup/validation and
external native read/permutation/validation are excluded. Cold scopes differ,
and exact dense FP64 versus thresholded TRS2 is not equal accuracy.
Later cuSOLVERMp releases have untested algorithm/Blackwell improvements;
these are measurements of the installed version.

The prior native warm four-GPU profile breaks183/323 ms identity/banded
into approximately117/146 ms iteration,25/29 ms export/GPU cleanup,
31/90 ms final transform/cleanup,4/50 ms working H and the remaining
small setup stages. Independent stage maxima need not sum to the API.
Eight-GPU rank-zero kernel time is80.9% NCCL send/receive plus15.1%
scalar allreduce; it includes peer waits and is not network bandwidth.
Raw outputs, all failures, binary hashes, stage provenance, CSV and replayable
analysis are in outer validation/rikyu-solver-scaling-20261111.
