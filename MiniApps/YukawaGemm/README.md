# YukawaGemm: a sparse GEMM mini application

A standalone driver for benchmarking NTPoly's distributed sparse
matrix-matrix multiply (`MatrixMultiply` from `PSMatrixAlgebraModule`) — the
kernel that dominates density matrix purification, and the first thing to
benchmark while porting NTPoly to GPU.

The workflow is intentionally split:

* **Python (via [uv](https://docs.astral.sh/uv/)) only generates data**:
  synthetic Hamiltonian-like matrices built from a screened Coulomb (Yukawa)
  potential on realistic atomic geometries.
* **A Fortran mini application does all the timing**, linked against the
  library exactly the way a solver would. It is part of the CMake build, so
  it compiles everywhere NTPoly does and will exercise the GPU backend the
  moment the backend is behind `MatrixMultiply`.

## Build

From the repository root (Fortran only, no Python bindings needed):

```sh
cd Build
cmake -DCMAKE_TOOLCHAIN_FILE=../Targets/mymachine.cmake \
      -DCMAKE_BUILD_TYPE=Release -DFORTRAN_ONLY=Yes ..
make -j
```

The executable lands in `Build/bin/yukawa_gemm`.

## Generate test matrices

`generate_matrices.py` has inline script metadata, so `uv run` resolves its
environment automatically (numpy + scipy):

```sh
cd MiniApps/YukawaGemm
uv run generate_matrices.py --list-presets
uv run generate_matrices.py --preset small          # writes yukawa_sc_32768.mtx (+ .json sidecar)
```

| preset    | geometry | atoms   | ~nnz/row | notes                    |
| --------- | -------- | ------- | -------- | ------------------------ |
| micro     | sc 8^3   | 512     | ~46      | smoke test               |
| tiny      | sc 16^3  | 4,096   | ~46      |                          |
| small     | sc 32^3  | 32,768  | ~46      | laptop-friendly          |
| small_md  | sc 32^3  | 32,768  | ~46      | jittered (amorphous)     |
| medium    | fcc 25^3 | 62,500  | ~77      |                          |
| large     | sc 96^3  | 884,736 | ~46      | cluster-sized            |
| xlarge    | fcc 63^3 | 1.0 M   | ~77      | cluster-sized            |

Geometry: simple-cubic (a = 2.7 A) or FCC (a = 3.6 A), one orbital per site,
Yukawa interaction `exp(-r/lambda)/r` with screening length `lambda = 1.0 A`
cut off at `Rcut = 6.0 A`, constant diagonal. Every parameter can be
overridden on the command line (`--geometry --side --lattice-const
--screening --rcut --jitter --seed --diag`).

For very large matrices, convert once to NTPoly's binary format (much faster
to read than MatrixMarket text):

```sh
mpirun -np 8 ../../../Build/bin/yukawa_gemm --input yukawa_sc_884736.mtx --convert
# produces yukawa_sc_884736.bin
```

## Run the benchmark

```sh
mpirun -np 4 Build/bin/yukawa_gemm \
    --input MiniApps/YukawaGemm/yukawa_sc_32768.mtx \
    --mode square --rounds 10 --warmup 3 --threshold 1e-4 \
    --out results/small_np4.json
```

| flag                | meaning                                                              |
| ------------------- | -------------------------------------------------------------------- |
| `--input FILE`      | `.mtx` (MatrixMarket) or `.bin` (NTPoly binary) input matrix.        |
| `--mode`            | `fixed` (default): steady-state `C <- A*A` — constant sparsity and bounded values, the mode to use for timing comparisons. `square`: `B <- B*B` each round — fill evolves like a purification chain, but values grow unboundedly (raw squaring has no polynomial rescaling), so runs saturate dense; use it to study fill dynamics, not for timing tables. |
| `--rounds N`        | Timed rounds after warmup (default 10).                              |
| `--warmup N`        | Untimed rounds (default 3) — page in memory, warm the pool.          |
| `--threshold T`     | Sparsity threshold passed to the multiply (default 1e-4).            |
| `--prows/--pcols/--pslices` | Explicit process grid (default: automatic square grid).      |
| `--out FILE`        | JSON summary path (default `yukawa_gemm_results.json`).              |
| `--direct`          | Bypass distributed wrapper (np=1, fixed mode only) — call local GEMM directly. |
| `--cusparse`        | Use cuSPARSE SpGEMM + CPU threshold (np=1, fixed mode only). Falls back to `--direct`/distributed on failure. |
| `--cusparse-mpi`    | Multi-GPU cuSPARSE with 1D row decomposition (fixed mode only). Each rank gets its own GPU via `OMPI_COMM_WORLD_LOCAL_RANK`. |
| `--cusparse-2d`     | Multi-GPU cuSPARSE with 2D SUMMA decomposition (fixed mode only). Requires square process grid (`--prows` = `--pcols`, `--pslices 1`). Set `OMP_NUM_THREADS=1`. |
| `--symmetric`       | Exploit A symmetry in 2D SUMMA: only compute upper triangle, transpose and send C^T to lower-triangle partner. Eliminates cuSPARSE SpGEMM pathology on slow off-diagonal blocks. |
| `--load-balance`    | Replace A with P·A·P^T (fixed-seed random permutation; NTPoly `LoadBalancerModule`). Blocks become uniform — use WITHOUT `--symmetric`. Do NOT set `UCX_TLS=sm,sysv,posix` on single-node runs with this flag. |

The multiply uses a memory pool that is constructed once and reused across
rounds, matching how the purification solvers call it.

## Benchmark validity: two things the benchmarks surfaced

1. **Fill is a function of the process grid.** At a fixed threshold, the
   number of nonzeros in the product grows with the MPI rank count
   (e.g. `tiny`, fixed mode, threshold 1e-4: 68 nnz/row at 1 rank, 122 at
   2, 172 at 4). Partial products are thresholded within the distributed
   algorithm, so a different process grid keeps a different set of
   near-threshold entries. Consequences:
   * per-rank-count results are perfectly deterministic (reruns are
     bit-identical), so **same-np comparisons (CPU vs GPU backend) are
     exact**;
   * **cross-np timings are only comparable at equal nnz/row** — always
     check `rounds_nnz` before comparing medians.
2. **Perfect crystals sit on the threshold.** In a perfectly ordered
   lattice, Yukawa products take discrete values and whole populations pile
   up at the cut, amplifying the effect above and breaking symmetry
   (`asymmetry` fingerprint becomes nonzero). The presets therefore default
   to a small position jitter (5% of the lattice constant); pass
   `--jitter 0` to generate a perfect crystal deliberately.

## Output

Per round: GEMM time (max over ranks), the untimed-but-reported rotate copy
(square mode), and the nonzero count. The JSON summary adds min/median/mean/
max GEMM times and a **fingerprint** (trace, grand sum, asymmetry, final
nnz). When comparing backends or processor counts:

* timings are comparable only at similar `nnz/row` (check `rounds_nnz`);
* the **fingerprint must match** — it is the correctness gate for the GPU
  port (asymmetry should stay ~0; a nonzero asymmetry means the backend is
  filling the wrong blocks).

## RIKYU full-node CPU baseline

The reference CPU numbers for the GPU port, run on one GB200 node (144 Grace
cores, 1.6 TB RAM) with `nvhpc-hpcx` (nvfortran + HPC-X MPI). Single-node jobs
need `UCX_TLS=sm,sysv,posix` and `UCX_NET_DEVICES=""` to disable InfiniBand
(mlx5 UAR exhaustion at high MPI rank counts) and `mpirun` not `srun` (srun
fails `MPI_Init` on RIKYU).

**Matrix:** sc 40³ = 64,000 atoms. This is the largest matrix that fits on a
single node at all process-grid configurations. The local sparse GEMM memory
pool allocates dense `columns × rows` arrays per rank (≈ 24 × (N/√np)² bytes),
which limits matrix size:

| matrix        | atoms    | fits? | bottleneck                          |
|---------------|----------|-------|-------------------------------------|
| sc 40³ (64k)  | 64,000   | yes   | all configs fit (≤ 2.7 GB/rank)     |
| sc 56³ (175k) | 175,616  | no    | OOM at np=144 (5.1 GB pool, 11 GB cgroup) |
| sc 96³ (885k) | 884,736  | no    | OOM at all single-node configs       |

Full results and analysis in `results/rikyu-fullnode/SUMMARY.md`.

| config        | median (ms) | min (ms) | nnz/row |
|---------------|-------------|----------|---------|
| np4_omp1      | 550         | 539      | 222.7   |
| np4_omp36     | 331         | 269      | 172.2   |
| np8_omp1      | 300         | 296      | 198.5   |
| np8_omp18     | 217         | 180      | 198.5   |
| np16_omp1     | 280         | 275      | 198.5   |
| np16_omp9     | 191         | 165      | 198.5   |
| np36_omp1     | 201         | 200      | 181.6   |
| **np36_omp4** | **128**     | **105**  | 181.6   |
| np72_omp1     | 153         | 148      | 152.1   |
| np144_omp1    | 200         | 169      | 152.0   |

**Best: np=36 omp=4 at 128 ms median (105 ms min).** The GPU port target is
to beat this on a single GB200 GPU while reproducing the trace fingerprint
(64287.048350011…).

## The memory bottleneck

The local sparse GEMM kernel (`SMatrixAlgebraModule`: `GemmMatrix_lsr` →
`SparseBranch_lsr` → `MultiplyBlock_lsr`) uses a memory pool
(`ConstructMatrixMemoryPool.f90`) that allocates **five dense arrays** sized
`columns × rows` per rank:

```fortran
ALLOCATE(this%pruned_list(columns * rows))      ! INTEGER
ALLOCATE(this%value_array(columns, rows))       ! REAL(NTREAL)
ALLOCATE(this%dirty_array(columns, rows))       ! LOGICAL
ALLOCATE(this%hash_index(columns, rows))         ! INTEGER
ALLOCATE(this%inserted_per_bucket(columns, rows)) ! INTEGER
```

At ~24 bytes/element this dominates memory and is the bottleneck for both
large matrices and the GPU port. Replacing this dense accumulator with a
sparse accumulator or GPU-friendly hash table is the first porting task.

## GPU benchmark: cuSPARSE SpGEMM

A standalone CUDA benchmark (`bench_cusparse.cu`) measures what a single
GB200 GPU can do on the same Yukawa matrices using cuSPARSE SpGEMM (no
NTPoly involvement — just raw `cusparseSpGEMM` on CSR A, computing C=A*A).

**Build on RIKYU:**
```sh
module load nvhpc-hpcx
CUDA_HOME=/shared/software/hpc_sdk/Linux_aarch64/26.3/cuda
nvcc -O2 -std=c++17 -I$CUDA_HOME/include -o bench_cusparse bench_cusparse.cu \
    -lcusparse -L/shared/software/hpc_sdk/Linux_aarch64/26.3/math_libs/lib64
```

**Results (single GB200 GPU):**

| matrix        | atoms   | GPU median  | CPU best   | speedup | CPU can run? |
|---------------|---------|-------------|------------|---------|--------------|
| sc 40³ (64k)  | 64,000  | **13.9 ms** | 128 ms     | **9.2×** | yes          |
| sc 56³ (175k) | 175,616 | **36.6 ms** | OOM        | —       | no           |
| sc 96³ (885k) | 884,736 | **185 ms**  | OOM        | —       | no           |

The GPU eliminates the memory bottleneck entirely — all three matrices run
with no issues. The 885k matrix (which would need 12+ CPU nodes) runs in 185
ms on one GPU. Trace matches the CPU baseline (64287.048350011 at 64k).

Full comparison in `results/rikyu-fullnode/GPU_VS_CPU.md`.

**Porting conclusion:** Replace NTPoly's local GEMM kernel (`MultiplyBlock` +
`PruneList` + the dense accumulator) with cuSPARSE SpGEMM + a post-hoc
thresholding kernel. This is both faster and removes the memory ceiling.

## Fill sweep: does the GPU advantage hold at realistic DFT density?

The benchmark above uses ~0.06% fill — far sparser than real DFT density
matrices (~1-5%). A fill sweep on sc 32³ = 32,768 atoms (varying rcut to
sweep fill from 0.12% to 4%) shows the GPU advantage collapses:

| Input fill | Input nnz/row | CPU median | GPU (cuSPARSE) | Speedup |
|-----------:|--------------:|-----------:|---------------:|--------:|
| 0.12%      | 40            | 64.6 ms    | 7.0 ms         | **9.2×** |
| 0.58%      | 189           | 234 ms     | 118.6 ms       | **2.0×** |
| 1.83%      | 599           | 553 ms     | **CRASH**      | —       |
| 4.0%       | 1314          | 1747 ms    | **CRASH**      | —       |

cuSPARSE crashes at 1.8%+ fill because it computes the **full product**
(no thresholding during multiply), which becomes nearly dense at realistic
fill. All four cuSPARSE algorithms (DEFAULT, ALG1, ALG2, ALG3) return
`CUSPARSE_STATUS_INSUFFICIENT_RESOURCES`.

NTPoly's CPU code keeps working because it thresholds **during** the multiply
— the product stays at ~336 nnz/row regardless of input fill. This
thresholding is the key algorithmic feature that cuSPARSE lacks.

**Revised porting conclusion:** cuSPARSE SpGEMM is NOT viable at realistic
DFT fill. The GPU port must replicate NTPoly's thresholding-during-multiply
algorithm (the hash accumulator in `MultiplyBlock` + `PruneList`), likely via
OpenACC directives on the existing Fortran kernel.

Full analysis in `results/rikyu-fullnode/FILL_SWEEP.md`.

## GPU port: OpenACC

The existing Fortran GEMM kernel (`MultiplyBlock_lsr` + `PruneList_lsr`) was
ported to GPU with OpenACC directives (`!$acc parallel loop independent`).
Built with nvfortran `-acc -gpu=mem:managed -gpu=cc100`. A `--direct` flag
bypasses the distributed wrapper at np=1 (eliminating per-round
gather/compose/transpose overhead).

Key fixes along the way:
- **Bubble sort → insertion sort** in `SortTripletList.f90` (was 97% of CPU
  overhead — 5020 ms on 64k).
- **INT32 overflow** in sort dispatch condition (caused 48 GB dense buffer
  allocation at 64k).
- **Direct path** calls `GemmMatrix_lsr` directly at np=1.

**Results (64k, np=1, `--direct`):**

| Backend              | Median   | vs CPU   |
|----------------------|----------|----------|
| CPU np=36 omp=4      | 128 ms   | 1.0×     |
| OpenACC np=1 direct  | 536 ms   | 0.24×    |

GPU compute is ~371 ms (MultiplyBlock 336 ms + PruneList 35 ms); the remaining
~165 ms is CPU↔GPU managed-memory page migration and CSR construction. The GPU
is **4.2× slower** than the CPU — the random-access hash accumulator pattern
is not GPU-friendly.

Full analysis in `FINDINGS.md` (Finding 7).

## GPU port: cuSPARSE hybrid (current best)

A hybrid approach: compute the **full product** with cuSPARSE SpGEMM on GPU,
then **threshold on GPU** (mark + count + compact kernels). NTPoly's CSC
arrays are passed directly as CSR to cuSPARSE (computes C^T in CSR = C in CSC
for symmetric C = A*A — no format conversion). Implemented as a C wrapper
(`cusparse_gemm.cu`) callable from Fortran via ISO_C_BINDING, compiled with
nvcc at CMake configure time.

The `--cusparse` flag enables the cuSPARSE backend. Combined with `--direct`,
it automatically falls back to the OpenACC path when cuSPARSE fails at high
fill.

### Optimization progression

| Version | Key change | 64k median | vs CPU |
|---------|-----------|-----------:|------:|
| v1 | Per-round cudaMalloc/cudaFree everywhere | 61 ms | 2.1× |
| v3 | Persistent A on GPU, skip per-round C→Fortran copy | 18.2 ms | 7.0× |
| **v4** | **Persistent C/work buffers + warm-call (skip size queries)** | **7.58 ms** | **16.9×** |

### v4 key optimizations

- **Persistent buffers** — A, C (full product), work buffers, threshold
  pipeline buffers, and host output buffers are all allocated once and reused.
  Eliminates ~10 ms of GPU idle time caused by per-round `cudaMalloc` of 180 MB.
- **Warm-call optimization** — after the first SpGEMM, skip the cuSPARSE
  size-query calls and reuse cached work buffer sizes (saves ~1.5 ms/round).
- **GPU-side thresholding** — mark (above threshold) + count (per row) + compact,
  all on GPU. Only the n-element prefix scan runs on CPU.
- **Algorithm selection** — `CUSPARSE_ALG` env var (0=DEFAULT, 1=ALG1,
  2=ALG2, 3=ALG3, -2=benchmark all). ALG1 is fastest; ALG2/ALG3 silently fail.

### Results (job 281613)

**64k matrix (0.06% fill):**

| Backend                      | Median   | Min    | vs CPU   | Trace              | Asymmetry |
|-----------------------------|----------|--------|----------|--------------------|-----------|
| CPU np=36 omp=4 (baseline)  | 128 ms   | 105 ms | 1.0×     | 64287.048350011    | —         |
| OpenACC np=1 direct          | 536 ms   | 536 ms | 0.24×    | 64287.048350011 ✓  | 0         |
| cuSPARSE hybrid v1           | 61 ms    | 61 ms  | 2.1×     | 64287.048350011 ✓  | 1.6e-16   |
| **cuSPARSE hybrid v4**       | **7.58 ms** | **7.56 ms** | **16.9×** | **64287.048350011 ✓** | **1.6e-16** |

**32k fill sweep:**

| rcut | Fill  | Backend      | GPU median | CPU (np=36) | Speedup | nnz       | Trace                |
|------|-------|-------------|-----------:|------------:|--------:|-----------|----------------------|
| 6    | 0.12% | cuSPARSE    | 4.05 ms    | 64.6 ms     | **16.0×** | 2,444,530 | 32913.930651629  ✓ |
| 10   | 0.58% | cuSPARSE    | 50.4 ms    | 234 ms      | **4.6×**  | 3,747,210 | 32914.015272987  ✓ |
| 15   | 1.83% | OpenACC     | 850 ms     | 553 ms      | 0.65×  | 3,751,584 | 32914.015291414  ✓ |
| 20   | 4.0%  | OpenACC     | 3203 ms    | 1747 ms     | 0.55×  | 3,751,584 | 32914.015291415  ✓ |

cuSPARSE achieves **4.6–16.9× speedup at low fill (≤0.58%)**. At high fill
(≥1.83%), cuSPARSE fails (silent C_nnz=0 for ALG2/ALG3, INSUFFICIENT_RESOURCES
for ALG1) and falls back to OpenACC, which is slower than CPU.

### Per-phase GPU timing (64k steady-state)

| Phase                        | GPU time | % of total |
|------------------------------|---------:|-----------:|
| cuSPARSE compute (warm)      | 6.3 ms   | 84%        |
| cuSPARSE copy                | 0.11 ms  | 1.5%       |
| GPU threshold mark + count   | 0.15 ms  | 2%         |
| GPU compact                  | 0.50 ms  | 6.7%       |
| D2H download (thresholded)   | 0.32 ms  | 4.3%       |
| Other (CPU scan + overhead)  | 0.20 ms  | 2.5%       |
| **Total**                    | **7.58 ms** | 100%    |

The cuSPARSE compute kernel (6.3 ms) is 84% of the total — it is the
theoretical floor for this approach.

### Bandwidth analysis

SpGEMM is **not a streaming operation** — it cannot hit peak HBM bandwidth:

| Metric                          | Value       |
|--------------------------------|------------:|
| GB200 HBM bandwidth (peak)     | 8 TB/s      |
| Estimated memory traffic        | ~241 MB     |
| Streaming BW floor              | 0.030 ms    |
| **Actual cuSPARSE compute**     | **6.3 ms**  |
| Achieved bandwidth               | 38 GB/s     |
| **% of peak HBM**               | **0.48%**   |

SpGEMM is latency-bound by hash table operations and random memory access,
not bandwidth-bound. Dense GEMM (cuBLAS) can hit >90% of peak, but at 0.06%
fill it would do 99.94% wasted work.

Full analysis in `FINDINGS.md` (Finding 9).

## Multi-GPU scaling: cuSPARSE with 1D row decomposition

The `--cusparse-mpi` mode scales the cuSPARSE hybrid across multiple GPUs
using a 1D row decomposition: each MPI rank reads the full matrix, extracts
its row-block, and calls cuSPARSE on its own GPU. No communication during
the multiply — only `MPI_Allreduce` for the fingerprint at the end.

### Results (64k matrix, 0.06% fill)

| np | GPUs | Nodes | median (ms) | vs CPU (128ms) | efficiency |
|----|------|-------|------------:|---------------:|-----------:|
| 1  | 1    | 1     | 7.76        | 16.5×          | 100%       |
| 2  | 2    | 1     | 4.33        | 29.6×          | 89.7%      |
| 4  | 4    | 1     | 2.35        | 54.5×          | 82.6%      |
| 8  | 8    | 2     | 1.39        | **92.1×**      | 69.8%      |

All traces match the CPU baseline (64287.048350011). The scaling is
compute-bound (communication is <1% of total time). Efficiency drops from
90% (np=2) to 70% (np=8) due to cuSPARSE per-call overhead and the fixed
threshold pipeline cost, which don't shrink with fewer rows per rank.

## Multi-GPU scaling: 2D SUMMA + symmetry exploitation

The 1D row decomposition above is fast but has two problems for production:
(1) it replicates the full matrix A on every rank (OOM risk at scale), and
(2) the output is in 1D row layout, requiring redistribution each round
during purification. The 2D SUMMA approach solves both: A is stored in 2D
block layout (no full replication), and the output stays in 2D block layout
(no redistribution for purification).

### Design

Process grid: pr × pc × 1 (square: pr = pc). Each rank stores block
(my_row, my_col) of A in CSC (NTPoly's native format). Per round:
1. AllGather row-block along `row_comm` → A_rows (m_loc × N)
2. Bcast B from rank my_col along `col_comm`
3. Transpose A_rows CSC → CSR → A_local (m_loc × N)
4. B = A[my_col,:] in CSC, aliased as CSR of N×m_loc (symmetric A trick)
5. cuSPARSE SpGEMM: C = A_local × B, then threshold

### cuSPARSE SpGEMM pathology

At np=4, two off-diagonal blocks with identical flop_count (1.49M) showed a
**120× throughput gap** (0.49 ms vs 58 ms). All cuSPARSE algorithms give the
same result — this is a structural cuSPARSE issue, not a workload issue.

### Symmetry exploitation (`--symmetric`)

For symmetric A, C = A² is symmetric: C[my_row,my_col] = C[my_col,my_row]^T.
Upper-triangle ranks compute SpGEMM, transpose, and send C^T to their
lower-triangle partner. Lower-triangle ranks skip SpGEMM and receive C^T.
This eliminates the slow off-diagonal block entirely.

**Result: np=4 from 58 ms → 4.0 ms (14.4× improvement).**

### Results — 64k matrix (0.06% fill), no-LB + symmetry

| np | Grid | GPUs | Nodes | Median (ms) | vs CPU (128ms) | vs np=1 |
|----|------|------|-------|------------:|---------------:|--------:|
| 1  | 1×1  | 1    | 1     | 7.57        | 16.9×          | 1.0×    |
| 4  | 2×2  | 4    | 1     | 4.13        | 31.0×          | 1.83×   |
| 9  | 3×3  | 12   | 3     | 4.18        | 30.6×          | 1.81×   |
| 16 | 4×4  | 16   | 4     | 2.73        | 46.9×          | 2.77×   |

Note: np=9 is not faster than np=4 — the bottleneck (diagonal SpGEMM) has
highly variable cuSPARSE throughput (8.7M–15.3M flops/ms) that doesn't
correlate with problem size, and one diagonal "pathological block" per grid
runs ~1.4× slower than its siblings with identical flops (the MAX-over-ranks
wall always pays the slowest).

### Results — 64k matrix, load-balanced (`--load-balance`, no `--symmetric`)

NTPoly's native `LoadBalancerModule` replaces A with P·A·Pᵀ (fixed-seed
random permutation; trace/grand-sum/nnz of A² invariant). Blocks become
statistically uniform, which (a) eliminates the pathological-block lottery
and (b) makes the symmetry exchange unnecessary — every rank computes its
own C block with **zero per-round communication**:

| np | Grid | GPUs | Nodes | Median (ms) | vs CPU (128ms) | vs np=1 (7.57) |
|----|------|------|-------|------------:|---------------:|--------:|
| 1  | 1×1  | 1    | 1     | 7.92        | 16.2×          | 0.96×   |
| 4  | 2×2  | 4    | 1     | 3.08        | 41.6×          | 2.46×   |
| 9  | 3×3  | 12   | 3     | 1.78        | 71.9×          | 4.25×   |
| 16 | 4×4  | 16   | 4     | **1.27**    | **100.8×**     | **5.96×** |

All traces/grand-sums match. LB + `--symmetric` is a net LOSS (exchange
messages become large once blocks are uniform — 5.36 ms at np=16); LB
WITHOUT `--symmetric` is the production config. Strong-scaling plot:
`strong_scaling.png` (regenerate with `uv run plot_strong_scaling.py`).

> **RIKYU transport caveat:** single-node runs that touch NTPoly's
> distributed path (`PermuteMatrix`) must NOT set
> `UCX_TLS=sm,sysv,posix` — the UCX shared-memory transport deterministically
> corrupts its non-blocking collectives at np≥2 on one node. Our blocking-
> collective cuSPARSE 2D SUMMA path is unaffected. For LB runs on a single
> node, omit `UCX_TLS` (IB loopback works at np=4).

### Results — 175k matrix (OOM on CPU)

The 175k matrix (175,616 × 175,616) OOM'd on CPU at all single-node configs.
With 2D SUMMA on 16 GPUs across 4 nodes:

| np | Grid | GPUs | Nodes | Median (ms) | vs np=1 | Trace |
|----|------|------|-------|------------:|--------:|------:|
| 1  | 1×1  | 1    | 1     | 20.60       | 1.0×    | 176409.895 ✓ |
| 4  | 2×2  | 4    | 1     | 10.65       | 1.93×   | 176409.895 ✓ |
| 16 | 4×4  | 16   | 4     | 6.18 (no-LB sym) | 3.33× | 176409.895 ✓ |
| 16 | 4×4  | 16   | 4     | **2.85 (LB no-sym)** | **7.69×** | 176409.895 ✓ |

### Results — 885k matrix (OOM on CPU)

The 885k matrix (884,736 × 884,736, 37M nnz) — the largest test matrix:

| np | Grid | GPUs | Nodes | Median (ms) | vs np=1 | Trace |
|----|------|------|-------|------------:|--------:|------:|
| 1  | 1×1  | 1    | 1     | 104.88      | 1.0×    | 888767.956 ✓ |
| 4  | 2×2  | 4    | 1     | 53.70       | 1.95×   | 888767.956 ✓ |
| 16 | 4×4  | 16   | 4     | 28.00 (no-LB sym) | 3.75× | 888767.956 ✓ |
| 16 | 4×4  | 16   | 4     | **10.46 (LB no-sym, pre-permuted)** | **10.03×** | 888767.956 ✓ |

> **885k load-balancing caveat:** runtime `--load-balance` OOMs at 885k —
> `PermuteMatrix` builds the permutation through NTPoly's distributed
> multiply, whose dense local accumulator (Finding 1) needs ~1.1 TB for the
> 221184² local blocks at the (4,4) grid. Workaround: pre-permute the
> matrix file once with `uv run permute_binary.py
> matrices/yukawa_sc_884736.bin matrices/yukawa_sc_884736_perm.bin`
> (fixed-seed random symmetric permutation; trace/grand-sum invariant) and
> run WITHOUT `--load-balance`. 175k works with runtime LB; at 885k the
> pre-permuted file gives **10.46 ms** vs 28.0 ms no-LB symmetric (2.68×).

Multi-GPU scaling improves with matrix size (np=16: 2.77× for 64k, 3.33×
for 175k, 3.75× for 885k) — larger matrices have more work per GPU.

### Comparison: 1D vs 2D SUMMA

| Approach | np=4 | np=16 | OOM risk | Purification-ready |
|----------|-----:|------:|----------|-------------------|
| 1D row | 2.35 ms | — | High (full A replicated) | No (redistribution each round) |
| 2D SUMMA + sym | 4.13 ms | 2.73 ms | Low (block layout) | **Yes** |

The 1D approach is faster at np=4 (no communication, simpler cuSPARSE calls),
but the 2D SUMMA + symmetry approach is the production choice: it handles
large matrices without OOM and integrates with purification without
redistribution.

### 2.5D k-split (`--ksplit`) — tested and rejected

The 2.5D variant splits the shared dimension across process slices
(grid pr×pc×ps; slice s computes the partial C^(s) over its k-slab, partials
summed across slices and thresholded). Implemented as `--ksplit` (with
`--cusparse-2d --load-balance`; partials merged on GPU via a dense
scatter-add, see Finding 13). Numerically correct (trace/grand-sum match to
13 significant figures at (2,2,2) and (2,2,4)) but far slower:

| Config | Round wall | vs equivalent 2D |
|--------|-----------:|------------------|
| np=8 (2,2,2) k-split | 84.6 ms | 47× slower than 1.78 ms |
| np=16 (2,2,4) k-split | 3040 ms | 2400× slower than 1.27 ms |

At 0.06% fill the unthresholded partial products are ~40–90× denser than A,
so each round moves and merges ~100× the input data to save at most ~2× on
latency-bound SpGEMM compute. Closed for this regime — 2D + load balancing
is the right algorithmic point. (Side note: per-block C density differs by
block type even under LB — diagonal blocks ~1.5× denser than off-diagonal;
see Finding 13's correction.)

### RIKYU multi-node job submission

Keep the default UCX transports for NTPoly solver and load-balancing paths.
Do not force `UCX_TLS=sm,sysv,posix`: later checks found collective corruption
with that setting (Finding 12). This also applies to `--nccl` runs, since
the embedded NTPoly reference still uses MPI.

Multi-node (8+ GPUs): do NOT set these — use IB defaults. RIKYU allows
GPU counts [1, 2, 3, 4, 8, 12, 16] per job (4 GPUs/node, max 4 nodes).

Full analysis in `FINDINGS.md` (Finding 11).

## Current state of the GPU port

The table below contains historical fixed-GEMM measurements. Current solver
correctness, CPU comparisons and NCCL results are in Findings 18–22.

| Approach              | 64k vs CPU | Multi-GPU? | High fill? | Status |
|-----------------------|------------|------------|------------|--------|
| CPU baseline          | 1.0×       | —          | yes        | reference |
| OpenACC               | 0.24×      | no         | yes        | working |
| cuSPARSE hybrid v4    | 16.9×      | no         | fallback   | working |
| cuSPARSE 1D (8 GPU)   | 92.1×      | yes        | fallback   | working |
| cuSPARSE 2D+sym (16 GPU) | 46.9×   | yes (4 nodes) | fallback   | working |
| cuSPARSE 2.5D k-split | 0.04×      | yes (4 nodes) | fallback   | correct, rejected |
| **cuSPARSE 2D + LB no-sym (16 GPU)** | **100.8×** | **yes (4 nodes)** | fallback | **working** |

The 2D SUMMA + load balancing (no symmetric exchange) is the production
choice: **100.8× over CPU** on 64k at np=16, **2.85 ms at 175k** and
**10.46 ms at 885k** (both OOM on CPU; 885k uses the pre-permuted matrix
file, see the caveat above), and integrates with purification without
per-round redistribution. All results are bit-identical to the CPU baseline
(trace = 64287.048350011 at 64k). The 2.5D k-split variant (Finding 13) is
numerically correct but 47–2400× slower at this fill — rejected.

## NCCL solver communication (single node, 1–4 GPUs)

Build with NCCL headers/library available (CMake checks `NCCL_HOME` and the
RIKYU SDK installation). Select `--nccl` to use NCCL device-buffer gathers
and scalar reductions in the TRS2 loop. MPI remains responsible for
communicator bootstrap and profiling reductions. The NTPoly reference is
unchanged. NCCL absence is an explicit error; there is no silent fallback.

```sh
OMP_NUM_THREADS=1 mpirun -np 4 ./yukawa_gemm \
  --input yukawa_sc_64000_gap.bin --mode fixed --solver trs2 \
  --nelec 32000 --threshold 1e-4 --conv-diff 1e-6 --max-iter 60 \
  --load-balance --prows 2 --pcols 2 --pslices 1 --nccl
```

Use grids 1×1, 1×2 and 2×2 for 1, 2 and 4 GPUs. Rectangular layouts now
merge/split NTPoly's local sub-blocks correctly and use true CSR/CSC
conversion. The NCCL path rejects multi-node runs or more than four ranks.
Omit `--nccl` to compare the same solver with MPI communication.

The default path builds CSR panels on the CPU, with bulk
host/device staging around each NCCL exchange. Variable-count payloads use
grouped NCCL send/receive; fixed-count metadata uses AllGather. First-use
communicator costs are included in iteration 1, so report full solves and
later-iteration timings separately. See Finding 20 for measured results.

### GPU accumulation experiments (Finding 21)

The solver also accepts these optional paths:

- `--device-partials`: retain cuSPARSE's raw chunk products on the GPU and
  accumulate immediately. Only the final thresholded product is downloaded.
  Thresholding still occurs after all chunks have contributed.
- `--gpu-dense`: multiply sparse panels directly into the GPU dense
  accumulator, bypassing chunk packing, cuSPARSE symbolic work, and raw
  partial arrays. This is an experimental alternative to chunked cuSPARSE.
- `--nccl-warmup`: exercise collective and P2P paths before the loop.
  `nccl_setup_seconds` includes communicator creation and warmup;
  `nccl_loop_plus_setup_seconds` reports setup plus loop. Spectral bounds,
  X0 construction and upload remain outside these numbers.
- `YUKAWA_NCCL_SHARED_COMM=1`: use the world NCCL communicator for mapped
  row/column P2P exchanges, reducing redundant communicator setup. MPI
  bootstraps the peer maps. Default behavior retains separate communicators.

For the validated four-GPU direct-accumulator experiment, append
`--gpu-dense --nccl-warmup` to the command above. Set
`YUKAWA_NCCL_SHARED_COMM=1` on every rank to test shared communicator setup.
For chunked cuSPARSE, use `--device-partials`; `--chunk-work 64000000`
was faster on the gapped 64k case, while the conservative default remains
4000000. Changing this budget changes workspace and intermediate fill.

Both accumulation paths require a dense local `m × k` buffer, roughly
8.192 GB/GPU for the 64k, 2×2 grid. They are not general sparse-memory
algorithms. Without `--resident-panels`, panel construction/conversion
remains on the CPU; the accumulation flags alone do not imply residency. The direct path passed GPU memory
checking, 1/2/4-GPU tiny cases and full 24-iteration matrix parity.

The profile CSV now has a companion `.ranks.csv`, gathered and written
after timing, retaining coherent phase vectors for each rank. Independent
phase maxima must not be subtracted from wall time to predict speedups.

When profiling only rank 0 with Nsight's MPI tracing, set
`NSYS_MPI_STORE_TEAMS_PER_RANK=1` as required by NVIDIA's subset-rank
profiling instructions. Preserve successful unprofiled measurements and
keep profiling failures separate from application failures.

### Resident panels and bitmap accumulation (Finding 22)

Add `--resident-panels` to `--gpu-dense --nccl` to keep CSR payloads,
row/column panel assembly, multiplication and the product slot on GPU.
Only counts/scalars are staged for communication. The CPU CSC matrix is
materialized after each solve for validation/output. Tested grids are
1×1, 1×2, 2×1 and 2×2, with nonsymmetric and empty-block CUDA tests.

`YUKAWA_BITMAP_ACCUM=1` marks contributed entries with a touched bitmap.
Count/compaction skip untouched values; later iterations clear previously
touched positions, including cancellations. This adds one bit per dense
element (128 MB/GPU at 64k/2×2). The numerical dense allocation remains;
this option saves work, not O(N²/GPUs) storage. It is off by default and
affects the resident direct multiply.

For a dense diagonalization comparison, the standalone MPI/CUDA
[`cusolvermp_bench.cu`](cusolvermp_bench.cu) forms the occupied density using
cuSOLVERMp and cuBLASMp. See [build instructions and timing scopes](CUSOLVERMP.md).

`--solver-reps N` restores cached GPU X0 before N independently converged
solves, reusing communicators and allocations. Every result receives full
reference parity. This repeats the same H; changed H needs new bounds/X0.
Profile paths gain `.repN` suffixes when N>1.

```sh
OMP_NUM_THREADS=1 YUKAWA_NCCL_SHARED_COMM=1 YUKAWA_BITMAP_ACCUM=1 \
mpirun -np 4 ./yukawa_gemm \
  --input yukawa_sc_64000_gap.bin --mode fixed --solver trs2 \
  --nelec 32000 --threshold 1e-4 --conv-diff 1e-6 --max-iter 60 \
  --load-balance --prows 2 --pcols 2 --pslices 1 \
  --nccl --nccl-warmup --gpu-dense --resident-panels --solver-reps 5
```

`gpu_bind_seconds` measures explicit runtime binding;
`gpu_numerical_setup_seconds` includes communicator/warmup, spectral
bounds, X0 construction/uploads and the optional snapshot.
`gpu_chain_seconds` sums loops and state resets;
`gpu_chain_with_setup_seconds` adds measured binding/numerical setup.
MPI startup, file I/O/permutation, input validation, embedded reference,
parity checks and CPU output materialization are excluded. For this 64k
case one active GPU wins a five-solve numerical batch; four GPUs gives
the fastest later solve. These are not total executable wall times.

`resident_panels_test` accepts an optional matrix size (e.g. `66` checks
32-column bitmap tails). It tests nonsymmetric products, cancellation,
empty peers, reset and rectangular layouts. OpenMPI/UCX's handled CUDA
pointer/context probes generate API-error reports under Compute Sanitizer.
Use `--tool memcheck --report-api-errors no` for this MPI test: memory
instrumentation remains active and backend CUDA statuses are checked.
Preserve original probe-error logs; this mode does not test every API return.

Multi-node NCCL is an experimental opt-in: set
`YUKAWA_NCCL_ALLOW_MULTINODE=1` in addition to the environment above.
Without it the solver retains the single-node/four-rank restriction.
Use `mpirun` and explicit placement, e.g. `-np 8 --map-by ppr:4:node`
with a 2×4 or 4×2 process grid inside a two-node allocation. GPU binding
uses the node-local MPI rank. Capture NCCL transport logs when comparing
nodes; a successful run alone does not establish that GPU Direct RDMA was
selected. The distributed tests accept eight ranks; use matrix size 68
for `resident_panels_test` so both grid dimensions divide the matrix.

### Normal library TRS2 and small-matrix coalescing (Finding 24)

The optional library backend runs the existing `DensityMatrixSolversModule::TRS2`
iteration with owned GPU matrices. Calling code and the `Matrix_ps` layout are
unchanged. This is a first library integration: other solvers and ordinary
matrix operations outside TRS2 continue to use their existing paths. Build
with CMake >=3.20 and CUDA/cuSPARSE/NCCL available:

```sh
cmake -S . -B build-native -DCMAKE_TOOLCHAIN_FILE=Targets/rikyu-cpu.cmake \
  -DCMAKE_BUILD_TYPE=Release -DFORTRAN_ONLY=ON -DNOSWIG=ON \
  -DNTPOLY_ENABLE_GPU_TRS2=ON -DYUKAWA_CUDA_PROFILE=OFF
cmake --build build-native --target trs2_reference native_gpu_test -j 16
NTPOLY_GPU_TRS2=1 YUKAWA_BITMAP_ACCUM=1 YUKAWA_NCCL_SHARED_COMM=1 \
  mpirun -np 4 --map-by ppr:4:node ./build-native/bin/trs2_reference \
  matrix.bin 2 2 5000 1e-4 1e-6 60 3
```

`NTPOLY_ENABLE_GPU_TRS2` defaults OFF. In an enabled build,
`NTPOLY_GPU_TRS2=1` selects the GPU iteration; otherwise TRS2 uses CPU.
Real nonidentity working-H and final basis transforms now use resident A×B
when admission succeeds, with host imports/exports at each transform boundary.
Bounds are computed from resident WH with O(local columns) host MPI reductions.
X0 is copied, shifted and scaled on GPU; only the sparse identity is imported.
Load balancing and the public matrix interface still use CPU paths.
GPU results are exported to the
normal distributed host matrix before TRS2 returns. The GPU matrix type owns
independent dynamic handles with deep-copy assignment and finalization;
legacy mini-app slots are reserved separately. Multiply/merge scratch is
still a process-wide cache, and one active GPU solver context/process is
supported. Sequential calls now borrow a cached communication context owned
by their process grid: full-grid and coalesced-root contexts are separate,
and each is initialized/warmed only on its first use. Matrix/panel payloads
are still released at the solve boundary. Changed H, dimensions, overlap and
CPU calls between GPU solves do not invalidate communication-only state.
`NTPoly GPU context reused: T/F` records reuse per GPU call.

`DestructProcessGrid` releases its cached NCCL communicators, streams, events
and staging buffers before freeing MPI handles; call it before MPI_Finalize.
`CopyProcessGrid` creates an independent grid/context lifetime. As with its
MPI handles, intrinsic assignment of owning grids is unsupported. Several
grids can retain independent idle contexts; nested/concurrent GPU calls are
unsupported. Communication settings are read when each context is created.
Cached library-internal NCCL storage remains until grid destruction; payload
memory admission still uses live free HBM. This does not cache numerical
bounds/X0, matrices or the quadratic product workspace per grid.

`NTPOLY_GPU_RANK_POLICY` accepts:

- `auto` (default): rank one-GPU and full-grid candidates by estimated work
  and communication, then admit only layouts that pass memory/index checks
  (Finding 25). A one-rank application avoids identity redistribution.
- `all`: request all application ranks; unsafe layouts fall back to CPU.
- `root`: request one active GPU; unsafe layouts fall back to CPU.

Coalescing gathers WH and the sparse identity, leaves other ranks outside the iteration,
and redistributes the result once. It broadcasts the iteration count, sigma
history and energy for consistent output/chemical potentials on every rank.
The application's input/output matrix distribution is preserved. The original
10k cutoff in Finding 24 has been replaced by the model in Finding 25.
It is tested at 1/4/8/16 ranks, including two independent application
communicators. The all-rank path requires one rank per GPU; more than four
ranks/multiple nodes also requires `YUKAWA_NCCL_ALLOW_MULTINODE=1`.
The current backend supports real matrices and one process-grid slice.

Finding 24 recorded, before the model change: at N=10000, two fresh processes ×3 solves/configuration give later normal
API times 0.138/0.083/0.075/0.078 s at 1/4/8/16 application ranks with
`auto`, all using one active GPU. Forced-all at 16 ranks is 2.131 s;
CPU (one OpenMP thread/rank) is 0.355 s. These include internal CPU setup,
CUDA/NCCL setup, import, iteration, export, teardown and final transform;
exclude MPI launch, input read/permutation and validation/output file writes.
All twelve saved matrices pass full comparison (<6.4e-16 relative norm).

`ntpoly_gpu_trs2_{setup,loop,export}_seconds` are per-call maximum-rank
phase counters. Setup now covers GPU runtime/communication/import, bounds and
X0 construction; it excludes preceding working-H construction and balancing.
Coalesced loop counters can include idle-rank arrival skew. Use the enclosing
API timing for latency comparisons; do not sum independent phase maxima.
Contexts are reused per process grid (Finding 27).

`NTPOLY_TRS2_PROFILE=1` adds maximum-rank stage counters in CPU and GPU
builds: allocate, working_h, balance, backend, bounds, initial_x, iteration,
export_cleanup and final. Allocation includes matrix/identity construction;
export_cleanup includes GPU matrix destruction. Final includes undo balancing,
basis transform, chemical potential and cleanup. Independent maxima need not
sum to the API maximum, and automatic host finalization can occur afterward.

The reference driver's optional arguments after repetitions are an output
filename (or `-`) and `identity|banded`. Banded constructs an SPD tridiagonal
inverse-overlap factor with diagonal1 and nearest neighbors0.01 outside
the API, in the already permuted input ordering. It tests nonidentity
transforms; the physical density's ordinary trace need not equal electron
count in a nonorthogonal basis. Compare full CPU/GPU matrices and energy.

The general resident product accepts independent real square A/B matrices
with matching distributed layouts; local blocks can be rectangular.
Aliasing the output with either input is rejected. The square ABI remains.
Identity construction now visits only the owned diagonal intersection,
including padded/empty blocks, for real and complex matrices.

### Reusing accessible input storage

The caller can construct or update H on CPU or GPU. Native imports now borrow
existing CSC arrays when CUDA can access them: managed allocations, allocations
on the current GPU, mapped host memory, and ordinary system allocations on
devices reporting `cudaDevAttrPageableMemoryAccess`. The source is read-only,
the conversion completes before return, and no borrowed pointer is retained.
The backend still needs an independently owned CSR representation; cuSPARSE
performs that required conversion on GPU without CPU payload staging.

For one local block, import skips CPU block composition entirely. Multiple
local blocks retain the existing CPU merge; coalescing retains the host MPI
gather. Unsupported memory falls back to the host conversion. Numerical data
is not cached by address, so changing values or patterns cannot reuse stale
payloads. CPU/GPU synchronization is explicit at the import boundary.

With identity overlap and no balancing, GPU TRS2 also borrows H directly as
the working-H input. It materializes a separate host WH only if CPU fallback
is selected. Balancing retains its independent working matrix. Public host
output and caller placement choices are preserved.

`NTPOLY_GPU_IMPORT=auto` (default) enables this path; `host` selects CPU
layout conversion as a control. `NTPOLY_GPU_IMPORT_REPORT=1` prints borrowed
CSC imports for diagnostics; leave it off for latency measurement.
`NTPOLY_TRS2_PROFILE=1` also reports whether identity H was borrowed.
The `import_csc_test` target exercises managed/device/system/mapped arrays,
device-only value changes on a separate stream, sorted output, changing
dimensions and empty matrices. On ATS systems, Compute Sanitizer may need
`NV_COMPUTE_SANITIZER_SHARED_ADDRESSING_SUPPORT=force` for system allocations.


### Work estimates and GPU memory admission (Finding 25)

`auto` now uses the working matrix's actual dimension, padded grid layout,
WH/X0 nonzero counts, estimated sparse-product work and native communication
startup cost. It ranks the existing single-GPU/full-grid choices; it does
not yet choose intermediate subgrids or rank CPU performance. Coefficients
are initial fits to the 10k/64k measurements, isolated in
`GPUPlanningModule.F90` for testing/calibration. They are estimates, not
promised timings. The sparse 12k case now selects one GPU too; the old 10k
boundary is gone. Higher fill at the same N can favor the full grid.

Memory admission includes the dense accumulator, padded bitmap, full-fill
bounds for changing CSR matrices, panel copies/gathers, merge buffers,
fixed WH, row/scan storage and a runtime allowance. Intermediate CSR counts
must also fit their current 32-bit representation. All ranks agree on the
plan, using the smallest full-grid memory budget; root mode checks only
its active GPU. The selector checks physical GPU UUIDs before all-rank
NCCL initialization, including per-rank device visibility. Duplicate
mappings or an unavailable multi-node opt-in cannot select that layout.

The default budget is **80% of currently free GPU memory**. Set
`NTPOLY_GPU_MEMORY_MB=<nonnegative integer MiB>` to reduce it further;
zero disables admission. The cap never increases the live-memory budget.
Neither `root` nor `all` bypasses safety checks. If the requested/admissible
GPU layouts do not fit, the same TRS2 call executes its CPU iteration with
its already-constructed WH/X0. No GPU iteration state needs reconstruction.
Memory is checked again collectively after CUDA/NCCL setup and before
matrix payload allocation; a rejected bootstrap is torn down before CPU
execution. CUDA memory-query/selection time is included in API time, but
is outside the narrower GPU-setup phase counter.

This is conservative admission for native TRS2, not a hard allocator quota
or a VRAM reservation. Full-fill bounds may reject sparse cases that would
fit in practice. Other processes can consume memory after a query, and
opaque runtime workspace is covered by an allowance/headroom rather than
a formal bound. Unexpected CUDA/NCCL errors retain the existing error path;
there is no mid-iteration repartition/restart. CPU setup/fallback and raw
mini-app helpers have their own memory requirements.

GPU validation covers normal 10k/12k selection, a 6 GiB cap that changes
coalescing to distribution, CPU fallback for zero/tiny budgets or unsafe
forced-root selection, duplicate GPU mappings, hidden GPUs, and absent
multi-node permission. All nineteen saved small output matrices agree
with CPU within 2.4e-15 relative norm; the 64k path retains its scalar
fingerprints. The pure planning unit suite also checks high-fill ranking,
32-bit panel limits and arithmetic overflow without allocating huge matrices.
Replay: `validation/rikyu-memory-policy-20261111/analyze.py`.

### Host regression checkpoint

The host suite checks the actual operation bodies under serial CUDA emulation,
Fortran CSR/planning helpers, memory/index limits and MPI 1/2/4-rank solver
parity. It does not establish CUDA runtime correctness or GPU performance.
With a CPU NTPoly build available:

```sh
python3 MiniApps/YukawaGemm/tests/run_host_checks.py --ntpoly-build /path/to/cpu-build
```

The default build directory is `BuildCheck`. The checkpoint was also tested
against a fresh CPU-only build made with:

```sh
cmake -S . -B /tmp/ntg-checkpoint-cpu -DFORTRAN_ONLY=ON -DNOSWIG=ON \
  -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_FLAGS=-ffree-line-length-none
cmake --build /tmp/ntg-checkpoint-cpu --target NTPoly -j 4
```

GPU runtime/transport and native-API validation remain separate RIKYU jobs,
as recorded in Findings 22–25. Generated matrices and transfer fragments,
local builds and machine-specific toolchain settings are ignored by Git.
