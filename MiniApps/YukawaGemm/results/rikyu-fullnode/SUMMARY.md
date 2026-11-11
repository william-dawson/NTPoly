# Full-Node CPU Baseline — RIKYU GB200

**Matrix:** sc 40³ = 64,000 atoms (Yukawa, lattice=2.7, screening=1.0, rcut=6.0, jitter=0.05, seed=42, diag=1.0)
**Input nnz:** 2,591,820 (40.5/row) → Output nnz varies by process grid (152–223/row)
**Hardware:** 1× GB200 node, 144 Grace cores, 1.6 TB RAM
**Toolchain:** nvhpc-hpcx (nvfortran + HPC-X MPI), `UCX_TLS=sm,sysv,posix` (single-node, no IB)
**Mode:** `fixed` (C ← A·A), 3 warmup + 10 timed rounds
**Date:** 2026-11-04

## Results

| Config        | np  | omp | Median (ms) | Min (ms) | nnz/row | Trace          |
|---------------|-----|-----|-------------|----------|---------|----------------|
| np4_omp1      | 4   | 1   | 550         | 539      | 222.7   | 64287.04835…   |
| np4_omp36     | 4   | 36  | 331         | 269      | 172.2   | 64287.04835…   |
| np8_omp1      | 8   | 1   | 300         | 296      | 198.5   | 64287.04835…   |
| np8_omp18     | 8   | 18  | 217         | 180      | 198.5   | 64287.04835…   |
| np16_omp1     | 16  | 1   | 280         | 275      | 198.5   | 64287.04835…   |
| np16_omp9     | 16  | 9   | 191         | 165      | 198.5   | 64287.04835…   |
| np36_omp1     | 36  | 1   | 201         | 200      | 181.6   | 64287.04835…   |
| **np36_omp4** | **36**  | **4**   | **128**         | **105**      | **181.6**   | **64287.04835…**   |
| np72_omp1     | 72  | 1   | 153         | 148      | 152.1   | 64287.04835…   |
| np144_omp1    | 144 | 1   | 200         | 169      | 152.0   | 64287.04835…   |

## Key Findings

### Correctness
- **Trace is bit-identical** across all 10 configs (64287.048350011…), confirming the distributed
  multiply produces the same result regardless of process grid shape. This is the
  correctness anchor for the GPU port — any GPU implementation must reproduce this trace.
- Grand sum and asymmetry vary slightly (98291–98361, 1.4e-3–2.0e-3) because the pruning
  threshold falls on different product values depending on how the matrix is block-distributed.
  This is expected and does not indicate errors.

### Performance
- **Best overall: np=36 omp=4 at 128 ms median (105 ms min).** Hybrid MPI×OpenMP wins because
  OMP threads parallelize the local block multiply without MPI communication overhead.
- **Best pure-MPI: np=72 at 153 ms median (148 ms min).** np=144 is slower (200 ms) due to
  communication overhead at this matrix size — more ranks = more Allgather/Reduce traffic
  for the same work.
- **Scaling is sublinear:** 4→72 ranks (18×) gives only 550→153 ms (3.6× speedup). The local
  GEMM's dense accumulator creates memory pressure and cache thrashing that limits scaling.

### Fill (nnz/row) varies with process grid
- np=4 → 223, np=16 → 199, np=36 → 182, np=72/144 → 152
- Fewer ranks = larger local blocks = more products above threshold = more fill
- Timings are only directly comparable at equal nnz/row:
  - np=72 (152.1) vs np=144 (152.0): np=72 wins, 153 vs 200 ms
  - np=8/16 (198.5): hybrid helps — np=16 omp=9 (191 ms) beats np=16 omp=1 (280 ms)
  - np=36 (181.6): hybrid helps enormously — omp=4 (128 ms) vs omp=1 (201 ms)

### Memory bottleneck
- The local sparse GEMM memory pool (`ConstructMatrixMemoryPool.f90`) allocates **dense**
  arrays sized `columns × rows` per rank (5 arrays: pruned_list, value_array, dirty_array,
  hash_index, inserted_per_bucket). Memory per rank ≈ 24 × (N/√np)² bytes.
- This limits matrix size on a single node:
  - 64k atoms: fits at all configs (largest: 2.7 GB/rank at np=4)
  - 175k atoms: OOMs at np=144 (5.1 GB pool, but 11 GB cgroup limit + MPI overhead)
  - 885k atoms: OOMs at all single-node configs (needs 12+ nodes, 1700+ ranks)
- **This dense accumulator is the primary target for the GPU port** — replacing it with a
  sparse accumulator or GPU-friendly hash table will unlock both larger matrices and GPU acceleration.

## Reference for GPU Port

The GPU port must reproduce the **trace = 64287.048350011…** at 64k atoms.
The performance target is to beat **128 ms median** (np=36 omp=4) on a single GB200 GPU.
