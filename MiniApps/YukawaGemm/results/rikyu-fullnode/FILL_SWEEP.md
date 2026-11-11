# Fill Sweep: GPU vs CPU at Realistic DFT Density

**Date:** 2026-11-04
**Hardware:** NVIDIA GB200, 144 Grace cores, 1× Blackwell GPU (198 GB HBM)
**Matrix:** sc 32³ = 32,768 atoms, Yukawa potential, varying rcut to sweep fill
**Mode:** C ← A·A (fixed), threshold=1e-4, np=36 omp=4 (CPU best config)

## The concern

Real DFT density matrices have ~1-5% fill (nonzero fraction). Our original
benchmark used ~0.06% fill — far sparser than realistic. This sweep tests
whether the GPU advantage holds at realistic fill levels.

## Results

| rcut | Input fill | Input nnz/row | CPU median | GPU (cuSPARSE) | Speedup | CPU nnz/row (thresh.) | GPU nnz/row (full) |
|------|-----------:|---------------:|-----------:|---------------:|--------:|----------------------:|-------------------:|
| 6    | 0.12%      | 40             | 64.6 ms    | **7.0 ms**     | **9.2×** | 171                  | 228                |
| 10   | 0.58%      | 189            | 234 ms     | **118.6 ms**   | **2.0×** | 340                  | 1143               |
| 15   | 1.83%      | 599            | 553 ms     | **CRASH**      | —       | 321                  | ~dense (cuSPARSE refuses) |
| 20   | 4.0%       | 1314           | 1747 ms    | **CRASH**      | —       | 336                  | ~dense (cuSPARSE refuses) |

## What happened

### The GPU speedup collapses with fill

At 0.12% fill (our original regime), the GPU is 9.2× faster. At 0.58%,
it's only 2× faster. At 1.83% and above, **cuSPARSE crashes** — all four
SpGEMM algorithms (DEFAULT, ALG1, ALG2, ALG3) return
`CUSPARSE_STATUS_INSUFFICIENT_RESOURCES`.

### Why cuSPARSE fails at high fill

cuSPARSE SpGEMM computes the **full product** — every nonzero in C=A*A,
including entries below the threshold. The full product size scales as
~z² per row (where z = input nnz/row):

- z=40: product ~228 nnz/row (fine)
- z=189: product ~1143 nnz/row (3× more than thresholded)
- z=599: product ~359k nnz/row → **exceeds N=32k** → product is nearly dense
- z=1314: product ~1.7M nnz/row → **way over N** → product IS dense

cuSPARSE's internal data structures can't handle a product that approaches
dense. The "insufficient resources" error occurs in the work estimation
step, before any computation begins.

### Why the CPU keeps working

NTPoly's local GEMM applies the threshold **during** the multiply — it
prunes each partial product as it accumulates. This keeps the output at
**~336 nnz/row regardless of input fill**:

```
rcut=6:  input 40 nnz/row  → product 171 nnz/row (thresholded)
rcut=10: input 189 nnz/row → product 340 nnz/row (thresholded)
rcut=15: input 599 nnz/row → product 321 nnz/row (thresholded)
rcut=20: input 1314 nnz/row → product 336 nnz/row (thresholded)
```

The threshold (1e-4) cuts off the same fraction of products regardless
of input density. The CPU only stores and processes entries above
threshold — the rest are discarded immediately during accumulation.

### The unfair comparison

At rcut=10, the GPU computes 1143 nnz/row (full product) while the CPU
computes 340 nnz/row (thresholded) — 3.4× more work for the GPU. The GPU's
per-nnz throughput is 6.7× higher than the CPU's (2.0× speedup despite 3.4×
more work), but the extra work erodes the advantage.

At rcut=15+, the full product exceeds the matrix dimension — it's
essentially dense — and cuSPARSE can't handle it at all.

## Implication for the GPU port

**cuSPARSE SpGEMM is NOT viable for realistic DFT fill levels.**

The porting strategy must change:

1. **Cannot replace NTPoly's local GEMM with cuSPARSE.** cuSPARSE has no
   thresholding — it computes the full product, which becomes dense at
   realistic fill. All algorithms crash.

2. **Must port NTPoly's algorithm to the GPU.** The hash-based accumulator
   with per-element thresholding (MultiplyBlock + PruneList) is what keeps
   the product small. This algorithm must be implemented on the GPU.

3. **The dense accumulator is still a problem.** NTPoly's current
   implementation uses dense `columns × rows` arrays per rank. On GPU, this
   limits matrix size (though GPU HBM is 198 GB vs 400 GB CPU cgroup).
   A GPU-friendly sparse accumulator (e.g., a hash table in GPU memory)
   would be better.

4. **OpenACC is the right approach.** Parallelizing the existing Fortran
   algorithm with OpenACC directives (parallelizing over columns, keeping
   the thresholding logic) is the path forward. nvfortran supports OpenACC
   with `-acc -gpu=managed`.

### Alternative: cuSPARSE with structure reuse

One possible hybrid: use cuSPARSE to compute the product STRUCTURE (which
entries are nonzero), apply thresholding to the structure, then compute
only the surviving values. cuSPARSE doesn't support this directly, but
it could be done with:
1. cuSPARSE SpGEMM with `beta=0` to get structure only
2. A custom GPU kernel to threshold and compact
3. A second cuSPARSE pass (or custom kernel) to compute surviving values

This is more complex but might work at moderate fill. At high fill (5%+),
even the structure is nearly dense, so this won't help.

## Verdict

The original 9.2× GPU speedup at 0.06% fill was misleading. At realistic
DFT fill (1-5%), cuSPARSE can't run at all. The GPU port must replicate
NTPoly's thresholding-during-multiply algorithm, not use cuSPARSE.
