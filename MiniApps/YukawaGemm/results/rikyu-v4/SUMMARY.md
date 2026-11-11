# v4 cuSPARSE Hybrid Benchmark Results

**Date:** 2026-11-04
**Jobs:** 281609 (alg sweep), 281610 (per-phase), 281613 (final benchmark)
**RIKYU account:** rkp00012
**Hardware:** GB200 (Blackwell B200, SM 10.0, 198 GB HBM, 8 TB/s HBM bandwidth)

## Optimization progression

| Version | Key change | 64k median | vs CPU (128 ms) |
|---------|-----------|-----------:|------:|
| v1 | Per-round cudaMalloc/cudaFree | 61 ms | 2.1× |
| v3 | Persistent A, skip per-round copy | 18.2 ms | 7.0× |
| **v4** | **Persistent C/work buffers + warm calls** | **7.58 ms** | **16.9×** |

## Final results (job 281613)

### 64k matrix (sc 40³, 64,000 atoms, 0.06% fill)

| Metric | Value |
|--------|------:|
| Input nnz | 2,591,820 |
| Full product nnz | 15,040,832 |
| Thresholded nnz | 4,865,120 |
| **GEMM median** | **7.583 ms** |
| GEMM min | 7.561 ms |
| GEMM max | 7.610 ms |
| GEMM mean | 7.587 ms |
| Trace | 64287.048350011 |
| Grand sum | 98232.413499784 |
| Asymmetry | 1.56e-16 |
| **Speedup vs CPU** | **16.9×** |

### 32k fill sweep (sc 32³, 32,768 atoms)

| rcut | Fill | Input nnz | Backend | GPU median | CPU (np=36) | Speedup | Output nnz | Trace | Asymmetry |
|-----:|-----:|----------:|---------|----------:|------------:|--------:|-----------:|------:|----------:|
| 6 | 0.12% | 1,308,276 | cuSPARSE | 4.054 ms | 64.6 ms | **16.0×** | 2,444,530 | 32913.930651629 | 1.42e-16 |
| 10 | 0.58% | 6,175,328 | cuSPARSE | 50.42 ms | 234 ms | **4.6×** | 3,747,210 | 32914.015272987 | 1.80e-16 |
| 15 | 1.83% | 19,618,918 | OpenACC | 850.3 ms | 553 ms | 0.65× | 3,751,584 | 32914.015291414 | 0 |
| 20 | 4.0% | 43,056,662 | OpenACC | 3203.2 ms | 1747 ms | 0.55× | 3,751,584 | 32914.015291415 | 0 |

## Algorithm benchmark (64k, job 281609)

| Algorithm | First call (cold) | Steady-state | Status |
|-----------|------------------:|-------------:|--------|
| DEFAULT | 30.8 ms | 6.5 ms | ✓ (maps to ALG1) |
| ALG1 | 30.7 ms | **6.5 ms** | ✓ Best |
| ALG2 | 20.2 ms | 0.13 ms (C_nnz=0) | ✗ Silent failure |
| ALG3 | 17.6 ms | 0.13 ms (C_nnz=0) | ✗ Silent failure |

## Per-phase GPU timing (64k steady-state, job 281610)

| Phase | GPU time | % of total |
|-------|---------:|-----------:|
| cuSPARSE compute (warm) | 6.3 ms | 84% |
| cuSPARSE copy | 0.11 ms | 1.5% |
| cuSPARSE workEstimation | ~0.1 ms | 1.3% |
| GPU threshold mark + count | 0.15 ms | 2% |
| CPU prefix scan | 0.04 ms | 0.5% |
| GPU compact | 0.50 ms | 6.7% |
| D2H download | 0.32 ms | 4.3% |
| **Total** | **7.58 ms** | 100% |

## Bandwidth model (64k)

| Metric | Value |
|--------|------:|
| GB200 HBM bandwidth (peak) | 8 TB/s |
| Estimated memory traffic | ~241 MB |
| Streaming BW floor | 0.030 ms |
| Compute floor (FP64) | 0.013 ms |
| Actual cuSPARSE compute | 6.3 ms |
| Achieved bandwidth | 38 GB/s |
| % of peak HBM | 0.48% |

SpGEMM is latency-bound (hash table ops, random access), not bandwidth-bound.
