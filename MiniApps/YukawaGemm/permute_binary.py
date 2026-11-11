#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy"]
# ///
"""Apply a fixed-seed random symmetric permutation P*A*P^T to an NTPoly
binary matrix, writing a new binary file.

This achieves the same block-homogenizing effect as --load-balance without
running NTPoly's PermuteMatrix at runtime.  PermuteMatrix computes the
permutation via NTPoly's distributed multiply, whose dense local accumulator
(~24 B per element of the local block) OOMs for very large matrices on the
(4,4) grid at 885k (FINDINGS F1).  Pre-permuting the file sidesteps that:
runs then use the permuted file WITHOUT --load-balance.

Trace/grand-sum/nnz are invariant under the symmetric permutation
(trace(P A^2 P^T) = trace(A^2)), up to floating-point summation order.

Binary format (NTPoly): int32 rows, int32 cols, int32 complex_flag,
int64 nnz, then nnz COO triplets (int32 col, int32 row, double value),
1-based indices.
"""

import sys

import numpy as np


def main(src, dst, seed=424242):
    triplet_dtype = np.dtype([("c", "<i4"), ("r", "<i4"), ("v", "<f8")])
    with open(src, "rb") as f:
        rows, cols, cflag = np.fromfile(f, dtype=np.int32, count=3)
        nnz = int(np.fromfile(f, dtype=np.int64, count=1)[0])
        data = np.fromfile(f, dtype=triplet_dtype, count=nnz)
    assert data.shape[0] == nnz, (
        f"expected {nnz} triplets, file has {data.shape[0]}")

    rng = np.random.default_rng(seed)
    # perm[old_0based] -> new_0based
    perm = rng.permutation(int(rows)).astype(np.int64)
    r_new = perm[data["r"].astype(np.int64) - 1] + 1
    c_new = perm[data["c"].astype(np.int64) - 1] + 1

    # Sort by (row, col) so the file has canonical CSR-like ordering.
    order = np.lexsort((c_new, r_new))
    out = np.empty(nnz, dtype=triplet_dtype)
    out["c"] = c_new[order].astype(np.int32)
    out["r"] = r_new[order].astype(np.int32)
    out["v"] = data["v"][order]

    with open(dst, "wb") as f:
        np.int32(rows).tofile(f)
        np.int32(cols).tofile(f)
        np.int32(cflag).tofile(f)
        np.int64(nnz).tofile(f)
        out.tofile(f)
    print(f"permuted {src} -> {dst}: n={rows} nnz={nnz} seed={seed}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
