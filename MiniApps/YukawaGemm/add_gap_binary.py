#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy"]
# ///
"""Open a spectral gap at half filling in an NTPoly binary matrix.

Modifies the DIAGONAL entries in place: diag(i) += delta * (+1/2) for even
sites, -= delta/2 for odd sites (ionic-Hubbard style alternating on-site
energies). On a simple-cubic grid with even dimensions the index parity is
the checkerboard sublattice, so the result is an insulator with a gap at
half filling (trace target = n/2) — the regime where TRS2 purification
converges cleanly. The unmodified matrix (metallic, no gap) remains
available for stress tests.

Binary format (NTPoly): int32 rows, int32 cols, int32 complex_flag,
int64 nnz, then nnz COO triplets (int32 col, int32 row, double value),
1-based indices.
"""

import sys

import numpy as np


def main(src, dst, delta=1.0):
    triplet_dtype = np.dtype([("c", "<i4"), ("r", "<i4"), ("v", "<f8")])
    with open(src, "rb") as f:
        rows, cols, cflag = np.fromfile(f, dtype=np.int32, count=3)
        nnz = int(np.fromfile(f, dtype=np.int64, count=1)[0])
        data = np.fromfile(f, dtype=triplet_dtype, count=nnz)
    assert data.shape[0] == nnz

    diag = data["r"] == data["c"]
    n_diag = int(diag.sum())
    shift = np.where(data["r"][diag] % 2 == 0, +delta / 2.0, -delta / 2.0)
    data["v"][diag] += shift

    with open(dst, "wb") as f:
        np.int32(rows).tofile(f)
        np.int32(cols).tofile(f)
        np.int32(cflag).tofile(f)
        np.int64(nnz).tofile(f)
        data.tofile(f)
    lo = data["v"][diag].min()
    hi = data["v"][diag].max()
    print(f"gapped {src} -> {dst}: n={rows} nnz={nnz} "
          f"diag entries={n_diag} diag range=[{lo:.3f},{hi:.3f}] delta={delta}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], float(sys.argv[3]) if len(sys.argv) > 3 else 1.0)
