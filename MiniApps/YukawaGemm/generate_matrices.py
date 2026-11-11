#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy", "scipy"]
# ///
"""Generate synthetic sparse matrices from a Yukawa potential on realistic
atomic geometries, for benchmarking NTPoly's sparse GEMM kernel.

The matrix mimics a single-orbital Hamiltonian: atoms are placed on a
simple-cubic or FCC lattice (with realistic lattice constants, optionally
jittered to mimic thermal disorder / amorphous structure), and

    V_ij = exp(-r_ij / lambda) / r_ij    for r_ij <= Rcut
         = 0                             otherwise

with a constant diagonal (on-site energy). The result is a symmetric matrix
with exponentially decaying off-diagonals and a physically sensible
neighbor count (~30-80 per row at the default cutoff), which is the sparsity
regime density-matrix purification operates in.

Run with uv (it creates an isolated environment from this file's metadata):

    uv run generate_matrices.py --preset small
    uv run generate_matrices.py --preset small --out matrices/small.mtx

or override the geometry directly:

    uv run generate_matrices.py --geometry fcc --side 25 --lattice-const 3.6 \
        --screening 1.0 --rcut 6.0 --jitter 0.1 --seed 7 --out my.mtx
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.io import mmwrite
from scipy.sparse import coo_matrix
from scipy.spatial import cKDTree

# Presets: "realistic" geometries in the spirit of metals/semiconductors.
#  - sc, a=2.7 A : simple-cubic, one orbital per site, nearest-neighbor
#                  spacing 2.7 A (generic metal).
#  - fcc, a=3.6 A: face-centered cubic (Cu-like), nearest-neighbor 2.55 A.
#  - jitter=0.05 : Gaussian displacement (fraction of the lattice constant,
#                  ~0.14 A rms - thermal disorder). This also smears the
#                  degenerate Yukawa product values of a perfect crystal,
#                  which otherwise pile up at the sparsity threshold and
#                  make fill (and fingerprints) rank-count dependent.
#  - lambda=1.0 A: Thomas-Fermi-like screening length.
#  - rcut=6.0 A  : interaction cutoff -> ~46 nnz/row (sc) / ~77 (fcc).
PRESETS: dict[str, dict] = {
    "micro":     dict(geometry="sc",  side=8,  lattice_const=2.7, screening=1.0, rcut=6.0, jitter=0.0),   # 512, smoke test
    "tiny":      dict(geometry="sc",  side=16, lattice_const=2.7, screening=1.0, rcut=6.0, jitter=0.05),  # 4,096
    "small":     dict(geometry="sc",  side=32, lattice_const=2.7, screening=1.0, rcut=6.0, jitter=0.05),  # 32,768
    "small_md":  dict(geometry="sc",  side=32, lattice_const=2.7, screening=1.0, rcut=6.0, jitter=0.12), # 32,768 strongly disordered
    "medium":    dict(geometry="fcc", side=25, lattice_const=3.6, screening=1.0, rcut=6.0, jitter=0.05),  # 62,500
    "large":     dict(geometry="sc",  side=96, lattice_const=2.7, screening=1.0, rcut=6.0, jitter=0.05),  # 884,736
    "xlarge":    dict(geometry="fcc", side=63, lattice_const=3.6, screening=1.0, rcut=6.0, jitter=0.05),  # 1,000,188
}


def build_positions(
    geometry: str, side: int, a: float, jitter: float, rng: np.random.Generator
) -> np.ndarray:
    """Return (n_atoms, 3) positions in Angstrom."""
    cells = np.indices((side, side, side), dtype=float).reshape(3, -1).T
    if geometry == "sc":
        pos = cells
    elif geometry == "fcc":
        basis = np.array(
            [[0.0, 0.0, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]]
        )
        pos = (cells[:, None, :] + basis[None, :, :]).reshape(-1, 3)
    else:
        raise ValueError(f"unknown geometry '{geometry}' (use 'sc' or 'fcc')")
    pos *= a
    if jitter > 0.0:
        pos += rng.normal(0.0, jitter * a, pos.shape)
    return pos


def yukawa_matrix(
    pos: np.ndarray, lam: float, rcut: float, diag: float
) -> coo_matrix:
    """Build the sparse symmetric Yukawa matrix for the given positions."""
    n = pos.shape[0]
    tree = cKDTree(pos)
    pairs = tree.query_pairs(rcut, output_type="ndarray")
    d = np.linalg.norm(pos[pairs[:, 0]] - pos[pairs[:, 1]], axis=1)
    d = np.maximum(d, 1e-12)  # guard against coincident points after jitter
    v = np.exp(-d / lam) / d

    rows = np.concatenate([pairs[:, 0], pairs[:, 1], np.arange(n)])
    cols = np.concatenate([pairs[:, 1], pairs[:, 0], np.arange(n)])
    vals = np.concatenate([v, v, np.full(n, diag, dtype=v.dtype)])
    return coo_matrix((vals, (rows, cols)), shape=(n, n))


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate Yukawa-potential sparse matrices for GEMM benchmarking."
    )
    parser.add_argument("--preset", choices=sorted(PRESETS), help="use a named preset")
    parser.add_argument("--list-presets", action="store_true", help="list presets and exit")
    parser.add_argument("--geometry", choices=["sc", "fcc"], help="lattice geometry")
    parser.add_argument("--side", type=int, help="lattice cells per dimension")
    parser.add_argument("--lattice-const", type=float, help="lattice constant in Angstrom")
    parser.add_argument("--screening", type=float, help="Yukawa screening length lambda in Angstrom")
    parser.add_argument("--rcut", type=float, help="interaction cutoff radius in Angstrom")
    parser.add_argument("--jitter", type=float, help="position jitter as a fraction of the lattice constant")
    parser.add_argument("--diag", type=float, default=1.0, help="diagonal (on-site) value; 1.0 keeps repeated squaring alive (density-matrix-like), 0.0 lets the threshold annihilate the matrix")
    parser.add_argument("--seed", type=int, default=42, help="random seed (jitter)")
    parser.add_argument("--out", help="output .mtx path (default: yukawa_<geometry>_<n>.mtx)")
    args = parser.parse_args()

    if args.list_presets:
        print(f"{'name':<10} {'geom':<5} {'side':>5} {'atoms':>10} {'~nnz/row':>9}")
        for name, p in PRESETS.items():
            n = p["side"] ** 3 * (4 if p["geometry"] == "fcc" else 1)
            per_row = sphere_density_nnz(p)
            print(f"{name:<10} {p['geometry']:<5} {p['side']:>5} {n:>10,} {per_row:>9.0f}")
        return

    if not args.preset:
        parser.error("either --preset or explicit geometry options are required")
    params = dict(PRESETS[args.preset])
    for key, flag in [
        ("geometry", args.geometry),
        ("side", args.side),
        ("lattice_const", args.lattice_const),
        ("screening", args.screening),
        ("rcut", args.rcut),
        ("jitter", args.jitter),
    ]:
        if flag is not None:
            params[key] = flag

    rng = np.random.default_rng(args.seed)
    pos = build_positions(params["geometry"], params["side"], params["lattice_const"], params["jitter"], rng)
    n = pos.shape[0]
    print(f"geometry   : {params['geometry']} {params['side']}^"
          f"{'3' if params['geometry'] == 'sc' else '3 x 4 atoms/cell'}")
    print(f"atoms      : {n:,}")
    print("building neighbor pairs ...", flush=True)
    mat = yukawa_matrix(pos, params["screening"], params["rcut"], args.diag)
    nnz = mat.nnz

    out = args.out or f"yukawa_{params['geometry']}_{n}.mtx"
    mmwrite(out, mat, symmetry="general", precision=17)
    file_bytes = Path(out).stat().st_size

    sidecar = {
        "generator": "generate_matrices.py",
        "preset": args.preset,
        **params,
        "diag": args.diag,
        "seed": args.seed,
        "n_atoms": n,
        "nnz": nnz,
        "nnz_per_row": nnz / n,
        "file": out,
        "file_bytes": file_bytes,
        "note": "V_ij = exp(-r/lambda)/r for r<=rcut, else 0; constant diagonal",
    }
    with open(str(out) + ".json", "w") as fh:
        json.dump(sidecar, fh, indent=2)

    print(f"matrix dim : {n:,}")
    print(f"nnz        : {nnz:,} ({nnz / n:.1f} per row)")
    print(f"wrote      : {out} ({file_bytes / 1e6:.1f} MB)")
    print(f"sidecar    : {out}.json")


def sphere_density_nnz(p: dict) -> float:
    """Estimated neighbors within rcut for an infinite lattice of this geometry."""
    atoms_per_a3 = (4.0 if p["geometry"] == "fcc" else 1.0) / p["lattice_const"] ** 3
    return 4.0 / 3.0 * np.pi * p["rcut"] ** 3 * atoms_per_a3


if __name__ == "__main__":
    main()
