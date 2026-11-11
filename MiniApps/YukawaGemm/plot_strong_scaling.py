# /// script
# requires-python = ">=3.10"
# dependencies = ["matplotlib>=3.8", "numpy>=1.24"]
# ///
"""Strong scaling plot for NTPoly GPU GEMM benchmarks on RIKYU.

Two configurations on the 64k Yukawa matrix:
  - 2D SUMMA + symmetry exploitation (no load balance): diagonal-block
    cuSPARSE pathology caps scaling.
  - 2D SUMMA + load balance (random permutation, no symmetric exchange):
    uniform blocks, no exchange communication.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ── Data (64k matrix, median GEMM ms per round) ─────────────────────────────
sym     = {1: 7.57, 4: 4.13, 9: 4.18, 16: 2.73}   # no-LB + --symmetric
lbnosym = {1: 7.92, 4: 3.08, 9: 1.78, 16: 1.27}  # LB, no --symmetric (IB)

cpu_64k = 128.0

fig, ax = plt.subplots(figsize=(9, 6), dpi=150, squeeze=True)

nps_sym = np.array(sorted(sym))
nps_lb  = np.array(sorted(lbnosym))
t_sym = np.array([sym[n] for n in nps_sym])
t_lb  = np.array([lbnosym[n] for n in nps_lb])

ax.loglog(nps_sym, t_sym, "o-", color="#1f77b4", lw=2, ms=9,
          label="2D SUMMA + symmetry (no LB)")
ax.loglog(nps_lb, t_lb, "D-", color="#2ca02c", lw=2, ms=8,
          label="2D SUMMA + load balance (no sym exchange)")

# Ideal lines from each config's np=1
xx = np.linspace(1, 16, 100)
ax.loglog(xx, sym[1] / xx, ":", color="#1f77b4", lw=0.9, alpha=0.6)
ax.loglog(xx, 7.57 / xx, "--", color="k", lw=0.8, alpha=0.4,
          label="ideal (t(1)/np)")

# CPU baseline
ax.axhline(cpu_64k, color="#7f7f7f", ls="-.", lw=1.2, alpha=0.7)
ax.text(4.4, cpu_64k * 1.04, "CPU baseline (64k): 128 ms",
        fontsize=8.5, color="#7f7f7f")

# Value annotations
for n, t in sym.items():
    ax.annotate(f"{t:.2f}", xy=(n, t), xytext=(n * 1.08, t * 1.13),
                fontsize=8, color="#1f77b4")
for n, t in lbnosym.items():
    if n == 9:
        ax.annotate(f"{t:.2f}", xy=(n, t), xytext=(n * 0.82, t * 1.30),
                    fontsize=8, color="#2ca02c")
    else:
        ax.annotate(f"{t:.2f}", xy=(n, t), xytext=(n * 1.08, t * 0.82),
                    fontsize=8, color="#2ca02c")

# Speedup box (inside axes, near bottom-center)
ax.text(4.5, 1.12,
        f"np=16:  {sym[1]/sym[16]:.1f}×  (sym)\nnp=16:  {sym[1]/lbnosym[16]:.1f}×  (LB no-sym)",
        fontsize=10, fontfamily="monospace",
        bbox=dict(boxstyle="round,pad=0.4", fc="lightyellow", ec="gray", alpha=0.9))

ax.set_xlabel("MPI ranks (= GPUs)", fontsize=12, fontweight="bold")
ax.set_ylabel("GEMM time (ms)", fontsize=12, fontweight="bold")
ax.set_title("Strong Scaling — 64k Yukawa matrix on RIKYU (GB200)\n"
             "cuSPARSE 2D SUMMA: symmetry exploit vs load-balanced no-exchange",
             fontsize=13, fontweight="bold")
ax.set_xlim(0.75, 20)
ax.set_ylim(1.0, 160)
ax.set_xticks([1, 2, 4, 8, 16])
ax.set_xticklabels(["1", "2", "4", "8", "16"])
ax.set_yticks([1, 2, 4, 8, 16, 32, 64, 128])
ax.set_yticklabels(["1", "2", "4", "8", "16", "32", "64", "128"])
ax.grid(True, which="both", ls="-", lw=0.3, alpha=0.4)
ax.legend(loc="lower left", fontsize=9.5)

fig.tight_layout()
fig.savefig("strong_scaling.png", dpi=150, bbox_inches="tight")
print("Saved strong_scaling.png")
