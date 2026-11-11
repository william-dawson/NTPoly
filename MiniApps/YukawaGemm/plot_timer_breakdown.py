# /// script
# requires-python = ">=3.10"
# dependencies = ["matplotlib>=3.8", "numpy>=1.24"]
# ///
"""Timer breakdown + cuSPARSE per-block throughput analysis (64k matrix)."""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ── Parsed from RIKYU job output (warm rounds) ─────────────────────────────
# Diagonal-block SpGEMM *compute* (warm) times and flop counts per block.
#   np: (list of compute ms per diagonal block, list of flops per block (millions))
diag_blocks = {
    1:  ([6.31], [106.8]),
    4:  ([3.25, 3.26], [51.9, 51.9]),
    9:  ([2.20, 3.00, 3.02], [33.5, 33.8, 33.5]),
    16: ([1.65, 1.65, 1.66, 2.28], [24.3, 24.6, 24.6, 24.2]),
}

# Wall-time decomposition per np (slowest diagonal rank)
#   np: (slowest SpGEMM total ms, threshold pipeline ms on that rank, measured wall ms)
wall = {
    1:  (6.55, 1.02, 7.57),
    4:  (3.46, 0.63, 4.13),
    9:  (3.94, 0.43, 4.18),
    16: (2.42, 0.28, 2.73),
}

nps = [1, 4, 9, 16]
t1 = wall[1][2]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5), dpi=150)

# ── Panel 1: wall-time decomposition ──────────────────────────────────────
x = np.arange(len(nps))
w = 0.5
spgemm = [wall[n][0] for n in nps]
thresh = [wall[n][1] for n in nps]
resid  = [wall[n][2] - wall[n][0] - wall[n][1] for n in nps]

b1 = ax1.bar(x, spgemm, w, label="cuSPARSE SpGEMM (slowest diag)", color="#1f77b4")
b2 = ax1.bar(x, thresh, w, bottom=spgemm, label="GPU threshold pipeline", color="#ff7f0e")
b3 = ax1.bar(x, resid,  w, bottom=[spgemm[i]+thresh[i] for i in range(4)],
             label="MPI sync / residual", color="#2ca02c")

# Ideal-scaling reference curve (t(1)/np)
xx = np.linspace(1, 16, 200)
ideal = t1 / xx
ax1.plot(np.interp(xx, nps, x), ideal, "k--", lw=1.2, alpha=0.7,
         label="ideal strong scaling t(1)/np")

for i, n in enumerate(nps):
    eff = t1 / (wall[n][2] * n) * 100
    ax1.annotate(f"{wall[n][2]:.2f} ms\n({eff:.0f}% eff)",
                 xy=(x[i], wall[n][2]), xytext=(x[i], wall[n][2] + 0.35),
                 ha="center", fontsize=9, fontweight="bold")

ax1.set_xticks(x)
ax1.set_xticklabels([f"np={n}" for n in nps], fontsize=11)
ax1.set_ylabel("GEMM time per round (ms)", fontsize=12, fontweight="bold")
ax1.set_title("Panel 1 — Where the wall time goes (64k matrix)\n"
              "wall = slowest diagonal block's SpGEMM + its threshold",
              fontsize=12, fontweight="bold")
ax1.set_ylim(0, 9)
ax1.legend(fontsize=9, loc="upper right")
ax1.grid(True, axis="y", lw=0.3, alpha=0.4)

# ── Panel 2: cuSPARSE per-block throughput ────────────────────────────────
for j, n in enumerate(nps):
    for comp, fl in zip(diag_blocks[n][0], diag_blocks[n][1]):
        gflops = fl / comp  # Mflop / ms = Gflop/s
        color = "#d62728" if gflops < 13 else "#1f77b4"
        ax2.scatter(n, gflops, s=90, color=color, zorder=5)

# median line
for a, b in [(1, 4), (4, 9), (9, 16)]:
    ma = np.median([f / c for c, f in zip(diag_blocks[a][0], diag_blocks[a][1])])
    mb = np.median([f / c for c, f in zip(diag_blocks[b][0], diag_blocks[b][1])])
    ax2.plot([a, b], [ma, mb], color="#1f77b4", lw=2, alpha=0.6)

ax2.axhline(y=16.9, color="k", ls="--", lw=1.2, alpha=0.7,
            label="np=1 throughput (16.9 GF/s)")
ax2.annotate("pathological block\n(same flops, ~1.4× slower,\nsets the wall time)",
             xy=(16, 10.6), xytext=(10.5, 12.2), fontsize=9, color="#d62728",
             arrowprops=dict(arrowstyle="->", color="#d62728", lw=1.2))
ax2.annotate("np=9 diag blocks:\n21k×21k SLOWER than\nnp=4's 32k×32k",
             xy=(9, 11.1), xytext=(5.0, 8.6), fontsize=9, color="#d62728",
             arrowprops=dict(arrowstyle="->", color="#d62728", lw=1.2))

ax2.set_xlabel("MPI ranks", fontsize=12, fontweight="bold")
ax2.set_ylabel("cuSPARSE throughput per diagonal block (GF/s)", fontsize=12, fontweight="bold")
ax2.set_title("Panel 2 — cuSPARSE SpGEMM throughput on diagonal blocks\n"
              "each dot = one block; red = pathological slow regime",
              fontsize=12, fontweight="bold")
ax2.set_xscale("log", base=2)
ax2.set_xticks(nps)
ax2.set_xticklabels([f"np={n}\n({int(64/np.sqrt(n)) if n>1 else 64}k block)" for n in nps], fontsize=10)
ax2.set_ylim(8, 20)
ax2.legend(fontsize=9)
ax2.grid(True, lw=0.3, alpha=0.4)

fig.suptitle("NTPoly 2D SUMMA + symmetry — 64k Yukawa matrix on GB200", fontsize=13, fontweight="bold", y=1.02)
fig.tight_layout()
out = "timer_breakdown.png"
fig.savefig(out, dpi=150, bbox_inches="tight")
print(f"Saved {out}")

# ── Text summary ──────────────────────────────────────────────────────────
print("\n=== Wall-time decomposition (64k) ===")
print(f"{'np':>4} {'SpGEMM':>8} {'thresh':>8} {'resid':>7} {'wall':>7} {'ideal':>7} {'eff':>6}")
for n in nps:
    sp, th, wa = wall[n]
    print(f"{n:>4} {sp:>8.2f} {th:>8.2f} {wa-sp-th:>7.2f} {wa:>7.2f} {t1/n:>7.2f} {t1/(wa*n)*100:>5.0f}%")

print("\n=== cuSPARSE per-block throughput ===")
for n in nps:
    for comp, fl in zip(diag_blocks[n][0], diag_blocks[n][1]):
        print(f"np={n:>2}  block {fl:5.1f} Mflop  {comp:5.2f} ms  ->  {fl/comp:5.1f} GF/s")
