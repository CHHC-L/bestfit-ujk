"""Figures for the write-up (needs matplotlib): figures/classification.png and figures/boundary.png.

classification.png: the U{j,k} grid for k <= 25 from the verified results/classification.csv.
boundary.png: the linear-waste region in the (k, j) plane: boundary points K_j certified (rows j <= 12)
and simulated (results/boundary_sweep.csv, not a proof), the fit j* ~ 2.52 sqrt(k), and the proven
bounded region k >= j(j+3)/2.
usage: make_figures.py [--kmax 25]
"""
import argparse
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from make_tables import status  # noqa: E402

CLASSES = [("", "#ffffff", None),
           (".", "#eeeeee", "bounded, trivial (j ≤ 2)"),
           ("b", "#b9d3ee", "bounded, known before"),
           ("b*", "#7fa7d6", "bounded, known, re-proved here"),
           ("B", "#1f4e9c", "bounded, new"),
           ("l", "#f6c9a8", "linear, known before"),
           ("l*", "#f0a070", "linear, known (1993), re-proved here"),
           ("L", "#c0392b", "linear, new"),
           ("S", "#999999", "Θ(√n) (j = k−1)"),
           ("?", "#ffff66", "open")]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kmax", type=int, default=25)
    a = ap.parse_args()
    cls = {(int(r["j"]), int(r["k"])): r for r in csv.DictReader(open(os.path.join(ROOT, "results", "classification.csv")))}
    os.makedirs(os.path.join(ROOT, "figures"), exist_ok=True)

    # grid
    code = {c[0]: i for i, c in enumerate(CLASSES)}
    ks = list(range(3, a.kmax + 1))
    js = list(range(1, a.kmax))
    Z = [[code[status(cls, j, k) if j < k else ""] for j in js] for k in ks]
    fig, ax = plt.subplots(figsize=(12.5, 7.5))
    ax.imshow(Z, cmap=ListedColormap([c[1] for c in CLASSES]), vmin=0, vmax=len(CLASSES) - 1,
              origin="lower", aspect="equal")
    ax.set_xticks(range(len(js)), js, fontsize=8)
    ax.set_yticks(range(len(ks)), ks, fontsize=8)
    ax.set_xticks([x - 0.5 for x in range(len(js) + 1)], minor=True)
    ax.set_yticks([y - 0.5 for y in range(len(ks) + 1)], minor=True)
    ax.grid(which="minor", color="white", linewidth=0.6)
    ax.tick_params(which="minor", length=0)
    ax.set_xlabel("j (items uniform on {1..j})")
    ax.set_ylabel("k (bin capacity)")
    ax.set_title("Best Fit under U{j,k}: expected waste O(1) (blue) or Θ(n) (red)")
    used = {status(cls, j, k) for k in ks for j in js if j < k}
    ax.legend(handles=[Patch(color=c[1], label=c[2]) for c in CLASSES if c[2] and c[0] in used],
              loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=9, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(ROOT, "figures", "classification.png"), dpi=150)

    # boundary K_j
    sim = [(int(r["j"]), int(r["K_j"])) for r in csv.DictReader(open(os.path.join(ROOT, "results", "boundary_sweep.csv")))]
    cert = []
    for j in range(8, 13):
        lin = [k for (jj, k), r in cls.items() if jj == j and r["result"] == "linear"]
        if lin:
            cert.append((j, max(lin)))
    fig, ax = plt.subplots(figsize=(8, 5))
    kk = [x / 4 for x in range(40, 1041)]
    lower = [2.52 * k ** 0.5 - 0.31 for k in kk]
    ax.fill_between(kk, lower, [min(k - 3, 48) for k in kk], where=[k - 3 > l for k, l in zip(kk, lower)],
                    color="#c0392b", alpha=0.12, lw=0, label="linear waste (j* ≤ j ≤ k−3)")
    ax.plot(kk, lower, "-", color="#c0392b", lw=1, label="fit j* ≈ 2.52 √k − 0.31")
    ax.plot([K for _, K in sim], [j for j, _ in sim], "o", ms=4, color="#c0392b",
            label="boundary (largest linear k per j), simulation, not a proof")
    ax.plot([K for _, K in cert], [j for j, _ in cert], "s", ms=9, mfc="none", mec="black", mew=1.5,
            label="certified boundary (rows j ≤ 12)")
    ax.plot(kk, [(2 * k + 2.25) ** 0.5 - 1.5 for k in kk], "--", color="#1f4e9c", lw=1.2,
            label="proved bounded below this: j ≤ √(2k) − O(1) (1991)")
    ax.set_xlim(10, 260)
    ax.set_ylim(0, 45)
    ax.set_xlabel("k (bin capacity)")
    ax.set_ylabel("j (largest item size)")
    ax.set_title("Best Fit has linear waste iff 2.52√k ≲ j ≤ k−3 (empirically)")
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(ROOT, "figures", "boundary.png"), dpi=150)
    print("wrote figures/classification.png, figures/boundary.png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
