"""Empirical (simulation-only) boundary of linear waste: for each j, the largest k with clearly
linear Best Fit waste under U{j,k}, located by bisection on k (assumes monotonicity in k).
Writes results/boundary_sweep.csv with columns j, K_j, j^2/K_j.

usage: boundary_sweep.py JMIN JMAX [--n 40000000] [--thr 2e-4]
"""
import argparse
import csv
import os
import subprocess

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def slope(j, k, n):
    out = subprocess.run([os.path.join(ROOT, "bfsim"), str(j), str(k), str(n), "2"], capture_output=True,
                         text=True, check=True).stdout
    return float(out.split("slope=")[1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("jmin", type=int)
    ap.add_argument("jmax", type=int)
    ap.add_argument("--n", type=int, default=40_000_000)
    ap.add_argument("--thr", type=float, default=2e-4)
    a = ap.parse_args()
    path = os.path.join(ROOT, "results", "boundary_sweep.csv")
    new = not os.path.exists(path)
    with open(path, "a", newline="") as f:
        w = csv.writer(f)
        if new:
            w.writerow(["j", "K_j", "j2_over_K", "slope_at_K", "slope_at_K+1"])
        for j in range(a.jmin, a.jmax + 1):
            lo, hi = j + 3, max(j + 4, int(j * j / 3.5))   # lo linear (if any), hi bounded
            if slope(j, lo, a.n) <= a.thr:
                w.writerow([j, "none", "", "", ""])
                f.flush()
                continue
            while slope(j, hi, a.n) > a.thr:
                hi *= 2
            while hi - lo > 1:
                mid = (lo + hi) // 2
                if slope(j, mid, a.n) > a.thr:
                    lo = mid
                else:
                    hi = mid
            s1, s2 = slope(j, lo, a.n), slope(j, hi, a.n)
            w.writerow([j, lo, round(j * j / lo, 3), s1, s2])
            f.flush()
            print(j, lo, round(j * j / lo, 3), s1, s2, flush=True)


if __name__ == "__main__":
    main()
