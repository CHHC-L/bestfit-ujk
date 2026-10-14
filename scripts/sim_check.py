"""Simulation sanity check of the verified classification (not part of any proof).

For every pair in results/classification.csv, simulate Best Fit (bfsim, n items, 2 runs) and compare
the waste growth per item over the second half with the certificate:
  linear : Theorem A gives liminf E W_n / n >= gamma / C (C = max_g c_g^+ / g over small gaps); the
           simulated slope must not be clearly below this certified lower bound;
  bounded: the simulated slope must be ~ 0.
Writes results/sim_check.csv.   usage: sim_check.py [--n 20000000] [--jobs 4]
"""
import argparse
import csv
import gzip
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def certified_rate(r):
    """gamma / C from the certificate (linear mode)."""
    j, k = int(r["j"]), int(r["k"])
    d0 = len(r["M"].split("-"))
    tok = gzip.open(os.path.join(ROOT, "certs", r["certificate"]), "rt").read().split()
    c = [int(v) for v in tok[5:5 + k - 1]]
    C = max(Fraction(max(c[g - 1], 0), g) for g in range(1, d0 + 1)) / int(tok[3])
    return Fraction(r["drift"]) / C


def simulate(j, k, n):
    out = subprocess.run([os.path.join(ROOT, "bfsim"), str(j), str(k), str(n), "2"], capture_output=True,
                         text=True, check=True).stdout
    return float(out.split("slope=")[1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=20_000_000)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--noise", type=float, default=5e-5, help="allowed |slope| for bounded / shortfall for linear")
    a = ap.parse_args()
    rows = list(csv.DictReader(open(os.path.join(ROOT, "results", "classification.csv"))))
    with ThreadPoolExecutor(a.jobs) as ex:
        slopes = list(ex.map(lambda r: simulate(int(r["j"]), int(r["k"]), a.n), rows))
    bad = []
    with open(os.path.join(ROOT, "results", "sim_check.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["j", "k", "result", "certified_rate_lower", "simulated_slope", "consistent"])
        for r, s in zip(rows, slopes):
            if r["result"] == "linear":
                lo = float(certified_rate(r))
                ok = s >= lo - a.noise and s > a.noise
            else:
                lo = ""
                ok = abs(s) <= a.noise
            w.writerow([r["j"], r["k"], r["result"], f"{lo:.6g}" if lo != "" else "", f"{s:.6g}", ok])
            if not ok:
                bad.append((r["j"], r["k"], r["result"], lo, s))
    lin = [s for r, s in zip(rows, slopes) if r["result"] == "linear"]
    bnd = [s for r, s in zip(rows, slopes) if r["result"] == "bounded"]
    print(f"{len(rows)} pairs simulated (n = {a.n:.0e}, 2 runs): linear slopes in [{min(lin):.2e}, {max(lin):.2e}], "
          f"bounded |slope| <= {max(abs(s) for s in bnd):.1e}; inconsistent: {len(bad)}")
    for b in bad:
        print("INCONSISTENT", b)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
