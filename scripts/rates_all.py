"""Certified waste-rate brackets (Theorem D) via the decoupled Y-chain.

For each pair: ylo / yhi certificates bounding the gap-1 creation rate of the Y-chain; with a linear and
a bounded2 certificate and `bfcert ycover`, liminf/limsup of the waste rate lie in [lo - 1/j, hi - 1/j]
(notes/method.md, Theorem D).  Attempts are appended to results/search/rates_log.csv.
usage: rates_all.py j:k:M [j:k:M ...] [--threads T]
(M given with dashes for gaps 2..d0, e.g. 9:13:5-5-4-4-4; M_1 = 0 is prepended)
"""
import argparse
import csv
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.makedirs(os.path.join(ROOT, "results", "search", "logs"), exist_ok=True)
OUT = os.path.join(ROOT, "results", "search", "rates_log.csv")


def run(j, k, M, mode, threads, extra):
    env = dict(os.environ, REACH="reach", SKIP_INDEP="1")
    out = subprocess.run([os.path.join(ROOT, "scripts", "run_vi_case.sh"), str(j), str(k), M, mode,
                          "--threads", str(threads)] + extra, capture_output=True, text=True, env=env).stdout
    val = None
    log = os.path.join(ROOT, "results", "search", "logs", f"{j}_{k}_{mode}_vi_M{M.replace(',', '-')}_R.log")
    for line in open(log):
        if line.startswith("creation_lower") or line.startswith("creation_upper"):
            val = float(line.split("=")[1].split()[0])
    return "VERDICT: CERTIFIED" in out, val


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("specs", nargs="+")
    ap.add_argument("--threads", type=int, default=1)
    a = ap.parse_args()
    for spec in a.specs:
        j, k, Mrest = spec.split(":")
        j, k = int(j), int(k)
        M = "0," + Mrest.replace("-", ",")
        okb = subprocess.run([os.path.join(ROOT, "bfcert"), "ycover", str(j), str(k), M]).returncode == 0
        okl, lo = run(j, k, M, "ylo", a.threads, ["--opt", "--reltol", "0.0005", "--vi", "6000"])
        okh, hi = run(j, k, M, "yhi", a.threads, ["--opt", "--reltol", "0.0005", "--vi", "6000"])
        ok = okb and okl and okh
        row = [j, k, M.replace(",", "-"), okb, lo, hi,
               (lo - 1 / j) if lo is not None else None, (hi - 1 / j) if hi is not None else None]
        print(row, flush=True)
        new = not os.path.exists(OUT)
        with open(OUT, "a", newline="") as f:
            w = csv.writer(f)
            if new:
                w.writerow(["j", "k", "M", "ycover", "pi_a_lo", "pi_a_hi", "rate_lo", "rate_hi"])
            w.writerow(row)
    return 0


if __name__ == "__main__":
    sys.exit(main())
