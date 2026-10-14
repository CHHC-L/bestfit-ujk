"""Batch classification of Best Fit under U{j,k}: search and verify drift certificates.

For each pair: simulation picks the expected regime and per-gap truncation levels (simulation only
guides the search; every claim rests on the exactly verified certificate).  We then run the
Kelley/value-iteration search and the exact C++ verifier, escalating the truncation on failure.
Attempts are appended to results/search/classify_log.csv; the verified table
results/classification.csv is produced by scripts/verify_all.py.

usage: classify.py --pairs 9:15,10:15 | --k 15 [--jmin 3 --jmax k-3] [--budget 3e6] [--threads 8]
"""
import argparse
import csv
import gzip
import os
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.makedirs(os.path.join(ROOT, "results", "search", "logs"), exist_ok=True)
CSV = os.path.join(ROOT, "results", "search", "classify_log.csv")  # search log (not the verified table)


def known(j, k):
    if j == k - 1:
        return "sqrt(n) [CCGJMSWY91]"
    if j == k - 2:
        return "bounded [KRS98]"
    if 2 * k >= j * (j + 3):
        return "bounded [CCGJMSWY91 Thm]"
    return None


def stats(j, k, n):
    out = subprocess.run([os.path.join(ROOT, "bfstats"), str(j), str(k), str(n)], capture_output=True,
                         text=True, check=True).stdout
    slope, tails = None, {}
    for line in out.splitlines():
        t = line.split()
        if t[0] == "slope":
            slope = float(t[1])
        elif t[0] == "tail":
            tails[int(t[1])] = [float(v) for v in t[2:]]
    return slope, tails


def choose_M(k, tails, eps, mode, mmax, d0=None):
    d0 = d0 if d0 else (k - 1) // 2
    M = []
    for g in range(1, d0 + 1):
        tg = tails[g]
        if mode == "linear" and g == 1:
            M.append(1)
            continue
        m = next((i for i, p in enumerate(tg) if p <= eps), len(tg))  # tg[i] = P(x >= i+1)
        M.append(max(1, min(mmax, m + 1)))
    return M


def count(j, k, M):
    out = subprocess.run([os.path.join(ROOT, "bfcert"), "count", str(j), str(k), ",".join(map(str, M)),
                          "reach"], capture_output=True, text=True, check=True).stdout.split()
    return int(out[1])


def attempt(j, k, M, mode, threads, vi, c0=None):
    Ms = ",".join(map(str, M))
    env = dict(os.environ, REACH="reach", SKIP_INDEP="1")
    t0 = time.time()
    out = subprocess.run([os.path.join(ROOT, "scripts", "run_vi_case.sh"), str(j), str(k), Ms, mode,
                          "--threads", str(threads), "--vi", str(vi), "--it", "60"] +
                         (["--c0", ",".join(f"{v:.8g}" for v in c0)] if c0 is not None else []),
                         capture_output=True, text=True, env=env).stdout
    ok = "VERDICT: CERTIFIED" in out
    bestc = None
    log = os.path.join(ROOT, "results", "search", "logs", f"{j}_{k}_{mode}_vi_M{'-'.join(map(str, M))}_R.log")
    if os.path.exists(log):
        for line in open(log):
            if line.startswith("best c ="):
                bestc = [float(v) for v in line.split("=")[1].split()]
    drift = None
    for line in out.splitlines():
        if "extreme j*D*drift" in line:
            drift = float(line.split("(drift = ")[1].split(")")[0])
    return ok, drift, time.time() - t0, out, bestc


def family_c0(j, k, mode, d0=None):
    """Linear part of the certified (j-1, k-1) case of the same mode, resized to d0 = (k-1)//2."""
    if not os.path.exists(CSV):
        return None
    for r in csv.DictReader(open(CSV)):
        if int(r["j"]) == j - 1 and int(r["k"]) == k - 1 and r["result"] == mode:
            cert = os.path.join(ROOT, "certs", f"{j - 1}_{k - 1}_{mode}_vi_M{r['M']}_R.cert.gz")
            if not os.path.exists(cert):
                continue
            tok = gzip.open(cert, "rt").read().split()
            D = int(tok[3])
            c = [int(v) / D for v in tok[5:5 + (k - 2)]]
            d0 = d0 if d0 else (k - 1) // 2
            c = (c + [0.0] * d0)[:d0]
            if mode == "bounded":
                c = [max(v, 1.0) for v in c]
            return c
    return None


def record(row):
    new = not os.path.exists(CSV)
    with open(CSV, "a", newline="") as f:
        w = csv.writer(f)
        if new:
            w.writerow(["j", "k", "result", "mode", "M", "patterns", "drift", "seconds", "slope_sim"])
        w.writerow(row)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", default=None)
    ap.add_argument("--k", type=int, default=None)
    ap.add_argument("--jmin", type=int, default=3)
    ap.add_argument("--jmax", type=int, default=None)
    ap.add_argument("--budget", type=float, default=3e6)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--vi", type=int, default=2500)
    ap.add_argument("--nsim", type=int, default=20_000_000)
    ap.add_argument("--mmax", type=int, default=12)
    ap.add_argument("--redo", action="store_true", help="also prove pairs settled by known theorems")
    ap.add_argument("--small-d0", action="store_true", help="use d0 = j small gaps when j < (k-1)/2")
    ap.add_argument("--eps", default="0.3,0.1,0.03,0.01,0.003,0.001,0.0003,0.0001",
                    help="occupation-tail thresholds tried in order (larger = coarser truncation)")
    a = ap.parse_args()
    a.eps_list = [float(v) for v in a.eps.split(",")]
    if a.pairs:
        pairs = [tuple(int(v) for v in p.split(":")) for p in a.pairs.split(",")]
    else:
        jmax = a.jmax if a.jmax else a.k - 2
        pairs = [(j, a.k) for j in range(a.jmin, jmax + 1)]
    for j, k in pairs:
        kn = known(j, k)
        if kn and not a.redo:
            print(f"U{{{j},{k}}}: {kn}", flush=True)
            record([j, k, kn, "", "", "", "", "", ""])
            continue
        if j == k - 1:
            print(f"U{{{j},{k}}}: {kn} (not a bounded/linear case)", flush=True)
            continue
        slope, tails = stats(j, k, a.nsim)
        modes = ["linear", "bounded"] if slope > 2e-4 else ["bounded", "linear"]
        if kn and "bounded" in kn:      # re-proving a known bounded result: only that mode
            modes = ["bounded"]
        if a.redo and k <= 14 and (j, k) in {(8, 11), (9, 12)}:
            modes = ["linear"]
        done = False
        for mode in modes:
            d0 = min(j, (k - 1) // 2) if a.small_d0 else (k - 1) // 2
            c0 = family_c0(j, k, mode, d0)
            prevM = None
            for eps in a.eps_list:
                M = choose_M(k, tails, eps, mode, a.mmax, d0)
                if M == prevM:
                    continue
                prevM = M
                n = count(j, k, M)
                if n > a.budget:
                    print(f"U{{{j},{k}}} {mode} eps={eps} M={M}: {n} patterns > budget", flush=True)
                    break
                ok, drift, sec, out, bestc = attempt(j, k, M, mode, a.threads, a.vi, c0)
                if bestc is not None:
                    c0 = bestc[:d0]
                print(f"U{{{j},{k}}} {mode} eps={eps} M={M} patterns={n}: "
                      f"{'CERTIFIED' if ok else 'failed'} drift={drift} ({sec:.0f}s)", flush=True)
                if ok:
                    record([j, k, mode, mode, "-".join(map(str, M)), n, drift, round(sec), slope])
                    done = True
                    break
            if done:
                break
        if not done:
            record([j, k, "OPEN", "", "", "", "", "", slope])
    return 0


if __name__ == "__main__":
    sys.exit(main())
