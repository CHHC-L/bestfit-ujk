"""For every certified linear case in results/classification.csv, search a "bounded2" certificate
(Theorem C: all bins with gap >= 2 stay O(1), so the linear waste is carried by gap-1 bins).

Starts from the truncation of the linear certificate and enlarges it on failure.
Attempts are appended to results/search/structure_log.csv.   usage: structure_all.py [--kmax K] [--threads T]
"""
import argparse
import csv
import gzip
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.makedirs(os.path.join(ROOT, "results", "search", "logs"), exist_ok=True)
OUT = os.path.join(ROOT, "results", "search", "structure_log.csv")  # search log


def family_c0(j, k):
    """bounded2 linear part of the (j-1, k-1) structure certificate, resized; default (0,1,...,1)."""
    d0 = (k - 1) // 2
    c = [0.0] + [1.0] * (d0 - 1)
    if os.path.exists(OUT):
        for r in csv.DictReader(open(OUT)):
            if int(r["j"]) == j - 1 and int(r["k"]) == k - 1 and r["result"] == "bounded2":
                cert = os.path.join(ROOT, "certs", f"{j - 1}_{k - 1}_bounded2_vi_M{r['M']}_R.cert.gz")
                if os.path.exists(cert):
                    tok = gzip.open(cert, "rt").read().split()
                    D = int(tok[3])
                    old = [int(v) / D for v in tok[5:5 + (k - 2)]]
                    c = [0.0] + [max(1.0, v) for v in (old + [1.0] * d0)[1:d0]]
    return c


def run(j, k, M, threads):
    Ms = ",".join(map(str, M))
    c0 = family_c0(j, k)
    env = dict(os.environ, REACH="reach", SKIP_INDEP="1")
    out = subprocess.run([os.path.join(ROOT, "scripts", "run_vi_case.sh"), str(j), str(k), Ms, "bounded2",
                          "--threads", str(threads), "--vi", "2500", "--it", "40",
                          "--c0", ",".join(f"{v:.8g}" for v in c0)],
                         capture_output=True, text=True, env=env).stdout
    drift = None
    for line in out.splitlines():
        if "(drift = " in line:
            drift = float(line.split("(drift = ")[1].split(")")[0])
    return "VERDICT: CERTIFIED" in out, drift


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kmax", type=int, default=99)
    ap.add_argument("--threads", type=int, default=9)
    a = ap.parse_args()
    done = set()
    if os.path.exists(OUT):
        done = {(int(r[0]), int(r[1])) for r in csv.reader(open(OUT)) if r and r[0] != "j"}
    rows = list(csv.DictReader(open(os.path.join(ROOT, "results", "classification.csv"))))
    seen = set()
    for r in rows:
        j, k = int(r["j"]), int(r["k"])
        if r["result"] != "linear" or (j, k) in done or (j, k) in seen or k > a.kmax:
            continue
        seen.add((j, k))
        M = [int(v) for v in r["M"].split("-")]
        M[0] = max(M[0], 2)
        ok, drift = False, None
        for extra in (-2, 0, 1, 2):
            Mt = [M[0]] + [max(1, m + extra) for m in M[1:]]
            ok, drift = run(j, k, Mt, a.threads)
            print(f"U{{{j},{k}}} bounded2 M={Mt}: {'CERTIFIED' if ok else 'failed'} drift={drift}", flush=True)
            if ok:
                break
        new = not os.path.exists(OUT)
        with open(OUT, "a", newline="") as f:
            w = csv.writer(f)
            if new:
                w.writerow(["j", "k", "result", "M", "drift"])
            w.writerow([j, k, "bounded2" if ok else "OPEN", "-".join(map(str, Mt)), drift])
    return 0


if __name__ == "__main__":
    sys.exit(main())
