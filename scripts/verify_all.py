"""Verify every certificate in certs/ and derive the result tables from the verified certificates.

For each certs/*.cert.gz (file name {j}_{k}_{mode}[_vi]_M{M_1-...-M_d0}[_R].cert.gz):
  1. `bfcert states`  writes the psi index table (lp/psi/...);
  2. `bfcert verify`  exact check (__int128);
  3. `verify2`        independent exact check (own Best Fit, own reachability BFS);
  4. both must say CERTIFIED and agree on the size of the reachable closure and on the extreme drift
     (an exact integer j*D*drift);
  5. ylo/yhi (Theorem D) additionally need `bfcert ycover`;
  6. optionally (--python N) scripts/verify_indep.py on certificates with at most N patterns.
Writes results/verification.csv (one row per certificate) and, from the passing certificates only,
results/classification.csv (one row per (j,k)), results/structure.csv and results/rates.csv.
Fails (exit 1) on any verification failure or if some pair has both a linear and a bounded certificate.

usage: python3 scripts/verify_all.py [--jobs 8] [--python 300000] [--only 9_13] [--max-bytes 1000000]
"""
import argparse
import csv
import glob
import os
import re
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NAME = re.compile(r"^(\d+)_(\d+)_([a-z0-9]+)(?:_vi)?_M([\d-]+)(_R)?\.cert\.gz$")
BIG_BYTES = 1_000_000          # certificates this large are verified one at a time (verify2 memory)
big_lock = threading.Lock()


def run(cmd):
    p = subprocess.run(cmd, capture_output=True, text=True, cwd=ROOT)
    return p.returncode, p.stdout + p.stderr


def check(path, python_max):
    name = os.path.basename(path)
    m = NAME.match(name)
    if not m:
        return {"cert": name, "status": "BAD NAME"}
    j, k, mode, Mdash, reach = int(m[1]), int(m[2]), m[3], m[4], m[5]
    M = Mdash.replace("-", ",")
    row = {"cert": name, "j": j, "k": k, "mode": mode, "M": Mdash}
    if not reach:
        row["status"] = "NOT REACH"   # every certificate here is a reach certificate
        return row
    t0 = time.time()
    states = os.path.join("lp", "psi", name[:-len(".cert.gz")])  # per certificate: no sharing between jobs
    rc, out = run(["./bfcert", "states", str(j), str(k), M, states])
    if rc:
        row["status"] = "STATES FAILED: " + out.strip()[-200:]
        return row
    rc1, out1 = run(["./bfcert", "verify", str(j), str(k), M, os.path.join("certs", name), "reach"])
    rc2, out2 = run(["./verify2", str(j), str(k), M, states, os.path.join("certs", name), "reach"])
    v1 = "VERDICT: CERTIFIED" in out1 and rc1 == 0
    v2 = "VERDICT: CERTIFIED" in out2 and rc2 == 0
    a = re.search(r"patterns=(\d+)", out1)
    b = re.search(r"verify2: mode=\S+ patterns=(\d+) checked=(\d+) extreme j\*D\*drift=(-?\d+) \(D=(\d+)", out2)
    e1 = re.search(r"extreme j\*D\*drift = (-?\d+)", out1)
    row["bfcert"] = "CERTIFIED" if v1 else "FAILED"
    row["verify2"] = "CERTIFIED" if v2 else "FAILED"
    agree = bool(a and b and e1) and a[1] == b[1] and e1[1] == b[3]
    row["agree"] = agree
    if b:
        row["patterns"], row["checked"], num, D = int(b[1]), int(b[2]), int(b[3]), int(b[4])
        row["drift"] = Fraction(num, j * D)
    ok = v1 and v2 and agree
    if mode in ("ylo", "yhi"):
        rc3, out3 = run(["./bfcert", "ycover", str(j), str(k), M])
        row["ycover"] = "COVERED" if rc3 == 0 and "VERDICT: COVERED" in out3 else "NOT COVERED"
        ok = ok and row["ycover"] == "COVERED"
    if python_max and row.get("patterns", 1 << 60) <= python_max:
        rc4, out4 = run([sys.executable, "scripts/verify_indep.py", str(j), str(k), M, states,
                         os.path.join("certs", name), "reach"])
        row["python"] = "CERTIFIED" if rc4 == 0 and "VERDICT: CERTIFIED" in out4 else "FAILED"
        ok = ok and row["python"] == "CERTIFIED"
    row["seconds"] = round(time.time() - t0, 1)
    row["status"] = "OK" if ok else "FAILED"
    return row


def check_sched(path, python_max):
    if os.path.getsize(path) >= BIG_BYTES:
        with big_lock:
            return check(path, python_max)
    return check(path, python_max)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    ap.add_argument("--python", type=int, default=0, help="also run verify_indep.py up to this many patterns")
    ap.add_argument("--only", default="", help="restrict to certificates whose name starts with this")
    ap.add_argument("--max-bytes", type=int, default=0, help="restrict to certificate files smaller than this")
    a = ap.parse_args()
    for exe in ("bfcert", "verify2"):
        if not os.path.exists(os.path.join(ROOT, exe)):
            sys.exit(f"missing ./{exe}: run `make` first")
    paths = sorted(glob.glob(os.path.join(ROOT, "certs", a.only + "*.cert.gz")),
                   key=lambda p: -os.path.getsize(p))
    if a.max_bytes:
        paths = [p for p in paths if os.path.getsize(p) < a.max_bytes]
    print(f"verifying {len(paths)} certificates with {a.jobs} jobs", flush=True)
    rows = []
    with ThreadPoolExecutor(a.jobs) as ex:
        for r in ex.map(lambda p: check_sched(p, a.python), paths):
            rows.append(r)
            drift = f"{float(r['drift']):.6g}" if "drift" in r else "-"
            print(f"[{len(rows)}/{len(paths)}] {r['status']:6s} {r['cert']}  "
                  f"patterns={r.get('patterns', '-')} drift={drift} {r.get('seconds', '')}s"
                  + (f" python={r['python']}" if "python" in r else ""), flush=True)
    rows.sort(key=lambda r: (r.get("k", 0), r.get("j", 0), r["cert"]))
    if a.only or a.max_bytes:  # partial run: report only, do not rewrite the tables
        bad = [r for r in rows if r["status"] != "OK"]
        if a.python:
            print(f"also checked by verify_indep.py: {sum(r.get('python') == 'CERTIFIED' for r in rows)}")
        print("FAILED:" if bad else "all OK", *[r["cert"] for r in bad])
        return 1 if bad else 0

    cols = ["cert", "j", "k", "mode", "M", "patterns", "checked", "drift", "bfcert", "verify2", "agree",
            "ycover", "python", "seconds", "status"]
    with open(os.path.join(ROOT, "results", "verification.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, cols, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({**r, "drift": f"{float(r['drift']):.12g}" if "drift" in r else ""})

    good = [r for r in rows if r["status"] == "OK"]
    bad = [r for r in rows if r["status"] != "OK"]

    # classification: per pair, the smallest passing linear/bounded certificate
    cls, conflicts = {}, []
    for r in good:
        if r["mode"] not in ("linear", "bounded"):
            continue
        key = (r["j"], r["k"])
        if key in cls and cls[key]["mode"] != r["mode"]:
            conflicts.append(key)
        if key not in cls or r["patterns"] < cls[key]["patterns"]:
            cls[key] = r
    with open(os.path.join(ROOT, "results", "classification.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["j", "k", "result", "certificate", "M", "patterns", "drift"])
        for (j, k) in sorted(cls, key=lambda t: (t[1], t[0])):
            r = cls[(j, k)]
            w.writerow([j, k, r["mode"], r["cert"], r["M"], r["patterns"], f"{float(r['drift']):.12g}"])

    # Theorem C: per pair, a passing bounded2 certificate
    struct = {}
    for r in good:
        if r["mode"] == "bounded2" and ((r["j"], r["k"]) not in struct or
                                        r["patterns"] < struct[(r["j"], r["k"])]["patterns"]):
            struct[(r["j"], r["k"])] = r
    with open(os.path.join(ROOT, "results", "structure.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["j", "k", "certificate", "M", "patterns", "drift", "linear_certified"])
        for (j, k) in sorted(struct, key=lambda t: (t[1], t[0])):
            r = struct[(j, k)]
            lin = cls.get((j, k), {}).get("mode") == "linear"
            w.writerow([j, k, r["cert"], r["M"], r["patterns"], f"{float(r['drift']):.12g}", lin])

    # Theorem D: best lower (ylo) and upper (yhi) bracket per pair; needs linear + bounded2 too
    lo, hi = {}, {}
    for r in good:
        key = (r["j"], r["k"])
        if r["mode"] == "ylo" and (key not in lo or r["drift"] > lo[key]["drift"]):
            lo[key] = r
        if r["mode"] == "yhi" and (key not in hi or r["drift"] < hi[key]["drift"]):
            hi[key] = r
    with open(os.path.join(ROOT, "results", "rates.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["j", "k", "creation_lo", "creation_hi", "rate_lo", "rate_hi", "ylo_cert", "yhi_cert",
                    "theorem_D_applies"])
        for key in sorted(set(lo) & set(hi), key=lambda t: (t[1], t[0])):
            j, k = key
            applies = cls.get(key, {}).get("mode") == "linear" and key in struct
            l, h = lo[key]["drift"], hi[key]["drift"]
            w.writerow([j, k, f"{float(l):.12g}", f"{float(h):.12g}", f"{float(l - Fraction(1, j)):.12g}",
                        f"{float(h - Fraction(1, j)):.12g}", lo[key]["cert"], hi[key]["cert"], applies])
            if not applies:
                conflicts.append(("theorem D prerequisites missing", key))

    print(f"\n{len(rows)} certificates: {len(good)} OK, {len(bad)} failed")
    print(f"classified pairs: {len(cls)} ({sum(r['mode'] == 'linear' for r in cls.values())} linear, "
          f"{sum(r['mode'] == 'bounded' for r in cls.values())} bounded); structure (C): {len(struct)}; "
          f"rate brackets (D): {len(set(lo) & set(hi))}")
    if a.python:
        print(f"also checked by verify_indep.py: {sum(r.get('python') == 'CERTIFIED' for r in rows)}")
    for r in bad:
        print("FAILED:", r["cert"], r.get("status"))
    for c in conflicts:
        print("CONFLICT:", c)
    return 1 if bad or conflicts else 0


if __name__ == "__main__":
    sys.exit(main())
