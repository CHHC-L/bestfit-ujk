"""End-to-end sanity check of a certificate on states that Best Fit actually visits.

Unlike the verifiers, this does not use truncation patterns or a reachability closure. It runs Best Fit
on an explicit vector of open-bin gaps for n random items starting from no bins, and at every visited
state x it evaluates the exact one-step drift
    DV(x) = (1/j) sum_s V(BF(x, s)) - V(x),   V(x) = sum_g c_g x_g + psi(min(x, M), big gap)
and checks the inequality that the certificate claims for x:
    linear   : DV(x) >= certified minimum (> 0)
    bounded  : DV(x) <= certified maximum (< 0) whenever some small gap has x_g > M_g
    bounded2 : the same with gap 1 excluded
    ylo/yhi  : at states with x_1 >= 1, the Y-chain drift (size-1 items are no-ops, gap-1 bins are
               counted when created) is >= / <= the certified bound (this is step (iii) of Theorem D).
It also checks the two invariants used by the chain (pairwise gaps < k; at most one gap > j).
The certified extreme drift is taken (exactly) from `bfcert verify`.

usage: e2e_check.py certs/<file>.cert.gz [...] [--steps 20000] [--seed 1]
"""
import argparse
import gzip
import os
import random
import re
import subprocess
import sys
from fractions import Fraction

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NAME = re.compile(r"^(\d+)_(\d+)_([a-z0-9]+)(?:_vi)?_M([\d-]+)_R\.cert\.gz$")


def best_fit(x, s, k):
    """Best Fit on gap counts x[0..k-1] (x[0] unused): fullest bin that fits, else a new bin."""
    y = list(x)
    for g in range(s, k):
        if y[g] > 0:
            y[g] -= 1
            if g > s:
                y[g - s] += 1
            return y
    y[k - s] += 1
    return y


def check(cert, steps, seed):
    name = os.path.basename(cert)
    m = NAME.match(name)
    j, k, mode, Mdash = int(m[1]), int(m[2]), m[3], m[4]
    Ms = [int(v) for v in Mdash.split("-")]
    d0 = len(Ms)
    M = [0] + Ms

    states = os.path.join(ROOT, "lp", "psi", f"{j}_{k}_M{Mdash}")
    if not os.path.exists(os.path.join(states, "psi_states.txt")):
        subprocess.run([os.path.join(ROOT, "bfcert"), "states", str(j), str(k), Mdash.replace("-", ","), states],
                       check=True, capture_output=True)
    psi_index = {}
    for line in open(os.path.join(states, "psi_states.txt")):
        left, big = line.split("|")
        parts = left.split()
        psi_index[(tuple(int(v) for v in parts[1:]), int(big))] = int(parts[0])
    tok = gzip.open(cert, "rt").read().split()
    D = int(tok[3])
    c = [0] + [int(v) for v in tok[5:5 + k - 1]]
    psi = [int(v) for v in tok[5 + k:]]

    # the certified extreme drift, exactly, as reported by the verifier
    out = subprocess.run([os.path.join(ROOT, "bfcert"), "verify", str(j), str(k), Mdash.replace("-", ","),
                          cert, "reach"], capture_output=True, text=True).stdout
    if "VERDICT: CERTIFIED" not in out:
        sys.exit(f"{name}: bfcert verify does not certify it")
    bound = Fraction(int(re.search(r"extreme j\*D\*drift = (-?\d+)", out)[1]), j * D)

    ymode = mode in ("ylo", "yhi")

    def V(x, created=0):
        q = tuple(0 if (ymode and g == 1) else min(x[g], M[g]) for g in range(1, d0 + 1))
        bigs = [g for g in range(d0 + 1, k) if x[g] > 0]
        b = bigs[0] if bigs else 0
        # Y-chain: gap-1 bins enter V only through c_1 * (number created); X-chain: through c_1 * x_1
        lin = sum(c[g] * x[g] for g in range(2 if ymode else 1, k)) + (c[1] * created if ymode else 0)
        return lin + psi[psi_index[(q, b)]]

    def drift(x):
        v0 = V(x)
        tot = 0
        for s in range(1, j + 1):
            if ymode and s == 1:
                continue  # Y-chain: a size-1 item only fills a gap-1 bin
            y = best_fit(x, s, k)
            created = 1 if y[1] > x[1] else 0
            tot += V(y, created) - v0
        return Fraction(tot, j * D)

    rng = random.Random(seed)
    x = [0] * k
    worst = None
    checked = 0
    for t in range(steps):
        # invariants of the chain
        gaps = sorted((g for g in range(1, k) for _ in range(min(x[g], 2))), reverse=True)
        assert len(gaps) < 2 or gaps[0] + gaps[1] < k, ("pairwise invariant", t, x)
        assert sum(x[g] for g in range(j + 1, k)) <= 1, ("at most one gap > j", t, x)
        sat = [g for g in range(1, d0 + 1) if x[g] > M[g]]
        test = {"linear": True, "bounded": bool(sat), "bounded2": any(g >= 2 for g in sat),
                "ylo": x[1] >= 1, "yhi": x[1] >= 1}[mode]
        if test:
            d = drift(x)
            checked += 1
            if mode in ("linear", "ylo"):
                ok = d >= bound
                worst = d if worst is None or d < worst else worst
            else:
                ok = d <= bound
                worst = d if worst is None or d > worst else worst
            if not ok:
                print(f"FAIL {name}: step {t} state {x[1:]} drift {float(d)} vs certified {float(bound)}")
                return 1
        x = best_fit(x, rng.randint(1, j), k)
    waste = sum(g * x[g] for g in range(1, k))
    print(f"e2e {name}: {steps} steps, {checked} states checked, extreme visited drift "
          f"{float(worst) if worst is not None else float('nan'):.6g} vs certified {float(bound):.6g}; "
          f"final waste {waste}, x_1 = {x[1]}  OK")
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("certs", nargs="+")
    ap.add_argument("--steps", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    return max(check(c, a.steps, a.seed) for c in a.certs)


if __name__ == "__main__":
    sys.exit(main())
