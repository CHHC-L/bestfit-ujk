"""Cutting-plane search over the linear part c, with value iteration (bfvi) for the correction psi.

linear : maximise  G(c) = min_policies c.d_pi   (c_1 = 1, |c_g| <= B)       need G > 0
bounded: minimise  H(c) = max_policies c.d_pi   (1 <= c_g <= B)              need H < 0
G is concave and H convex in c; bfvi returns a value and a subgradient d_pi.  The final (c, psi)
is rationalised into a certificate for `bfcert verify` and checked by both verifiers.

usage: kelley.py j k M mode [--reach] [--B 50] [--it 40] [--vi 3000] [--threads 8]
"""
import argparse
import os
import subprocess
import sys

import highspy
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def run_vi(a, c, prefix, iters, warm=False):
    cf = prefix + ".c.txt"
    with open(cf, "w") as f:
        f.write(" ".join(f"{v:.17g}" for v in c) + "\n")
    cmd = [os.path.join(ROOT, "bfvi"), str(a.j), str(a.k), a.M, a.mode, cf, prefix, str(iters),
           str(a.threads)] + (["reach"] if a.reach else []) + (["warm"] if warm else [])
    out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
    res = {}
    for line in out.splitlines():
        t = line.split()
        if t[0] in ("gain_lo", "gain_hi"):
            res[t[0]] = float(t[1])
        elif t[0] == "subgrad":
            res["d"] = np.array([float(v) for v in t[1:]])
        elif t[0] == "iters":
            res["iters"] = int(t[1])
    return res


def next_c(a, cuts, d0, center=None, radius=None):
    """Kelley master LP: variables c_1..c_d0, t.  With a center, the box is intersected with a
    trust region |c - center| <= radius (bundle-style stabilisation); without, the value is a valid
    global bound."""
    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    inf = highspy.kHighsInf
    lin = a.mode in ("linear", "ylo")
    for g in range(d0):
        if a.mode in ("ylo", "yhi"):   # c_1 = 1 counts gap-1 creations; others free
            lo, hi = (1.0, 1.0) if g == 0 else (-a.B, a.B)
        elif lin:
            lo, hi = (1.0, 1.0) if g == 0 else (-a.B, (g + 1.0) if a.rate else a.B)
        elif a.mode == "upper":   # V >= W - const  needs  c_g >= g
            lo, hi = g + 1.0, a.B
        elif a.mode in ("bounded2", "ybounded") and g == 0:   # gap 1 is left out of V
            lo, hi = 0.0, 0.0
        else:
            lo, hi = 1.0, a.B
        if center is not None and lo < hi:
            lo, hi = max(lo, center[g] - radius), min(hi, center[g] + radius)
        h.addVar(lo, hi)
    h.addVar(-inf, inf)
    t = d0
    h.changeColCost(t, -1.0 if lin else 1.0)
    for d in cuts:
        idx = list(range(d0)) + [t]
        if lin:   # t <= c.d  ->  c.d - t >= 0
            h.addRow(0.0, inf, d0 + 1, idx, list(d) + [-1.0])
        else:     # t >= c.d  ->  c.d - t <= 0
            h.addRow(-inf, 0.0, d0 + 1, idx, list(d) + [-1.0])
    if lin:
        h.changeColBounds(t, -1e3, 1e3)
    else:
        h.changeColBounds(t, -1e5, 1e5)
    h.run()
    x = np.array(h.getSolution().col_value)
    return x[:d0], x[t]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("j", type=int)
    ap.add_argument("k", type=int)
    ap.add_argument("M")
    ap.add_argument("mode", choices=["linear", "bounded", "upper", "bounded2", "ylo", "yhi", "ybounded"])
    ap.add_argument("--rate", action="store_true", help="linear mode: require c_g <= g so that the "
                    "certified drift is itself a lower bound on the waste rate")
    ap.add_argument("--reach", action="store_true")
    ap.add_argument("--B", type=float, default=50.0)
    ap.add_argument("--it", type=int, default=40)
    ap.add_argument("--vi", type=int, default=3000)
    ap.add_argument("--vifinal", type=int, default=20000)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--D", type=int, default=10**9)
    ap.add_argument("--reltol", type=float, default=0.05)
    ap.add_argument("--c0", default=None, help="comma separated initial c (small gaps)")
    ap.add_argument("--opt", action="store_true", help="optimise the margin instead of stopping at the "
                    "first certificate of the right sign")
    a = ap.parse_args()
    d0 = len(a.M.split(",")) if "," in a.M else (a.k - 1) // 2   # (k-1)//2, or j (see Chain)
    lin = a.mode in ("linear", "ylo")
    yrate = a.mode in ("ylo", "yhi")
    tag = f"{a.j}_{a.k}_{a.mode}{'_rate' if a.rate else ''}_vi_M{a.M.replace(',', '-')}" + ("_R" if a.reach else "")
    os.makedirs(os.path.join(ROOT, "lp"), exist_ok=True)
    prefix = os.path.join(ROOT, "lp", tag)
    if a.c0:
        c = np.array([float(v) for v in a.c0.split(",")])
    else:
        c = np.ones(d0)
        if lin:
            c[1:] = -0.5
        elif a.mode == "upper":
            c = np.arange(1, d0 + 1, dtype=float) * 1.5
        elif a.mode in ("bounded2", "ybounded"):
            c[0] = 0.0
        if yrate:
            c = np.zeros(d0)
            c[0] = 1.0
    cuts, best, bestc = [], None, None
    radius = 0.5 if a.c0 else a.B
    for it in range(a.it):
        r = run_vi(a, c, prefix, a.vi, warm=it > 0)
        val = r["gain_lo"] if lin else r["gain_hi"]
        improved = best is None or (val > best if lin else val < best)
        if improved:
            best, bestc = val, c.copy()
            radius = min(2 * radius, a.B) if it > 0 else radius
        else:
            radius = max(radius / 2, 1e-3)
        cuts.append(r["d"])
        _, bound = next_c(a, cuts, d0)                     # global bound (for stopping)
        c, _ = next_c(a, cuts, d0, center=bestc, radius=radius)
        print(f"it {it:3d} value {val:+.8f} best {best:+.8f} master {bound:+.8f} radius {radius:.3g} "
              f"vi_iters {r['iters']}", flush=True)
        if (lin and bound - best < 1e-6) or (not lin and best - bound < 1e-6):
            break
        if not a.opt and a.mode not in ("upper", "ylo", "yhi") and \
                ((lin and best > 1e-5) or (not lin and best < -1e-5)):
            break
        if a.mode not in ("upper", "ylo", "yhi") and ((lin and bound < -1e-4) or (not lin and bound > 1e-4)):
            break  # no certificate in this class
        if (lin and (best > 0 or yrate) and bound - best < a.reltol * abs(best)) or \
           (not lin and (best < 0 or a.mode in ("upper", "yhi")) and best - bound < a.reltol * abs(best)):
            break
    print("best c =", " ".join(f"{v:.6f}" for v in bestc))
    # final long value iteration with the best c, then write the certificate
    r = run_vi(a, bestc, prefix, a.vifinal, warm=True)
    print(f"final VI: gain in [{r['gain_lo']:.10f}, {r['gain_hi']:.10f}] after {r['iters']} iterations")
    psi = np.fromfile(prefix + ".psi", dtype=np.float64)
    D = a.D
    cfull = np.zeros(a.k - 1)
    cfull[:d0] = bestc
    cn = np.rint(cfull * D).astype(np.int64)
    pn = np.rint((psi - psi.min()) * D).astype(np.int64)
    cert = os.path.join(ROOT, "certs", tag + ".cert")
    with open(cert, "w") as f:
        f.write(f"mode {a.mode}\nD {D}\nc " + " ".join(map(str, cn)) + "\npsi\n")
        for i in range(0, len(pn), 64):
            f.write(" ".join(map(str, pn[i:i + 64])) + "\n")
    print("certificate:", cert)
    return 0


if __name__ == "__main__":
    sys.exit(main())
