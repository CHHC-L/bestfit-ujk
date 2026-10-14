"""Independent exact verifier for Best Fit U{j,k} drift certificates (no shared code with bfcert).

For every truncation pattern it builds concrete gap-count vectors x (several representatives per
saturated coordinate), applies Best Fit to the explicit counts for every item size, and evaluates
V(x) = sum_g c_g x_g + psi(min(x, M)) with exact rationals.  psi is looked up through
psi_states.txt (the index -> truncated-state table written by `bfcert gen`), so only that
bijection is shared with the C++ tool.

usage: verify_indep.py j k M_1,..,M_d0 lpdir certfile [reach]
With "reach", the pattern set is the closure of the empty state under abstract Best Fit moves
(computed here from explicit Best Fit on representative states); otherwise all patterns whose
gaps satisfy the pairwise invariant g + h < k.
"""
import gzip
import itertools
import sys
from fractions import Fraction


def best_fit(x, s, k):
    """x: dict gap -> count (positive entries only).  Returns new dict."""
    y = dict(x)
    fits = [g for g, n in y.items() if n > 0 and g >= s]
    if fits:
        g = min(fits)
        y[g] -= 1
        if y[g] == 0:
            del y[g]
        if g - s > 0:
            y[g - s] = y.get(g - s, 0) + 1
    else:
        y[k - s] = y.get(k - s, 0) + 1
    return y


def reachable_superset(x, k):
    gaps = sorted(g for g, n in x.items() for _ in range(min(n, 2)))
    return len(gaps) < 2 or gaps[-1] + gaps[-2] < k


def main():
    j, k = int(sys.argv[1]), int(sys.argv[2])
    d0 = (k - 1) // 2
    Ms = [int(t) for t in sys.argv[3].split(",")]
    if len(Ms) == 1:
        Ms = Ms * d0
    elif len(Ms) != d0:
        # d0 = j < (k-1)/2 is allowed: at most one open bin has a gap > j (checked below via V())
        assert len(Ms) == j and j < d0, "invalid number of truncation levels"
        d0 = j
    M = {g: Ms[g - 1] for g in range(1, d0 + 1)}
    lpdir, certfile = sys.argv[4], sys.argv[5]

    psi_index = {}
    with open(f"{lpdir}/psi_states.txt") as f:
        for line in f:
            left, big = line.split("|")
            parts = left.split()
            psi_index[(tuple(int(v) for v in parts[1:]), int(big))] = int(parts[0])

    opener = gzip.open if certfile.endswith(".gz") else open
    with opener(certfile, "rt") as fh:
        tok = fh.read().split()
    assert tok[0] == "mode"
    mode = tok[1]
    D = int(tok[3])
    assert tok[4] == "c"
    c = {g: int(tok[5 + g - 1]) for g in range(1, k)}
    assert tok[5 + k - 1] == "psi"
    psi = [int(v) for v in tok[5 + k:]]
    assert len(psi) == len(psi_index), (len(psi), len(psi_index))

    def V(x):
        small = tuple(min(x.get(g, 0), M[g]) for g in range(1, d0 + 1))
        bigs = [g for g in range(d0 + 1, k) if x.get(g, 0) > 0]
        assert len(bigs) <= 1 and all(x[g] == 1 for g in bigs)
        b = bigs[0] if bigs else 0
        # exact integer numerator of V(x) (V = numerator / D)
        return sum(c[g] * n for g, n in x.items()) + psi[psi_index[(small, b)]]

    big_opts = [0] + list(range(d0 + 1, k))

    def pattern_of(x):
        p = tuple(min(x.get(g, 0), M[g] + 1) for g in range(1, d0 + 1))
        bigs = [g for g in range(d0 + 1, k) if x.get(g, 0) > 0]
        return p, (bigs[0] if bigs else 0)

    def all_patterns():
        for b in big_opts:
            for p in itertools.product(*[range(M[g] + 2) for g in range(1, d0 + 1)]):
                yield p, b

    def reach_patterns():
        start = (tuple([1 if (ymode and g == 1) else 0 for g in range(1, d0 + 1)]), 0)
        seen = {start}
        stack = [start]
        while stack:
            p, b = stack.pop()
            yield p, b
            x = {g: p[g - 1] for g in range(1, d0 + 1) if p[g - 1] > 0}
            if b:
                x[b] = 1
            for s in range(1, j + 1):
                y = step(x, s, k)
                succ = [pattern_of(y)]
                # a saturated level that was decremented may still be saturated
                for g in range(1, d0 + 1):
                    if p[g - 1] == M[g] + 1 and y.get(g, 0) < x.get(g, 0):
                        q = list(succ[0][0])
                        q[g - 1] = M[g] + 1
                        succ.append((tuple(q), succ[0][1]))
                for nxt in succ:
                    if nxt not in seen:
                        seen.add(nxt)
                        stack.append(nxt)

    reach = len(sys.argv) > 6 and sys.argv[6] == "reach"
    ymode = mode in ("ylo", "yhi", "ybounded")
    if ymode:
        assert M[1] == 0, "Y-chain certificates need M_1 = 0"
        base_fit = best_fit

        def step(x, s, k):  # gap-1 bins form an infinite reservoir: size-1 items change nothing
            return dict(x) if s == 1 else base_fit(x, s, k)
    else:
        step = best_fit
    extreme = None
    count = 0
    for p, b in (reach_patterns() if reach else all_patterns()):
        if True:
            x0 = {g: p[g - 1] for g in range(1, d0 + 1) if p[g - 1] > 0}
            if b:
                x0[b] = 1
            if not reachable_superset(x0, k):
                if reach:
                    print("abstractly reachable pattern violates pairwise invariant", p, b)
                    return 1
                continue
            sat = [g for g in range(1, d0 + 1) if p[g - 1] == M[g] + 1]
            if mode in ("bounded2", "ybounded"):  # only gaps >= 2 define the exceptional set
                sat_rest = [g for g in sat if g >= 2]
                if not sat_rest:
                    continue
            if mode == "bounded" and not sat:  # linear and upper check every pattern
                continue
            count += 1
            drifts = set()
            # representatives: saturated coordinates at M+1, M+2, M+7 (independently for the first)
            reps = [0, 1, 6]
            for shift in reps:
                for which in ([None] if not sat else sat[:1] + [None]):
                    x = dict(x0)
                    for g in sat:
                        if which is None or g == which:
                            x[g] = M[g] + 1 + shift
                    vx = V(x)
                    tot = 0  # j * D * drift, exact integer
                    for s in range(1, j + 1):
                        tot += V(step(x, s, k)) - vx
                    drifts.add(tot)
            if len(drifts) != 1:
                print("drift not constant on pattern", p, b, drifts)
                return 1
            d = drifts.pop()
            if extreme is None or (d < extreme if mode in ("linear", "ylo") else d > extreme):
                extreme = d
    extreme = Fraction(extreme, j * D)
    if mode == "linear":
        ok = extreme > 0
    elif mode == "upper":
        ok = all(c[g] >= g * D for g in range(1, d0 + 1))
    elif mode in ("ylo", "yhi"):
        ok = c[1] == D
        print(f"{'creation_lower' if mode == 'ylo' else 'creation_upper'} = {float(extreme):.12g}; "
              f"waste rate {'>=' if mode == 'ylo' else '<='} {float(extreme - Fraction(1, j)):.12g}")
    elif mode in ("bounded2", "ybounded"):
        ok = extreme < 0 and c[1] == 0 and all(c[g] > 0 for g in range(2, d0 + 1))
    else:
        ok = extreme < 0 and all(c[g] > 0 for g in range(1, d0 + 1))
    print(f"independent check: mode={mode} patterns={count} extreme drift={float(extreme):.12g} "
          f"({extreme})")
    print("VERDICT:", "CERTIFIED" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
