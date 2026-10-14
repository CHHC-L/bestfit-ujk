// Best Fit bin packing under U{j,k}: gap-count Markov chain and truncation patterns.
//
// State x = (x_1..x_{k-1}), x_g = number of open bins with gap g.  Item s ~ U{1..j}.
// Best Fit: g* = min{g >= s : x_g > 0}; if it exists x_{g*}--, x_{g*-s}++ (if g* > s),
// otherwise a new bin is opened: x_{k-s}++.
//
// Invariant (Coffman-Johnson-Shor-Weber 1993): any two open bins have gaps g,h with g+h < k.
// Hence "small" gaps g with 2g < k (g <= d0 = (k-1)/2) can have any multiplicity, while at most
// one bin with a "big" gap G (2G >= k) exists, and it excludes small gaps h >= k-G.
//
// A pattern is p = (p_1..p_d0, b) with p_g in {0..M_g+1} (M_g+1 means ">= M_g+1") and
// b in {0} U {d0+1..k-1} (the big bin's gap, 0 = none).  The one-step drift of
// V(x) = sum_g c_g x_g + psi(min(x,M)) depends on x only through its pattern.
#pragma once
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <unordered_set>

struct Move {            // effect of one item on a pattern
    int dec = 0;         // gap whose count decreases (0 = none)
    int inc = 0;         // gap whose count increases (0 = none)
};

struct Chain {
    int j, k, d0;
    std::vector<int> M;          // M[g] for g = 1..d0 (index 0 unused)
    std::vector<int64_t> pstride, qstride;  // mixed radix strides for patterns / psi states
    int64_t nsmallP = 1, nsmallQ = 1;       // number of small-part combinations
    int nb;                      // number of big options = 1 + (k-1-d0)
    std::vector<int64_t> qmap;   // dense psi index -> compressed (-1 if incompatible)
    int64_t nq = 0;              // number of compressed psi states
    int satFrom = 1;             // saturation is checked on gaps satFrom..d0 ("bounded2": 2)
    // Decoupled "Y-chain": gap-1 bins are an infinite reservoir, size-1 items are no-ops, and
    // gap-1 bins are only counted when created (requires M_1 = 0; the BFS starts with p_1 = 1).
    bool ychain = false;

    Chain(int j_, int k_, const std::vector<int>& Mg) : j(j_), k(k_) {
        // d0 = number of "small" gaps with unbounded counts.  Default (k-1)/2 (pairwise invariant);
        // when j < (k-1)/2 one may take d0 = j, since at most one open bin has a gap > j (a newer
        // bin is opened by an item of size <= j, which would have fitted into such a bin).
        d0 = (int)Mg.size() - 1;
        if (!(d0 == (k - 1) / 2 || (d0 == j && j < (k - 1) / 2))) {
            fprintf(stderr, "invalid number of truncation levels %d (need %d or j=%d)\n", d0, (k - 1) / 2, j);
            exit(1);
        }
        M.assign(d0 + 1, 0);
        for (int g = 1; g <= d0; ++g) M[g] = Mg[g];
        pstride.assign(d0 + 2, 1);
        qstride.assign(d0 + 2, 1);
        for (int g = 1; g <= d0; ++g) {
            pstride[g] = nsmallP; nsmallP *= (M[g] + 2);
            qstride[g] = nsmallQ; nsmallQ *= (M[g] + 1);
        }
        nb = 1 + (k - 1 - d0);
        qmap.assign(nsmallQ * nb, -1);
        std::vector<int> q(d0 + 1, 0);
        for (int bi = 0; bi < nb; ++bi) {
            int b = bi == 0 ? 0 : d0 + bi;
            for (int64_t t = 0; t < nsmallQ; ++t) {
                decodeQ(t, q);
                if (compatible(q, b)) qmap[bi * nsmallQ + t] = nq++;
            }
        }
    }
    static int bigIndex(int b, int d0) { return b == 0 ? 0 : b - d0; }
    void decodeP(int64_t t, std::vector<int>& p) const {
        for (int g = 1; g <= d0; ++g) { p[g] = t % (M[g] + 2); t /= (M[g] + 2); }
    }
    void decodeQ(int64_t t, std::vector<int>& q) const {
        for (int g = 1; g <= d0; ++g) { q[g] = t % (M[g] + 1); t /= (M[g] + 1); }
    }
    bool compatible(const std::vector<int>& p, int b) const {
        if (b == 0) return true;
        for (int h = k - b; h <= d0; ++h) if (h >= 1 && p[h] > 0) return false;
        return true;
    }
    // Best Fit move for item s from a state with small levels p (only zero/nonzero matters) and big gap b.
    Move move(const std::vector<int>& p, int b, int s) const {
        Move mv;
        if (ychain && s == 1) return mv;
        int target = 0;
        for (int g = s; g <= d0; ++g) if (p[g] > 0) { target = g; break; }
        if (!target && b && b >= s) target = b;
        if (target) {
            mv.dec = target;
            if (target > s) mv.inc = target - s;
        } else {
            mv.inc = k - s;
        }
        return mv;
    }
    // psi index (compressed) of the truncated state min(p + move, M), with new big gap.
    int64_t qindexAfter(const std::vector<int>& p, int b, const Move& mv) const {
        int64_t t = 0;
        int nbig = b;
        for (int g = 1; g <= d0; ++g) {
            int v = p[g];
            if (g == mv.dec) v -= 1;
            if (g == mv.inc) v += 1;
            if (v > M[g]) v = M[g];
            if (v < 0) { fprintf(stderr, "negative level\n"); exit(1); }
            t += (int64_t)v * qstride[g];
        }
        if (mv.dec > d0) nbig = 0;
        if (mv.inc > d0) {
            if (nbig != 0) { fprintf(stderr, "two big bins\n"); exit(1); }
            nbig = mv.inc;
        }
        int64_t r = qmap[bigIndex(nbig, d0) * nsmallQ + t];
        if (r < 0) { fprintf(stderr, "incompatible successor\n"); exit(1); }
        return r;
    }
    int64_t qindex(const std::vector<int>& p, int b) const {
        Move none;
        return qindexAfter(p, b, none);
    }
    bool saturated(const std::vector<int>& p) const {
        for (int g = satFrom; g <= d0; ++g) if (p[g] == M[g] + 1) return true;
        return false;
    }
    int64_t pindex(const std::vector<int>& p, int b) const {
        int64_t t = 0;
        for (int g = 1; g <= d0; ++g) t += (int64_t)p[g] * pstride[g];
        return (int64_t)bigIndex(b, d0) * nsmallP + t;
    }
    int decodePattern(int64_t idx, std::vector<int>& p) const {
        int bi = (int)(idx / nsmallP);
        decodeP(idx % nsmallP, p);
        return bi == 0 ? 0 : d0 + bi;
    }
    // Abstract successors of pattern (p,b) under item s.  The concrete successor's pattern is
    // min(x + Delta, M+1); when a saturated level (M+1, i.e. x >= M+1) is decremented the
    // successor level is M or M+1, so both are produced.  This over-approximates the set of
    // patterns of reachable states.
    template <class F> void abstractSucc(const std::vector<int>& p, int b, int s, F&& f) const {
        Move mv = move(p, b, s);
        std::vector<int> p2 = p;
        int b2 = b;
        if (mv.dec) {
            if (mv.dec > d0) b2 = 0;
            else p2[mv.dec] -= 1;
        }
        if (mv.inc) {
            if (mv.inc > d0) b2 = mv.inc;
            else p2[mv.inc] = std::min(p2[mv.inc] + 1, M[mv.inc] + 1);
        }
        f(p2, b2);
        if (mv.dec && mv.dec <= d0 && p[mv.dec] == M[mv.dec] + 1) {
            p2[mv.dec] += 1;
            f(p2, b2);
        }
    }
    // Patterns reachable from the empty state in the abstract pattern graph (sorted indices).
    std::vector<int64_t> reachable() const {
        int64_t total = nsmallP * nb;
        bool dense = total <= ((int64_t)1 << 33);  // bitmap up to 1 GiB, else hash set
        std::vector<bool> seenBits(dense ? total : 0, false);
        std::unordered_set<int64_t> seenSet;
        auto visit = [&](int64_t id) {  // returns true if newly seen
            if (dense) {
                if (seenBits[id]) return false;
                seenBits[id] = true;
                return true;
            }
            return seenSet.insert(id).second;
        };
        std::vector<int64_t> stack, out;
        std::vector<int> p(d0 + 1, 0);
        if (ychain) p[1] = 1;
        int64_t start = pindex(p, 0);
        visit(start);
        stack.push_back(start);
        while (!stack.empty()) {
            int64_t cur = stack.back();
            stack.pop_back();
            out.push_back(cur);
            int b = decodePattern(cur, p);
            for (int s = 1; s <= j; ++s)
                abstractSucc(p, b, s, [&](const std::vector<int>& p2, int b2) {
                    int64_t id = pindex(p2, b2);
                    if (visit(id)) stack.push_back(id);
                });
        }
        std::sort(out.begin(), out.end());
        return out;
    }
    // Enumerate patterns: all compatible ones, or only the abstractly reachable ones.
    template <class F> void forPatterns(F&& f, bool reachOnly = false) const {
        std::vector<int> p(d0 + 1, 0);
        if (reachOnly) {
            for (int64_t id : reachable()) {
                int b = decodePattern(id, p);
                if (!compatible(p, b)) { fprintf(stderr, "reachable but incompatible\n"); exit(1); }
                f(p, b);
            }
            return;
        }
        for (int bi = 0; bi < nb; ++bi) {
            int b = bi == 0 ? 0 : d0 + bi;
            for (int64_t t = 0; t < nsmallP; ++t) {
                decodeP(t, p);
                if (!compatible(p, b)) continue;
                f(p, b);
            }
        }
    }
};
