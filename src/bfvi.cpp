// bfvi: value iteration for the adversarial truncated Best Fit chain.
//
// For a fixed linear part c (small gaps only), the best correction psi in
//     V(x) = sum_g c_g x_g + psi(min(x, M))
// is the relative value function of an average-reward MDP on truncated states q, in which an
// adversary picks the pattern p in pre(q) (which saturated coordinates are "really" large):
//   linear : (T psi)(q) = min_{p in pre(q)}        (1/j) sum_s [ c.Delta_s(p) + psi(q'_s(p)) ]
//   bounded: (T psi)(q) = max_{p in pre_sat(q)}    (same)       over q with a saturated preimage
// Any psi with min_q (T psi - psi)(q) > 0 (linear) / max_q (T psi - psi)(q) < 0 (bounded) is a
// drift certificate.  We also propagate the vector-valued values of the greedy policy to get a
// subgradient of the optimal gain with respect to c (used by the outer cutting-plane loop).
//
//   bfvi j k M_1,..,M_d0 mode cfile outprefix iters threads [reach]
// cfile: d0 numbers c_1..c_d0.  Writes outprefix.psi (binary doubles over compressed psi index)
// and prints: lower/upper gain bounds and the subgradient at the extreme state.
#include "bfchain.hpp"
#include <algorithm>
#include <atomic>
#include <cmath>
#include <fstream>
#include <thread>

struct Pat {
    int64_t q;          // compressed psi index of min(p, M)
    bool sat;
};

int main(int argc, char** argv) {
    if (argc < 9) { fprintf(stderr, "usage: bfvi j k M mode cfile outprefix iters threads [reach]\n"); return 2; }
    int j = atoi(argv[1]), k = atoi(argv[2]);
    std::vector<int> M;
    int d0;
    {
        std::vector<int> v;
        std::string s = argv[3], tok;
        size_t pos = 0;
        while (true) {
            size_t e = s.find(',', pos);
            v.push_back(std::stoi(s.substr(pos, e - pos)));
            if (e == std::string::npos) break;
            pos = e + 1;
        }
        d0 = v.size() == 1 ? (k - 1) / 2 : (int)v.size();  // see Chain for the allowed values
        M.assign(d0 + 1, 0);
        for (int g = 1; g <= d0; ++g) M[g] = v.size() == 1 ? v[0] : v[g - 1];
    }
    // linear: min-adversary over all patterns; bounded: max-adversary over saturated patterns;
    // upper: max-adversary over all patterns (rate upper bounds with c_g >= g)
    // Y-chain modes: ylo (min, all patterns), yhi (max, all patterns), ybounded (as bounded2)
    std::string md = argv[4];
    bool ymode = md == "ylo" || md == "yhi" || md == "ybounded";
    bool linear = md == "linear" || md == "ylo";
    // bounded2: like bounded but only gaps >= 2 count as saturated (gap 1 may grow; c_1 = 0)
    bool bounded2 = md == "bounded2" || md == "ybounded";
    bool satOnly = md == "bounded" || bounded2;
    std::vector<double> c(d0 + 1, 0.0);
    {
        std::ifstream in(argv[5]);
        for (int g = 1; g <= d0; ++g) in >> c[g];
    }
    std::string outp = argv[6];
    long iters = atol(argv[7]);
    int nth = atoi(argv[8]);
    bool reach = false, warm = false;
    for (int a = 9; a < argc; ++a) {
        if (std::string(argv[a]) == "reach") reach = true;
        if (std::string(argv[a]) == "warm") warm = true;  // start from outprefix.psi if present
    }

    Chain ch(j, k, M);
    if (bounded2) ch.satFrom = 2;
    if (ymode) {
        if (M[1] != 0) { fprintf(stderr, "Y-chain modes need M_1 = 0\n"); return 2; }
        ch.ychain = true;
    }
    // collect patterns, grouped by q
    std::vector<int64_t> pidx;   // pattern index (Chain::pindex)
    std::vector<Pat> pats;
    ch.forPatterns([&](const std::vector<int>& p, int b) {
        bool sat = ch.saturated(p);
        if (satOnly && !sat) return;
        pidx.push_back(ch.pindex(p, b));
        pats.push_back({ch.qindex(p, b), sat});
    }, reach);
    int64_t np = (int64_t)pats.size();
    std::vector<int64_t> order(np);
    for (int64_t i = 0; i < np; ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&](int64_t a, int64_t b) { return pats[a].q < pats[b].q; });
    // successor table
    std::vector<int32_t> succ(np * j);
    std::vector<int8_t> dec(np * j), inc(np * j);
    std::vector<int64_t> qOf(np);
    {
        std::vector<int> p(d0 + 1);
        for (int64_t r = 0; r < np; ++r) {
            int64_t i = order[r];
            int b = ch.decodePattern(pidx[i], p);
            qOf[r] = pats[i].q;
            for (int s = 1; s <= j; ++s) {
                Move mv = ch.move(p, b, s);
                succ[r * j + s - 1] = (int32_t)ch.qindexAfter(p, b, mv);
                dec[r * j + s - 1] = (int8_t)(mv.dec <= d0 ? mv.dec : 0);
                inc[r * j + s - 1] = (int8_t)(mv.inc <= d0 ? mv.inc : 0);
            }
        }
    }
    // group starts
    std::vector<int64_t> gstart;
    std::vector<int64_t> gq;
    for (int64_t r = 0; r < np; ++r)
        if (r == 0 || qOf[r] != qOf[r - 1]) { gstart.push_back(r); gq.push_back(qOf[r]); }
    gstart.push_back(np);
    int64_t ng = (int64_t)gq.size();
    int64_t nq = ch.nq;
    fprintf(stderr, "bfvi j=%d k=%d mode=%s patterns=%lld groups=%lld nq=%lld\n", j, k, argv[4],
            (long long)np, (long long)ng, (long long)nq);

    std::vector<double> psi(nq, 0.0), npsi(nq, 0.0);
    if (warm) {
        FILE* wf = fopen((outp + ".psi").c_str(), "rb");
        if (wf) {
            size_t got = fread(psi.data(), sizeof(double), psi.size(), wf);
            fclose(wf);
            if (got != psi.size()) std::fill(psi.begin(), psi.end(), 0.0);
        }
    }
    std::vector<double> W((size_t)nq * d0, 0.0), nW((size_t)nq * d0, 0.0);
    std::vector<double> gain(ng);
    const double tau = 0.9;  // aperiodicity transform
    double lo = 0, hi = 0;
    double trackBest = linear ? -INFINITY : INFINITY;  // monotone bound we wait on
    long trackIt = 0;
    int64_t extremeG = 0;
    for (long it = 0; it < iters; ++it) {
        auto work = [&](int t) {
            std::vector<double> wv(d0 + 1);
            for (int64_t gi = t; gi < ng; gi += nth) {
                double best = linear ? INFINITY : -INFINITY;
                int64_t bestr = -1;
                for (int64_t r = gstart[gi]; r < gstart[gi + 1]; ++r) {
                    double v = 0;
                    for (int s = 0; s < j; ++s) {
                        int64_t e = r * j + s;
                        if (dec[e]) v -= c[dec[e]];
                        if (inc[e]) v += c[inc[e]];
                        v += psi[succ[e]];
                    }
                    v /= j;
                    if (linear ? v < best : v > best) { best = v; bestr = r; }
                }
                int64_t q = gq[gi];
                npsi[q] = (1 - tau) * psi[q] + tau * best;
                gain[gi] = npsi[q] - psi[q];
                // vector values under the greedy choice
                for (int g = 1; g <= d0; ++g) wv[g] = 0;
                for (int s = 0; s < j; ++s) {
                    int64_t e = bestr * j + s;
                    if (dec[e]) wv[dec[e]] -= 1;
                    if (inc[e]) wv[inc[e]] += 1;
                    for (int g = 1; g <= d0; ++g) wv[g] += W[(size_t)succ[e] * d0 + g - 1];
                }
                for (int g = 1; g <= d0; ++g)
                    nW[(size_t)q * d0 + g - 1] = (1 - tau) * W[(size_t)q * d0 + g - 1] + tau * wv[g] / j;
            }
        };
        std::vector<std::thread> th;
        for (int t = 0; t < nth; ++t) th.emplace_back(work, t);
        for (auto& x : th) x.join();
        lo = INFINITY; hi = -INFINITY;
        for (int64_t gi = 0; gi < ng; ++gi) {
            if (gain[gi] < lo) { lo = gain[gi]; if (linear) extremeG = gi; }
            if (gain[gi] > hi) { hi = gain[gi]; if (!linear) extremeG = gi; }
        }
        if (it % 200 == 0) fprintf(stderr, "  it=%ld gain in [%.8g, %.8g]\n", it, lo / tau, hi / tau);
        // the lower (upper) gain bound is monotone under value iteration; stop when it stagnates
        double cur = linear ? lo : hi;
        if (linear ? cur > trackBest + 1e-10 : cur < trackBest - 1e-10) { trackBest = cur; trackIt = it; }
        bool stagnant = it > 400 && it - trackIt > 300;
        bool last = it == iters - 1 || hi - lo < 1e-10 * tau || stagnant;
        if (last) {
            // psi (not yet replaced) is the function whose one-step gains were just computed.
            int64_t qe = gq[extremeG];
            printf("iters %ld\ngain_lo %.12g\ngain_hi %.12g\nsubgrad", it + 1, lo / tau, hi / tau);
            for (int g = 0; g < d0; ++g)
                printf(" %.12g", (nW[(size_t)qe * d0 + g] - W[(size_t)qe * d0 + g]) / tau);
            printf("\n");
            break;
        }
        // relative value iteration: renormalise and swap
        double ref = npsi[gq[0]];
        for (int64_t gi = 0; gi < ng; ++gi) npsi[gq[gi]] -= ref;
        std::vector<double> wref(d0);
        for (int g = 0; g < d0; ++g) wref[g] = nW[(size_t)gq[0] * d0 + g];
        for (int64_t gi = 0; gi < ng; ++gi)
            for (int g = 0; g < d0; ++g) nW[(size_t)gq[gi] * d0 + g] -= wref[g];
        std::swap(psi, npsi);
        std::swap(W, nW);
    }
    FILE* f = fopen((outp + ".psi").c_str(), "wb");
    fwrite(psi.data(), sizeof(double), psi.size(), f);
    fclose(f);
    return 0;
}
