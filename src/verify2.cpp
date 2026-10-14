// verify2: a second, independent exact verifier for Best Fit U{j,k} drift certificates.
//
// Re-implements scripts/verify_indep.py in C++ (for speed) without using bfchain.hpp: Best Fit is
// applied to explicit gap-count vectors, the abstract reachability closure is computed by its own BFS
// over representative states, and psi is looked up through psi_states.txt (the only shared artefact).
// For every pattern the drift is evaluated on several concrete representatives (saturated coordinates
// at M+1, M+2, M+7) and must agree.
//
//   verify2 j k M_1,..,M_d0 lpdir certfile [reach]
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <unordered_set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

typedef __int128 i128;
static int J, K, D0;
static std::vector<int> Mv;   // Mv[g], g = 1..D0

static std::string slurp(const std::string& path) {
    std::string out;
    if (path.size() > 3 && path.compare(path.size() - 3, 3, ".gz") == 0) {
        FILE* p = popen(("gzip -dc '" + path + "'").c_str(), "r");
        char buf[1 << 16];
        size_t n;
        while ((n = fread(buf, 1, sizeof buf, p)) > 0) out.append(buf, n);
        pclose(p);
        return out;
    }
    std::ifstream f(path);
    std::stringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

// explicit Best Fit on counts x[1..K-1]
static void bestFit(std::vector<long>& x, int s) {
    for (int g = s; g < K; ++g)
        if (x[g] > 0) {
            x[g]--;
            if (g > s) x[g - s]++;
            return;
        }
    x[K - s]++;
}

static bool pairwiseOK(const std::vector<long>& x) {
    int a = 0, b = 0;  // two largest gaps (with multiplicity)
    for (int g = K - 1; g >= 1 && !b; --g)
        for (long c = 0; c < x[g] && !b; ++c) {
            if (!a) a = g;
            else b = g;
        }
    return !(a && b && a + b >= K);
}

struct Pat {
    std::vector<int> p;  // levels for gaps 1..D0 (0..M+1)
    int big;
};

// This verifier's own encodings (independent of bfchain.hpp): pattern levels p_g in 0..M_g+1 and
// truncated levels q_g in 0..M_g, mixed radix in g, then the big gap (0 or D0+1..K-1).
static std::vector<int64_t> PST, QST;
static int64_t NSP = 1, NSQ = 1;
static int bigIdx(int big) { return big == 0 ? 0 : big - D0; }
static int64_t encodeP(const Pat& P) {
    int64_t t = 0;
    for (int g = 1; g <= D0; ++g) t += (int64_t)P.p[g] * PST[g];
    return (int64_t)bigIdx(P.big) * NSP + t;
}
static Pat decodeP(int64_t code) {
    Pat P;
    P.p.assign(D0 + 1, 0);
    int bi = (int)(code / NSP);
    P.big = bi == 0 ? 0 : D0 + bi;
    int64_t t = code % NSP;
    for (int g = 1; g <= D0; ++g) { P.p[g] = (int)(t % (Mv[g] + 2)); t /= (Mv[g] + 2); }
    return P;
}
static int64_t encodeQ(const std::vector<int>& q, int big) {
    int64_t t = 0;
    for (int g = 1; g <= D0; ++g) t += (int64_t)q[g] * QST[g];
    return (int64_t)bigIdx(big) * NSQ + t;
}

static Pat patternOf(const std::vector<long>& x, bool& ok) {
    Pat P;
    P.p.assign(D0 + 1, 0);
    P.big = 0;
    ok = true;
    for (int g = 1; g <= D0; ++g) P.p[g] = (int)std::min<long>(x[g], Mv[g] + 1);
    for (int g = D0 + 1; g < K; ++g)
        if (x[g] > 0) {
            if (x[g] > 1 || P.big) ok = false;  // the big-gap invariant would be violated
            P.big = g;
        }
    return P;
}

static std::vector<long> representative(const Pat& P, int shift, int only) {
    std::vector<long> x(K, 0);
    for (int g = 1; g <= D0; ++g) {
        x[g] = P.p[g];
        if (P.p[g] == Mv[g] + 1 && (only == 0 || only == g)) x[g] = Mv[g] + 1 + shift;
    }
    if (P.big) x[P.big] = 1;
    return x;
}

int main(int argc, char** argv) {
    if (argc < 6) { fprintf(stderr, "usage: verify2 j k M lpdir cert [reach]\n"); return 2; }
    J = atoi(argv[1]);
    K = atoi(argv[2]);
    std::vector<int> v;
    {
        std::stringstream ss(argv[3]);
        std::string t;
        while (std::getline(ss, t, ',')) v.push_back(atoi(t.c_str()));
    }
    D0 = v.size() == 1 ? (K - 1) / 2 : (int)v.size();
    if (!(D0 == (K - 1) / 2 || (D0 == J && J < (K - 1) / 2))) { fprintf(stderr, "bad M\n"); return 2; }
    Mv.assign(D0 + 1, 0);
    for (int g = 1; g <= D0; ++g) Mv[g] = v.size() == 1 ? v[0] : v[g - 1];
    PST.assign(D0 + 1, 0);
    QST.assign(D0 + 1, 0);
    for (int g = 1; g <= D0; ++g) {
        PST[g] = NSP; NSP *= Mv[g] + 2;
        QST[g] = NSQ; NSQ *= Mv[g] + 1;
    }
    std::string lpdir = argv[4], certfile = argv[5];
    bool reach = argc > 6 && std::string(argv[6]) == "reach";

    // psi index table: "idx q_1 .. q_D0 | big"
    std::unordered_map<int64_t, long> psiIndex;
    {
        std::ifstream f(lpdir + "/psi_states.txt");
        std::string line;
        while (std::getline(f, line)) {
            std::stringstream ss(line);
            long idx;
            ss >> idx;
            std::vector<int> q(D0 + 1, 0);
            for (int g = 1; g <= D0; ++g) ss >> q[g];
            std::string bar;
            int big;
            ss >> bar >> big;
            if (!psiIndex.emplace(encodeQ(q, big), idx).second) { fprintf(stderr, "duplicate psi state\n"); return 2; }
        }
    }
    // certificate
    std::string text = slurp(certfile);
    std::stringstream cs(text);
    std::string tag, mode;
    long long Dden;
    cs >> tag >> mode >> tag >> Dden >> tag;
    std::vector<long long> c(K, 0);
    for (int g = 1; g <= K - 1; ++g) cs >> c[g];
    cs >> tag;
    std::vector<long long> psi;
    long long val;
    while (cs >> val) psi.push_back(val);
    if (psi.size() != psiIndex.size()) { fprintf(stderr, "psi size mismatch %zu vs %zu\n", psi.size(), psiIndex.size()); return 2; }
    bool ymode = mode == "ylo" || mode == "yhi" || mode == "ybounded";
    bool minMode = mode == "linear" || mode == "ylo";
    int satFrom = (mode == "bounded2" || mode == "ybounded") ? 2 : 1;
    bool allPatterns = mode == "linear" || mode == "upper" || mode == "ylo" || mode == "yhi";
    if (ymode && Mv[1] != 0) { fprintf(stderr, "Y-chain needs M_1 = 0\n"); return 2; }

    auto V = [&](const std::vector<long>& x) -> i128 {
        std::vector<int> q(D0 + 1, 0);
        int big = 0;
        for (int g = 1; g <= D0; ++g) q[g] = (int)std::min<long>(x[g], Mv[g]);
        for (int g = D0 + 1; g < K; ++g) if (x[g] > 0) big = g;
        auto it = psiIndex.find(encodeQ(q, big));
        if (it == psiIndex.end()) { fprintf(stderr, "state missing from psi table\n"); exit(3); }
        i128 r = psi[it->second];
        for (int g = 1; g < K; ++g) r += (i128)c[g] * x[g];
        return r;
    };
    auto step = [&](std::vector<long>& x, int s) {
        if (ymode && s == 1) return;  // gap-1 reservoir absorbs size-1 items
        bestFit(x, s);
    };

    // pattern set
    std::vector<int64_t> pats;  // pattern codes of the closure
    if (reach) {
        std::unordered_set<int64_t> seen;
        std::vector<int64_t> stack;
        Pat start;
        start.p.assign(D0 + 1, 0);
        start.big = 0;
        if (ymode) start.p[1] = 1;
        seen.insert(encodeP(start));
        stack.push_back(encodeP(start));
        while (!stack.empty()) {
            int64_t code = stack.back();
            stack.pop_back();
            pats.push_back(code);
            Pat P = decodeP(code);
            std::vector<long> x = representative(P, 0, 0);
            for (int s = 1; s <= J; ++s) {
                std::vector<long> y = x;
                step(y, s);
                bool ok;
                Pat Q = patternOf(y, ok);
                if (!ok || !pairwiseOK(y)) { fprintf(stderr, "invariant violated in BFS\n"); return 3; }
                std::vector<Pat> succ{Q};
                for (int g = 1; g <= D0; ++g)
                    if (P.p[g] == Mv[g] + 1 && y[g] < x[g]) {  // decremented saturated level may stay saturated
                        Pat R = Q;
                        R.p[g] = Mv[g] + 1;
                        succ.push_back(R);
                    }
                for (auto& S : succ) {
                    int64_t c2 = encodeP(S);
                    if (seen.insert(c2).second) stack.push_back(c2);
                }
            }
        }
    } else {
        fprintf(stderr, "verify2 supports reach mode only\n");
        return 2;
    }
    bool have = false;
    i128 extreme = 0;
    long checked = 0;
    for (int64_t code : pats) {
        Pat P = decodeP(code);
        std::vector<int> sat;
        for (int g = 1; g <= D0; ++g) if (P.p[g] == Mv[g] + 1) sat.push_back(g);
        bool satSel = false;
        for (int g : sat) if (g >= satFrom) satSel = true;
        if (!allPatterns && !satSel) continue;
        ++checked;
        bool first = true;
        i128 d0v = 0;
        for (int shift : {0, 1, 6}) {
            std::vector<int> onlys{0};
            if (!sat.empty()) onlys.push_back(sat[0]);
            for (int only : onlys) {
                std::vector<long> x = representative(P, shift, only);
                i128 vx = V(x), tot = 0;
                for (int s = 1; s <= J; ++s) {
                    std::vector<long> y = x;
                    step(y, s);
                    tot += V(y) - vx;
                }
                if (first) { d0v = tot; first = false; }
                else if (tot != d0v) { fprintf(stderr, "drift not constant on a pattern\n"); return 3; }
            }
        }
        if (!have || (minMode ? d0v < extreme : d0v > extreme)) { extreme = d0v; have = true; }
    }
    bool ok;
    if (mode == "linear") ok = extreme > 0;
    else if (mode == "ylo" || mode == "yhi") ok = c[1] == Dden;
    else if (mode == "upper") {
        ok = true;
        for (int g = 1; g <= D0; ++g) if (c[g] < (long long)g * Dden) ok = false;
    } else {
        ok = extreme < 0;
        for (int g = satFrom; g <= D0; ++g) if (c[g] <= 0) ok = false;
        if (satFrom == 2 && c[1] != 0) ok = false;
    }
    printf("verify2: mode=%s patterns=%zu checked=%ld extreme j*D*drift=%lld (D=%lld, j=%d) drift=%.12Lg\n",
           mode.c_str(), pats.size(), checked, (long long)extreme, Dden, J,
           (long double)extreme / ((long double)J * Dden));
    printf("VERDICT: %s\n", ok ? "CERTIFIED" : "FAILED");
    return ok ? 0 : 1;
}
