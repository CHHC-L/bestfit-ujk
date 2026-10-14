// bfcert: generate drift LPs and exactly verify Lyapunov certificates for Best Fit under U{j,k}.
//
//   bfcert gen    j k M_1,...,M_d0 outdir [reach]    write LP rows (one per pattern) in binary CSR
//   bfcert verify j k M_1,...,M_d0 certfile [reach]  exact integer check of a certificate
// With "reach", only patterns reachable from the empty state in the abstract pattern graph are used.
//   bfcert test   j k nsteps [smalld0]          cross-check pattern moves against explicit Best Fit
//   bfcert states j k M_1,...,M_d0 outdir         write only psi_states.txt (psi index table)
//   bfcert count  j k M_1,...,M_d0 [reach]        number of patterns and psi states
//   bfcert ycover j k 0,M_2,...,M_d0              Theorem D coverage check (X closure vs Y closure)
//
// Row of pattern p (values scaled by j):  j*drift(p) = sum_s [ c.Delta_s + psi(q'_s) - psi(q) ].
// Columns: 0..k-2 -> c_1..c_{k-1};  k-1.. -> psi (compressed truncated states).
#include "bfchain.hpp"
#include <algorithm>
#include <map>
#include <random>
#include <fstream>
#include <unistd.h>
#include <sstream>

static std::string readMaybeGz(const std::string& path) {
    std::string out;
    if (path.size() > 3 && path.compare(path.size() - 3, 3, ".gz") == 0) {
        std::string cmd = "gzip -dc '" + path + "'";
        FILE* p = popen(cmd.c_str(), "r");
        if (!p) return out;
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

static std::vector<int> parseM(int k, const std::string& s) {
    // a single value means "same for all (k-1)/2 small gaps"; otherwise one level per small gap
    // (either (k-1)/2 of them, or j of them when j < (k-1)/2 -- see Chain).
    std::vector<int> v;
    std::stringstream ss(s);
    std::string tok;
    while (std::getline(ss, tok, ',')) v.push_back(std::stoi(tok));
    int d0 = v.size() == 1 ? (k - 1) / 2 : (int)v.size();
    std::vector<int> M(d0 + 1, 0);
    for (int g = 1; g <= d0; ++g) M[g] = v.size() == 1 ? v[0] : v[g - 1];
    return M;
}

// psi state table (compressed index, truncated small levels, big gap): the only artefact the
// independent verifiers take from this program.  It depends only on (j, k, M).
static void writePsiStates(const Chain& ch, const std::string& out) {
    std::string cmd = "mkdir -p '" + out + "'";
    if (system(cmd.c_str()) != 0) exit(1);
    // write to a temporary file and rename, so a concurrent reader never sees a partial table
    std::string tmp = out + "/psi_states.txt.tmp" + std::to_string((long long)getpid());
    FILE* g = fopen(tmp.c_str(), "w");
    if (!g) { fprintf(stderr, "cannot write %s\n", tmp.c_str()); exit(1); }
    std::vector<int> q(ch.d0 + 1);
    for (int bi = 0; bi < ch.nb; ++bi)
        for (int64_t t = 0; t < ch.nsmallQ; ++t) {
            int64_t r = ch.qmap[bi * ch.nsmallQ + t];
            if (r < 0) continue;
            ch.decodeQ(t, q);
            fprintf(g, "%lld", (long long)r);
            for (int h = 1; h <= ch.d0; ++h) fprintf(g, " %d", q[h]);
            fprintf(g, " | %d\n", bi == 0 ? 0 : ch.d0 + bi);
        }
    if (fclose(g) != 0 || rename(tmp.c_str(), (out + "/psi_states.txt").c_str()) != 0) {
        fprintf(stderr, "cannot write %s/psi_states.txt\n", out.c_str());
        exit(1);
    }
}

static int cmd_gen(int j, int k, const std::vector<int>& M, const std::string& out, bool reach) {
    Chain ch(j, k, M);
    std::vector<int64_t> rowptr{0};
    std::vector<int32_t> col;
    std::vector<int32_t> val;
    std::vector<uint8_t> flag;
    std::vector<std::pair<int64_t, int>> acc;
    ch.forPatterns([&](const std::vector<int>& p, int b) {
        std::map<int64_t, int> row;  // column -> coefficient
        int64_t q0 = ch.qindex(p, b);
        for (int s = 1; s <= j; ++s) {
            Move mv = ch.move(p, b, s);
            if (mv.dec) row[mv.dec - 1] -= 1;
            if (mv.inc) row[mv.inc - 1] += 1;
            int64_t q1 = ch.qindexAfter(p, b, mv);
            row[k - 1 + q1] += 1;
            row[k - 1 + q0] -= 1;
        }
        for (auto& [c, v] : row) if (v != 0) { col.push_back((int32_t)c); val.push_back(v); }
        rowptr.push_back((int64_t)col.size());
        flag.push_back(ch.saturated(p) ? 1 : 0);
    }, reach);
    int64_t nrows = (int64_t)flag.size(), ncols = (k - 1) + ch.nq, nnz = (int64_t)col.size();
    std::string cmd = "mkdir -p " + out;
    if (system(cmd.c_str()) != 0) return 1;
    FILE* f = fopen((out + "/lp.bin").c_str(), "wb");
    int64_t hdr[4] = {nrows, ncols, nnz, (int64_t)ch.nq};
    fwrite(hdr, sizeof(int64_t), 4, f);
    fwrite(rowptr.data(), sizeof(int64_t), rowptr.size(), f);
    fwrite(col.data(), sizeof(int32_t), col.size(), f);
    fwrite(val.data(), sizeof(int32_t), val.size(), f);
    fwrite(flag.data(), 1, flag.size(), f);
    fclose(f);
    writePsiStates(ch, out);
    printf("j=%d k=%d d0=%d patterns=%lld psi=%lld nnz=%lld\n", j, k, ch.d0, (long long)nrows,
           (long long)ch.nq, (long long)nnz);
    return 0;
}

// Certificate file (text):
//   mode linear|bounded
//   D <denominator>
//   c <k-1 integers>            numerators of c_1..c_{k-1}
//   psi <nq integers>           numerators of psi over compressed truncated states
// Linear: verify min over ALL patterns of j*D*drift > 0.
// Bounded: verify c_g > 0 for all g and max over saturated patterns of j*D*drift < 0.
static int cmd_verify(int j, int k, const std::vector<int>& M, const std::string& cf, bool reach) {
    Chain ch(j, k, M);
    const std::string text = readMaybeGz(cf);   // certificates may be gzip-compressed
    std::istringstream in(text);
    std::string mode0;
    {
        std::istringstream peek(text);
        std::string t;
        peek >> t >> mode0;
    }
    if (mode0 == "bounded2" || mode0 == "ybounded") ch.satFrom = 2;  // gap 1 excluded from V
    bool ymode = mode0 == "ylo" || mode0 == "yhi" || mode0 == "ybounded";
    if (ymode) {
        if (M[1] != 0) { fprintf(stderr, "Y-chain certificates need M_1 = 0\n"); return 2; }
        ch.ychain = true;
    }
    std::string tag, mode;
    long long D;
    in >> tag >> mode;
    in >> tag >> D;
    std::vector<long long> c(k), psi(ch.nq);
    in >> tag;
    for (int g = 1; g <= k - 1; ++g) in >> c[g];
    in >> tag;
    for (int64_t i = 0; i < ch.nq; ++i) in >> psi[i];
    if (!in) { fprintf(stderr, "bad certificate file\n"); return 2; }
    {
        std::string extra;
        if (in >> extra) { fprintf(stderr, "certificate has more psi values than psi states\n"); return 2; }
    }
    bool linear = mode == "linear" || mode == "ylo";  // minimum drift over all patterns
    bool upper = mode == "upper" || mode == "yhi";     // maximum drift over all patterns
    __int128 best = linear ? ((__int128)1 << 100) : -((__int128)1 << 100);
    long long npat = 0, nsat = 0;
    std::vector<int> worstP;
    int worstB = 0;
    ch.forPatterns([&](const std::vector<int>& p, int b) {
        ++npat;
        bool sat = ch.saturated(p);
        if (!linear && !upper && !sat) return;
        nsat += sat;
        int64_t q0 = ch.qindex(p, b);
        __int128 num = 0;
        for (int s = 1; s <= j; ++s) {
            Move mv = ch.move(p, b, s);
            if (mv.dec) num -= c[mv.dec];
            if (mv.inc) num += c[mv.inc];
            num += psi[ch.qindexAfter(p, b, mv)] - psi[q0];
        }
        if (linear ? num < best : num > best) { best = num; worstP = p; worstB = b; }
    }, reach);
    long double drift = (long double)best / ((long double)j * (long double)D);
    bool ok;
    // small-gap ratios c_g / g (big gaps hold at most one bin and only shift V by a constant)
    long double cratio_max = 0, cratio_min = 1e30L;
    for (int g = 1; g <= ch.d0; ++g) {
        long double r = (long double)c[g] / ((long double)g * D);
        cratio_max = std::max(cratio_max, r);
        cratio_min = std::min(cratio_min, r);
    }
    if (mode == "ylo" || mode == "yhi") {
        // bound on the stationary gap-1 creation rate pi(a) of the Y-chain (needs c_1 = 1)
        ok = c[1] == D;
        if (ok) printf("%s = %.12Lg   [pi(a) %s this; waste rate = pi(a) - 1/j %s %.12Lg]\n",
                       mode == "ylo" ? "creation_lower" : "creation_upper", drift,
                       mode == "ylo" ? ">=" : "<=", mode == "ylo" ? ">=" : "<=",
                       drift - 1.0L / j);
    } else if (linear) {
        ok = best > 0;
        if (ok && cratio_max > 0)
            printf("rate_lower = drift / max_g(c_g^+/g) = %.12Lg   [liminf E W_n / n]\n", drift / cratio_max);
    } else if (upper) {
        ok = cratio_min >= 1;  // V >= W - const
        if (ok) printf("rate_upper = %.12Lg   [limsup E W_n / n]\n", drift);
    } else {
        ok = best < 0;
        for (int g = ch.satFrom; g <= ch.d0; ++g) if (c[g] <= 0) ok = false;  // big gaps hold <= 1 bin
        if (ch.satFrom == 2 && c[1] != 0) ok = false;
    }
    long long cmax = 0, psimax = 0;
    for (int g = 1; g <= k - 1; ++g) cmax = std::max(cmax, std::llabs(c[g]));
    for (auto v : psi) psimax = std::max(psimax, std::llabs(v));
    printf("mode=%s j=%d k=%d patterns=%lld checked_saturated=%lld\n", mode.c_str(), j, k, npat, nsat);
    printf("extreme j*D*drift = %lld  (drift = %.12Lg)  D=%lld  max|c|=%lld max|psi|=%lld\n",
           (long long)best, drift, D, cmax, psimax);
    printf("worst pattern:");
    for (int g = 1; g <= ch.d0; ++g) printf(" %d", worstP.empty() ? -1 : worstP[g]);
    printf(" | big %d\n", worstB);
    printf("VERDICT: %s\n", ok ? "CERTIFIED" : "FAILED");
    return ok ? 0 : 1;
}

// Cross-check: simulate explicit Best Fit with an explicit list of bins and compare each step's
// effect on gap counts with Chain::move applied to the pattern of the current state.
static int cmd_test(int j, int k, long nsteps, bool smalld0) {
    std::vector<int> Mg(((smalld0 && j < (k - 1) / 2) ? j : (k - 1) / 2) + 1, 2);
    Chain ch(j, k, Mg);
    std::mt19937_64 rng(12345 + 1000 * j + k);
    std::vector<int> bins;  // gaps of open bins (gap>0), in opening order
    long checked = 0;
    for (long t = 0; t < nsteps; ++t) {
        int s = 1 + (int)(rng() % j);
        // explicit best fit
        int bi = -1;
        for (int i = 0; i < (int)bins.size(); ++i)
            if (bins[i] >= s && (bi < 0 || bins[i] < bins[bi])) bi = i;
        std::vector<long> before(k, 0), after(k, 0);
        for (int g : bins) before[g]++;
        if (bi >= 0) {
            bins[bi] -= s;
            if (bins[bi] == 0) bins.erase(bins.begin() + bi);
        } else bins.push_back(k - s);
        for (int g : bins) after[g]++;
        // invariant: pairwise gap sums < k
        std::vector<int> sorted = bins;
        std::sort(sorted.begin(), sorted.end());
        if (sorted.size() >= 2 && sorted[sorted.size() - 1] + sorted[sorted.size() - 2] >= k) {
            fprintf(stderr, "invariant violated\n");
            return 1;
        }
        // pattern move
        std::vector<int> p(ch.d0 + 1, 0);
        int b = 0;
        for (int g = 1; g < k; ++g) {
            if (g <= ch.d0) p[g] = before[g] > 0 ? 1 : 0;
            else if (before[g] > 0) {
                if (before[g] > 1 || b) { fprintf(stderr, "big invariant violated\n"); return 1; }
                b = g;
            }
        }
        Move mv = ch.move(p, b, s);
        std::vector<long> pred = before;
        if (mv.dec) pred[mv.dec]--;
        if (mv.inc) pred[mv.inc]++;
        for (int g = 1; g < k; ++g)
            if (pred[g] != after[g]) { fprintf(stderr, "mismatch at step %ld\n", t); return 1; }
        ++checked;
        if (bins.size() > 2000) bins.clear();  // keep it cheap: restart from empty
    }
    printf("test j=%d k=%d: %ld steps agree with explicit Best Fit\n", j, k, checked);
    return 0;
}

// Theorem D needs the Y-chain certificates to hold at every state the full chain X visits while
// x_1 >= 1 (then gaps >= 2 move exactly as in Y).  Check that every abstractly reachable pattern of X
// with p_1 = 1 (M_1 = 0, so p_1 = 1 means x_1 >= 1) lies in the abstract reachable closure of Y.
static int cmd_ycover(int j, int k, const std::vector<int>& M) {
    if (M[1] != 0) { fprintf(stderr, "ycover needs M_1 = 0\n"); return 2; }
    Chain cx(j, k, M), cy(j, k, M);
    cy.ychain = true;
    std::vector<int64_t> rx = cx.reachable(), ry = cy.reachable();
    std::vector<int> p(cx.d0 + 1);
    long long withGap1 = 0, missing = 0;
    for (int64_t id : rx) {
        cx.decodePattern(id, p);
        if (p[1] != 1) continue;
        ++withGap1;
        if (!std::binary_search(ry.begin(), ry.end(), id)) ++missing;
    }
    printf("ycover j=%d k=%d: X patterns with x_1>=1: %lld, Y patterns: %zu, not covered: %lld\n", j, k,
           withGap1, ry.size(), missing);
    printf("VERDICT: %s\n", missing == 0 ? "COVERED" : "NOT COVERED");
    return missing == 0 ? 0 : 1;
}

int main(int argc, char** argv) {
    if (argc < 4) { fprintf(stderr, "usage: see header\n"); return 2; }
    std::string cmd = argv[1];
    int j = atoi(argv[2]), k = atoi(argv[3]);
    if (j < 1 || j > k - 1) { fprintf(stderr, "need 1 <= j <= k-1\n"); return 2; }
    if (cmd == "test") return cmd_test(j, k, atol(argv[4]), argc > 5 && std::string(argv[5]) == "smalld0");
    std::vector<int> M = parseM(k, argv[4]);
    if (cmd == "count") {  // bfcert count j k M [reach]: number of patterns and psi states
        bool r = argc > 5 && std::string(argv[5]) == "reach";
        Chain ch(j, k, M);
        long long n = 0;
        ch.forPatterns([&](const std::vector<int>&, int) { ++n; }, r);
        printf("patterns %lld psi %lld\n", n, (long long)ch.nq);
        return 0;
    }
    if (cmd == "ycover") return cmd_ycover(j, k, M);
    if (cmd == "states") {  // bfcert states j k M outdir: write only psi_states.txt
        if (argc < 6) { fprintf(stderr, "usage: bfcert states j k M outdir\n"); return 2; }
        Chain ch(j, k, M);
        writePsiStates(ch, argv[5]);
        printf("psi states %lld\n", (long long)ch.nq);
        return 0;
    }
    bool reach = argc > 6 && std::string(argv[6]) == "reach";
    if (cmd == "gen") return cmd_gen(j, k, M, argv[5], reach);
    if (cmd == "verify") return cmd_verify(j, k, M, argv[5], reach);
    fprintf(stderr, "unknown command\n");
    return 2;
}
