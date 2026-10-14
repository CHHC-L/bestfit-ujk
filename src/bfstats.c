// Simulation statistics for Best Fit under U{j,k} (used only to choose proof parameters).
// Output (machine readable):
//   slope <waste growth per item over the second half>
//   tail g P(x_g>=1) P(x_g>=2) ... P(x_g>=16)     for g = 1..k-1 (second half of the run)
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
static uint64_t s[2] = {0x9E3779B97F4A7C15ULL, 0xD1B54A32D192ED03ULL};
static inline uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static inline uint64_t nxt(void) {
    uint64_t s0 = s[0], s1 = s[1], r = s0 + s1;
    s1 ^= s0;
    s[0] = rotl(s0, 55) ^ s1 ^ (s1 << 14);
    s[1] = rotl(s1, 36);
    return r;
}
#define MT 16
int main(int argc, char** argv) {
    if (argc < 4) { fprintf(stderr, "usage: bfstats j k n\n"); return 2; }
    int j = atoi(argv[1]), k = atoi(argv[2]);
    long n = atol(argv[3]);
    if (j < 1 || j >= k || k > 511) { fprintf(stderr, "need 1 <= j < k <= 511\n"); return 2; }
    s[0] ^= 1000003ULL * j + k;
    s[1] ^= 7919ULL * k + j;
    long N[512] = {0};
    static double ge[512][MT + 1];
    long waste = 0, whalf = 0;
    for (long t = 1; t <= n; t++) {
        int sz = 1 + (int)((nxt() >> 11) % (uint64_t)j);
        int g = sz;
        while (g < k && N[g] == 0) g++;
        if (g < k) { N[g]--; if (g - sz > 0) N[g - sz]++; waste -= sz; }
        else { if (k - sz > 0) N[k - sz]++; waste += k - sz; }
        if (t == n / 2) whalf = waste;
        if (t > n / 2)
            for (int h = 1; h < k; h++)
                for (int m = 1; m <= MT && N[h] >= m; m++) ge[h][m] += 1.0;
    }
    double T = (double)(n - n / 2);
    printf("slope %.8g\n", (double)(waste - whalf) / T);
    for (int h = 1; h < k; h++) {
        printf("tail %d", h);
        for (int m = 1; m <= MT; m++) printf(" %.4g", ge[h][m] / T);
        printf("\n");
    }
    return 0;
}
