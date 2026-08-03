// src/montecarlo.c — Haar-measure Monte-Carlo over Z_p and the zeta integral.
#include "padic/montecarlo.h"

#include <assert.h>
#include <math.h>

// ── exact closed forms ───────────────────────────────────────────────────────
double padic_zeta_integral(int p, double s) {
    assert(p >= 2 && s > -1.0);                 // integral diverges for s <= -1
    return (1.0 - 1.0 / (double)p) / (1.0 - pow((double)p, -(s + 1.0)));
}

double padic_zeta_variance(int p, double s) {
    assert(s > -0.5);                           // I(2s) needs 2s > -1, else infinite variance
    double i1 = padic_zeta_integral(p, s);
    double i2 = padic_zeta_integral(p, 2.0 * s);
    return i2 - i1 * i1;
}

// ── splitmix64 (each TU keeps its own tiny deterministic RNG) ─────────────────
static unsigned long long sm64(unsigned long long *st) {
    unsigned long long z = (*st += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

// ── Haar sampling and the estimator ──────────────────────────────────────────
void padic_haar_sample(padic *out, int p, int k, unsigned long long *rng) {
    assert(p >= 2 && k >= 1 && k <= PADIC_MAXK);
    *out = padic_zero(p, k);
    for (int i = 0; i < k; i++)                  // each digit uniform on {0..p-1} = Haar on Z/p^k
        out->a[i] = (int)(sm64(rng) % (unsigned long long)p);
}

double padic_mc_integral(int p, int k, double s, long N, unsigned long long seed, double *se) {
    // s > -1 keeps I(s) finite (the mean estimates it); the reported se is a valid
    // error bar only for s > -1/2, where the variance I(2s)-I(s)^2 is finite.
    assert(p >= 2 && k >= 1 && k <= PADIC_MAXK && N > 0 && s > -1.0);
    unsigned long long rng = seed;
    double mean = 0.0, m2 = 0.0;                 // Welford running mean + M2
    for (long i = 1; i <= N; i++) {
        padic x;
        padic_haar_sample(&x, p, k, &rng);
        int v = padic_valuation(&x);             // capped at k for the x == 0 sample
        double f = pow((double)p, -s * (double)v);   // |x|_p^s = p^{-s v}; exact 1 when s == 0
        double d = f - mean;
        mean += d / (double)i;
        m2 += d * (f - mean);
    }
    if (se) *se = (N > 1) ? sqrt(m2 / (double)(N - 1) / (double)N) : 0.0;  // sample sd / sqrt(N)
    return mean;
}
