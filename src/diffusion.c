// src/diffusion.c — ensemble diffusion of the ultrametric walk.
#include "padic/diffusion.h"

#include "padic/tree.h"     // pwalk_*
#include "padic/padic.h"    // padic_valuation, padic_abs

#include <assert.h>
#include <math.h>

// A decorrelated seed for walker k. NOTE: seed + k*GOLDEN alone is wrong here —
// splitmix64 is a counter that steps by GOLDEN, so consecutive seeds produce the
// SAME stream shifted by one call, making the "ensemble" one correlated
// trajectory. Run it through the splitmix finalizer so the per-walker starting
// states are scattered off that lattice and the streams are effectively independent.
static unsigned long long subseed(unsigned long long seed, long k) {
    unsigned long long z = seed + (unsigned long long)k * 0x9E3779B97F4A7C15ULL;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

// ── exact closed forms ───────────────────────────────────────────────────────
double padic_escape_prob(int p, int n, double alpha, int m) {
    // P(jump separation s > m) with s ~ p^{-alpha s} / Z over s = 1..n. Weigh
    // RELATIVE to the largest term (s = 1): g(s) = p^{-alpha(s-1)} <= 1. The ratio
    // is unchanged, but nothing underflows to 0/0 for large alpha — so escape from
    // level 0 stays exactly 1, as it must (any jump leaves the single start leaf).
    double z = 0.0, tail = 0.0;
    for (int s = 1; s <= n; s++) {
        double w = pow((double)p, -alpha * (double)(s - 1));
        z += w;
        if (s > m) tail += w;
    }
    return tail / z;
}

double padic_escape_time(int p, int n, double alpha, int m) {
    double pe = padic_escape_prob(p, n, alpha, m);
    return pe > 0.0 ? 1.0 / pe : INFINITY;   // m >= n: the ball is everything, no escape
}

// Equilibrium is the uniform law over the p^n leaves. The FRACTION of leaves at
// separation s is (p-1)p^{s-1}/p^n = (p-1)p^{s-1-n} — a bounded double for any n,
// so folding the 1/p^n into the exponent avoids overflowing p^n or the u64 shell
// count (which wraps for p>=5 well before n=40). Such a leaf has |X|_p = p^{s-n};
// the origin (s = 0) has |0|_p = 0.
double padic_eq_sep(int p, int n) {
    double sum = 0.0;
    for (int s = 1; s <= n; s++)
        sum += (double)(p - 1) * pow((double)p, (double)(s - 1 - n)) * (double)s;
    return sum;
}
double padic_eq_abs(int p, int n) {
    double sum = 0.0;
    for (int s = 1; s <= n; s++)
        sum += (double)(p - 1) * pow((double)p, (double)(s - 1 - n)) * pow((double)p, (double)(s - n));
    return sum;
}
double padic_eq_msd(int p, int n) {
    double sum = 0.0;
    for (int s = 1; s <= n; s++)
        sum += (double)(p - 1) * pow((double)p, (double)(s - 1 - n)) * pow((double)p, 2.0 * (double)(s - n));
    return sum;
}

// ── Monte-Carlo measurements ─────────────────────────────────────────────────
double padic_mean_escape_steps(int p, int n, double alpha, int m,
                               long walkers, unsigned long long seed) {
    assert(m >= 0 && m < n && walkers > 0);    // m == n would never escape

    // Escape must be Monte-Carlo-samplable. When p_esc(m) falls below double
    // resolution the shell sampler's cumulative rounds to 1.0 and can never draw
    // s > m, so the walk would never escape (an ~1/p_esc-step event). Large alpha,
    // or m near n on a deep tree (e.g. p=2,n=60,m=55), land here — return the exact
    // closed form instead of looping forever.
    double pe = padic_escape_prob(p, n, alpha, m);
    if (pe < 1e-9) return padic_escape_time(p, n, alpha, m);

    double total = 0.0;
    for (long k = 0; k < walkers; k++) {
        pwalk w;
        pwalk_init(&w, p, n, alpha, subseed(seed, k), 0);
        long steps = 0;
        for (;;) {
            pwalk_step(&w);
            steps++;
            if (n - padic_valuation(&w.pos) > m) break;   // separation from start 0 exceeds m
        }
        total += (double)steps;
    }
    return total / (double)walkers;
}

void padic_diffuse(int p, int n, double alpha, long walkers, int T, int m0,
                   unsigned long long seed,
                   double *sep, double *absp, double *msd, double *surv) {
    assert(walkers > 0 && T >= 0 && m0 >= 0 && m0 <= n);
    for (int t = 0; t <= T; t++) {
        if (sep)  sep[t]  = 0.0;
        if (absp) absp[t] = 0.0;
        if (msd)  msd[t]  = 0.0;
        if (surv) surv[t] = 0.0;
    }
    for (long k = 0; k < walkers; k++) {
        pwalk w;
        pwalk_init(&w, p, n, alpha, subseed(seed, k), 0);
        if (surv) surv[0] += 1.0;              // t = 0: at the origin, s = 0 <= m0
        for (int t = 1; t <= T; t++) {
            pwalk_step(&w);
            int s = n - padic_valuation(&w.pos);
            double a = padic_abs(&w.pos);      // = p^{-(n-s)}
            if (sep)  sep[t]  += (double)s;
            if (absp) absp[t] += a;
            if (msd)  msd[t]  += a * a;
            if (surv && s <= m0) surv[t] += 1.0;
        }
    }
    double inv = 1.0 / (double)walkers;
    for (int t = 0; t <= T; t++) {
        if (sep)  sep[t]  *= inv;
        if (absp) absp[t] *= inv;
        if (msd)  msd[t]  *= inv;
        if (surv) surv[t] *= inv;
    }
}
