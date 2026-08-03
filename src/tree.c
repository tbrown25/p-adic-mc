// src/tree.c — the p-adic tree metric and the ultrametric random walk.
#include "padic/tree.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

// ── tree metric ──────────────────────────────────────────────────────────────
int padic_tree_ancestor_depth(const padic *x, const padic *y) {
    padic d = padic_sub(x, y);          // v_p(x-y) = # matching low digits; caps at k
    return padic_valuation(&d);
}

int padic_tree_sep(const padic *x, const padic *y) {
    return x->k - padic_tree_ancestor_depth(x, y);   // n - c
}

unsigned long long padic_shell_count(int p, int n, int s) {
    if (s < 1 || s > n) return 0;       // (p-1) choices at the split digit, p each above it
    unsigned long long c = (unsigned long long)(p - 1);
    for (int i = 0; i < s - 1; i++) c *= (unsigned long long)p;  // u64: (p-1)p^{s-1}
    return c;                           // overflows int well within n<=64, so accumulate wide
}

// ── deterministic RNG: splitmix64 (tiny, dependency-free, reproducible) ───────
static unsigned long long sm64(unsigned long long *s) {
    unsigned long long z = (*s += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}
static double sm64_unit(unsigned long long *s) {       // uniform in [0, 1)
    return (double)(sm64(s) >> 11) * (1.0 / 9007199254740992.0);
}
static int sm64_below(unsigned long long *s, int m) {  // uniform in [0, m), m small
    return (int)(sm64(s) % (unsigned long long)m);
}

// ── the walk ─────────────────────────────────────────────────────────────────
void pwalk_init(pwalk *w, int p, int n, double alpha, unsigned long long seed,
                unsigned long long start) {
    assert(p >= 2 && n >= 1 && n <= PADIC_MAXK && alpha > 0.0);
    w->p = p;
    w->n = n;
    w->alpha = alpha;
    w->rng = seed;
    w->pos = padic_from_ull(p, n, start);
    w->last_sep = 0;

    // Shell weights are p^{-alpha·s} (s = 1..n). Work RELATIVE to the largest
    // term (s = 1): g(s) = p^{-alpha·(s-1)} <= 1. Identical normalized
    // probabilities, but no underflow-to-zero / 0-over-0 for large alpha — there
    // the mass correctly collapses onto shell 1 (the local limit) instead of a
    // NaN table that would default every draw to the farthest shell.
    double z = 0.0;
    for (int s = 1; s <= n; s++) z += pow((double)p, -alpha * (double)(s - 1));
    double acc = 0.0;
    w->shell_cum[0] = 0.0;
    for (int s = 1; s <= n; s++) {
        acc += pow((double)p, -alpha * (double)(s - 1)) / z;
        w->shell_cum[s] = acc;
    }
    w->shell_cum[n] = 1.0;              // pin the top against floating-point drift
}

unsigned long long pwalk_step(pwalk *w) {
    // 1) choose a shell s in {1..n} with probability p^{-alpha·s} / Z.
    double u = sm64_unit(&w->rng);
    int s = w->n;
    for (int t = 1; t <= w->n; t++) {
        if (u <= w->shell_cum[t]) { s = t; break; }
    }
    w->last_sep = s;

    // 2) a uniform leaf at separation s: agree in the low n-s digits, differ at
    //    digit d = n-s, free above it. Then v_p(diff) = d = c, so sep = n-c = s.
    int d = w->n - s;
    padic y = w->pos;
    int delta = 1 + sm64_below(&w->rng, w->p - 1);          // 1..p-1: guarantees a change
    y.a[d] = (w->pos.a[d] + delta) % w->p;
    for (int i = d + 1; i < w->n; i++) y.a[i] = sm64_below(&w->rng, w->p);
    // digits 0..d-1 are left equal to w->pos.

    w->pos = y;
    return padic_to_ull(&w->pos);
}

int pwalk_last_sep(const pwalk *w) { return w->last_sep; }

void pwalk_run(pwalk *w, int steps, void *out, int header) {
    FILE *f = (FILE *)out;
    padic start = w->pos;
    padic prev = w->pos;
    if (f && header)
        fprintf(f, "step,site,sep_prev,absp_prev,sep_start,absp_start\n");
    for (int i = 1; i <= steps; i++) {
        unsigned long long site = pwalk_step(w);
        if (f) {
            fprintf(f, "%d,%llu,%d,%.6g,%d,%.6g\n",
                    i, site,
                    w->last_sep, padic_dist(&prev, &w->pos),
                    padic_tree_sep(&start, &w->pos), padic_dist(&start, &w->pos));
        }
        prev = w->pos;
    }
}
