// tests/test_tree.c — the tree metric and the ultrametric walk. No framework.
#include "padic/tree.h"

#include <math.h>
#include <stdio.h>

static int fails = 0;
static void check(const char *name, int ok) {
    printf("  [%s] %s\n", ok ? "PASS" : "FAIL", name);
    if (!ok) fails++;
}

// ── the shell MARGINAL: P(sep == s) should equal p^{-alpha·s} / Z. ────────────
// Tolerance is ~6 binomial standard errors, so a mis-set tail exponent (e.g. a
// sampler that ignores alpha, or is off by 10%) is caught, while the correct
// sampler passes deterministically for the fixed seed.
static int shell_marginal_ok(int p, int n, double alpha, unsigned long long seed, long N) {
    pwalk w;
    pwalk_init(&w, p, n, alpha, seed, 0);
    double theo[PADIC_MAXK + 1] = {0}, z = 0.0;
    for (int s = 1; s <= n; s++) z += pow((double)p, -alpha * s);
    for (int s = 1; s <= n; s++) theo[s] = pow((double)p, -alpha * s) / z;

    long cnt[PADIC_MAXK + 1] = {0};
    for (long i = 0; i < N; i++) {
        pwalk_step(&w);
        cnt[pwalk_last_sep(&w)]++;
    }
    for (int s = 1; s <= n; s++) {
        double emp = (double)cnt[s] / (double)N;
        double se = sqrt(theo[s] * (1.0 - theo[s]) / (double)N);
        if (fabs(emp - theo[s]) > 6.0 * se + 1e-9) return 0;
    }
    return 1;
}

// ── WITHIN a shell the leaf must be uniform. From a fixed origin 0, a jump to
// separation s sets the split digit (position n-s) uniformly over the p-1
// nonzero values and every higher "free" digit uniformly over 0..p-1. Resetting
// to 0 each step samples that conditional directly. This is what the shell
// marginal alone cannot see — it catches a deterministic split digit or
// un-randomized free digits (both of which keep the separation correct).
static int within_shell_ok(int p, int n, double alpha, unsigned long long seed, long N) {
    pwalk w;
    pwalk_init(&w, p, n, alpha, seed, 0);
    long splitc[8] = {0}, freec[8] = {0};      // p <= 7 for these tests
    long split_total = 0, free_total = 0;
    for (long i = 0; i < N; i++) {
        w.pos = padic_from_ull(p, n, 0);       // resample a jump FROM the origin
        pwalk_step(&w);
        int s = pwalk_last_sep(&w);
        splitc[w.pos.a[n - s]]++;              // the split digit (must be in 1..p-1)
        split_total++;
        if (s >= 2) { freec[w.pos.a[n - 1]]++; free_total++; }  // a genuinely free digit
    }
    int ok = (splitc[0] == 0);                 // the split digit is never equal (never 0 vs origin)
    double e_split = (double)split_total / (double)(p - 1);
    for (int v = 1; v < p; v++) {              // uniform over the p-1 nonzero values
        double se = sqrt(e_split * (1.0 - 1.0 / (double)(p - 1)));
        if (fabs((double)splitc[v] - e_split) > 6.0 * se + 1e-9) ok = 0;
    }
    if (free_total > 0) {                      // uniform over all p values
        double e_free = (double)free_total / (double)p;
        for (int v = 0; v < p; v++) {
            double se = sqrt(e_free * (1.0 - 1.0 / (double)p));
            if (fabs((double)freec[v] - e_free) > 6.0 * se + 1e-9) ok = 0;
        }
    }
    return ok;
}

int main(void) {
    printf("p-adic tree — metric, shells, and the ultrametric walk\n\n");

    const int p = 3, n = 6;

    // ── separation identities ────────────────────────────────────────────────
    {
        padic x = padic_from_ull(p, n, 42);
        check("sep(x, x) == 0", padic_tree_sep(&x, &x) == 0);

        int ok_sep = 1, ok_dist = 1;
        for (int d = 0; d < n; d++) {
            padic y = x;
            y.a[d] = (x.a[d] + 1) % p;          // differ only at digit d
            if (padic_tree_sep(&x, &y) != n - d) ok_sep = 0;
            if (fabs(padic_dist(&x, &y) - pow(p, -d)) > 1e-12) ok_dist = 0;
        }
        check("differ at digit d  =>  sep == n-d", ok_sep);
        check("differ at digit d  =>  |x-y|_p == p^-d", ok_dist);
    }

    // ── shell counts (p-1)p^{s-1}, and they tile all p^m - 1 other leaves ─────
    {
        int m = 3;
        long total = 1;
        for (int i = 0; i < m; i++) total *= p;
        long counted[PADIC_MAXK + 1] = {0};
        padic x = padic_from_ull(p, m, 0);
        for (long y = 0; y < total; y++) {
            padic yy = padic_from_ull(p, m, (unsigned long long)y);
            counted[padic_tree_sep(&x, &yy)]++;
        }
        int ok = (counted[0] == 1);
        long sum = 0;
        for (int s = 1; s <= m; s++) {
            if ((unsigned long long)counted[s] != padic_shell_count(p, m, s)) ok = 0;
            sum += counted[s];
        }
        check("shell counts == (p-1)p^{s-1}", ok);
        check("shells tile all p^m - 1 other leaves", sum == total - 1);
    }

    // ── the tree is vertex-transitive: the shell counts hold from EVERY leaf, ─
    //    not just 0 (guards a base-point-dependent off-by-one in separation).
    {
        int m = 3, total = 1;
        for (int i = 0; i < m; i++) total *= p;
        int ok = 1;
        for (int ai = 0; ai < total && ok; ai++) {
            long cc[PADIC_MAXK + 1] = {0};
            padic a = padic_from_ull(p, m, (unsigned)ai);
            for (int bi = 0; bi < total; bi++) {
                padic b = padic_from_ull(p, m, (unsigned)bi);
                cc[padic_tree_sep(&a, &b)]++;
            }
            if (cc[0] != 1) ok = 0;
            for (int s = 1; s <= m; s++)
                if ((unsigned long long)cc[s] != padic_shell_count(p, m, s)) ok = 0;
        }
        check("shell counts hold from every base leaf (vertex-transitive)", ok);
    }

    // ── the tree separation is an ultrametric (strong triangle) ──────────────
    {
        int m = 3, total = 1;
        for (int i = 0; i < m; i++) total *= p;  // 27 leaves -> 19683 triples
        int ok = 1;
        for (int ai = 0; ai < total && ok; ai++)
            for (int bi = 0; bi < total && ok; bi++)
                for (int ci = 0; ci < total && ok; ci++) {
                    padic a = padic_from_ull(p, m, (unsigned)ai);
                    padic b = padic_from_ull(p, m, (unsigned)bi);
                    padic c = padic_from_ull(p, m, (unsigned)ci);
                    int sac = padic_tree_sep(&a, &c);
                    int sab = padic_tree_sep(&a, &b);
                    int sbc = padic_tree_sep(&b, &c);
                    int mx = sab > sbc ? sab : sbc;
                    if (sac > mx) ok = 0;
                }
        check("ultrametric: sep(a,c) <= max(sep(a,b), sep(b,c))", ok);
    }

    // ── the sampler's shell marginal matches p^{-alpha·s}/Z across p and alpha ─
    {
        check("P(shell s) ~ p^{-1.0 s}  (p=3)", shell_marginal_ok(3, 6, 1.0, 0xC0FFEEULL, 300000));
        check("P(shell s) ~ p^{-0.5 s}  (p=3)", shell_marginal_ok(3, 6, 0.5, 0xBEEF01ULL, 400000));
        check("P(shell s) ~ p^{-2.0 s}  (p=3)", shell_marginal_ok(3, 6, 2.0, 0x5EED09ULL, 300000));
        check("P(shell s) ~ p^{-1.0 s}  (p=2)", shell_marginal_ok(2, 10, 1.0, 0xA11CE0ULL, 300000));
    }

    // ── within-shell uniformity across p and alpha ───────────────────────────
    {
        check("within-shell uniform  (p=3, alpha=1)", within_shell_ok(3, 6, 1.0, 0xD00D22ULL, 200000));
        check("within-shell uniform  (p=3, alpha=2)", within_shell_ok(3, 6, 2.0, 0x1234ABULL, 300000));
        check("within-shell uniform  (p=2, alpha=1)", within_shell_ok(2, 10, 1.0, 0x99AA33ULL, 200000));
    }

    // ── the walk stays in range and every jump actually moves ────────────────
    {
        pwalk w;
        pwalk_init(&w, p, n, 1.0, 0x777ULL, 0);
        long cap = 1;
        for (int i = 0; i < n; i++) cap *= p;     // p^n
        int in_range = 1, no_self = 1;
        unsigned long long prev = padic_to_ull(&w.pos);
        for (long i = 0; i < 100000; i++) {
            unsigned long long site = pwalk_step(&w);
            if (site >= (unsigned long long)cap) in_range = 0;
            if (site == prev) no_self = 0;
            prev = site;
        }
        check("every site stays in [0, p^n)", in_range);
        check("every jump actually moves (no self-jump)", no_self);
    }

    // ── determinism: same seed -> identical trajectory; differ otherwise ─────
    {
        pwalk a, b, c;
        pwalk_init(&a, p, n, 1.0, 12345ULL, 4);
        pwalk_init(&b, p, n, 1.0, 12345ULL, 4);
        pwalk_init(&c, p, n, 1.0, 99999ULL, 4);
        int same = 1, differs = 0;
        for (int i = 0; i < 64; i++) {
            unsigned long long sa = pwalk_step(&a);
            unsigned long long sb = pwalk_step(&b);
            unsigned long long sc = pwalk_step(&c);
            if (sa != sb) same = 0;
            if (sa != sc) differs = 1;
        }
        check("same seed => identical trajectory", same);
        check("different seed => trajectory diverges", differs);
    }

    printf(fails ? "\n%d check(s) FAILED\n" : "\nAll checks passed.\n", fails);
    return fails ? 1 : 0;
}
