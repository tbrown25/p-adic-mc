// tests/test_montecarlo.c — Haar sampling on Z_p and the zeta integral. No framework.
#include "padic/montecarlo.h"
#include "padic/padic.h"

#include <math.h>
#include <stdio.h>

static int fails = 0;
static void check(const char *name, int ok) {
    printf("  [%s] %s\n", ok ? "PASS" : "FAIL", name);
    if (!ok) fails++;
}
static int close_rel(double a, double b, double rel) { return fabs(a - b) <= rel * fabs(b) + 1e-12; }

// splitmix64, used here only to mint independent run-seeds for the N^-1/2 test
static unsigned long long sm(unsigned long long *s) {
    unsigned long long z = (*s += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

int main(void) {
    printf("p-adic Monte-Carlo — Haar sampling and the zeta integral\n\n");

    // ── closed forms (I(1), I(2) also equal milestone 3's equilibrium moments) ─
    check("I(p,0) == 1",
          close_rel(padic_zeta_integral(3, 0.0), 1.0, 1e-12)
          && close_rel(padic_zeta_integral(2, 0.0), 1.0, 1e-12));
    check("I(3,1) == 3/4  (= milestone 3 <|X|_p>)",   close_rel(padic_zeta_integral(3, 1.0), 0.75, 1e-12));
    check("I(3,2) == 9/13 (= milestone 3 <|X|_p^2>)", close_rel(padic_zeta_integral(3, 2.0), 9.0 / 13.0, 1e-12));
    check("I(2,1) == 2/3",                            close_rel(padic_zeta_integral(2, 1.0), 2.0 / 3.0, 1e-12));

    // ── the Haar sampler: valuation law P(v=j) = (1-1/p)p^{-j}, uniform digits ─
    {
        int p = 3, k = 20;
        long N = 1000000;
        unsigned long long rng = 0xC0FFEEC0FFEEULL;
        long vc[PADIC_MAXK + 1] = {0};
        long dig[8] = {0}, digtot = 0;              // per-VALUE digit histogram (p <= 7 here)
        for (long i = 0; i < N; i++) {
            padic x;
            padic_haar_sample(&x, p, k, &rng);
            vc[padic_valuation(&x)]++;
            for (int d = 0; d < 5; d++) { dig[x.a[d]]++; digtot++; }
        }
        int ok = 1;
        for (int j = 0; j <= 4; j++) {
            double emp = (double)vc[j] / (double)N;
            double theo = (1.0 - 1.0 / p) * pow((double)p, -(double)j);
            if (fabs(emp - theo) > 6.0 * sqrt(theo * (1.0 - theo) / (double)N)) ok = 0;
        }
        check("Haar valuation law  P(v=j) = (1-1/p) p^{-j}", ok);
        // EVERY residue 0..p-1 must appear ~1/p (not just digit 0) — a sampler that
        // skips a nonzero value passes the valuation law but fails here.
        int ud = 1;
        double e = (double)digtot / (double)p;
        for (int v = 0; v < p; v++)
            if (fabs((double)dig[v] - e) > 6.0 * sqrt(e * (1.0 - 1.0 / p))) ud = 0;
        check("Haar digits uniform on {0..p-1}  P(digit=v) = 1/p", ud);
    }

    // ── MC estimate converges to the closed form (s = 0 exact) ───────────────
    {
        int p = 3, k = 25;
        long N = 1000000;
        double ss[] = {-0.4, 0.0, 1.0, 2.0, 3.0};   // -0.4 is the heavy-tail regime (still s > -1/2)
        int ok = 1;
        for (int j = 0; j < 5; j++) {
            double s = ss[j], se = -1.0;
            double est = padic_mc_integral(p, k, s, N, 0xC0DE01u + (unsigned)j * 911u, &se);
            double I = padic_zeta_integral(p, s);
            if (s == 0.0) { if (!(est == 1.0 && se == 0.0)) ok = 0; }
            else if (fabs(est - I) > 5.0 * se) ok = 0;
        }
        check("MC within 5*se of I(s), s in {-0.4,0,1,2,3}; s=0 exact", ok);
    }

    // ── the truncation bias is real and matched: E_k = I(s) + p^{-k(s+1)}(1-I(s)) ─
    // A small k makes the O(p^{-k(s+1)}) offset resolvable, pinning the valuation
    // cap / x=0 handling that k=25 hides.
    {
        int p = 3, k = 2;
        double s = 1.0;
        long N = 4000000;
        double se = -1.0;
        double mc = padic_mc_integral(p, k, s, N, 0xB1A5E5u, &se);
        double I = padic_zeta_integral(p, s);
        double bias = pow((double)p, -(double)k * (s + 1.0)) * (1.0 - I);
        check("small-k truncation bias matches I(s)+p^{-k(s+1)}(1-I(s))",
              fabs(mc - (I + bias)) < 5.0 * se && bias > 10.0 * se);
    }

    // ── the reported standard error matches theory sqrt(Var/N) ───────────────
    {
        int p = 3, k = 25;
        long N = 1000000;
        int ok = 1;
        for (double s = 1.0; s <= 2.0; s += 1.0) {
            double se = -1.0;
            padic_mc_integral(p, k, s, N, 0xABCDEFu, &se);
            if (!close_rel(se, sqrt(padic_zeta_variance(p, s) / (double)N), 0.10)) ok = 0;
        }
        check("reported standard error ~ sqrt(Var/N)", ok);
    }

    // ── the N^{-1/2} law: run-to-run scatter over R runs == sqrt(Var/N) ───────
    {
        int p = 3, k = 25;
        double s = 1.0;
        int R = 200;
        long Nr = 20000;
        unsigned long long master = 0x5EED1234ULL;
        double est[256], sum = 0.0;
        for (int r = 0; r < R; r++) {
            est[r] = padic_mc_integral(p, k, s, Nr, sm(&master), NULL);   // independent sub-seeds
            sum += est[r];
        }
        double mean = sum / R, var = 0.0;
        for (int r = 0; r < R; r++) { double d = est[r] - mean; var += d * d; }
        var /= (double)(R - 1);
        double emp_std = sqrt(var);
        double theo_std = sqrt(padic_zeta_variance(p, s) / (double)Nr);
        check("run-to-run std ~ sqrt(Var/N)  (the N^-1/2 law)", close_rel(emp_std, theo_std, 0.15));
        check("grand mean over R runs ~ I(s)",
              fabs(mean - padic_zeta_integral(p, s)) < 5.0 * theo_std / sqrt((double)R));
    }

    printf(fails ? "\n%d check(s) FAILED\n" : "\nAll checks passed.\n", fails);
    return fails ? 1 : 0;
}
