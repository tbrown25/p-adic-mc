// tests/test_diffusion.c — ensemble diffusion: escape-time hierarchy, the
// log-periodic scale p^alpha, and relaxation to the equilibrium closed forms.
#include "padic/diffusion.h"
#include "padic/tree.h"     // padic_shell_count

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int fails = 0;
static void check(const char *name, int ok) {
    printf("  [%s] %s\n", ok ? "PASS" : "FAIL", name);
    if (!ok) fails++;
}
static int close_rel(double a, double b, double rel) {   // |a-b| <= rel*|b| (+tiny abs)
    return fabs(a - b) <= rel * fabs(b) + 1e-9;
}

int main(void) {
    printf("p-adic diffusion — escape hierarchy, log-periodic scale, equilibrium\n\n");

    // ── first-passage escape time is exactly geometric: MC mean == 1/p_esc(m) ─
    {
        int p = 3, n = 8;
        double alpha = 1.0;
        // m = 0: any jump escapes the single start leaf, so it always takes 1 step.
        check("escape(m=0) == 1 step exactly",
              close_rel(padic_mean_escape_steps(p, n, alpha, 0, 20000, 0x1), 1.0, 1e-9)
              && close_rel(padic_escape_time(p, n, alpha, 0), 1.0, 1e-12));
        // m = 1, 2: MC mean matches the closed form 1/p_esc(m) within ~1/sqrt(M).
        double mc1 = padic_mean_escape_steps(p, n, alpha, 1, 60000, 0xA1);
        double mc2 = padic_mean_escape_steps(p, n, alpha, 2, 60000, 0xB2);
        check("escape(m=1) MC matches 1/p_esc(1)", close_rel(mc1, padic_escape_time(p, n, alpha, 1), 0.04));
        check("escape(m=2) MC matches 1/p_esc(2)", close_rel(mc2, padic_escape_time(p, n, alpha, 2), 0.05));
    }

    // ── edge cases: no underflow, no hang, no overflow ───────────────────────
    {
        // huge alpha: escape from level 0 is always 1 step (relative weights, no 0/0 NaN)
        check("escape_time(m=0) == 1 even for huge alpha",
              close_rel(padic_escape_time(3, 6, 700.0, 0), 1.0, 1e-9)
              && close_rel(padic_escape_prob(3, 6, 700.0, 0), 1.0, 1e-9));
        // infeasible escape (p_esc below fp resolution) returns the closed form
        // instead of looping forever: p=2,n=60,m=55 is an ORDINARY-alpha wedge.
        double inf_mc = padic_mean_escape_steps(2, 60, 1.0, 55, 100, 0x9);
        check("infeasible escape falls back to closed form (no hang)",
              isfinite(inf_mc) && close_rel(inf_mc, padic_escape_time(2, 60, 1.0, 55), 1e-9));
        // equilibrium closed forms stay correct where the u64 shell count would overflow
        check("eq_abs(p=5,n=30) ~ (p-1)/p / (1-p^-2) (no u64 overflow)",
              close_rel(padic_eq_abs(5, 30), 0.8 / 0.96, 1e-4));
    }

    // ── the escape times march by p^alpha (discrete scale invariance) ─────────
    // Away from the top of a finite tree the ratio converges to p^alpha; the
    // last few levels carry an O(p^{-alpha(n-m)}) boundary correction, so we
    // check the deep levels (gap n-m >= 6, where the correction is < a few %).
    {
        int n = 12;
        struct { int p; double alpha; } cs[] = {{3, 1.0}, {2, 1.0}, {3, 0.5}, {2, 1.5}};
        int ok = 1;
        for (int ci = 0; ci < 4; ci++) {
            int p = cs[ci].p; double a = cs[ci].alpha, scale = pow((double)p, a);
            for (int m = 1; m <= n - 6; m++) {
                double r = padic_escape_time(p, n, a, m + 1) / padic_escape_time(p, n, a, m);
                if (!close_rel(r, scale, 0.06)) ok = 0;
            }
        }
        check("tau(m+1)/tau(m) -> p^alpha  (p in {2,3}, alpha in {0.5,1,1.5})", ok);
    }

    // ── relaxation to the equilibrium closed forms (uniform over p^n leaves) ──
    {
        int p = 3, n = 6, m0 = 3;
        double alpha = 1.0;
        int T = 5000;
        long M = 20000;
        double *sep = calloc(T + 1, sizeof(double));
        double *absp = calloc(T + 1, sizeof(double));
        double *msd = calloc(T + 1, sizeof(double));
        double *surv = calloc(T + 1, sizeof(double));
        if (!sep || !absp || !msd || !surv) { printf("  [FAIL] allocation\n"); return 1; }
        padic_diffuse(p, n, alpha, M, T, m0, 0xD1FF, sep, absp, msd, surv);

        // late-time average (second half, well past tau_{n-1} ~ p^{alpha(n-1)});
        // each timepoint is an average over M independent walkers.
        double es = padic_eq_sep(p, n);
        double a_sep = 0, a_abs = 0, a_msd = 0, a_surv = 0, max_dev = 0;
        int lo = T / 2, cnt = 0;
        for (int t = lo; t <= T; t++) {
            a_sep += sep[t]; a_abs += absp[t]; a_msd += msd[t]; a_surv += surv[t]; cnt++;
            double d = fabs(sep[t] - es);
            if (d > max_dev) max_dev = d;
        }
        a_sep /= cnt; a_abs /= cnt; a_msd /= cnt; a_surv /= cnt;

        // Independence guard: with M independent walkers each sep[t] sits within
        // ~sigma/sqrt(M) (a few x 0.01) of equilibrium. A correlated ensemble (the
        // subseed-lattice bug) behaves like ~one walker and swings by ~sigma ~ 0.5+.
        check("equilibrated ensemble is independent (sep[t] never swings from eq)", max_dev < 0.15);

        // equilibrium survival closed form: leaves within separation m0 over p^n
        double leaves_le_m0 = 1.0;   // the origin itself (s = 0)
        for (int s = 1; s <= m0; s++) leaves_le_m0 += (double)padic_shell_count(p, n, s);
        double eq_surv = leaves_le_m0 / pow((double)p, (double)n);

        check("<s(t)> relaxes to <s>_eq",       close_rel(a_sep, padic_eq_sep(p, n), 0.03));
        check("<|X|_p> relaxes to <|X|_p>_eq",  close_rel(a_abs, padic_eq_abs(p, n), 0.03));
        check("MSD relaxes to <|X|_p^2>_eq",    close_rel(a_msd, padic_eq_msd(p, n), 0.04));
        check("ball occupation relaxes to eq",  close_rel(a_surv, eq_surv, 0.06));

        // and it genuinely relaxes FROM the origin (not already at equilibrium)
        int relaxed = (sep[0] == 0.0) && (absp[0] == 0.0) && (surv[0] == 1.0)
                      && (sep[T] > sep[1]) && (absp[T] > absp[1]) && (surv[T] < surv[1]);
        check("starts at the origin and spreads toward equilibrium", relaxed);

        free(sep); free(absp); free(msd); free(surv);
    }

    printf(fails ? "\n%d check(s) FAILED\n" : "\nAll checks passed.\n", fails);
    return fails ? 1 : 0;
}
