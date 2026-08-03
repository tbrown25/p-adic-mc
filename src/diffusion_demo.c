// src/diffusion_demo.c — ensemble diffusion: the escape-time hierarchy (the
// log-periodic backbone) plus a CSV relaxation time series.
//
//   ./build/diffusion_demo [p] [n] [alpha] [walkers] [T] [m0] [seed]
//
// stderr: parameters, the exact escape-time hierarchy tau(m) with tau(m+1)/tau(m)
//         approaching p^alpha, and the equilibrium closed forms.
// stdout: CSV  t, sep, absp, msd, surv  (ensemble means at each step) — redirect
//         to a .csv to plot the log-periodic staircase (milestone 5).
#include "padic/diffusion.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
    int p                 = argc > 1 ? atoi(argv[1]) : 3;
    int n                 = argc > 2 ? atoi(argv[2]) : 8;
    double alpha          = argc > 3 ? atof(argv[3]) : 1.0;
    long walkers          = argc > 4 ? atol(argv[4]) : 20000;
    int T                 = argc > 5 ? atoi(argv[5]) : 6000;
    int m0                = argc > 6 ? atoi(argv[6]) : 3;
    unsigned long long sd = argc > 7 ? strtoull(argv[7], NULL, 10) : 1618033989ULL;

    if (p < 2 || n < 2 || n > 40 || alpha <= 0.0 || walkers < 1 || T < 1 || m0 < 0 || m0 > n) {
        fprintf(stderr, "usage: diffusion_demo [p>=2] [2<=n<=40] [alpha>0] [walkers] [T>=1] [0<=m0<=n] [seed]\n");
        return 2;
    }

    fprintf(stderr, "# p-adic diffusion: p=%d n=%d alpha=%g walkers=%ld T=%d m0=%d seed=%llu\n",
            p, n, alpha, walkers, T, m0, sd);
    fprintf(stderr, "#\n# escape-time hierarchy (exact): tau(m)=1/P(jump sep>m); ratio -> p^alpha=%.4f\n",
            pow((double)p, alpha));
    fprintf(stderr, "#   m      tau(m)      tau(m)/tau(m-1)\n");
    double prev = 0.0;
    for (int m = 0; m < n; m++) {
        double tau = padic_escape_time(p, n, alpha, m);
        if (m == 0) fprintf(stderr, "#  %2d  %12.4f            -\n", m, tau);
        else        fprintf(stderr, "#  %2d  %12.4f       %8.4f\n", m, tau, tau / prev);
        prev = tau;
    }
    fprintf(stderr, "#\n# equilibrium (uniform over p^n leaves): <s>=%.4f  <|X|_p>=%.5f  <|X|_p^2>=%.5f\n",
            padic_eq_sep(p, n), padic_eq_abs(p, n), padic_eq_msd(p, n));

    double *sep  = calloc((size_t)T + 1, sizeof(double));
    double *absp = calloc((size_t)T + 1, sizeof(double));
    double *msd  = calloc((size_t)T + 1, sizeof(double));
    double *surv = calloc((size_t)T + 1, sizeof(double));
    if (!sep || !absp || !msd || !surv) { fprintf(stderr, "out of memory\n"); return 1; }

    padic_diffuse(p, n, alpha, walkers, T, m0, sd, sep, absp, msd, surv);

    printf("t,sep,absp,msd,surv\n");
    for (int t = 0; t <= T; t++)
        printf("%d,%.6g,%.6g,%.6g,%.6g\n", t, sep[t], absp[t], msd[t], surv[t]);

    free(sep); free(absp); free(msd); free(surv);
    return 0;
}
