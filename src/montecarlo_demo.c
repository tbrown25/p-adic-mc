// src/montecarlo_demo.c — estimate the p-adic zeta integral by Haar Monte-Carlo
// and show the N^{-1/2} convergence against the closed form.
//
//   ./build/montecarlo_demo [p] [k] [N] [seed]
//
// Prints, for several s, the closed form I(s) vs the Monte-Carlo estimate with
// its standard error; then a convergence table for s=1 where err*sqrt(N) should
// hover around sqrt(Var) (the Monte-Carlo rate).
#include "padic/montecarlo.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
    int p                 = argc > 1 ? atoi(argv[1]) : 3;
    int k                 = argc > 2 ? atoi(argv[2]) : 30;
    long N                = argc > 3 ? atol(argv[3]) : 2000000;
    unsigned long long sd = argc > 4 ? strtoull(argv[4], NULL, 10) : 271828ULL;

    if (p < 2 || k < 1 || k > PADIC_MAXK || N < 1) {
        fprintf(stderr, "usage: montecarlo_demo [p>=2] [1<=k<=%d] [N>=1] [seed]\n", PADIC_MAXK);
        return 2;
    }

    printf("# p-adic Haar Monte-Carlo of  I(s) = integral_{Z_p} |x|_p^s dx = (1-1/p)/(1-p^-(s+1))\n");
    printf("# p=%d  precision k=%d  N=%ld  seed=%llu\n#\n", p, k, N, sd);
    printf("#    s     I(s) closed     I(s) MC        std err     |error|   error/se\n");
    double slist[] = {0.0, 0.5, 1.0, 2.0, 3.0};
    for (int j = 0; j < 5; j++) {
        double s = slist[j], se = 0.0;
        double mc = padic_mc_integral(p, k, s, N, sd + (unsigned long long)j * 1009ULL, &se);
        double I = padic_zeta_integral(p, s);
        double err = fabs(mc - I);
        printf("  %5.2f   %12.7f   %12.7f   %10.3g   %9.3g   %6.2f\n",
               s, I, mc, se, err, se > 0 ? err / se : 0.0);
    }

    printf("#\n# convergence at s=1: err should fall like N^-1/2, so err*sqrt(N) ~ sqrt(Var)=%.4f\n",
           sqrt(padic_zeta_variance(p, 1.0)));
    printf("#        N        |error|     err*sqrt(N)\n");
    for (long n = 100; n <= N; n *= 10) {
        double mc = padic_mc_integral(p, k, 1.0, n, sd + 77ULL, NULL);
        double err = fabs(mc - padic_zeta_integral(p, 1.0));
        printf("  %10ld   %10.3g   %10.4f\n", n, err, err * sqrt((double)n));
    }
    return 0;
}
