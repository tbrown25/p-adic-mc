// src/walk_demo.c — dump one ultrametric random-walk trajectory as CSV.
//
//   ./build/walk_demo [p] [n] [alpha] [steps] [seed] [start]
//
// Columns: step, site, sep_prev, absp_prev, sep_start, absp_start.
// `sep_*` is the tree separation (levels up to the common ancestor); `absp_*` is
// the p-adic distance |Δ|_p. Redirect to a .csv for milestone 5's visualization.
#include "padic/tree.h"

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
    int p                 = argc > 1 ? atoi(argv[1]) : 3;
    int n                 = argc > 2 ? atoi(argv[2]) : 6;
    double alpha          = argc > 3 ? atof(argv[3]) : 1.0;
    int steps             = argc > 4 ? atoi(argv[4]) : 40;
    unsigned long long sd = argc > 5 ? strtoull(argv[5], NULL, 10) : 2718281828ULL;
    unsigned long long st = argc > 6 ? strtoull(argv[6], NULL, 10) : 0ULL;

    if (p < 2 || n < 1 || n > PADIC_MAXK || alpha <= 0.0) {
        fprintf(stderr, "usage: walk_demo [p>=2] [1<=n<=%d] [alpha>0] [steps] [seed] [start]\n",
                PADIC_MAXK);
        return 2;
    }

    fprintf(stderr, "# p-adic walk: p=%d n=%d alpha=%g steps=%d seed=%llu start=%llu\n",
            p, n, alpha, steps, sd, st);
    pwalk w;
    pwalk_init(&w, p, n, alpha, sd, st);
    pwalk_run(&w, steps, stdout, 1);
    return 0;
}
