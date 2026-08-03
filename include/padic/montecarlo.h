// padic/montecarlo.h — Monte-Carlo over Z_p with the Haar measure (uniform base-p
// digits), and the canonical p-adic integral it estimates.
//
// The normalized Haar measure on Z_p gives each residue ball mod p^k measure
// p^{-k}, so a Haar sample truncated to precision k is exactly k independent
// uniform digits — a uniform element of Z/p^k. The textbook local zeta integral
//
//     I(s) = integral_{Z_p} |x|_p^s dx = (1 - p^{-1}) / (1 - p^{-(s+1)}),   s > -1,
//
// is the p-adic analog of "sample uniformly, average": estimate it as the mean of
// |x_i|_p^s over N Haar samples. I(1) = p/(p+1) and I(2) recover milestone 3's
// equilibrium <|X|_p> and <|X|_p^2>. The estimator is unbiased up to an
// O(p^{-k(s+1)}) truncation and converges at the Monte-Carlo rate N^{-1/2} with
// variance I(2s) - I(s)^2 (finite for s > -1/2).
//
// Milestone 4 of p-adic-mc. Depends on padic.h; dependency-free C.
#ifndef PADIC_MONTECARLO_H
#define PADIC_MONTECARLO_H

#include "padic/padic.h"

// ── exact closed forms ───────────────────────────────────────────────────────
double padic_zeta_integral(int p, double s);   // I(s) = (1-1/p)/(1-p^{-(s+1)}),  s > -1
double padic_zeta_variance(int p, double s);   // Var(|x|_p^s) = I(2s) - I(s)^2,  finite for s > -1/2

// ── Haar sampling and the Monte-Carlo estimator ──────────────────────────────
// Draw a Haar-uniform element of Z/p^k (k independent uniform digits) into `out`,
// advancing the caller's splitmix64 state `*rng`.
void padic_haar_sample(padic *out, int p, int k, unsigned long long *rng);

// Monte-Carlo estimate of I(s) from N Haar samples at precision k (requires
// s > -1). If `se` is not NULL, *se receives the estimated standard error
// (sample sd / sqrt(N)) — a valid error bar only for s > -1/2, where the variance
// is finite. Deterministic given `seed`.
double padic_mc_integral(int p, int k, double s, long N, unsigned long long seed, double *se);

#endif // PADIC_MONTECARLO_H
