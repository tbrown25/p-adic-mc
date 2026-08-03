// padic/diffusion.h — ensemble diffusion of the ultrametric walk: the p-adic
// MSD(t) and the log-periodic, hierarchical relaxation that is the fingerprint
// of a tree geometry (rather than a line).
//
// An ensemble of independent walkers (tree.h) starts at the origin leaf 0. We
// watch, versus step-time t, the p-adic displacement d(t) = |X(t) - X(0)|_p =
// |X(t)|_p and the tree level s(t) = n - v_p(X(t)) it has reached.
//
// The exact backbone. A jump leaves the "level-m ball" {leaves within
// separation m of the start} IFF its own separation s' exceeds m — and that is
// independent of where the walker sits inside the ball and of its history. So
// the first-passage escape time from level m is exactly Geometric with success
// probability
//
//     p_esc(m) = sum_{s = m+1}^{n} P(s),   P(s) = p^{-alpha s} / Z,
//
// giving E[escape time] = 1 / p_esc(m) in closed form, and the geometric spacing
//
//     tau(m+1) / tau(m)  ->  p^alpha            (discrete scale invariance)
//
// which is the origin of the log-periodicity: timescales march by p^alpha, so
// the relaxation climbs a staircase that repeats every factor p^alpha in t. On
// a finite tree every observable relaxes to the uniform-over-p^n equilibrium,
// whose closed forms are below.
//
// Milestone 3 of p-adic-mc. Depends on tree.h / padic.h; dependency-free C.
#ifndef PADIC_DIFFUSION_H
#define PADIC_DIFFUSION_H

// ── exact closed forms ───────────────────────────────────────────────────────
double padic_escape_prob(int p, int n, double alpha, int m);  // P(jump separation > m); m in 0..n
double padic_escape_time(int p, int n, double alpha, int m);  // 1 / escape_prob = E[first-passage steps]

double padic_eq_sep(int p, int n);   // <s>_eq   over the uniform stationary law on p^n leaves
double padic_eq_abs(int p, int n);   // <|X|_p>_eq
double padic_eq_msd(int p, int n);   // <|X|_p^2>_eq  (the equilibrium of the p-adic MSD)

// ── Monte-Carlo measurements ─────────────────────────────────────────────────
// Mean first-passage steps for a walker from the origin until its separation
// from the start first exceeds m. (Should match padic_escape_time within ~1/sqrt(M).)
double padic_mean_escape_steps(int p, int n, double alpha, int m,
                               long walkers, unsigned long long seed);

// Run `walkers` independent walkers from the origin for T steps; fill each
// non-NULL array (length T+1, index t = 0..T) with the ensemble mean at step t:
//   sep[t]  = <s(t)>            mean tree level reached
//   absp[t] = <|X(t)|_p>        mean p-adic displacement
//   msd[t]  = <|X(t)|_p^2>      the p-adic MSD
//   surv[t] = P(s(t) <= m0)     occupation of the level-m0 ball
// Deterministic given `seed` (walker k draws a decorrelated sub-seed).
void padic_diffuse(int p, int n, double alpha, long walkers, int T, int m0,
                   unsigned long long seed,
                   double *sep, double *absp, double *msd, double *surv);

#endif // PADIC_DIFFUSION_H
