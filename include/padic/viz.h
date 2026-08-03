// padic/viz.h — the three p-adic figures, drawn into an svg canvas from C.
//
// Each scene draws into a caller-owned canvas (call svg_begin before / svg_end
// after), so a test can inspect the canvas's element count and bounding box.
//
// Milestone 5 of p-adic-mc. Depends on svg.h + tree.h + diffusion.h + montecarlo.h.
#ifndef PADIC_VIZ_H
#define PADIC_VIZ_H

#include "padic/svg.h"

// The depth-n p-ary tree of Z/p^n balls with one ultrametric walk overlaid: each
// jump is an arc up to the two leaves' common ancestor (height = separation) and
// back down, colored by separation. (Keep p^n modest — e.g. p=2,n=6 — to stay legible.)
void padic_viz_tree(svg *s, int p, int n, double alpha, int steps, unsigned long long seed);

// The log-periodic staircase: mean tree level <s(t)> vs t on a log-t axis, with
// vertical gridlines at t = p^k and the equilibrium <s>_eq marked.
void padic_viz_relaxation(svg *s, int p, int n, double alpha, long walkers, int T);

// Monte-Carlo convergence: RMS error of the zeta integral I(s) vs N on log-log
// axes, against a slope -1/2 reference line.
void padic_viz_convergence(svg *s, int p, int k, double sexp, long Nmax, unsigned long long seed);

#endif // PADIC_VIZ_H
