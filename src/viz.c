// src/viz.c — the three p-adic figures, drawn from C into an svg canvas.
#include "padic/viz.h"

#include "padic/tree.h"          // pwalk_*
#include "padic/diffusion.h"     // padic_diffuse, padic_eq_sep
#include "padic/montecarlo.h"    // padic_mc_integral, padic_zeta_integral
#include "padic/padic.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const char *INK = "#c9d4e3", *FAINT = "#39404d", *AXIS = "#6b7686";

static long ipow(int p, int e) { long r = 1; while (e-- > 0) r *= p; return r; }

// prefix_val of the first d LOW digits of x (a[0] most significant in the layout):
// sum_{i<d} a_i p^{d-1-i}. At d=n this is the leaf's left-to-right layout index.
static long prefix_val(long x, int p, int d) {
    long j = 0, t = x;
    for (int i = 0; i < d; i++) { int a = (int)(t % p); t /= p; j += (long)a * ipow(p, d - 1 - i); }
    return j;
}

// blue (separation 1, local) -> red (separation n, root-deep)
static void jump_color(int sj, int n, char *out) {
    double t = (n > 1) ? (double)(sj - 1) / (double)(n - 1) : 0.0;
    int r = (int)(70 + t * 175), g = (int)(130 - t * 70), b = (int)(210 - t * 150);
    snprintf(out, 8, "#%02x%02x%02x", r, g, b);
}

// ── figure 1: the tree of Z/p^n balls + one ultrametric walk ─────────────────
void padic_viz_tree(svg *s, int p, int n, double alpha, int steps, unsigned long long seed) {
    double mx = 44, top = 56, bot = 30;
    double plotw = s->w - 2 * mx, ploth = s->h - top - bot;
    long leaves = ipow(p, n);
    double leafsp = plotw / (double)leaves, levelh = ploth / (double)n;
    #define NX(d, j) (mx + ((double)(j) + 0.5) * (double)ipow(p, n - (d)) * leafsp)
    #define NY(d)    (top + (double)(d) * levelh)

    svg_text(s, mx, 30, 15, INK, "start",
             "p-adic tree of ℤ/pⁿ balls  —  one ultrametric walk (jumps up to the common ancestor)");

    // tree edges (each node to its parent), then faint leaf dots
    for (int d = 1; d <= n; d++)
        for (long j = 0; j < ipow(p, d); j++)
            svg_line(s, NX(d, j), NY(d), NX(d - 1, j / p), NY(d - 1), FAINT, 1.0);
    for (long j = 0; j < leaves; j++) svg_circle(s, NX(n, j), NY(n), 1.6, FAINT, NULL, 0);

    // the walk
    long *val = malloc((size_t)(steps + 1) * sizeof(long));
    int *sepv = malloc((size_t)(steps + 1) * sizeof(int));
    if (!val || !sepv) { free(val); free(sepv); return; }
    pwalk w;
    pwalk_init(&w, p, n, alpha, seed, 0);
    val[0] = 0;
    for (int t = 1; t <= steps; t++) { val[t] = (long)pwalk_step(&w); sepv[t] = pwalk_last_sep(&w); }

    char col[8];
    double px[2 * PADIC_MAXK + 2], py[2 * PADIC_MAXK + 2];
    for (int t = 1; t <= steps; t++) {
        long a = val[t - 1], b = val[t];
        int c = n - sepv[t];                       // common-ancestor depth
        int m = 0;
        for (int d = n; d >= c; d--) { px[m] = NX(d, prefix_val(a, p, d)); py[m] = NY(d); m++; }
        for (int d = c + 1; d <= n; d++) { px[m] = NX(d, prefix_val(b, p, d)); py[m] = NY(d); m++; }
        jump_color(sepv[t], n, col);
        svg_polyline(s, px, py, m, col, 1.7, 0.55);
    }
    svg_circle(s, NX(n, prefix_val(0, p, n)), NY(n), 4, "#ffd54a", "#20242c", 0.8);  // start leaf

    free(val); free(sepv);
    #undef NX
    #undef NY
}

// ── figure 2: the log-periodic staircase <s(t)> ──────────────────────────────
void padic_viz_relaxation(svg *s, int p, int n, double alpha, long walkers, int T) {
    double mx = 56, top = 56, bot = 48;
    double plotw = s->w - mx - 24, ploth = s->h - top - bot;
    double *sep = calloc((size_t)T + 1, sizeof(double));
    double *px = malloc((size_t)T * sizeof(double)), *py = malloc((size_t)T * sizeof(double));
    if (!sep || !px || !py) { free(sep); free(px); free(py); return; }
    padic_diffuse(p, n, alpha, walkers, T, n / 2, 0xF16A5EEDULL, sep, NULL, NULL, NULL);

    double lT = log((double)T);
    #define XT(t) (mx + (log((double)(t)) / lT) * plotw)
    #define YS(v) (top + ploth - ((v) / (double)n) * ploth)

    svg_text(s, mx, 30, 15, INK, "start",
             "log-periodic staircase: mean tree level ⟨s(t)⟩ climbs ~1 level per factor pᵅ");
    svg_line(s, mx, top, mx, top + ploth, AXIS, 1.2);            // y axis
    svg_line(s, mx, top + ploth, mx + plotw, top + ploth, AXIS, 1.2);  // x axis
    svg_text(s, mx - 8, YS(0) + 4, 11, AXIS, "end", "0");
    svg_text(s, mx - 8, YS(n) + 4, 11, AXIS, "end", "n");

    char lab[32];
    for (long t = 1; t <= T; t *= p) {                          // gridlines at t = p^k
        svg_line(s, XT(t), top, XT(t), top + ploth, FAINT, 0.8);
        snprintf(lab, sizeof lab, "%ld", t);
        svg_text(s, XT(t), top + ploth + 16, 10, AXIS, "middle", lab);
    }
    double eqs = padic_eq_sep(p, n);                            // equilibrium level
    svg_line(s, mx, YS(eqs), mx + plotw, YS(eqs), "#5cc28d", 1.0);
    svg_text(s, mx + plotw, YS(eqs) - 5, 10, "#5cc28d", "end", "⟨s⟩ (eq)");

    int m = 0;
    for (int t = 1; t <= T; t++) { px[m] = XT(t); py[m] = YS(sep[t]); m++; }
    svg_polyline(s, px, py, m, "#5aa9ff", 1.8, 1.0);

    free(sep); free(px); free(py);
    #undef XT
    #undef YS
}

// ── figure 3: Monte-Carlo N^-1/2 convergence (log-log) ───────────────────────
void padic_viz_convergence(svg *s, int p, int k, double sexp, long Nmax, unsigned long long seed) {
    double mx = 60, top = 56, bot = 48;
    double plotw = s->w - mx - 24, ploth = s->h - top - bot;

    long Ns[64];
    double err[64];
    int cnt = 0;
    for (long N = 10; N <= Nmax && cnt < 64; ) {               // ~half-decade steps
        Ns[cnt++] = N;
        long next = (long)llround((double)N * 3.1622776601683795);
        N = (next > N) ? next : N + 1;
    }
    double I = padic_zeta_integral(p, sexp);
    int R = 24;
    unsigned long long master = seed;
    for (int i = 0; i < cnt; i++) {
        double se2 = 0.0;
        for (int r = 0; r < R; r++) {
            unsigned long long z = (master += 0x9E3779B97F4A7C15ULL);   // splitmix64 finalizer -> independent seeds
            z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
            z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
            z ^= z >> 31;
            double est = padic_mc_integral(p, k, sexp, Ns[i], z, NULL);
            se2 += (est - I) * (est - I);
        }
        err[i] = sqrt(se2 / R);
    }

    if (cnt < 2) {   // too few N samples to span a log-log axis (Nmax < ~32)
        svg_text(s, mx, 30, 15, INK, "start", "Monte-Carlo convergence — need a larger N range");
        svg_line(s, mx, top, mx, top + ploth, AXIS, 1.2);
        svg_line(s, mx, top + ploth, mx + plotw, top + ploth, AXIS, 1.2);
        return;
    }

    double xlo = log10((double)Ns[0]), xhi = log10((double)Ns[cnt - 1]);
    double ylo = log10(err[cnt - 1]), yhi = log10(err[0]);      // err shrinks with N
    ylo -= 0.15; yhi += 0.15;
    #define XN(N) (mx + (log10((double)(N)) - xlo) / (xhi - xlo) * plotw)
    #define YE(e) (top + (yhi - log10(e)) / (yhi - ylo) * ploth)

    svg_text(s, mx, 30, 15, INK, "start",
             "Monte-Carlo of the ζ integral: RMS error ∝ 1/√N  (slope −½, log-log)");
    svg_line(s, mx, top, mx, top + ploth, AXIS, 1.2);
    svg_line(s, mx, top + ploth, mx + plotw, top + ploth, AXIS, 1.2);
    svg_text(s, mx + plotw, top + ploth + 18, 11, AXIS, "end", "N (log)");
    svg_text(s, mx, top - 8, 11, AXIS, "start", "RMS error (log)");

    // slope -1/2 reference through the first point
    double rx[2] = {XN(Ns[0]), XN(Ns[cnt - 1])};
    double ry[2] = {YE(err[0]), YE(err[0] * sqrt((double)Ns[0] / (double)Ns[cnt - 1]))};
    svg_polyline(s, rx, ry, 2, "#5cc28d", 1.4, 1.0);

    double px[64], py[64];
    for (int i = 0; i < cnt; i++) { px[i] = XN(Ns[i]); py[i] = YE(err[i]); }
    svg_polyline(s, px, py, cnt, "#5aa9ff", 1.6, 1.0);
    for (int i = 0; i < cnt; i++) svg_circle(s, px[i], py[i], 3, "#ffd54a", "#20242c", 0.6);

    #undef XN
    #undef YE
}
