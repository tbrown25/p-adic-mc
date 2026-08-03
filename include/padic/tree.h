// padic/tree.h — the p-adic tree over Z/p^n and one ultrametric random walk.
//
// The p^n residues mod p^n are the leaves of a rooted p-ary tree of depth n
// (the depth-n truncation of the tree of Z_p-balls). Two leaves x, y share a
// common ancestor at depth
//
//     c(x, y) = v_p(x - y)         (the number of matching LOW-order digits),
//
// and their separation in the tree is
//
//     s(x, y) = n - c(x, y)        (0 iff x == y, at most n at the root).
//
// The ultrametric distance is |x - y|_p = p^{-c} = p^{-(n-s)}, so a metric that
// decays in s is, up to a constant, a power of |x - y|_p.
//
// The walk is the embedded jump chain of the p-adic (Vladimirov) diffusion: from
// x it jumps to y != x with weight
//
//     w(x, y)  ∝  p^{-(alpha+1) · s(x,y)}  =  const · |x - y|_p^{-(alpha+1)},
//
// the discrete Vladimirov / alpha-stable Lévy kernel (alpha > 0 sets the tail).
// Because exactly (p-1)p^{s-1} leaves sit at separation s, the probability of
// landing in *shell* s collapses to the closed form
//
//     P(shell s)  ∝  p^{-alpha·s},        s = 1 .. n,
//
// which milestone 2 samples against. That structure also makes a step O(n)
// instead of O(p^n): pick a shell, then a uniform leaf inside it.
//
// Milestone 2 of p-adic-mc. Depends only on padic.h; dependency-free C.
#ifndef PADIC_TREE_H
#define PADIC_TREE_H

#include "padic/padic.h"

// ── the tree metric (x, y are leaves at the same precision n = their k) ───────
int padic_tree_ancestor_depth(const padic *x, const padic *y);  // c = v_p(x-y), capped at n
int padic_tree_sep(const padic *x, const padic *y);             // s = n - c  (0 iff x == y)
unsigned long long padic_shell_count(int p, int n, int s);      // #leaves at separation s: (p-1)p^{s-1}

// ── the ultrametric random walk ──────────────────────────────────────────────
typedef struct {
    int p;                          // the prime
    int n;                          // walk depth: sites live in Z/p^n (n <= PADIC_MAXK)
    double alpha;                   // kernel exponent (> 0): larger = lighter tails, more local
    unsigned long long rng;         // splitmix64 state (deterministic given the seed)
    padic pos;                      // current leaf (a padic of precision n)
    int last_sep;                   // separation of the most recent jump (1..n)
    double shell_cum[PADIC_MAXK + 1]; // cumulative P(shell <= s), s = 1..n; index 0 unused
} pwalk;

// Initialize a walker at `start` (reduced mod p^n) with a deterministic seed.
void pwalk_init(pwalk *w, int p, int n, double alpha, unsigned long long seed,
                unsigned long long start);

// Take one jump; returns the new site as an integer in [0, p^n). Never
// self-jumps. The integer id is faithful only when p^n <= 2^64; deeper trees
// still walk correctly on the digit representation, but the returned id wraps.
unsigned long long pwalk_step(pwalk *w);

// The separation of the last jump (set by pwalk_step): s in {1..n}. Undefined
// before the first step.
int pwalk_last_sep(const pwalk *w);

// Run `steps` jumps, writing one CSV row per step to `out` (with a header when
// `header` is nonzero): step, site, sep_from_prev, absp_from_prev,
// sep_from_start, absp_from_start. Pass out = NULL to walk silently.
void pwalk_run(pwalk *w, int steps, void *out, int header);

#endif // PADIC_TREE_H
