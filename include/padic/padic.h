// padic/padic.h — p-adic integers in Z_p, to fixed precision.
//
// A p-adic integer truncated to precision k is a base-p digit string, held
// little-endian (a[0] is the units digit):
//
//     x = a[0] + a[1]*p + a[2]*p^2 + ... + a[k-1]*p^(k-1)   (mod p^k),
//     0 <= a[i] < p.
//
// This is the ring Z/p^k, the depth-k truncation of Z_p. Arithmetic carries to
// the LEFT (toward higher powers of p) and the top carry is dropped — that is
// what "mod p^k" means here. Digits give the p-adic valuation for free (the
// index of the first nonzero digit), and hence |x|_p and the ultrametric.
//
// Milestone 1 of p-adic-mc: the number type + arithmetic the random walk and
// Monte-Carlo layers will stand on. Dependency-free C; see the Makefile.
#ifndef PADIC_PADIC_H
#define PADIC_PADIC_H

#define PADIC_MAXK 64   // max precision (digits); keeps everything allocation-free

typedef struct {
    int p;                 // the prime
    int k;                 // precision (number of digits, 1..PADIC_MAXK)
    int a[PADIC_MAXK];     // digits, little-endian: a[0] = units
} padic;

// ── construction ─────────────────────────────────────────────────────────────
padic padic_zero(int p, int k);
padic padic_from_ull(int p, int k, unsigned long long n);  // n reduced mod p^k
unsigned long long padic_to_ull(const padic *x);           // value in [0, p^k) (needs p^k to fit in u64)

// ── arithmetic in Z/p^k  (operands must share p and k) ───────────────────────
padic padic_add(const padic *x, const padic *y);
padic padic_neg(const padic *x);                           // additive inverse: p^k - x
padic padic_sub(const padic *x, const padic *y);
padic padic_mul(const padic *x, const padic *y);

// ── p-adic structure ─────────────────────────────────────────────────────────
int    padic_valuation(const padic *x);                    // v_p: index of first nonzero digit; k if x == 0
double padic_abs(const padic *x);                          // |x|_p = p^(-v_p(x)); 0 when x == 0 to precision
double padic_dist(const padic *x, const padic *y);         // |x - y|_p  (the ultrametric)

// ── helpers ──────────────────────────────────────────────────────────────────
int padic_is_zero(const padic *x);
int padic_equal(const padic *x, const padic *y);

#endif // PADIC_PADIC_H
