// src/padic.c — Z/p^k arithmetic on base-p digit strings.
#include "padic/padic.h"

#include <math.h>

padic padic_zero(int p, int k) {
    padic r;
    r.p = p;
    r.k = k;
    for (int i = 0; i < k; i++) r.a[i] = 0;
    return r;
}

padic padic_from_ull(int p, int k, unsigned long long n) {
    padic r = padic_zero(p, k);
    for (int i = 0; i < k; i++) {          // repeated division by p peels off digits;
        r.a[i] = (int)(n % (unsigned)p);   // stopping at k digits reduces mod p^k for free
        n /= (unsigned)p;
    }
    return r;
}

unsigned long long padic_to_ull(const padic *x) {
    unsigned long long v = 0;              // Horner from the most significant digit down
    for (int i = x->k - 1; i >= 0; i--)
        v = v * (unsigned long long)x->p + (unsigned long long)x->a[i];
    return v;
}

padic padic_add(const padic *x, const padic *y) {
    padic r = padic_zero(x->p, x->k);
    int carry = 0;
    for (int i = 0; i < x->k; i++) {
        int s = x->a[i] + y->a[i] + carry;
        r.a[i] = s % x->p;
        carry = s / x->p;                  // carry LEFT; the final carry-out is dropped (mod p^k)
    }
    return r;
}

padic padic_neg(const padic *x) {
    // -x = p^k - x. Complement each digit to (p-1), then add 1 — the base-p
    // analog of two's complement. Works for x == 0 too (gives p^k ≡ 0).
    padic c = padic_zero(x->p, x->k);
    for (int i = 0; i < x->k; i++) c.a[i] = (x->p - 1) - x->a[i];
    padic one = padic_from_ull(x->p, x->k, 1);
    return padic_add(&c, &one);
}

padic padic_sub(const padic *x, const padic *y) {
    padic ny = padic_neg(y);
    return padic_add(x, &ny);
}

padic padic_mul(const padic *x, const padic *y) {
    long acc[PADIC_MAXK];
    for (int i = 0; i < x->k; i++) acc[i] = 0;
    for (int i = 0; i < x->k; i++)         // schoolbook convolution; terms of degree >= k drop
        for (int j = 0; i + j < x->k; j++)
            acc[i + j] += (long)x->a[i] * (long)y->a[j];
    padic r = padic_zero(x->p, x->k);
    long carry = 0;
    for (int i = 0; i < x->k; i++) {
        long s = acc[i] + carry;
        r.a[i] = (int)(s % x->p);
        carry = s / x->p;
    }
    return r;
}

int padic_valuation(const padic *x) {
    for (int i = 0; i < x->k; i++)
        if (x->a[i] != 0) return i;
    return x->k;                            // zero to this precision: valuation is "at least k"
}

double padic_abs(const padic *x) {
    int v = padic_valuation(x);
    if (v >= x->k) return 0.0;              // |0|_p = 0
    return pow((double)x->p, -v);
}

double padic_dist(const padic *x, const padic *y) {
    padic d = padic_sub(x, y);
    return padic_abs(&d);
}

int padic_is_zero(const padic *x) {
    for (int i = 0; i < x->k; i++)
        if (x->a[i] != 0) return 0;
    return 1;
}

int padic_equal(const padic *x, const padic *y) {
    if (x->p != y->p || x->k != y->k) return 0;
    for (int i = 0; i < x->k; i++)
        if (x->a[i] != y->a[i]) return 0;
    return 1;
}
