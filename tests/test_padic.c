// tests/test_padic.c — no framework: synthetic checks + an exit code.
#include "padic/padic.h"

#include <math.h>
#include <stdio.h>

static int fails = 0;
static void check(const char *name, int ok) {
    printf("  [%s] %s\n", ok ? "PASS" : "FAIL", name);
    if (!ok) fails++;
}

int main(void) {
    printf("p-adic (Z/p^k) — arithmetic, valuation, ultrametric\n\n");

    const int p = 7, k = 12;

    // round-trip through an integer
    {
        unsigned long long n = 1234567ULL;
        padic x = padic_from_ull(p, k, n);
        check("from_ull / to_ull round-trip", padic_to_ull(&x) == n);
    }

    // valuation and |.|_p of p^3, a unit, and 0
    {
        padic p3 = padic_from_ull(p, k, (unsigned long long)(p * p * p)); // p^3
        padic unit = padic_from_ull(p, k, 10);                            // 10 not divisible by 7
        padic z = padic_zero(p, k);
        check("v_p(p^3) == 3", padic_valuation(&p3) == 3);
        check("v_p(unit) == 0", padic_valuation(&unit) == 0);
        check("v_p(0) == k", padic_valuation(&z) == k);
        check("|p^3|_p == p^-3", fabs(padic_abs(&p3) - pow(p, -3)) < 1e-12);
        check("|unit|_p == 1", fabs(padic_abs(&unit) - 1.0) < 1e-12);
        check("|0|_p == 0", padic_abs(&z) == 0.0);
    }

    // carry: (p-1) + 1 == p
    {
        padic a = padic_from_ull(p, k, (unsigned long long)(p - 1));
        padic one = padic_from_ull(p, k, 1);
        padic s = padic_add(&a, &one);
        check("(p-1) + 1 carries to p", padic_to_ull(&s) == (unsigned long long)p);
    }

    // subtraction: x - x == 0, and (x - y) + y == x
    {
        padic x = padic_from_ull(p, k, 98765);
        padic y = padic_from_ull(p, k, 4321);
        padic xmx = padic_sub(&x, &x);
        padic back = padic_sub(&x, &y);
        back = padic_add(&back, &y);
        check("x - x == 0", padic_is_zero(&xmx));
        check("(x - y) + y == x", padic_equal(&back, &x));
    }

    // multiplication agrees with integer arithmetic (within precision)
    {
        padic a = padic_from_ull(p, k, 123);
        padic b = padic_from_ull(p, k, 456);
        padic ab = padic_mul(&a, &b);
        check("123 * 456 == 56088", padic_to_ull(&ab) == 56088ULL);
    }

    // Hensel headline: 1/(1-p) = 1 + p + p^2 + ... , i.e. (1-p) * S == 1 mod p^k
    {
        padic S = padic_zero(p, k);
        for (int i = 0; i < k; i++) S.a[i] = 1;          // the geometric series, all digits 1
        padic one = padic_from_ull(p, k, 1);
        padic pp = padic_from_ull(p, k, (unsigned long long)p);
        padic one_minus_p = padic_sub(&one, &pp);        // 1 - p in Z/p^k
        padic prod = padic_mul(&one_minus_p, &S);
        check("(1 - p) * (1 + p + p^2 + ...) == 1", padic_equal(&prod, &one));
    }

    // distance + the strong triangle (ultrametric) inequality
    {
        padic x = padic_from_ull(p, k, 1);
        padic y = padic_from_ull(p, k, 1 + p * p);       // differ by p^2
        check("dist(x, x) == 0", padic_dist(&x, &x) == 0.0);
        check("dist(1, 1+p^2) == p^-2", fabs(padic_dist(&x, &y) - pow(p, -2)) < 1e-12);

        padic a = padic_from_ull(p, k, 6);
        padic b = padic_from_ull(p, k, 7 * 7 * 3);       // v=2
        padic apb = padic_add(&a, &b);
        double lhs = padic_abs(&apb);
        double rhs = padic_abs(&a) > padic_abs(&b) ? padic_abs(&a) : padic_abs(&b);
        check("ultrametric: |a+b|_p <= max(|a|,|b|)", lhs <= rhs + 1e-12);
    }

    printf(fails ? "\n%d check(s) FAILED\n" : "\nAll checks passed.\n", fails);
    return fails ? 1 : 0;
}
