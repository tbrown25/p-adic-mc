// tests/test_viz.c — the SVG writer and the three figures. No framework.
#include "padic/viz.h"
#include "padic/svg.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

static int fails = 0;
static void check(const char *name, int ok) {
    printf("  [%s] %s\n", ok ? "PASS" : "FAIL", name);
    if (!ok) fails++;
}
// all drawn content lies within the viewport (a real correctness property: a
// layout off-by-one would push coordinates out of frame)
static int inbox(const svg *s) {
    double e = 0.5;
    return s->minx >= -e && s->miny >= -e && s->maxx <= s->w + e && s->maxy <= s->h + e;
}

// render a scene into a text buffer so a test can inspect the emitted SVG
static void draw_tree(svg *s)  { padic_viz_tree(s, 2, 4, 1.0, 20, 12345ULL); }
static void draw_relax(svg *s) { padic_viz_relaxation(s, 3, 5, 1.0, 2000, 300); }
static void draw_conv(svg *s)  { padic_viz_convergence(s, 3, 15, 1.0, 3000, 999ULL); }

static size_t render(void (*draw)(svg *), char *buf, size_t cap) {
    FILE *f = tmpfile();
    svg s;
    svg_begin(&s, f, 900, 560, "#12151b");
    draw(&s);
    svg_end(&s);
    rewind(f);
    size_t n = fread(buf, 1, cap - 1, f);
    buf[n] = 0;
    fclose(f);
    return n;
}
static int count_sub(const char *b, const char *pat) {
    int c = 0;
    for (const char *q = b; (q = strstr(q, pat)); q++) c++;
    return c;
}

int main(void) {
    printf("p-adic viz — SVG writer and the three figures\n\n");

    // ── the SVG writer ───────────────────────────────────────────────────────
    {
        FILE *f = tmpfile();
        svg s;
        svg_begin(&s, f, 200, 100, "#000");     // bg rect is written directly, not counted
        svg_line(&s, 10, 10, 190, 90, "#fff", 1);
        svg_circle(&s, 100, 50, 5, "#f00", NULL, 0);
        double xs[3] = {20, 100, 180}, ys[3] = {80, 20, 80};
        svg_polyline(&s, xs, ys, 3, "#0f0", 1, 0.5);
        svg_rect(&s, 0, 0, 50, 50, "#123");
        svg_text(&s, 100, 50, 10, "#fff", "middle", "hi");
        long before = s.elements;
        check("five primitives counted", before == 5);
        check("no bad coordinates", s.bad_coord == 0);
        check("content within viewport", inbox(&s));

        svg_line(&s, NAN, 10, 20, 20, "#fff", 1);   // non-finite -> flagged + skipped
        check("non-finite coord flagged and skipped", s.bad_coord == 1 && s.elements == before);

        svg_end(&s);
        rewind(f);
        char buf[4096] = {0};
        size_t nb = fread(buf, 1, sizeof buf - 1, f);
        (void)nb;
        check("document starts with <svg", strncmp(buf, "<svg", 4) == 0);
        check("document contains </svg>", strstr(buf, "</svg>") != NULL);
        fclose(f);
    }

    // ── the three figures: well-formed, in-viewport, no bad coords ───────────
    {
        FILE *f = tmpfile();
        svg s;
        svg_begin(&s, f, 900, 560, "#12151b");
        padic_viz_tree(&s, 2, 4, 1.0, 20, 12345ULL);
        svg_end(&s);
        fclose(f);
        long edges = 0;
        for (int d = 1; d <= 4; d++) { long q = 1; for (int i = 0; i < d; i++) q *= 2; edges += q; }
        check("tree figure: >= tree edges, in-viewport, clean", s.elements >= edges && s.bad_coord == 0 && inbox(&s));
    }
    {
        FILE *f = tmpfile();
        svg s;
        svg_begin(&s, f, 900, 560, "#12151b");
        padic_viz_relaxation(&s, 3, 5, 1.0, 2000, 300);
        svg_end(&s);
        fclose(f);
        check("relaxation figure: nonempty, in-viewport, clean", s.elements > 0 && s.bad_coord == 0 && inbox(&s));
    }
    {
        FILE *f = tmpfile();
        svg s;
        svg_begin(&s, f, 900, 560, "#12151b");
        padic_viz_convergence(&s, 3, 15, 1.0, 3000, 999ULL);
        svg_end(&s);
        fclose(f);
        check("convergence figure: nonempty, in-viewport, clean", s.elements > 0 && s.bad_coord == 0 && inbox(&s));
    }

    // ── determinism: identical params/seed -> identical bytes ────────────────
    {
        static char a[1 << 17], b[1 << 17];
        FILE *f1 = tmpfile();
        svg s1; svg_begin(&s1, f1, 900, 560, "#12151b");
        padic_viz_tree(&s1, 2, 5, 1.0, 30, 77ULL); svg_end(&s1);
        rewind(f1); size_t na = fread(a, 1, sizeof a, f1); fclose(f1);

        FILE *f2 = tmpfile();
        svg s2; svg_begin(&s2, f2, 900, 560, "#12151b");
        padic_viz_tree(&s2, 2, 5, 1.0, 30, 77ULL); svg_end(&s2);
        rewind(f2); size_t nb = fread(b, 1, sizeof b, f2); fclose(f2);

        check("same seed => identical SVG bytes", na == nb && memcmp(a, b, na) == 0);
    }

    // ── tree geometry oracle: exact element mix + centered root ──────────────
    {
        static char buf[1 << 16];
        render(draw_tree, buf, sizeof buf);
        long edges = 0;
        for (int d = 1; d <= 4; d++) { long q = 1; for (int i = 0; i < d; i++) q *= 2; edges += q; }
        int okc = count_sub(buf, "<line") == edges
                  && count_sub(buf, "<circle") == 16 + 1     // 16 leaves + start marker
                  && count_sub(buf, "<polyline") == 20;      // one per jump
        check("tree: exact element counts (edges, leaves+1, jumps)", okc);
        double mx = 44, top = 56, plotw = 900 - 2 * mx;      // root at horizontal center, top row
        char want[64];
        snprintf(want, sizeof want, "x2=\"%.2f\" y2=\"%.2f\"", mx + plotw / 2, top);
        check("tree: edges converge on the centered root (layout scale/offset correct)", strstr(buf, want) != NULL);
    }

    // ── relaxation x-axis is LOG-t: gridlines at t=p^k are equally spaced ─────
    {
        static char buf[1 << 20];
        render(draw_relax, buf, sizeof buf);
        double gx[32];
        int ng = 0;
        for (const char *p = buf; (p = strstr(p, "<line")); p++) {
            double x1, y1, x2, y2;
            char strk[16] = {0};
            if (sscanf(p, "<line x1=\"%lf\" y1=\"%lf\" x2=\"%lf\" y2=\"%lf\" stroke=\"%15[^\"]\"",
                       &x1, &y1, &x2, &y2, strk) == 5
                && strcmp(strk, "#39404d") == 0 && fabs(x1 - x2) < 0.01 && ng < 32)
                gx[ng++] = x1;
        }
        for (int i = 0; i < ng; i++)
            for (int j = i + 1; j < ng; j++)
                if (gx[j] < gx[i]) { double t = gx[i]; gx[i] = gx[j]; gx[j] = t; }
        double dmin = 1e300, dmax = 0;
        for (int i = 1; i < ng; i++) { double d = gx[i] - gx[i - 1]; if (d < dmin) dmin = d; if (d > dmax) dmax = d; }
        // equal spacing => log-t axis; a linear axis would space t=1,p,p^2 ~p apart
        check("relaxation: gridlines equally spaced => log-t axis (not linear)", ng >= 4 && dmax / dmin < 1.25);
    }

    // ── the RNG-driven scenes are deterministic too ──────────────────────────
    {
        static char a2[1 << 20], b2[1 << 20];
        size_t na = render(draw_relax, a2, sizeof a2), nb = render(draw_relax, b2, sizeof b2);
        check("relaxation: deterministic (identical bytes)", na == nb && memcmp(a2, b2, na) == 0);
        size_t nc = render(draw_conv, a2, sizeof a2), nd = render(draw_conv, b2, sizeof b2);
        check("convergence: deterministic (identical bytes)", nc == nd && memcmp(a2, b2, nc) == 0);
    }

    printf(fails ? "\n%d check(s) FAILED\n" : "\nAll checks passed.\n", fails);
    return fails ? 1 : 0;
}
