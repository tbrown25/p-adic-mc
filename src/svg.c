// src/svg.c — the tiny SVG writer.
#include "padic/svg.h"

#include <math.h>

// Record a point in the bounding box; return 0 (and flag bad_coord) if non-finite.
static int track(svg *s, double x, double y) {
    if (!isfinite(x) || !isfinite(y)) { s->bad_coord = 1; return 0; }
    if (x < s->minx) s->minx = x;
    if (x > s->maxx) s->maxx = x;
    if (y < s->miny) s->miny = y;
    if (y > s->maxy) s->maxy = y;
    return 1;
}

void svg_begin(svg *s, FILE *f, double w, double h, const char *bg) {
    s->f = f;
    s->w = w;
    s->h = h;
    s->elements = 0;
    s->bad_coord = 0;
    s->minx = s->miny = 1e300;
    s->maxx = s->maxy = -1e300;
    fprintf(f, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%.0f\" height=\"%.0f\" "
               "viewBox=\"0 0 %.0f %.0f\">\n", w, h, w, h);
    if (bg) fprintf(f, "<rect x=\"0\" y=\"0\" width=\"%.0f\" height=\"%.0f\" fill=\"%s\"/>\n", w, h, bg);
}

void svg_end(svg *s) { fprintf(s->f, "</svg>\n"); }

void svg_line(svg *s, double x1, double y1, double x2, double y2, const char *stroke, double width) {
    if (!track(s, x1, y1) || !track(s, x2, y2)) return;
    fprintf(s->f, "<line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\" stroke=\"%s\" stroke-width=\"%.2f\"/>\n",
            x1, y1, x2, y2, stroke, width);
    s->elements++;
}

void svg_rect(svg *s, double x, double y, double w, double h, const char *fill) {
    if (!track(s, x, y) || !track(s, x + w, y + h)) return;
    fprintf(s->f, "<rect x=\"%.2f\" y=\"%.2f\" width=\"%.2f\" height=\"%.2f\" fill=\"%s\"/>\n", x, y, w, h, fill);
    s->elements++;
}

void svg_circle(svg *s, double cx, double cy, double r, const char *fill, const char *stroke, double sw) {
    if (!track(s, cx - r, cy - r) || !track(s, cx + r, cy + r)) return;  // bbox includes the disk extent
    fprintf(s->f, "<circle cx=\"%.2f\" cy=\"%.2f\" r=\"%.2f\" fill=\"%s\"", cx, cy, r, fill ? fill : "none");
    if (stroke) fprintf(s->f, " stroke=\"%s\" stroke-width=\"%.2f\"", stroke, sw);
    fprintf(s->f, "/>\n");
    s->elements++;
}

void svg_polyline(svg *s, const double *xs, const double *ys, int n, const char *stroke, double width, double opacity) {
    if (n < 2) return;
    for (int i = 0; i < n; i++)                 // all-or-nothing skip, like the other primitives
        if (!isfinite(xs[i]) || !isfinite(ys[i])) { s->bad_coord = 1; return; }
    fprintf(s->f, "<polyline fill=\"none\" stroke=\"%s\" stroke-width=\"%.2f\"", stroke, width);
    if (opacity < 1.0) fprintf(s->f, " stroke-opacity=\"%.2f\"", opacity);
    fprintf(s->f, " points=\"");
    for (int i = 0; i < n; i++) {
        track(s, xs[i], ys[i]);
        fprintf(s->f, "%s%.2f,%.2f", i ? " " : "", xs[i], ys[i]);
    }
    fprintf(s->f, "\"/>\n");
    s->elements++;
}

// escape the XML metacharacters so a label can't break the document (UTF-8 bytes pass through)
static void put_escaped(FILE *f, const char *t) {
    for (; *t; t++) {
        if (*t == '&') fputs("&amp;", f);
        else if (*t == '<') fputs("&lt;", f);
        else if (*t == '>') fputs("&gt;", f);
        else fputc(*t, f);
    }
}

void svg_text(svg *s, double x, double y, double size, const char *fill, const char *anchor, const char *txt) {
    if (!track(s, x, y)) return;
    fprintf(s->f, "<text x=\"%.2f\" y=\"%.2f\" font-size=\"%.1f\" fill=\"%s\" text-anchor=\"%s\" "
                  "font-family=\"sans-serif\">", x, y, size, fill, anchor ? anchor : "start");
    put_escaped(s->f, txt);
    fprintf(s->f, "</text>\n");
    s->elements++;
}
