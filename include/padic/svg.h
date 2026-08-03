// padic/svg.h — a tiny, dependency-free SVG writer.
//
// Just enough to draw the p-adic figures (the tree + walk, the log-periodic
// staircase, the Monte-Carlo convergence) straight from C. A canvas wraps an open
// FILE and, for testing, tracks how many primitives were drawn, the bounding box
// of the drawn content, and whether any non-finite coordinate was ever passed
// (such a primitive is skipped rather than written, so the SVG stays valid).
//
// Milestone 5 of p-adic-mc.
#ifndef PADIC_SVG_H
#define PADIC_SVG_H

#include <stdio.h>

typedef struct {
    FILE *f;
    double w, h;                    // viewport
    long elements;                  // primitives actually drawn
    int bad_coord;                  // set if any NaN/inf coordinate was passed (that primitive is skipped)
    double minx, miny, maxx, maxy;  // bounding box of drawn content
} svg;

void svg_begin(svg *s, FILE *f, double w, double h, const char *bg);  // bg NULL = transparent
void svg_end(svg *s);

void svg_line(svg *s, double x1, double y1, double x2, double y2, const char *stroke, double width);
void svg_rect(svg *s, double x, double y, double w, double h, const char *fill);
void svg_circle(svg *s, double cx, double cy, double r, const char *fill, const char *stroke, double sw);
void svg_polyline(svg *s, const double *xs, const double *ys, int n, const char *stroke, double width, double opacity);
void svg_text(svg *s, double x, double y, double size, const char *fill, const char *anchor, const char *txt);

#endif // PADIC_SVG_H
