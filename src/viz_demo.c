// src/viz_demo.c — render the three p-adic figures to SVG files.
//
//   ./build/viz_demo [out_dir]     (default ".")
//
// Writes tree.svg, relaxation.svg, convergence.svg. Deterministic (fixed seeds).
#include "padic/viz.h"
#include "padic/svg.h"

#include <stdio.h>

static long emit(const char *path, int which) {
    FILE *f = fopen(path, "w");
    if (!f) { perror(path); return -1; }
    svg s;
    svg_begin(&s, f, 900, 560, "#12151b");
    if (which == 0)      padic_viz_tree(&s, 2, 6, 1.0, 48, 20250803ULL);        // 64-leaf binary tree + walk
    else if (which == 1) padic_viz_relaxation(&s, 3, 8, 1.0, 20000, 6000);      // the staircase
    else                 padic_viz_convergence(&s, 3, 25, 1.0, 100000, 424242ULL); // N^-1/2
    svg_end(&s);
    fclose(f);
    return s.bad_coord ? -1 : s.elements;
}

int main(int argc, char **argv) {
    const char *dir = argc > 1 ? argv[1] : ".";
    const char *names[3] = {"tree.svg", "relaxation.svg", "convergence.svg"};
    char path[512];
    for (int i = 0; i < 3; i++) {
        snprintf(path, sizeof path, "%s/%s", dir, names[i]);
        long e = emit(path, i);
        if (e < 0) { fprintf(stderr, "failed to write %s\n", path); return 1; }
        printf("wrote %s  (%ld elements)\n", path, e);
    }
    return 0;
}
