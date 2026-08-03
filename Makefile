# p-adic-mc — build without CMake (just a C11 compiler).
#   make test     build + run all correctness checks (padic + tree + diffusion + montecarlo + viz)
#   make walk     single-trajectory demo (CSV to stdout)
#   make diffuse  ensemble-diffusion demo (escape hierarchy + CSV)
#   make mc       Haar Monte-Carlo demo (zeta integral + N^-1/2 convergence)
#   make viz      render the three figures to figures/*.svg
CC     ?= cc
CFLAGS ?= -std=c11 -O2 -Wall -Wextra -Iinclude
LDLIBS ?= -lm

CORE = src/padic.c src/tree.c src/diffusion.c src/montecarlo.c
HDRS = include/padic/padic.h include/padic/tree.h include/padic/diffusion.h include/padic/montecarlo.h

.PHONY: all test walk diffuse mc viz clean
all: test

build:
	@mkdir -p build

build/test_padic: src/padic.c tests/test_padic.c include/padic/padic.h | build
	$(CC) $(CFLAGS) src/padic.c tests/test_padic.c $(LDLIBS) -o $@

build/test_tree: src/padic.c src/tree.c tests/test_tree.c \
                 include/padic/padic.h include/padic/tree.h | build
	$(CC) $(CFLAGS) src/padic.c src/tree.c tests/test_tree.c $(LDLIBS) -o $@

build/test_diffusion: src/padic.c src/tree.c src/diffusion.c tests/test_diffusion.c \
                 include/padic/padic.h include/padic/tree.h include/padic/diffusion.h | build
	$(CC) $(CFLAGS) src/padic.c src/tree.c src/diffusion.c tests/test_diffusion.c $(LDLIBS) -o $@

build/test_montecarlo: src/padic.c src/montecarlo.c tests/test_montecarlo.c \
                 include/padic/padic.h include/padic/montecarlo.h | build
	$(CC) $(CFLAGS) src/padic.c src/montecarlo.c tests/test_montecarlo.c $(LDLIBS) -o $@

build/test_viz: $(CORE) src/svg.c src/viz.c tests/test_viz.c $(HDRS) \
                 include/padic/svg.h include/padic/viz.h | build
	$(CC) $(CFLAGS) $(CORE) src/svg.c src/viz.c tests/test_viz.c $(LDLIBS) -o $@

build/walk_demo: src/padic.c src/tree.c src/walk_demo.c \
                 include/padic/padic.h include/padic/tree.h | build
	$(CC) $(CFLAGS) src/padic.c src/tree.c src/walk_demo.c $(LDLIBS) -o $@

build/diffusion_demo: src/padic.c src/tree.c src/diffusion.c src/diffusion_demo.c \
                 include/padic/padic.h include/padic/tree.h include/padic/diffusion.h | build
	$(CC) $(CFLAGS) src/padic.c src/tree.c src/diffusion.c src/diffusion_demo.c $(LDLIBS) -o $@

build/montecarlo_demo: src/padic.c src/montecarlo.c src/montecarlo_demo.c \
                 include/padic/padic.h include/padic/montecarlo.h | build
	$(CC) $(CFLAGS) src/padic.c src/montecarlo.c src/montecarlo_demo.c $(LDLIBS) -o $@

build/viz_demo: $(CORE) src/svg.c src/viz.c src/viz_demo.c $(HDRS) \
                 include/padic/svg.h include/padic/viz.h | build
	$(CC) $(CFLAGS) $(CORE) src/svg.c src/viz.c src/viz_demo.c $(LDLIBS) -o $@

test: build/test_padic build/test_tree build/test_diffusion build/test_montecarlo build/test_viz  ## build + run every check
	./build/test_padic
	@echo
	./build/test_tree
	@echo
	./build/test_diffusion
	@echo
	./build/test_montecarlo
	@echo
	./build/test_viz

walk: build/walk_demo  ## single-trajectory demo
	@echo "built ./build/walk_demo — e.g. ./build/walk_demo 3 6 1.0 40 > walk.csv"

diffuse: build/diffusion_demo  ## ensemble-diffusion demo
	@echo "built ./build/diffusion_demo — e.g. ./build/diffusion_demo 3 8 1.0 20000 6000 3 > diffuse.csv"

mc: build/montecarlo_demo  ## Haar Monte-Carlo demo
	@echo "built ./build/montecarlo_demo — e.g. ./build/montecarlo_demo 3 30 2000000"

viz: build/viz_demo  ## render figures/*.svg
	@mkdir -p figures
	./build/viz_demo figures

clean:
	rm -rf build
