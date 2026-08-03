# p-adic-mc — build without CMake (just a C11 compiler).
#   make test     build + run all correctness checks (padic + tree + diffusion)
#   make walk     build the single-trajectory demo (writes CSV to stdout)
#   make diffuse  build the ensemble-diffusion demo (escape hierarchy + CSV)
CC     ?= cc
CFLAGS ?= -std=c11 -O2 -Wall -Wextra -Iinclude
LDLIBS ?= -lm

.PHONY: all test walk diffuse clean
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

build/walk_demo: src/padic.c src/tree.c src/walk_demo.c \
                 include/padic/padic.h include/padic/tree.h | build
	$(CC) $(CFLAGS) src/padic.c src/tree.c src/walk_demo.c $(LDLIBS) -o $@

build/diffusion_demo: src/padic.c src/tree.c src/diffusion.c src/diffusion_demo.c \
                 include/padic/padic.h include/padic/tree.h include/padic/diffusion.h | build
	$(CC) $(CFLAGS) src/padic.c src/tree.c src/diffusion.c src/diffusion_demo.c $(LDLIBS) -o $@

test: build/test_padic build/test_tree build/test_diffusion  ## build + run every check
	./build/test_padic
	@echo
	./build/test_tree
	@echo
	./build/test_diffusion

walk: build/walk_demo  ## single-trajectory demo
	@echo "built ./build/walk_demo — e.g. ./build/walk_demo 3 6 1.0 40 > walk.csv"

diffuse: build/diffusion_demo  ## ensemble-diffusion demo
	@echo "built ./build/diffusion_demo — e.g. ./build/diffusion_demo 3 8 1.0 20000 6000 3 > diffuse.csv"

clean:
	rm -rf build
