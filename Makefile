# p-adic-mc — build without CMake (just a C11 compiler).
#   make test   build + run all correctness checks (padic + tree)
#   make walk   build the trajectory-dumping demo (writes CSV to stdout)
CC     ?= cc
CFLAGS ?= -std=c11 -O2 -Wall -Wextra -Iinclude
LDLIBS ?= -lm

.PHONY: all test walk clean
all: test

build:
	@mkdir -p build

build/test_padic: src/padic.c tests/test_padic.c include/padic/padic.h | build
	$(CC) $(CFLAGS) src/padic.c tests/test_padic.c $(LDLIBS) -o $@

build/test_tree: src/padic.c src/tree.c tests/test_tree.c \
                 include/padic/padic.h include/padic/tree.h | build
	$(CC) $(CFLAGS) src/padic.c src/tree.c tests/test_tree.c $(LDLIBS) -o $@

build/walk_demo: src/padic.c src/tree.c src/walk_demo.c \
                 include/padic/padic.h include/padic/tree.h | build
	$(CC) $(CFLAGS) src/padic.c src/tree.c src/walk_demo.c $(LDLIBS) -o $@

test: build/test_padic build/test_tree  ## build + run every check
	./build/test_padic
	@echo
	./build/test_tree

walk: build/walk_demo  ## build the trajectory demo (run: ./build/walk_demo)
	@echo "built ./build/walk_demo — e.g. ./build/walk_demo 3 6 1.0 40 > walk.csv"

clean:
	rm -rf build
