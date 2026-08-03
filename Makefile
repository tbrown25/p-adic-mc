# p-adic-mc — build without CMake (just a C11 compiler). `make test` builds and
# runs the correctness checks for the p-adic number type.
CC     ?= cc
CFLAGS ?= -std=c11 -O2 -Wall -Wextra -Iinclude

.PHONY: all test clean
all: test

build/test_padic: src/padic.c tests/test_padic.c include/padic/padic.h
	@mkdir -p build
	$(CC) $(CFLAGS) src/padic.c tests/test_padic.c -lm -o $@

test: build/test_padic  ## build + run the p-adic checks
	./build/test_padic

clean:
	rm -rf build
