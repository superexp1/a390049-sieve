# OpenMP is deliberately NOT part of CXXFLAGS. CXXFLAGS is the variable users
# override, and `make CXXFLAGS=-O2` silently dropping -fopenmp produced a link
# failure that looks nothing like its cause. Override OMPFLAGS/OMPLIBS instead.
UNAME_S := $(shell uname -s)
UNAME_M := $(shell uname -m)

ifeq ($(UNAME_S),Darwin)
  # Apple's clang ships neither omp.h nor libomp and rejects -fopenmp outright;
  # Homebrew's libomp supplies both. -march=native is also not a valid argument
  # on arm64: "unsupported argument 'native' to option '-march='".
  ifeq ($(origin CXX),default)
    CXX := clang++
  endif
  BREW_OMP ?= $(shell brew --prefix libomp 2>/dev/null)
  ifeq ($(UNAME_M),arm64)
    ARCH   ?= -mcpu=native
  else
    ARCH   ?= -march=native
  endif
  OMPFLAGS ?= -Xpreprocessor -fopenmp -I$(BREW_OMP)/include
  OMPLIBS  ?= -L$(BREW_OMP)/lib -lomp
else
  CXX      ?= g++
  ARCH     ?= -march=native
  OMPFLAGS ?= -fopenmp
  OMPLIBS  ?=
endif

CXXFLAGS ?= -O3 $(ARCH) -Wall

# Homebrew GCC on macOS wants plain -fopenmp and no libomp:
#   make CXX=g++-14 OMPFLAGS=-fopenmp OMPLIBS=

all: sieve oddpart

sieve: src/sieve.cpp
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) $(LDFLAGS) -o $@ $< $(OMPLIBS)

oddpart: src/oddpart.cpp
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) $(LDFLAGS) -o $@ $< $(OMPLIBS)

# the v2 baseline, kept so the speedup claim in the README is checkable
bench/sieve-v2-baseline: bench/sieve-v2-baseline.cpp
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) $(LDFLAGS) -o $@ $< $(OMPLIBS)

# the baseline's BLOCK is a compile-time constant: it takes five arguments
bench: sieve bench/sieve-v2-baseline
	@echo "v2 baseline:"; ./bench/sieve-v2-baseline 17592000000000 17600000000000 $${T:-128} 262144 "" 2>&1 >/dev/null | tail -1
	@echo "current:   "; ./sieve            17592000000000 17600000000000 $${T:-128} 1048576 "" 8192 2>&1 >/dev/null | tail -1

KNOWN_TERMS = 40 208 928 3904 260608 1045504 16764928 268386304

# asserts, rather than printing and leaving the reader to count
check: sieve
	@echo "check: the eight known terms below 3e8"
	@got=$$(./sieve 2 300000000 $${T:-8} 1048576 "" 8192 2>/dev/null \
	        | awk '$$1=="TERM"{print $$2}' | sort -n | tr '\n' ' '); \
	 want=$$(printf '%s\n' $(KNOWN_TERMS) | sort -n | tr '\n' ' '); \
	 if [ "$$got" = "$$want" ]; then echo "check: PASS  $$got"; \
	 else echo "check: FAIL"; echo "  want: $$want"; echo "  got:  $$got"; exit 1; fi

# full regression suite; see AITESTING.md
test: sieve oddpart
	@tests/run-tests.sh

clean:
	rm -f sieve oddpart bench/sieve-v2-baseline

.PHONY: all bench check test clean
