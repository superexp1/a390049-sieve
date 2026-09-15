# Testing

## Environment

No virtualenv or package manager is needed. The programs are two
single-translation-unit C++ files; the two helper scripts use only the Python 3
standard library (`sys`, `struct`), deliberately, so that the reference
implementation cannot drift with a dependency.

| requirement | version used | notes |
|---|---|---|
| GCC or Clang | g++ 11.4.0, clang 20.1.8 | needs `__int128` and OpenMP |
| OpenMP runtime | libgomp (GCC) / libomp (Clang) | see README for macOS |
| GNU make | 4.3 | |
| Python 3 | 3.10 | stdlib only, no packages |

Build the release binaries before testing — everything below tests `-O3`
builds, which is what actually runs:

    make clean && make

## Running the tests

    make check     # the eight known terms below 3e8; asserts, exits non-zero on failure
    make test      # the full suite (tests/run-tests.sh), ~2 minutes
    make bench     # v2 baseline vs current, ~1 minute on 32 cores

`T=` sets the thread count for all three, e.g. `make test T=16`. The suite
works in a `mktemp -d` directory and leaves nothing in the tree.

## What the suite covers

`tests/run-tests.sh` — 29 cases, each one a defect that was live before it.

**Correctness**

- the eight known terms below 3e8
- agreement with `tests/reference.py`, an independent smallest-prime-factor
  implementation that factors each n and evaluates the identity directly. It
  shares no code with the sieve; its only job is to disagree if the sieve is
  wrong.
- `lo = 0`, where record 0 has `rest = 0` and `hit()` would otherwise spin
- 20 window x block combinations must produce identical term lists, including
  block > window and blocks far from the default

**Checkpointing** — the bitmap is positional, so applying one to the wrong run
silently skipped unswept windows:

- a resume of the *same* run still works (positive control, checked first)
- a state file from a different range is refused, exit 3, with a reason
- a state file with a different window size is refused
- a different *block* size still resumes: it does not move window boundaries
- a truncated file and a random file are both refused

**Argument validation** — each of these previously wrapped, divided by zero, or
asked the allocator for something absurd:

- `hi <= lo`, `hi == lo`, window 0, block 0, block > 2^40
- thread counts of 0 and -4, for both binaries
- `oddpart` with M below the smallest odd part
- a malformed `SIEVE_CHECKPOINT_SECS` warns and falls back

**Checkpoint cadence** — the interval is a time budget, so every resume test
above finishes inside one interval and only ever exercises the final save.
These force the periodic path with `SIEVE_CHECKPOINT_SECS=0`:

- a periodic checkpoint is actually written mid-run
- a killed run leaves an incomplete but *valid* checkpoint, which then resumes
  to completion

**oddpart** — the odd parts of every term with `2^m - 3 < 2000`; agreement
with `tests/oddpart-reference.py`, which derives each candidate from the
mathematics alone (`X = (P + omega(k)) / (4P - 3Q - R)`) and shares no code
with the C++; and the window boundaries, where the odd-only indexing
(`m = (lo|1) + 2i`) is most likely to be wrong.

## Checking a long run

Do not trust the log; it is opened append and survives restarts.

    tools/state-status.py run.state

prints the range the file belongs to and how many of its windows are complete,
and exits non-zero unless all of them are. Pre-header state files need
`--legacy`, which reports the bit count but cannot name a range.

## Things the suite does not cover

- **The 24-bit bucket slot ceiling.** `sieve` now refuses to start when there
  are more than 2^24 primes below the root (hi > 9.6e16), but reaching that
  limit needs a ~310 MB prime sieve and is not exercised here.
- **Real macOS.** The Darwin branch of the Makefile is verified only by
  expanding it (`make -n UNAME_S=Darwin UNAME_M=arm64`) and by confirming that
  `clang --target=arm64-apple-macos12 -march=native` is the error it claims to
  be. Nothing in this repo has been compiled or run on Apple hardware.
- **Multi-day runs.** Resume is tested across a kill, not across days.
