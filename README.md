# a390049-sieve

Superexponential AI Labs, 
Mike Krygier

A segmented sieve that computes sigma, psi, phi and omega simultaneously at
**2.5 billion integers per second** on 64 cores, and an enumerator that searches
the same problem by odd part instead of by size.

Written to attack [OEIS A390049](https://oeis.org/A390049) — the numbers k with

    sigma(k) = psi(k) + phi(k) + omega(k)

where sigma is the sum of divisors, psi the Dedekind psi function, phi Euler's
totient and omega the number of distinct prime factors. The entry listed eight
terms and asked whether any exist beyond the family 2^m(2^m - 3) with 2^m - 3
prime.

## What it found

    a(9)  = 1099508482048  = 2^20 (2^20 - 3)
    a(10) = 17592173461504 = 2^22 (2^22 - 3)

and, by exhausting every integer below 17592173461504, that these are
consecutive — nothing lies between them, and nothing else lies below.

The sieve was validated by recovering the eight previously known terms, and
found a(9) without being told it existed.

## Two searches that cover different escapes

`sieve` is complete below a bound on k and says nothing above it. `oddpart` is
complete below a bound on the *odd part* and says nothing above that, but places
no ceiling on k at all — because the power of 2 is not a free parameter. Writing
k = 2^a m with m odd and X = 2^(a-1), and P, Q, R for sigma, psi, phi of m:

    X = (P + omega(k)) / (4P - 3Q - R),   omega(k) = omega(m) + 1

so each odd part admits at most one candidate k, of any magnitude. One division,
no search.

`oddpart` has run to odd part 10^11: every term there is 2^m(2^m - 3) with
2^m - 3 prime, exactly twelve of them, m in {3,4,5,6,9,10,12,14,20,22,24,29}.
That excludes a term hiding far above any sieve's reach with a small odd part.

A term escaping both would need k >= 1.76e13 **and** odd part > 10^11.

## Performance

Measured on a 64-core EPYC 7B13, 128 threads, n ~ 1.76e13:

 **2509.8** | sigma/psi/phi accumulated as products over prime powers |

Reproduce with `make bench` (set `T=` for thread count).


**Every division in the inner loop is exact.** After applying p, phi becomes
(n/p)(p-1), and any later prime q | n divides n/p, so it divides phi — likewise
psi and the unfactored remainder. So x/p is x * inv(p) mod 2^64, one multiply.
AVX2 has no integer-divide instruction, which is why division was the one thing
in the hot loop that could never vectorise. Worth 1.27x.

**Accumulate multiplicatively.** Carrying sigma, psi and phi as products over
prime powers, rather than as running n-scaled values, makes the final
combination three multiplies instead of two divisions by a per-element runtime
value — divisions that ran for the ~69% of integers having a prime factor above
the sieve root. Worth a further 1.42x.


## Build and run

    make                      # builds sieve and oddpart
    make check                # asserts the eight known terms below 3e8
    make test                 # full regression suite; see AITESTING.md

    ./sieve <lo> <hi> [threads] [window] [statefile] [block]
    ./oddpart <M> [threads] [window] [block]

Defaults are window 2^20 and block 8192 (256 KB, half the 512 KB L2 of the
machine above); both were swept, and block is flat from 8K to 32K *on that
machine*. Apple Silicon has a 4-16 MB L2, so re-sweep block there rather than
trusting the default.

Needs GCC or Clang with OpenMP:
`oddpart` uses `__int128` and both use `__builtin_expect`.

`sieve` exits 0 when done, 2 on a bad argument, 3 on an unusable state file,
and 4 if the final checkpoint could not be written.

### macOS



    brew install libomp
    make

With Homebrew GCC instead:

    brew install gcc
    make CXX=g++-14 OMPFLAGS=-fopenmp OMPLIBS=

OpenMP lives in `OMPFLAGS`, not `CXXFLAGS`, so overriding `CXXFLAGS` cannot
silently drop it.

## Checkpointing and resume

`sieve` checkpoints to `statefile` and resumes from it. The checkpoint is a
**per-window bitmap, not a count** — with `schedule(dynamic)` windows finish out
of order, so a count is not a low-water mark and resuming at one would silently
skip unprocessed windows.

The bitmap is written on a **time budget**, not every N windows. It grows with
the range — 1.9 MB for the a(10) sweep, 31.5 MB for a(11) — and is fsync'd
inside the critical section, so a fixed window cadence means very different
things at different sizes. On the a(11) range it cost 2.7x: 2686.6 M n/s with
no state file against 995.1 M n/s with one.

| variable | default | what it bounds |
|---|---|---|
| `SIEVE_CHECKPOINT_SECS` | 60 | checkpoint I/O, and the work a crash discards |
| `SIEVE_PROGRESS_SECS` | 5 | progress lines; at the old cadence the a(10) log reached 11 MB on a single line |

Because a crash replays up to one checkpoint interval, `<statefile>.terms` can
re-append terms it had already recorded: `sort -u` it before use.

**Do not use the log to decide a run finished.** It is opened append and
survives restarts, so one crash leaves a completion marker in it forever. Ask
the state file instead:

    tools/state-status.py run.state

which prints the range the file belongs to and how many of its windows are
complete, and exits non-zero unless every one of them is.

State files written before the header existed are bare bitmaps. `sieve` will
not resume from one — it cannot tell which range it describes, which is the
whole point — but `tools/state-status.py` still reads one and reports the bit
count, so the earlier campaigns remain checkable:

    tools/state-status.py --legacy old-run.state

## Accuracy of the claims here

The terms are exhaustive results, not samples. Completeness of the second
campaign came from the bitmap above (15632795 of 15632795 windows), and was read
independently twice before being believed. `oddpart`'s null results are backed
by a positive control: the same machinery run over *all* odd parts recovers all
nine then-known terms from odd parts below 1.1e6, with no size search, and
independently produced 2^22-3, 2^24-3 and 2^29-3.

What is **not** claimed: that 2^24(2^24 - 3) and 2^29(2^29 - 3) are a(11) and
a(12). They are terms, but nothing has tested the range below them.

## License

MIT — see LICENSE.
