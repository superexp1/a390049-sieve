// A390049 sieve, cache-blocked.  sigma(k) = psi(k) + phi(k) + omega(k)
//
// v1 kept an array-of-structs so an element cost one cache line instead of five,
// which took 28 M n/s (numpy) to 128 M n/s. It was then bound by DRAM: with a
// 2M-entry window nothing stays resident, so every one of the ~3 touches per
// integer pulled a 64-byte line for a 32-byte record, about 192 B/integer.
//
// v2 processes the window in L2-sized BLOCKS, applying every prime to a block
// while it is still hot. That only pays off if a block can skip the primes that
// miss it: at n ~ 1.7e13 there are ~295k primes below sqrt(n) and a 16k-entry
// block is hit by almost none of them. So primes are split --
//
//   small (p < BLOCK): hit most blocks, walked directly;
//   large (p >= BLOCK): hit at most once per block, so their hits are bucketed
//                       by block once per window and then drained.
//
// Checkpointing added: v1 had none, and this runs for days.
//
// v3 moves init and the final test inside the block loop. In v2 both swept the
// whole window, so the buffer was written once and read once through DRAM per
// window no matter how well the middle was blocked -- which is why v2 got FASTER
// as the window got smaller, the opposite of what blocking should buy. With them
// inside, a block is initialised, sieved and tested while resident.
//
// v4 then attacks the arithmetic. Every division hit() performs is exact, so it
// becomes a multiply by a modular inverse; and accumulating sigma/psi/phi as
// products over prime powers rather than as running n-scaled values removes the
// last two divisions by a per-element runtime value. AVX2 (which -march=native
// enables, and which GCC already uses for record init) has no integer-divide
// instruction, so division was the one thing in the inner loop that could never
// vectorise. 20 div in the object code became 7, none of them in the hot path.
//
// Measured on a 64-core EPYC 7B13, 128 threads, n ~ 1.76e13:
//
//   v2, its own best window/block            271.6 M n/s
//   v3, per-block init+test                 1380.2
//   v4a, + modular-inverse exact division   1779.0
//   v4b, + multiplicative accumulation      2509.8      9.2x over v2
//
// v5 is a robustness pass, no arithmetic change. The checkpoint file now
// carries a header identifying the run: it was a bare bitmap, and a bitmap is
// POSITIONAL, so resuming a different range (or a different window size) with
// an old file skipped windows that had never been sieved and still exited
// "complete". Any mismatch is now a hard error. Arguments are validated
// (hi<=lo, zero window, zero block and negative thread counts each used to
// wrap, divide by zero or abort), the 24-bit bucket slot field is checked
// against the prime count instead of assumed, and the done-bitmap read is
// atomic rather than racing the writers.
//
// Block 8192 is 256 KB, half the 512 KB L2; window 2^20. Every version recovers
// all eight known terms below 3e8.
//
//   g++ -O3 -march=native -fopenmp -o sieve src/sieve.cpp
//   ./sieve <lo> <hi> [threads] [window] [statefile] [block]
//
// v6: checkpoint and progress cadence are on a TIME budget, not a window count.
// The bitmap grows with the range -- 1.9 MB for the a(10) sweep but 31.5 MB for
// a(11) -- and a fixed 256-window cadence writes and fsyncs the whole thing
// inside the critical section, which at full speed is every 0.1 s. Measured on
// the a(11) range: 2686.6 M n/s with no state file against 995.1 M n/s with
// one, a 2.7x loss that grows with the range. A time budget bounds both the
// I/O and the work a crash discards, independent of how big the range is.
//
// exit 0 done, 2 bad arguments, 3 unusable state file, 4 checkpoint write failed
//
// env: SIEVE_CHECKPOINT_SECS (default 60) seconds between checkpoints
//      SIEVE_PROGRESS_SECS   (default 5)  seconds between progress lines

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>
#include <unistd.h>
#include <fcntl.h>
#include <omp.h>

struct Rec { uint64_t sigma, psi, phi, rest; };

// Every division hit() performs is EXACT. phi starts at n and after applying p
// is (n/p)(p-1); any later prime q | n divides n/p, so it divides phi too, and
// likewise psi and rest. An exact division by p is x * inv(p) mod 2^64, one
// multiply. The "does p still divide rest" test is the same multiply against
// floor((2^64-1)/p). AVX2 has no integer-divide instruction, so this is also
// the only form of the inner loop that could ever be vectorised.
struct Div { uint64_t inv, lim, pm1, pp1; uint32_t p; };

static inline uint64_t inv64(uint64_t a) {      // Newton, a odd: a*inv == 1 mod 2^64
    uint64_t x = a;                             // 3 bits correct
    for (int i = 0; i < 5; ++i) x *= 2 - a * x; // doubles each step: 6,12,24,48,96
    return x;
}

static uint64_t BLOCK = 1u << 13;                // 8192 recs = 256 KB, half the 512 KB L2

static std::vector<uint32_t> primes_to(uint64_t limit) {
    std::vector<bool> comp(limit + 1, false);
    std::vector<uint32_t> out;
    for (uint64_t i = 2; i <= limit; ++i)
        if (!comp[i]) { out.push_back((uint32_t)i);
                        for (uint64_t j = i*i; j <= limit; j += i) comp[j] = true; }
    return out;
}

// apply one prime at one index: all four fields, one resident record.
//
// v4 accumulates sigma, psi and phi as PRODUCTS over prime powers rather than as
// running n-scaled values. Multiplicativity gives, for n = rest * prod p^e with
// rest prime or 1,
//     sigma = prod S(p,e) * (rest+1),  psi = prod p^(e-1)(p+1) * (rest+1),
//     phi   = prod p^(e-1)(p-1) * (rest-1),
// so the final combination is three multiplies instead of the two divisions by a
// per-element runtime rest that v3 needed -- and those ran for the ~69% of
// integers that have a prime factor above the sieve root. In hit(), e == 1 is the
// common case and costs one multiply per field by a precomputed p-1 / p+1, so
// phi and psi stop dividing as well.
static inline void hit(Rec& r, uint8_t& om, const Div& d) {
    const uint64_t p = d.p, inv = d.inv, lim = d.lim;
    om++;
    if (p == 2) {
        uint64_t pw = 2, s = 3, pe1 = 1;
        r.rest >>= 1;
        while ((r.rest & 1) == 0) { r.rest >>= 1; pw <<= 1; s += pw; pe1 <<= 1; }
        r.sigma *= s;
        r.phi   *= pe1;                         // 2^(e-1) * (2-1)
        r.psi   *= pe1 * 3;                     // 2^(e-1) * (2+1)
        return;
    }
    r.rest *= inv;                              // exact: p | rest
    if (__builtin_expect(r.rest * inv > lim, 1)) {      // e == 1, the common case
        r.sigma *= 1 + p;
        r.phi   *= d.pm1;
        r.psi   *= d.pp1;
        return;
    }
    uint64_t pw = p * p, s = 1 + p + pw, pe1 = p;       // e >= 2
    r.rest *= inv;
    while (r.rest * inv <= lim) { r.rest *= inv; pw *= p; s += pw; pe1 *= p; }
    r.sigma *= s;
    r.phi   *= pe1 * d.pm1;
    r.psi   *= pe1 * d.pp1;
}

// On macOS fsync() only hands the data to the drive; the drive is still free to
// hold it in a volatile cache. F_FULLFSYNC is the call that asks it to commit.
// The resume argument is "the term is on disk strictly before the window's bit
// is", so this has to be the strong one on every platform.
static int durable_sync(FILE* f) {
    const int fd = fileno(f);
#ifdef __APPLE__
    if (fcntl(fd, F_FULLFSYNC, 0) == 0) return 0;   // falls through on failure
#endif
    return fsync(fd);
}

// The checkpoint is a per-window bitmap, and a bitmap is POSITIONAL: bit w means
// "window w of THIS run". Without a header the same bits were silently accepted
// for any other run, so resuming a different range skipped windows that had
// never been sieved -- and the bitmap then read as complete. Identity is the
// range and the window size; block size only affects speed, so it is recorded
// for diagnostics but not enforced.
static const char STATE_MAGIC[8] = {'A','3','9','0','0','4','9','S'};
static const uint32_t STATE_VERSION = 1;

struct StateHdr {
    char     magic[8];
    uint32_t version, hdrsize;
    uint64_t lo0, hi0, window, block, nwin;
};

static bool same_run(const StateHdr& a, const StateHdr& b) {
    return a.lo0 == b.lo0 && a.hi0 == b.hi0 && a.window == b.window && a.nwin == b.nwin;
}

// Write to <state>.tmp, sync, then rename: the reader never sees a partial file.
static bool save_state(const std::string& state, const StateHdr& h,
                       const std::vector<uint8_t>& bits) {
    const std::string tmp = state + ".tmp";
    FILE* f = fopen(tmp.c_str(), "wb");
    if (!f) return false;
    bool ok = fwrite(&h, sizeof h, 1, f) == 1
           && (bits.empty() || fwrite(bits.data(), 1, bits.size(), f) == bits.size())
           && fflush(f) == 0
           && durable_sync(f) == 0;
    if (fclose(f) != 0) ok = false;
    if (!ok) { remove(tmp.c_str()); return false; }
    return rename(tmp.c_str(), state.c_str()) == 0;
}

// Seconds between checkpoints / progress lines. A checkpoint costs one write of
// the whole bitmap under the critical section, so the interval has to be set by
// time rather than by a window count that means different things at different
// range sizes. The cost of a crash is bounded by the same number.
static double env_secs(const char* name, double dflt) {
    const char* v = getenv(name);
    if (!v || !*v) return dflt;
    char* end = nullptr;
    const double d = strtod(v, &end);
    if (end == v || *end || !(d >= 0.0) || d > 86400.0) {
        fprintf(stderr, "warning: ignoring %s=%s, using %.0f\n", name, v, dflt);
        return dflt;
    }
    return d;
}

int main(int argc, char** argv) {
    if (argc < 3) { fprintf(stderr, "usage: sieve <lo> <hi> [threads] [window] [state] [block]\n"); return 2; }
    const uint64_t lo0 = strtoull(argv[1], nullptr, 10);
    const uint64_t hi0 = strtoull(argv[2], nullptr, 10);
    const long threadsL = argc > 3 ? strtol(argv[3], nullptr, 10) : omp_get_max_threads();
    const uint64_t W   = argc > 4 ? strtoull(argv[4], nullptr, 10) : (1u << 20);
    const std::string state = argc > 5 ? argv[5] : "";
    if (argc > 6) BLOCK = strtoull(argv[6], nullptr, 10);

    // Each of these used to be undefined rather than rejected: hi <= lo wrapped
    // the window count and reported a range "complete" without sieving it, a
    // zero window or block divided by zero, and a negative thread count reached
    // num_threads() and asked libgomp for a terabyte of stacks.
    if (hi0 <= lo0) {
        fprintf(stderr, "error: hi (%llu) must be greater than lo (%llu)\n",
                (unsigned long long)hi0, (unsigned long long)lo0);
        return 2;
    }
    if (W == 0)     { fprintf(stderr, "error: window must be at least 1\n"); return 2; }
    if (BLOCK == 0) { fprintf(stderr, "error: block must be at least 1\n");  return 2; }
    if (BLOCK > (1ull << 40)) {          // bucket words pack the offset above bit 24
        fprintf(stderr, "error: block must be at most 2^40\n"); return 2;
    }
    if (threadsL < 1 || threadsL > 65536) {
        fprintf(stderr, "error: threads must be between 1 and 65536 (got %ld)\n", threadsL);
        return 2;
    }
    const int threads = (int)threadsL;

    const uint64_t root = (uint64_t)std::sqrt((double)hi0) + 2;
    std::vector<uint32_t> primes = primes_to(root);

    // A bucket word is (block offset << 24) | prime index, so the index has to
    // fit in 24 bits. pi(x) reaches 2^24 at x = 310248241, i.e. hi = 9.63e16;
    // past that the field used to wrap silently and corrupt the answer.
    if (primes.size() > (1u << 24)) {
        fprintf(stderr, "error: %zu primes below the sieve root overflow the 24-bit\n"
                        "       bucket slot field; this build is limited to hi <= 9.6e16\n",
                primes.size());
        return 2;
    }
    size_t nsmall = 0;
    while (nsmall < primes.size() && primes[nsmall] < BLOCK) nsmall++;

    std::vector<Div> dv(primes.size());
    for (size_t i = 0; i < primes.size(); ++i) {
        uint64_t q = primes[i];
        dv[i].p = (uint32_t)q;
        dv[i].inv = (q & 1) ? inv64(q) : 0;     // p == 2 takes the shift path
        dv[i].lim = (q & 1) ? (~0ull) / q : 0;
        dv[i].pm1 = q - 1;
        dv[i].pp1 = q + 1;
    }

    const uint64_t nwin = (hi0 - lo0 + W - 1) / W;
    // Completed windows are tracked as a BITMAP, not a count. schedule(dynamic)
    // finishes windows out of order, so a count is not a low-water mark: resuming
    // at one would skip windows that were never processed, silently losing any
    // term inside them. One bit per window is ~1 MB for the whole a(10) sweep.
    std::vector<uint8_t> doneBits((nwin + 7) / 8, 0);

    StateHdr hdr;
    memset(&hdr, 0, sizeof hdr);
    memcpy(hdr.magic, STATE_MAGIC, sizeof hdr.magic);
    hdr.version = STATE_VERSION;  hdr.hdrsize = (uint32_t)sizeof(StateHdr);
    hdr.lo0 = lo0;  hdr.hi0 = hi0;  hdr.window = W;  hdr.block = BLOCK;  hdr.nwin = nwin;

    uint64_t already = 0;
    if (!state.empty()) {
        if (FILE* f = fopen(state.c_str(), "rb")) {     // absent = start from scratch
            StateHdr got;
            if (fread(&got, sizeof got, 1, f) != 1) {
                fprintf(stderr, "error: %s is too short to be a state file\n", state.c_str());
                fclose(f); return 3;
            }
            if (memcmp(got.magic, STATE_MAGIC, sizeof got.magic) != 0
                || got.version != STATE_VERSION || got.hdrsize != sizeof(StateHdr)) {
                fprintf(stderr, "error: %s is not a v%u sieve state file\n",
                        state.c_str(), (unsigned)STATE_VERSION);
                fclose(f); return 3;
            }
            if (!same_run(got, hdr)) {
                fprintf(stderr,
                    "error: %s belongs to a different run and cannot be resumed\n"
                    "  file: [%llu,%llu) window %llu -> %llu windows\n"
                    "  this: [%llu,%llu) window %llu -> %llu windows\n"
                    "  a bitmap is positional; applying it here would mark windows\n"
                    "  complete that were never sieved\n", state.c_str(),
                    (unsigned long long)got.lo0,(unsigned long long)got.hi0,
                    (unsigned long long)got.window,(unsigned long long)got.nwin,
                    (unsigned long long)hdr.lo0,(unsigned long long)hdr.hi0,
                    (unsigned long long)hdr.window,(unsigned long long)hdr.nwin);
                fclose(f); return 3;
            }
            const size_t nread = fread(doneBits.data(), 1, doneBits.size(), f);
            if (nread != doneBits.size()) {
                fprintf(stderr, "error: %s is truncated: %zu of %zu bitmap bytes\n",
                        state.c_str(), nread, doneBits.size());
                fclose(f); return 3;
            }
            fclose(f);
            if (got.block != hdr.block)
                fprintf(stderr, "note: state file was written with block %llu, now %llu"
                                " (block does not move window boundaries)\n",
                        (unsigned long long)got.block,(unsigned long long)hdr.block);
            for (uint64_t w = 0; w < nwin; ++w)
                if (doneBits[w >> 3] & (1u << (w & 7))) already++;
        }
    }
    fprintf(stderr, "range [%llu,%llu) threads %d window %llu block %llu\n"
                    "primes<=%llu: %zu (%zu small, %zu bucketed)  already complete: %llu/%llu\n",
            (unsigned long long)lo0,(unsigned long long)hi0,threads,
            (unsigned long long)W,(unsigned long long)BLOCK,
            (unsigned long long)root, primes.size(), nsmall, primes.size()-nsmall,
            (unsigned long long)already,(unsigned long long)nwin);

    const double ckpt_secs = env_secs("SIEVE_CHECKPOINT_SECS", 60.0);
    const double prog_secs  = env_secs("SIEVE_PROGRESS_SECS",    5.0);

    double t0 = omp_get_wtime();
    uint64_t done = 0;
    double last_ckpt = t0, last_prog = t0;   // shared, only touched under critical

    #pragma omp parallel num_threads(threads)
    {
        // v3: init and the final test moved INSIDE the block loop. In v2 they each
        // swept the whole window, so the buffer was written once and read once
        // through DRAM per window regardless of the blocking in between -- which
        // is why v2 got faster as the window got smaller, the opposite of what
        // blocking is supposed to buy. Here a block is initialised, sieved and
        // tested while resident, so the buffer is only ever BLOCK records and the
        // window is free to grow: its only remaining job is to amortise the
        // bucket build, one division per large prime rather than one per block.
        std::vector<Rec> buf(BLOCK);
        std::vector<uint8_t> om(BLOCK);
        std::vector<uint64_t> nextoff(nsmall ? nsmall : 1);
        const uint64_t nblk = (W + BLOCK - 1) / BLOCK;
        std::vector<std::vector<uint64_t>> bucket(nblk);   // packed (blockoff<<24)|slot

        #pragma omp for schedule(dynamic, 1)
        for (uint64_t w = 0; w < nwin; ++w) {
            uint8_t seen;                                       // atomic: the
            #pragma omp atomic read                             // writers below
            seen = doneBits[w >> 3];                            // run concurrently
            if (seen & (1u << (w & 7))) continue;               // already done
            const uint64_t lo = lo0 + w * W;
            const uint64_t hi = (lo + W < hi0) ? lo + W : hi0;
            const uint64_t n  = hi - lo;

            for (auto& b : bucket) b.clear();
            for (size_t k = nsmall; k < primes.size(); ++k) {   // bucket the large primes
                uint32_t p = primes[k];
                if ((uint64_t)p * p > hi) break;
                uint64_t start = ((lo + p - 1) / p) * p;
                if (start >= hi) continue;
                uint32_t slot = (uint32_t)k;    // 295k primes, fits the 24-bit field
                uint64_t i0 = start - lo;
                if (lo == 0 && i0 == 0) i0 = p;   // n=0 has rest=0; hit() would spin
                for (uint64_t i = i0; i < n; i += p) {
                    uint64_t bi = i / BLOCK;
                    bucket[bi].push_back(((i - bi * BLOCK) << 24) | slot);
                }
            }

            size_t nsm = 0;                                     // rolling offsets, one
            for (size_t k = 0; k < nsmall; ++k) {               // division per window
                uint32_t p = primes[k];
                if ((uint64_t)p * p > hi) break;
                uint64_t start = ((lo + p - 1) / p) * p;
                nextoff[k] = (start == 0) ? p : start - lo;
                nsm = k + 1;
            }

            for (uint64_t b = 0; b < nblk; ++b) {
                const uint64_t bs = b * BLOCK;
                if (bs >= n) break;
                const uint64_t be = (bs + BLOCK < n) ? bs + BLOCK : n;
                const uint64_t blen = be - bs, base = lo + bs;

                for (uint64_t j = 0; j < blen; ++j) {
                    buf[j] = Rec{1,1,1,base + j}; om[j] = 0;
                }
                if (base == 0) buf[0] = Rec{0,0,0,1};

                for (size_t k = 0; k < nsm; ++k) {              // small primes, walked
                    uint32_t p = primes[k];
                    uint64_t j = nextoff[k] - bs;
                    for (; j < blen; j += p) hit(buf[j], om[j], dv[k]);
                    nextoff[k] = j + bs;
                }
                for (uint64_t packed : bucket[b])               // large primes, drained
                    hit(buf[packed >> 24], om[packed >> 24], dv[packed & 0xffffff]);

                for (uint64_t j = 0; j < blen; ++j) {
                    Rec& r = buf[j];
                    uint64_t sg = r.sigma, ps = r.psi, ph = r.phi, o = om[j];
                    if (r.rest > 1) {           // rest is 1 or a single prime > root
                        sg *= r.rest + 1; ps *= r.rest + 1; ph *= r.rest - 1; o++;
                    }
                    if (sg == ps + ph + o && base + j > 1)
                        #pragma omp critical
                        { printf("TERM %llu\n", (unsigned long long)(base+j)); fflush(stdout);
                          if (!state.empty()) {
                              FILE* tf = fopen((state + ".terms").c_str(), "a");
                              bool ok = tf != nullptr;
                              if (ok) {
                                  ok = fprintf(tf, "%llu\n", (unsigned long long)(base+j)) > 0
                                    && fflush(tf) == 0 && durable_sync(tf) == 0;
                                  if (fclose(tf) != 0) ok = false;
                              }
                              if (!ok) fprintf(stderr, "warning: could not record term "
                                               "%llu in %s.terms\n",
                                               (unsigned long long)(base+j), state.c_str());
                          } }
                }
            }

            #pragma omp critical
            {
                done++;
                #pragma omp atomic update
                doneBits[w >> 3] |= (uint8_t)(1u << (w & 7));   // this window, specifically
                const double now = omp_get_wtime();
                if (!state.empty() && now - last_ckpt >= ckpt_secs) {
                    if (!save_state(state, hdr, doneBits))
                        fprintf(stderr, "\nwarning: could not write checkpoint %s\n", state.c_str());
                    last_ckpt = omp_get_wtime();       // measure from completion
                }
                if (now - last_prog >= prog_secs) {
                    last_prog = now;
                    const double el = now - t0;
                    fprintf(stderr, "\r%llu/%llu  %.1f M n/s  eta %.2f h    ",
                            (unsigned long long)(already+done),(unsigned long long)nwin,
                            done*(double)W/el/1e6,
                            (nwin-already-done)*el/done/3600.0);
                }
            }
        }
    }
    const bool saved = state.empty() || save_state(state, hdr, doneBits);
    double el = omp_get_wtime() - t0;
    fprintf(stderr, "\ncomplete: %.1f s, %.1f M n/s\n", el, (hi0-lo0)/el/1e6);
    if (!saved) {
        fprintf(stderr, "error: final checkpoint to %s failed; the run is NOT resumable\n",
                state.c_str());
        return 4;
    }
    return 0;
}
