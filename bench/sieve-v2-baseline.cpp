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
// FROZEN: kept only so the README's speedup number stays checkable. Its BLOCK
// is a compile-time constant, so it takes FIVE arguments, not six. The only
// changes since v2 are the two correctness fixes marked below.
//
//   g++ -O3 -march=native -fopenmp -o sieve2 sieve2.cpp
//   ./sieve2 <lo> <hi> [threads] [window] [statefile]

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <string>
#include <vector>
#include <unistd.h>
#include <omp.h>

struct Rec { uint64_t sigma, psi, phi, rest; };

static const uint64_t BLOCK = 1u << 14;          // 16k records = 512 KiB, fits L2

static std::vector<uint32_t> primes_to(uint64_t limit) {
    std::vector<bool> comp(limit + 1, false);
    std::vector<uint32_t> out;
    for (uint64_t i = 2; i <= limit; ++i)
        if (!comp[i]) { out.push_back((uint32_t)i);
                        for (uint64_t j = i*i; j <= limit; j += i) comp[j] = true; }
    return out;
}

// apply one prime at one index: all four fields, one resident record
static inline void hit(Rec& r, uint8_t& om, uint32_t p) {
    r.phi -= r.phi / p;
    r.psi += r.psi / p;
    om++;
    uint64_t pw = p, s = 1 + p;
    r.rest /= p;
    while (r.rest % p == 0) { r.rest /= p; pw *= p; s += pw; }
    r.sigma *= s;
}

int main(int argc, char** argv) {
    if (argc < 3) { fprintf(stderr, "usage: sieve2 <lo> <hi> [threads] [window] [state]\n"); return 2; }
    const uint64_t lo0 = strtoull(argv[1], nullptr, 10);
    const uint64_t hi0 = strtoull(argv[2], nullptr, 10);
    const int threads  = argc > 3 ? atoi(argv[3]) : omp_get_max_threads();
    const uint64_t W   = argc > 4 ? strtoull(argv[4], nullptr, 10) : (1u << 21);
    const std::string state = argc > 5 ? argv[5] : "";

    const uint64_t root = (uint64_t)std::sqrt((double)hi0) + 2;
    std::vector<uint32_t> primes = primes_to(root);
    size_t nsmall = 0;
    while (nsmall < primes.size() && primes[nsmall] < BLOCK) nsmall++;

    const uint64_t nwin = (hi0 - lo0 + W - 1) / W;
    // Completed windows are tracked as a BITMAP, not a count. schedule(dynamic)
    // finishes windows out of order, so a count is not a low-water mark: resuming
    // at one would skip windows that were never processed, silently losing any
    // term inside them. One bit per window is ~1 MB for the whole a(10) sweep.
    std::vector<uint8_t> doneBits((nwin + 7) / 8, 0);
    uint64_t already = 0;
    if (!state.empty()) {
        if (FILE* f = fopen(state.c_str(), "rb")) {
            // v5 fix: a short read was accepted silently.
            if (fread(doneBits.data(), 1, doneBits.size(), f) != doneBits.size())
                fprintf(stderr, "warning: short read from %s\n", state.c_str());
            fclose(f);
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

    double t0 = omp_get_wtime();
    uint64_t done = 0;

    #pragma omp parallel num_threads(threads)
    {
        std::vector<Rec> buf(W);
        std::vector<uint8_t> om(W);
        const uint64_t nblk = (W + BLOCK - 1) / BLOCK;
        std::vector<std::vector<uint64_t>> bucket(nblk);   // packed (idx<<24)|prime-index
        std::vector<uint32_t> bprime;

        #pragma omp for schedule(dynamic, 1)
        for (uint64_t w = 0; w < nwin; ++w) {
            if (doneBits[w >> 3] & (1u << (w & 7))) continue;   // already done
            const uint64_t lo = lo0 + w * W;
            const uint64_t hi = (lo + W < hi0) ? lo + W : hi0;
            const uint64_t n  = hi - lo;

            for (uint64_t i = 0; i < n; ++i) { uint64_t v = lo+i; buf[i] = Rec{1,v,v,v}; om[i]=0; }
            if (lo == 0) buf[0] = Rec{0,0,0,1};

            for (auto& b : bucket) b.clear();
            bprime.clear();
            for (size_t k = nsmall; k < primes.size(); ++k) {   // bucket the large primes
                uint32_t p = primes[k];
                if ((uint64_t)p * p > hi) break;
                uint64_t start = ((lo + p - 1) / p) * p;
                if (start >= hi) continue;
                uint32_t slot = (uint32_t)bprime.size();
                bprime.push_back(p);
                // v5 fix: n=0 has rest=0 and hit() would spin forever on it.
                uint64_t i0 = start - lo;
                if (lo == 0 && i0 == 0) i0 = p;
                for (uint64_t i = i0; i < n; i += p)
                    bucket[i / BLOCK].push_back((i << 24) | slot);
            }

            for (uint64_t b = 0; b < nblk; ++b) {
                const uint64_t bs = b * BLOCK, be = (bs + BLOCK < n) ? bs + BLOCK : n;
                if (bs >= n) break;
                for (size_t k = 0; k < nsmall; ++k) {           // small primes, walked
                    uint32_t p = primes[k];
                    if ((uint64_t)p * p > hi) break;
                    uint64_t start = ((lo + bs + p - 1) / p) * p;
                    if (start >= lo + be) continue;
                    uint64_t i = start - lo;
                    if (lo == 0 && i == 0) i = p;
                    for (; i < be; i += p) hit(buf[i], om[i], p);
                }
                for (uint64_t packed : bucket[b])               // large primes, drained
                    hit(buf[packed >> 24], om[packed >> 24], bprime[packed & 0xffffff]);
            }

            for (uint64_t i = 0; i < n; ++i) {
                Rec& r = buf[i];
                if (r.rest > 1) { r.phi -= r.phi/r.rest; r.psi += r.psi/r.rest;
                                  r.sigma *= r.rest + 1; om[i]++; }
                if (r.sigma == r.psi + r.phi + (uint64_t)om[i] && lo + i > 1)
                    #pragma omp critical
                    { printf("TERM %llu\n", (unsigned long long)(lo+i)); fflush(stdout);
                      if (!state.empty()) {
                          if (FILE* tf = fopen((state + ".terms").c_str(), "a")) {
                              fprintf(tf, "%llu\n", (unsigned long long)(lo+i));
                              fflush(tf); fsync(fileno(tf)); fclose(tf);
                          }
                      } }
            }

            #pragma omp critical
            {
                done++;
                doneBits[w >> 3] |= (uint8_t)(1u << (w & 7));   // this window, specifically
                if (!state.empty() && (done & 0xff) == 0) {
                    if (FILE* f = fopen((state + ".tmp").c_str(), "wb")) {
                        fwrite(doneBits.data(), 1, doneBits.size(), f);
                        fflush(f); fsync(fileno(f)); fclose(f);
                        rename((state + ".tmp").c_str(), state.c_str());
                    }
                }
                if ((done & 0x3f) == 0) {
                    double el = omp_get_wtime() - t0;
                    fprintf(stderr, "\r%llu/%llu  %.1f M n/s  eta %.2f h    ",
                            (unsigned long long)(already+done),(unsigned long long)nwin,
                            done*(double)W/el/1e6,
                            (nwin-already-done)*el/done/3600.0);
                }
            }
        }
    }
    if (!state.empty()) {
        if (FILE* f = fopen((state + ".tmp").c_str(), "wb")) {
            fwrite(doneBits.data(), 1, doneBits.size(), f);
            fflush(f); fsync(fileno(f)); fclose(f);
            rename((state + ".tmp").c_str(), state.c_str());
        }
    }
    double el = omp_get_wtime() - t0;
    fprintf(stderr, "\ncomplete: %.1f s, %.1f M n/s\n", el, (hi0-lo0)/el/1e6);
    return 0;
}
