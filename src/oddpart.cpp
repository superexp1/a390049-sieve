// A390049 by ODD PART rather than by size.
//
// k = 2^a * m with m odd and X = 2^(a-1) satisfies X = (P + omega)/D where
// P,Q,R are sigma,psi,phi of m and D = 4P - 3Q - R.  X is DETERMINED by m, so
// each odd m yields at most one candidate k -- of any magnitude.  Enumerating
// odd m <= M is therefore complete for every term whose ODD PART is at most M,
// with no ceiling on k at all.  That is complementary to the segmented sieve,
// which is complete below a bound on k but says nothing above it.
//
// Odd k is tested in the same pass by its own equation P = Q + R + omega.
//
// Odd m is enumerated over [0, M), i.e. m < M.
//
// The buffer holds only ODD values: index i is m = (lo|1) + 2i. Consecutive
// odd multiples of an odd p differ by 2p, so they step by p in index space.
// The previous version initialised and sieved every integer in the window and
// then threw the even half away at the test, doing twice the necessary work in
// twice the memory.
//
// The window is processed in cache-sized BLOCKS. It had none: a 2^20 window is
// 16 MB of records, so every record was written, scattered over, and read back
// through DRAM. Blocking initialises, sieves and tests a block while it is
// resident, which measured 21.5x. The window survives as the unit over which
// the per-prime start offset is amortised -- one division per prime per window,
// then a rolling offset per block -- so shrinking the window instead would
// trade the cache win against a start-offset scan that grows with M.
//
//   g++ -O3 -march=native -fopenmp -o oddpart src/oddpart.cpp
//   ./oddpart <M> [threads] [window] [block]
//
// exit 0 done, 2 bad arguments
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <vector>
#include <omp.h>

struct Rec { uint64_t sigma, psi, phi, rest; };

// Every division the sieve loop performs is EXACT, so it is a multiply by the
// modular inverse mod 2^64, and "does p still divide rest" is the same multiply
// against floor((2^64-1)/p). The main sieve has done this since v4; oddpart
// never got it. It was worth 3% while oddpart was DRAM-bound and much more once
// blocking made it resident.
struct Dv { uint64_t inv, lim; };
static inline uint64_t inv64(uint64_t a){      // Newton, a odd: a*inv==1 mod 2^64
    uint64_t x=a; for(int i=0;i<5;++i)x*=2-a*x; return x; }

static std::vector<uint32_t> primes_to(uint64_t n) {
    std::vector<bool> c(n+1,false); std::vector<uint32_t> o;
    for (uint64_t i=2;i<=n;++i) if(!c[i]){o.push_back((uint32_t)i);
        for(uint64_t j=i*i;j<=n;j+=i)c[j]=true;}
    return o;
}

int main(int argc,char**argv){
    if(argc<2){fprintf(stderr,"usage: oddpart <M> [threads] [window] [block]\n");return 2;}
    const uint64_t M=strtoull(argv[1],nullptr,10);
    // atoi used to accept a negative count, which reached num_threads() and had
    // libgomp try to allocate a terabyte of stacks before dying.
    const long thL=argc>2?strtol(argv[2],nullptr,10):omp_get_max_threads();
    if(thL<1||thL>65536){
        fprintf(stderr,"error: threads must be between 1 and 65536 (got %ld)\n",thL);
        return 2;}
    const int th=(int)thL;
    if(M<3){fprintf(stderr,"nothing to do: M must be at least 3 (smallest odd m)\n");return 0;}
    const uint64_t root=(uint64_t)std::sqrt((double)M)+2;
    std::vector<uint32_t> pr=primes_to(root);
    std::vector<uint32_t> opr;                // 2 has no odd multiples at all
    for(uint32_t q:pr) if(q!=2) opr.push_back(q);
    std::vector<Dv> dv(opr.size());
    for(size_t i=0;i<opr.size();++i){ const uint64_t q=opr[i];
        dv[i].inv=inv64(q); dv[i].lim=(~0ull)/q; }   // opr is all odd
    const uint64_t W=(argc>3)?strtoull(argv[3],nullptr,10):(1u<<20);
    // 4096 recs = 128 KB. The resident set is buf + om + nextoff, and nextoff
    // is 8 bytes per prime below the root -- 104 KB at M = 2e10 -- so the
    // records get about a quarter of a 512 KB L2, not half. Swept flat from
    // 2048 up on a 64-core EPYC 7B13; 1024 loses to the per-block prime scan.
    const uint64_t B=(argc>4)?strtoull(argv[4],nullptr,10):4096;
    if(W<2){fprintf(stderr,"error: window must be at least 2\n");return 2;}
    if(B<1){fprintf(stderr,"error: block must be at least 1\n");return 2;}
    const uint64_t nwin=(M+W-1)/W;
    double t0=omp_get_wtime(); uint64_t done=0;

    #pragma omp parallel num_threads(th)
    {
      std::vector<Rec> buf(B); std::vector<uint8_t> om(B);
      std::vector<uint64_t> nextoff(opr.size()?opr.size():1);
      #pragma omp for schedule(dynamic,1)
      for(uint64_t w=0;w<nwin;++w){
        const uint64_t lo=w*W, hi=(lo+W<M)?lo+W:M;
        const uint64_t lo1=lo|1;                  // first odd >= lo
        if(lo1>=hi)continue;
        const uint64_t nodd=(hi-lo1+1)/2;         // odd values in [lo,hi)

        size_t npr=0;                             // ONE division per prime per
        for(size_t k=0;k<opr.size();++k){         // window; blocks then roll
          const uint32_t p=opr[k];
          if((uint64_t)p*p>hi)break;
          uint64_t q=(lo1+p-1)/p; if((q&1)==0)++q;// first ODD multiplier of p
          nextoff[k]=(q*p-lo1)/2;                 // its odd-index in the window
          npr=k+1;
        }

        for(uint64_t bs=0;bs<nodd;bs+=B){
        const uint64_t be=(bs+B<nodd)?bs+B:nodd, blen=be-bs, base=lo1+2*bs;
        for(uint64_t i=0;i<blen;++i){buf[i]=Rec{1,1,1,base+2*i};om[i]=0;}
        for(size_t k=0;k<npr;++k){
          const uint32_t p=opr[k];
          uint64_t j=nextoff[k]-bs;               // >= 0: nextoff never trails bs
          const uint64_t inv=dv[k].inv, lim=dv[k].lim;
          for(;j<blen;j+=p){
            Rec&r=buf[j]; uint64_t pw=p,S=1+p,pe1=1;
            r.rest*=inv;                          // exact: p | rest
            while(r.rest*inv<=lim){r.rest*=inv;pw*=p;S+=pw;pe1*=p;}
            r.sigma*=S; r.phi*=pe1*(p-1); r.psi*=pe1*(p+1); om[j]++;
          }
          nextoff[k]=j+bs;
        }
        for(uint64_t i=0;i<blen;++i){
          const uint64_t m=base+2*i; if(m<3)continue;
          Rec&r=buf[i]; uint64_t P=r.sigma,Q=r.psi,R=r.phi; uint32_t o=om[i];
          if(r.rest>1){P*=r.rest+1;Q*=r.rest+1;R*=r.rest-1;o++;}
          // odd k = m itself
          if(P==Q+R+(uint64_t)o)
            #pragma omp critical
            {printf("ODD-TERM k=%llu omega=%u\n",(unsigned long long)m,o);fflush(stdout);}
          // even k = 2^a * m
          // X = num/D must be a power of two. Write num = 2^v*a, D = 2^w*b
          // with a,b odd; then num/D = 2^(v-w)*(a/b), and a/b is a ratio of odd
          // numbers, so it is a power of two only when it is 1. Hence a == b
          // and v >= w, with X = 2^(v-w). No division, which also makes this
          // half of the program portable to hardware with no integer divide.
          __int128 D=(__int128)4*P-(__int128)3*Q-(__int128)R;
          if(D>0){
            const uint64_t num=P+o+1, Du=(uint64_t)D;
            if((__int128)Du==D){
              const int vn=__builtin_ctzll(num), vd=__builtin_ctzll(Du);
              const uint64_t X=(vn>=vd)?(1ull<<(vn-vd)):0;
              if(vn>=vd&&(num>>vn)==(Du>>vd))
                #pragma omp critical
                {printf("TERM oddpart=%llu X=%llu omega=%u\n",
                   (unsigned long long)m,(unsigned long long)(uint64_t)X,o+1);fflush(stdout);}
            }
          }
        }
        }
        #pragma omp critical
        { if((++done&0x3f)==0){double e=omp_get_wtime()-t0;
            fprintf(stderr,"\r%llu/%llu %.0f M/s eta %.2fh   ",
              (unsigned long long)done,(unsigned long long)nwin,
              done*(double)W/e/1e6,(nwin-done)*e/done/3600.0);} }
      }
    }
    fprintf(stderr,"\ncomplete: %.1f s\n",omp_get_wtime()-t0);
    return 0;
}
