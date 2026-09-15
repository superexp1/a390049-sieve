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
//   g++ -O3 -march=native -fopenmp -o oddpart src/oddpart.cpp
//   ./oddpart <M> [threads]
//
// exit 0 done, 2 bad arguments
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <vector>
#include <omp.h>

struct Rec { uint64_t sigma, psi, phi, rest; };

static std::vector<uint32_t> primes_to(uint64_t n) {
    std::vector<bool> c(n+1,false); std::vector<uint32_t> o;
    for (uint64_t i=2;i<=n;++i) if(!c[i]){o.push_back((uint32_t)i);
        for(uint64_t j=i*i;j<=n;j+=i)c[j]=true;}
    return o;
}

int main(int argc,char**argv){
    if(argc<2){fprintf(stderr,"usage: oddpart <M> [threads]\n");return 2;}
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
    const uint64_t W=1u<<20;
    const uint64_t nwin=(M+W-1)/W;
    double t0=omp_get_wtime(); uint64_t done=0;

    #pragma omp parallel num_threads(th)
    {
      std::vector<Rec> buf(W/2+1); std::vector<uint8_t> om(W/2+1);
      #pragma omp for schedule(dynamic,1)
      for(uint64_t w=0;w<nwin;++w){
        const uint64_t lo=w*W, hi=(lo+W<M)?lo+W:M;
        const uint64_t lo1=lo|1;                  // first odd >= lo
        if(lo1>=hi)continue;
        const uint64_t nodd=(hi-lo1+1)/2;         // odd values in [lo,hi)
        for(uint64_t i=0;i<nodd;++i){buf[i]=Rec{1,1,1,lo1+2*i};om[i]=0;}
        for(uint32_t p:pr){
          if((uint64_t)p*p>hi)break;
          if(p==2)continue;                       // 2 has no odd multiples
          uint64_t q=(lo1+p-1)/p; if((q&1)==0)++q;// first ODD multiplier of p
          const uint64_t start=q*p;               // odd, and >= lo1 >= 1
          if(start>=hi)continue;
          for(uint64_t i=(start-lo1)/2;i<nodd;i+=p){
            Rec&r=buf[i]; uint64_t pw=p,S=1+p,pe1=1;
            r.rest/=p;
            while(r.rest%p==0){r.rest/=p;pw*=p;S+=pw;pe1*=p;}
            r.sigma*=S; r.phi*=pe1*(p-1); r.psi*=pe1*(p+1); om[i]++;
          }
        }
        for(uint64_t i=0;i<nodd;++i){
          const uint64_t m=lo1+2*i; if(m<3)continue;
          Rec&r=buf[i]; uint64_t P=r.sigma,Q=r.psi,R=r.phi; uint32_t o=om[i];
          if(r.rest>1){P*=r.rest+1;Q*=r.rest+1;R*=r.rest-1;o++;}
          // odd k = m itself
          if(P==Q+R+(uint64_t)o)
            #pragma omp critical
            {printf("ODD-TERM k=%llu omega=%u\n",(unsigned long long)m,o);fflush(stdout);}
          // even k = 2^a * m
          __int128 D=(__int128)4*P-(__int128)3*Q-(__int128)R;
          if(D>0){
            __int128 num=(__int128)P+o+1;
            if(num%D==0){
              __int128 X=num/D;
              if(X>=1&&(X&(X-1))==0)
                #pragma omp critical
                {printf("TERM oddpart=%llu X=%llu omega=%u\n",
                   (unsigned long long)m,(unsigned long long)(uint64_t)X,o+1);fflush(stdout);}
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
