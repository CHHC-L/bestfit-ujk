// Best Fit on U{j,k}: counts N[g] of open bins with gap g (1..k-1). Item s uniform in 1..j.
// Print waste at checkpoints, averaged over reps.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
static uint64_t s[2]={0x9E3779B97F4A7C15ULL,0xD1B54A32D192ED03ULL};
static inline uint64_t rotl(uint64_t x,int k){return (x<<k)|(x>>(64-k));}
static inline uint64_t nxt(void){uint64_t s0=s[0],s1=s[1],r=s0+s1;s1^=s0;s[0]=rotl(s0,55)^s1^(s1<<14);s[1]=rotl(s1,36);return r;}
int main(int argc,char**argv){
  if(argc<5){fprintf(stderr,"usage: bfsim j k n reps\n");return 2;}
  int j=atoi(argv[1]),k=atoi(argv[2]); long n=atol(argv[3]); int reps=atoi(argv[4]);
  if(j<1||j>=k||k>4095){fprintf(stderr,"need 1 <= j < k <= 4095\n");return 2;}
  s[0]^=j*1000003ULL+k; s[1]^=k*7919ULL;
  long N[4096]; double w1=0,w2=0;
  for(int r=0;r<reps;r++){
    for(int g=0;g<k;g++)N[g]=0;
    long waste=0; long half=n/2, wh=0;
    for(long t=1;t<=n;t++){
      int sz=1+(int)((nxt()>>11)%(uint64_t)j);
      int g=sz; while(g<k && N[g]==0) g++;
      if(g<k){N[g]--; if(g-sz>0)N[g-sz]++; waste-=sz;}
      else { if(k-sz>0)N[k-sz]++; waste+=k-sz; }
      if(t==half) wh=waste;
    }
    w1+=wh; w2+=waste;
  }
  printf("%d %d n=%ld waste(n/2)=%.1f waste(n)=%.1f ratio=%.3f slope=%.5f\n",j,k,n,w1/reps,w2/reps,(w2/(w1+1e-9)),(w2-w1)/reps/(n/2.0));
  return 0;
}
