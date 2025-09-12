     /* fxcm file compressor is based on paq8 model.

    Copyright (C) 2025 Kaido Orav

    LICENSE

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of
    the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details at
    Visit <http://www.gnu.org/copyleft/gpl.html>.


To compress:   fxcm c input output
To decompress: fxcm d input output

This compressor contains changed parts (mostly marked with comments) from fallowing open source programs:
paq8hp12
paq8px_v208
paq8pxv
cmix

*/

#define PLAINTEXT  // uncomment for plain text input, otherwise dictionray encoded input is expected
#ifdef PLAINTEXT
    #define VERSION 25
#else
    #define VERSION 26
#endif

#include <stdio.h>
#include <time.h>
#include <mem.h>
#define NDEBUG  // remove for debugging (turns on Array bound checks)
#include <assert.h>

#ifdef UNIX  // not tested
#include <stdio.h>
#include <sys/types.h>
#include <stdlib.h>
#include <memory.h>
#include <cstdio>
#include <ctype.h>
#include <sys/cdefs.h>
#include <dirent.h>
#include <errno.h>
#endif
// AVX2
#include <immintrin.h>
// 8, 16, 32 bit unsigned types (adjust as appropriate)
typedef unsigned char  U8;
typedef unsigned short U16;
typedef unsigned int   U32;
typedef unsigned long long int U64;
//
typedef __m128i XMM;
typedef __m256i YMM;

// min, max functions
#if  !defined(WINDOWS) || !defined (min)
inline int min(int a, int b) {return a<b?a:b;}
inline int max(int a, int b) {return a<b?b:a;}
#endif

#define ispowerof2(x) ((x&(x-1))==0)
#include <math.h>
#include <windows.h>

const int num_models = 560;// 16
short *model_predictions1,*model_predictions1_ptr;
unsigned int prediction_index=0;

void AddPrediction(int x) {
    model_predictions1[prediction_index++]=x;
    assert(prediction_index >= 0 && prediction_index < num_models);
}

void ResetPredictions() {
    assert(prediction_index >= 0 && prediction_index <= num_models);
    prediction_index = 0;
}

// Array
template <class T> void alloc(T*&ptr, int c) {
  ptr=(T*)calloc(c, sizeof(T));
  if (!ptr) exit(1);
}
 
// for aligned data
template <class T> void alloc1(T*&data, int c,T*&ptr,const int align=16) {
  ptr=(T*)calloc(c, sizeof(T));
  if (!ptr) exit(1);
  data=(T*)(((uintptr_t)ptr+(align-1)) & ~(uintptr_t)(align-1));
}

// Squash returns p = 1/(1 + exp(-d)), d scaled by 8 bits, p scaled by 12 bits
short sqt[4096];

int squashc(int d ) {
    if (d < -2047)return 1;
    if (d > 2047)return 4095;
    float p = 1.0f / (1.0f + exp(-d / 256.0));
    p *= 4096.0;
    U32 pi = (U32)round(p);
    if (pi > 4095)pi = 4095;
    if (pi < 1)pi = 1;
    return pi;
}

inline int squash(int d) {
  if (d < -2047)return 1;
  if (d > 2047)return 4095;
  return sqt[d + 2047];
}

// Stretch is inverse of squash. d = ln(p/(1-p)), d scaled by 8 bits, p by 12 bits.
// d has range -2047 to 2047 representing -8 to 8. p has range 0 to 4095.
short strt[4096];
inline short clp(int z) {
    if (z<-2047) {
        z=-2047;
    } else if (z>2047) {
        z=2047;
    }
    return z;
}
inline short clp1(int z){
    if (z<0) {
        z=0;
    } else if (z>4095) {
        z=4095;
    }
    return z;
}
int stretchc(int p) {
    assert(p>=0 && p<=4095);
    if (p==0) p=1;
    float f=p/4096.0f;
    float d=log(f/(1.0f-f))*256.0f;
    int di=(int)round(d);
    if (di>2047) di = 2047;
    if (di<-2047) di = -2047;
    return di;
}

inline short stretch(int p) { 
    assert(p>=0 && p<=4095);
    if (p<0) p=0;
    if (p>4095) p=4095;
    return strt[p];
}

template <const int S=256>
struct alignas(64) Inputs {
    short n[S];
    int ncount;     // mixer input count
    void add(int p){ 
        assert(p>-2048 && p<2048);// fixme, when enabled compression is different
        n[ncount++]=clp(p);
        assert(ncount>=0 && ncount<=S);
        AddPrediction((clp(p)));
    }
};

template <const int S>
struct BlockData {
    int y;        // Last bit, 0 or 1, set by encoder
    int c0;       // Last 0-7 bits of the partial byte with a leading 1 bit (1-255)
    U32 c4;       // Last 4 whole bytes, packed.
    int bpos;     // bits in c0 (0 to 7)
    int blpos;    // Relative position in block
    int bposshift;
    int c0shift_bpos;
    int cmBitState;
    Inputs<S> mxInputs1; // array of inputs, for two layers
    Inputs<64> mxInputs2;
    Inputs<64> mxInputs4;
    void Init() {
        y=0, c0=1, c4=0 ,bpos=0, blpos=0, bposshift=0, c0shift_bpos=0, cmBitState=0;
    }
};

BlockData<544> x; //maintains current global data block

// ilog(x) = round(log2(x) * 16), 0 <= x < 256
U8 ilog[256];
// Compute lookup table by numerical integration of 1/x
void InitIlog() {
    U32 x=14155776;
    for (int i=2; i<257; ++i) {
        x+=774541002/(i*2-1);  // numerator is 2^29/ln 2
        ilog[i-1]=x>>24;
    }
}

// State table
//   nex(state, 0) = next state if bit y is 0, 0 <= state < 256
//   nex(state, 1) = next state if bit y is 1
//   nex(state, 2) = number of zeros in bit history represented by state
//   nex(state, 3) = number of ones represented
//
// States represent a bit history within some context.
struct StateTable {
  int mdc; // maximum discount
  enum {B=5, N=64}; // sizes of b, t
  int b[6];  // x -> max y, y -> max x
  unsigned char ns[1024]; // state*4 -> next state if 0, if 1, n0, n1
  unsigned char t[N][N][2]={{{0}}};

int num_states(int x, int y) {
  if (x<y) return num_states(y, x);
  if (x<0 || y<0 || x>=N || y>=N || y>=B || x>=b[y]) return 0;
  return 1+(y>0 && x+y<b[5]);
}

// New value of count x if the opposite bit is observed
void discount(int& x) {
  int y=0;
  if (x>2){
    for (int i=1;i<mdc;i++) y+=x>=i;
    x=y;
  }
}

// compute next x,y (0 to N) given input b (0 or 1)
void next_state(int& x, int& y, int b) {
  if (x<y)
    next_state(y, x, 1-b);
  else {
    if (b) {
      ++y;
      discount(x);
    }
    else {
      ++x;
      discount(y);
    }
    while (!t[x][y][1]) {
      if (y<2) --x;
      else {
        x=(x*(y-1)+(y/2))/y;
        --y;
      }
    }
  }
}

// Initialize next state table ns[state*4] -> next if 0, next if 1, x, y
void generate() {
  memset(ns, 0, sizeof(ns));
  memset(t, 0, sizeof(t));
  // Assign states
  int state=0;
  for (int i=0; i<256; ++i) {
    for (int y=0; y<=i; ++y) {
      int x=i-y;
      int n=num_states(x, y);
      if (n) {
        t[x][y][0]=state;
        t[x][y][1]=n;
        state+=n;
      }
    }
  }

  // Print/generate next state table
  state=0;
  for (int i=0; i<N; ++i) {
    for (int y=0; y<=i; ++y) {
      int x=i-y;
      for (int k=0; k<t[x][y][1]; ++k) {
        int x0=x, y0=y, x1=x, y1=y;  // next x,y for input 0,1
        int ns0=0, ns1=0;
        next_state(x0, y0, 0);
        next_state(x1, y1, 1);
        ns[state*4]=ns0=t[x0][y0][0];
        ns[state*4+1]=ns1=t[x1][y1][0]+(t[x1][y1][1]>1);
        ns[state*4+2]=x;
        ns[state*4+3]=y;

          // uncomment to print table above
        //printf("{%3d,%3d,%2d,%2d},", ns[state*4], ns[state*4+1],
         // ns[state*4+2], ns[state*4+3]);
        //if (state%4==3) printf(" // %d-%d\n", state-3, state);
        if (state>0xff || t[x][y][1]==0 || t[x0][y0][1]==0 || t[x1][y1][1]==0) return;
        assert(state>=0 && state<256);
        assert(t[x][y][1]>0);
        assert(t[x][y][0]<=state);
        assert(t[x][y][0]+t[x][y][1]>state);
        assert(t[x][y][1]<=6);
        assert(t[x0][y0][1]>0);
        assert(t[x1][y1][1]>0);
        assert(ns0-t[x0][y0][0]<t[x0][y0][1]);
        assert(ns0-t[x0][y0][0]>=0);
        assert(ns1-t[x1][y1][0]<t[x1][y1][1]);
        assert(ns1-t[x1][y1][0]>=0);
        ++state;
        if (state>0xff) return;
      }
    }
  }
}

void __attribute__ ((noinline)) Init(int s0,int s1,int s2,int s3,int s4,int s5,int s6,U8 *table) {
    b[0]=s0;b[1]=s1;b[2]=s2;b[3]=s3;b[4]=s4;b[5]=s5;mdc=s6;
    generate(); 
    memcpy(table,  ns, 1024);
}
};
//State tables
U8 STA1[256][4];
U8 STA2[256][4];
U8 STA4[256][4];
U8 STA5[256][4];
U8 STA6[256][4];
U8 STA7[256][4];

// Dictionary reverse
int wfgets(char *str, int count, FILE  *fp) {
    int c, i = 0;
    while (i<count-1 && ((c=getc(fp))!=EOF)) {
        str[i++]=c; if (c=='\n')str[i-1]=0;
        if (c=='\n')
            break;
    }
    str[i]=0;
    return i;
}

char *s;
char *dictW[44515];
U8 dictWLen[44515];
char *dictM[472];
int codeword2sym[256]; 
int dict1size=80;
int dict2size=32;
int dict12size=dict1size*dict2size;
int sizeDict,sizeDict2;
// decoded codewords for wordmodel
U32 cwTEXT;           // text
U32 cwNOWIKI;         // nowiki
U32 cwMATH;           // math
U32 cwPRE;            // pre
U32 cwPAGE;           // page image category user wikipedia
U32 cwIMAGE;
U32 cwCATEGORY;
U32 cwUSER;
U32 cwWIKIPEDIA;
//
U32 cwTABLE;
U32 cwTD;
//
U32 cwSEE;
U32 cwALSO;
//
U32 cwEXTERNAL;
U32 cwLINKS;
//
U32 cwREFERENCES;     // References
U32 cwBIBLIOGRAPHY;   // Bibliography

U32 cwIS;
U32 cwWITH;
U32 cwTHE;
U32 cwON,cwIN;
U32 cwWWW,cwHTTP;
U32 cwISBN;

void loaddict(FILE  *file){
    int line_count=0,len=0;
    s=(char *)malloc(8192*8);
    while ((len=wfgets(s, 8192*8, file)) )  {
        dictW[line_count]=(char *)malloc(len);
        memcpy(dictW[line_count],  s, len);
        dictWLen[line_count]=U8(len-1);
        //printf("%d,%s\n",len,dictW[line_count]);
        if (cwTEXT==0 && strcmp(dictW[line_count], "text")==0) cwTEXT=line_count;
        else if (cwNOWIKI==0 && strcmp(dictW[line_count], "nowiki")==0) cwNOWIKI=line_count;
        else if (cwMATH==0 && strcmp(dictW[line_count], "math")==0) cwMATH=line_count;
        else if (cwPRE==0 && strcmp(dictW[line_count], "pre")==0) cwPRE=line_count;
        else if (cwPAGE==0 && strcmp(dictW[line_count], "page")==0) cwPAGE=line_count;
        else if (cwIMAGE==0 && strcmp(dictW[line_count], "image")==0) cwIMAGE=line_count;
        else if (cwCATEGORY==0 && strcmp(dictW[line_count], "category")==0) cwCATEGORY=line_count;
        else if (cwUSER==0 && strcmp(dictW[line_count], "user")==0) cwUSER=line_count;
        else if (cwWIKIPEDIA==0 && strcmp(dictW[line_count], "wikipedia")==0) cwWIKIPEDIA=line_count;
        else if (cwTABLE==0 && strcmp(dictW[line_count], "table")==0) cwTABLE=line_count;
        else if (cwTD==0 && strcmp(dictW[line_count], "td")==0) cwTD=line_count;
        
        else if (cwEXTERNAL==0 && strcmp(dictW[line_count], "external")==0) cwEXTERNAL=line_count;
        else if (cwLINKS==0 && strcmp(dictW[line_count], "links")==0) cwLINKS=line_count;
        else if (cwSEE==0 && strcmp(dictW[line_count], "see")==0) cwSEE=line_count;
        else if (cwALSO==0 && strcmp(dictW[line_count], "also")==0) cwALSO=line_count;
        else if (cwREFERENCES==0 && strcmp(dictW[line_count], "references")==0) cwREFERENCES=line_count;
        else if (cwBIBLIOGRAPHY==0 && strcmp(dictW[line_count], "bibliography")==0) cwBIBLIOGRAPHY=line_count;
        
        else if (cwIS==0 && strcmp(dictW[line_count], "is")==0) cwIS=line_count;
        else if (cwWITH==0 && strcmp(dictW[line_count], "with")==0) cwWITH=line_count;
        else if (cwTHE==0 && strcmp(dictW[line_count], "the")==0) cwTHE=line_count;
        
        else if (cwON==0 && strcmp(dictW[line_count], "on")==0) cwON=line_count;
        else if (cwIN==0 && strcmp(dictW[line_count], "in")==0) cwIN=line_count;
        
        else if (cwWWW==0 && strcmp(dictW[line_count], "www")==0) cwWWW=line_count;
        else if (cwHTTP==0 && strcmp(dictW[line_count], "http")==0) cwHTTP=line_count;
        
        else if (cwISBN==0 && strcmp(dictW[line_count], "isbn")==0) cwISBN=line_count;
        
        line_count++;
    }
    free(s);
    printf("Loaded %d words\n",line_count);
    //printf("Loaded %d %d %d %d %d %d %d %d %d, %d %d %d, %d %d,  %d %d, %d\n", cwTEXT,cwNOWIKI,cwMATH,cwPRE,cwPAGE,cwIMAGE,cwCATEGORY,cwUSER,cwWIKIPEDIA ,cwIS,cwWITH,cwTHE,cwON,cwIN,cwWWW ,cwHTTP,cwISBN);
    sizeDict=line_count;
}

inline int decodeCodeWord(int cw) {
    int i=0;
    int c=cw&255;
    if (codeword2sym[c]<dict1size) {
        i=codeword2sym[c];
        return i;
    }

    i=dict1size*(codeword2sym[c]-dict1size);
    c=(cw>>8)&255;

    if (codeword2sym[c]<dict1size) {
        i+=codeword2sym[c];
        return i+dict1size;
    }

    i=(i-dict12size)*dict2size;
    i+=dict1size*(codeword2sym[c]-dict1size);

    c=(cw>>16)&255;
    i+=codeword2sym[c];
    return i+80*49;
}

void dosym() {
    FILE *f=fopen("english.dic","rb");
    loaddict(f);
    fclose(f);
    for (int c=0; c<256; c++) {
        codeword2sym[c]=0;
    }
    int charsUsed=0;
    for (int c=128; c<256; c++) {
       codeword2sym[c]=charsUsed;
       charsUsed++;
    }
}

U32 lastCW=0;
void decodeWord(int c) {
        lastCW=decodeCodeWord(c);
        assert(lastCW>=0 && lastCW<44515);
}


struct Mix {
    int N;  // n
    int* wt;  // weights, scaled 24 bits
    int x1, x2;    // inputs, scaled 8 bits (-2047 to 2047)
    int cxt;  // last context (0..n-1)
    int pr;   // last output

    void Init(int n=256) {
       N=n, x1=0, x2=0, cxt=0, pr=0;
       alloc(wt, n*2);
       for (int i=0; i<N*2; ++i)
            wt[i]=1<<23;
    }
    int pp(int p1, int p2, int cx) {
        assert(cx>=0 && cx<N);
        cxt=cx*2;
        return pr=((x1=p1)*(wt[cxt]>>16)+(x2=p2)*(wt[cxt+1]>>16)+128)>>8;
    }
    void update(int y) {
        assert(y==0 || y==1);
        int err=((y<<12)-squash(pr));
        if ((wt[cxt]&3)<3)
            err*=4-(++wt[cxt]&3);
        err=(err+8)>>4;
        wt[cxt]+=x1*err&-4;
        wt[cxt+1]+=x2*err;
  }
};

// Mixer m(N, M, S=1, w=0) combines models using M neural networks with
//   N inputs each, of which up to S may be selected.  If S > 1 then
//   the outputs of these neural networks are combined using another
//   neural network (with parameters S, 1, 1).  If S = 1 then the
//   output is direct.  The weights are initially w (+-32K).
//   It is used as follows:
// m.update() trains the network where the expected output is the
//   last bit (in the global variable y).
// m.add(stretch(p)) inputs prediction from one of N models.  The
//   prediction should be positive to predict a 1 bit, negative for 0,
//   nominally +-256 to +-2K.  The maximum allowed value is +-32K but
//   using such large values may cause overflow if N is large.
// m.set(cxt, range) selects cxt as one of 'range' neural networks to
//   use.  0 <= cxt < range.  Should be called up to S times such
//   that the total of the ranges is <= M.
// m.p() returns the output prediction that the next bit is 1 as a
//   12 bit number (0 to 4095).

// Vector product a*b of n signed words, returning signed integer scaled down by 8 bits.
// n is rounded up to a multiple of 8.

//static int dot_product (const short* const t, const short* const w, int n);

// Train n neural network weights w[n] on inputs t[n] and err.
// w[i] += ((t[i]*2*err)+(1<<16))>>17 bounded to +- 32K.
// n is rounded up to a multiple of 8.
//int mix1SK=0,mix1TO=0;
struct Mixer1 { 
    int N, M;   // max inputs, max contexts, max context sets
    short *tx;  // N inputs from add()  
    short *wx ; // N*M weights
    short *ptr;
    int cxt;    // S contexts
    int pr;     // last result (scaled 12 bits)
    int shift1; 
    int elim;
    int uperr;
    int err;

    int dot_product (const short* const t, const short* const w, int n) {
        assert(n == ((n + 15) & -16));
        YMM sum = _mm256_setzero_si256 ();
        while ((n -= 16) >= 0) { // Each loop sums 16 products
            YMM tmp = _mm256_madd_epi16 (*(YMM *) &t[n], *(YMM *) &w[n]); // t[n] * w[n] + t[n+1] * w[n+1]
            tmp = _mm256_srai_epi32 (tmp, 8); //                                        (t[n] * w[n] + t[n+1] * w[n+1]) >> 8
            sum = _mm256_add_epi32 (sum, tmp); //                                sum += (t[n] * w[n] + t[n+1] * w[n+1]) >> 8
        } 
        sum =_mm256_hadd_epi32(sum,_mm256_setzero_si256 ());       //add [1]=[1]+[2], [2]=[3]+[4], [3]=0, [4]=0, [5]=[5]+[6], [6]=[7]+[8], [7]=0, [8]=0
        sum =_mm256_hadd_epi32(sum,_mm256_setzero_si256 ());       //add [1]=[1]+[2], [2]=0,       [3]=0, [4]=0, [5]=[5]+[6], [6]=0,       [7]=0, [8]=0
        XMM lo = _mm256_extractf128_si256(sum, 0);
        XMM hi = _mm256_extractf128_si256(sum, 1);
        XMM newsum = _mm_add_epi32(lo, hi);                    //sum last two
        return _mm_cvtsi128_si32(newsum);
    }

    void train (const short* const t, short* const w, int n, const int e) {
        assert(n == ((n + 15) & -16));
        if (e) {
            const YMM one = _mm256_set1_epi16 (1);
            const YMM err = _mm256_set1_epi16 (short(e));
            while ((n -= 16) >= 0) { // Each iteration adjusts 16 weights
                YMM tmp = _mm256_adds_epi16 (*(YMM *) &t[n], *(YMM *) &t[n]); // t[n] * 2
                tmp = _mm256_mulhi_epi16 (tmp, err); //                                     (t[n] * 2 * err) >> 16
                tmp = _mm256_adds_epi16 (tmp, one); //                                     ((t[n] * 2 * err) >> 16) + 1
                tmp = _mm256_srai_epi16 (tmp, 1); //                                      (((t[n] * 2 * err) >> 16) + 1) >> 1
                tmp = _mm256_adds_epi16 (tmp, *(YMM *) &w[n]); //                    ((((t[n] * 2 * err) >> 16) + 1) >> 1) + w[n]
                *(YMM *) &w[n] = tmp; //                                          save the new eight weights, bounded to +- 32K
            }
        }
    }

    /*
// dot_product returns dot product t*w of n elements.  n is rounded
// up to a multiple of 8.  Result is scaled down by 8 bits.
int dot_product(short *t, short *w, int n) {
int sum=0;
n=(n+15)&-16;
for (int i=0; i<n; i+=2)
    sum+=(t[i]*w[i]+t[i+1]*w[i+1]) >> 8;
return sum;
}

// Train neural network weights w[n] given inputs t[n] and err.
// w[i] += t[i]*err, i=0..n-1.  t, w, err are signed 16 bits (+- 32K).
// err is scaled 16 bits (representing +- 1/2).  w[i] is clamped to +- 32K
// and rounded.  n is rounded up to a multiple of 8.

void train(short *t, short *w, int n, int err) {
n=(n+15)&-16;
for (int i=0; i<n; ++i) {
    int wt=w[i]+(((t[i]*err*2>>16)+1)>>1);
    if (wt<-32768) wt=-32768;
    if (wt>32767) wt=32767;
    w[i]=wt;
}
}
*/

    // Adjust weights to minimize coding cost of last prediction
    void __attribute__ ((noinline)) update(int y) {
        err=((y<<12)-pr)*uperr/4;
        if (err>32767)
        err=32767;
        if (err<-32768)
        err=-32768;
        if(err>=-elim && err<=elim) err=0;
        train(&tx[0], &wx[cxt*N], N, err);
    }

    // predict next bit
    int __attribute__ ((noinline)) p() {
        assert(cxt>=0 && cxt<M);
        int dp=dot_product(&tx[0], &wx[cxt*N], N)*shift1>>11;
        return pr=squash(dp);
    }
    int __attribute__ ((noinline)) p1() {
        assert(cxt>=0 && cxt<M);
        int dp=dot_product(&tx[0], &wx[cxt*N], N)*shift1>>11;
        if (dp<-2047) {
            dp=-2047;
        }
        else if (dp>2047) {
            dp=2047;
        }
        pr=squash(dp);
        return dp;
    }
    void setTxWx(int n,short* mn) {
        N=n;
        alloc1(wx,(N*M)+32,ptr,32);
        tx=mn; 
        // Set bias
        for (int j=0; j<M*N; ++j) wx[j]=129;
    }
    void reset() {
        for (int j=0; j<M*N; ++j) wx[j]=wx[j]*2;
    }
    void Init(int m,  U32 s,U32 e,U32 ue) {
        M=m,  cxt=0, shift1=s,elim=e,uperr=ue;err=0;
        pr=2048; //initial p=0.5
    }
    void Free() {
        free(ptr);
    }
    void Print() {
        // print N weights averaged over context
        printf("Mixer(%d,%d): ", N, M);
        for (int i=0; i<N; ++i) {
            int w=0;
            for (int j=0; j<M; ++j)
            w+=wx[j*N+i];//,printf("%d ",wx[j*N+i]);;
            printf("%d ", w/M);
        }
        printf("\n");
    }
};

// A StateMap maps a context to a probability.

static int dt[1024];  // i -> 16K/(i+i+3)
// No update limit
struct StateMap {
    int N;        // Number of contexts
    int cxt;      // Context of last prediction
    U32 *t;       // cxt -> prediction in high 22 bits, count in low 10 bits
    int pr;
    const U8 *nn;
    int next(int i, int y){
        return nn[ y + i*4];
    }
    void __attribute__ ((noinline)) Init(int n, const U8 *nn1){
        nn=nn1;
        N=n, cxt=0, pr=2048;
        assert(ispowerof2(n));
        alloc(t,n);
        for (int i=0; i<N; ++i){
            U32 n0=next(i, 2)*3+1;
            U32 n1=next(i, 3)*3+1;
            t[i]=(((n1<<20) / (n0+n1)) << 12);
        }
    }
    void Free() {
        free(t);
    }
    inline void update() {    
        U32 *p=&t[cxt], p0=p[0];
        int pr1=p0>>14;  // count, prediction
        p0+=(x.y<<18)-pr1;
        p[0]=p0;
    }
    // update bit y (0..1), predict next bit in context cx
    void set(const int c) {  
        assert(cxt>=0 && cxt<N);
        update();
        pr=t[cxt=c]>>20;
    } 
}; 
// With update limit
struct StateMap1 {
    int N;        // Number of contexts
    int cxt;      // Context of last prediction
    U32 *t;       // cxt -> prediction in high 22 bits, count in low 10 bits
    int pr;
    int mask;
    int limit; 
    void __attribute__ ((noinline)) Init(int n, int lim) {
        N=n, cxt=0, pr=2048, mask=n-1,limit=lim;
        assert(ispowerof2(n));
        alloc(t,n);
        assert(limit>0 && limit<1024);
        for (int i=0; i<N; ++i)
        t[i]=1<<31;
    }
    void Free() {
        free(t);
    }
    inline void update() {    
        U32 *p=&t[cxt], p0=p[0];
        int n=p0&1023, pr1=p0>>12;  // count, prediction
        p0+=(n<limit);
        p0+=(((((x.y<<20)-pr1)))*dt[n]+512)&0xfffffc00;
        p[0]=p0;
    }
    // update bit y (0..1), predict next bit in context cx
    void set(const int c) {  
        assert(cxt>=0 && cxt<N);
        update();
        pr=t[cxt=(c&mask)]>>20;
    } 
}; 

// A RunContextMap maps a context into the next byte and a repeat
// count up to M.  Size should be a power of 2.  Memory usage is 3M/4.
struct RunContextMap {
    enum {B=4,M=4}; 
    U8 *t;   // hash t
    U8 *ptr;
    U8* cp;
    short rc[512];
    U8 tmp[B];
    U32 n;
    void Init(int m,int rcm_ml=8){ 
        alloc1(t,m,ptr,64);  
        n=(m/B-1);
        for (int r=0;r<B;r++) tmp[r]=0;
        cp=&t[0]+1;
        for (int r=0;r<256;r++) {
            int c=ilog[r]*8;
            if ((r&1)==0) c=c*rcm_ml/4;
            rc[r+256]=clp(c);
            rc[r]=clp(-c);
        }
    }
    void Free() {
        free(ptr);
    }
    void __attribute__ ((noinline)) set(U32 cx,U8 c1) {  // update count
        if (cp[0]==0) cp[0]=2, cp[1]=c1;
        else if (cp[1]!=c1) cp[0]=1, cp[1]=c1;
        else if (cp[0]<254) cp[0]=cp[0]+2;
        cp=find(cx)+1;
    }
    int p() {  // predict next bit
        int b=x.c0shift_bpos ^ (cp[1] >> x.bposshift);
        if (b<=1)
        return rc[b*256+cp[0]];
        else
        return 0;
    }
    int mix() {  // return run length
        x.mxInputs1.add(p());
        return cp[0]!=0;
    }

    inline  U8* find(U32 i) {
        U16 chk=(i>>16^i)&0xffff;
        i=i*M&n;
        U8 *p;
        U16 *cp1;
        int j;
        for (j=0; j<M; ++j) {
            p=&t[(i+j)*B];
            cp1=(U16*)p;
            if (p[2]==0) {*cp1=chk;break;}
            if (*cp1==chk) break;  // found
        }
        if (j==0) return p+1;  // front
        if (j==M) {
            --j;
            memset(&tmp, 0, B);
            memmove(&tmp, &chk, 2);
            if (M>2 && t[(i+j)*B+2]>t[(i+j-1)*B+2]) --j;
        }
        else memcpy(&tmp, cp1, B);
        memmove(&t[(i+1)*B], &t[i*B], j*B);
        memcpy(&t[i*B], &tmp, B);
        return &t[i*B+1];
    }
};

// Map for modelling contexts of (nearly-)stationary data.
// The context is looked up directly. For each bit modelled, a 16bit prediction is stored.
// The adaptation rate is controlled by the caller, see mix().

// - BitsOfContext: How many bits to use for each context. Higher bits are discarded.
// - InputBits: How many bits [1..8] of input are to be modelled for each context.
// New contexts must be set at those intervals.

// Uses (2^(BitsOfContext+1))*((2^InputBits)-1) bytes of memory.
int sscmrate=0;
struct SmallStationaryContextMap {
    U16 *Data;
    int Context, Mask, Stride, bCount, bTotal, B, N;
    U16 *cp;

    void __attribute__ ((noinline)) Init(int BitsOfContext,  int InputBits=8) {
        assert(InputBits>0 && InputBits<=8);
        Context=0, Mask=(1<<BitsOfContext)-1, 
        Stride=(1<<InputBits)-1, bCount=0, bTotal=InputBits, B=0;
        N=(1ull<<BitsOfContext)*((1ull<<InputBits)-1);
        alloc(Data,N);
        for (int i=0; i<N; ++i)
        Data[i]=0x7FFF;
        cp=&Data[0];
    }
    void Free() {
        free(Data);
    }
    void set(U32 ctx) {
        Context = (ctx&Mask)*Stride;
        bCount=B=0;
    }
    void __attribute__ ((noinline)) mix(int r) {
        int rate =r +7; const int Multiplier=1;const int Divisor=4;
        *cp+=((x.y<<16)-(*cp)+(1<<(rate-1)))>>rate;
        B+=(x.y && B>0);
        cp = &Data[Context+B];
        int Prediction = (*cp)>>4;
        x.mxInputs1.add((stretch(Prediction)*Multiplier)/Divisor);
        x.mxInputs1.add(((Prediction-2048)*Multiplier)/(Divisor*2));
        prediction_index--;
        bCount++; B+=B+1;
        if (bCount==bTotal)
        bCount=B=0;
    }
};



// Context map for large contexts.  Most modeling uses this type of context
// map.  It includes a built in RunContextMap to predict the last byte seen
// in the same context, and also bit-level contexts that map to a bit
// history state.
//
// Bit histories are stored in a hash table.  The table is organized into
// 64-byte buckets alinged on cache page boundaries.  Each bucket contains
// a hash chain of 7 elements, plus a 2 element queue (packed into 1 byte)
// of the last 2 elements accessed for LRU replacement.  Each element has
// a 2 byte checksum for detecting collisions, and an array of 7 bit history
// states indexed by the last 0 to 2 bits of context.  The buckets are indexed
// by a context ending after 0, 2, or 5 bits of the current byte.  Thus, each
// byte modeled results in 3 main memory accesses per context, with all other
// accesses to cache.
//
// On bits 0, 2 and 5, the context is updated and a new bucket is selected.
// The most recently accessed element is tried first, by comparing the
// 16 bit checksum, then the 7 elements are searched linearly.  If no match
// is found, then the element with the lowest priority among the 5 elements
// not in the LRU queue is replaced.  After a replacement, the queue is
// emptied (so that consecutive misses favor a LFU replacement policy).
// In all cases, the found/replaced element is put in the front of the queue.
//
// The priority is the state number of the first element (the one with 0
// additional bits of context).  The states are sorted by increasing n0+n1
// (number of bits seen), implementing a LFU replacement policy.
//
// When the context ends on a byte boundary (bit 0), only 3 of the 7 bit
// history states are used.  The remaining 4 bytes implement a run model
// as follows: <count:7,d:1> <b1> <unused> <unused> where <b1> is the last byte
// seen, possibly repeated.  <count:7,d:1> is a 7 bit count and a 1 bit
// flag (represented by count * 2 + d).  If d=0 then <count> = 1..127 is the
// number of repeats of <b1> and no other bytes have been seen.  If d is 1 then
// other byte values have been seen in this context prior to the last <count>
// copies of <b1>.
//
// As an optimization, the last two hash elements of each byte (representing
// contexts with 2-7 bits) are not updated until a context is seen for
// a second time.  This is indicated by <count,d> = <1,0> (2).  After update,
// <count,d> is updated to <2,0> or <1,1> (4 or 3).

inline int sc(int p){
    if (p>0) return p>>7;
    return (p+127)>>7;// p+((1<<s)-1);
}

// A BH maps a 32 bit hash to an array of B bytes (checksum and B-2 values)
//
// BH bh(N); creates N element table with B bytes each.
//   N must be a power of 2.  The first byte of each element is
//   reserved for a checksum to detect collisions.  The remaining
//   B-1 bytes are values, prioritized by the first value.  This
//   byte is 0 to mark an unused element.
//
// bh[i] returns a pointer to the i'th element, such that
//   bh[i][0] is a checksum of i, bh[i][1] is the priority, and
//   bh[i][2..B-1] are other values (0-255).
//   The low lg(n) bits as an index into the table.
//   If a collision is detected, up to M nearby locations in the same
//   cache line are tested and the first matching checksum or
//   empty element is returned.
//   If no match or empty element is found, then the lowest priority
//   element is replaced.

// 2 byte checksum with LRU replacement (except last 2 by priority)

template <const int A, const int B> // Warning: values 3, 7 for A are the only valid parameters
union  E {  // hash element, 64 bytes
    struct{ // this is bad uc
        U16 chk[A];  // byte context checksums
        U8 last;     // last 2 accesses (0-6) in low, high nibble
        U8 bh[A][7]; // byte context, 3-bit context -> bit history state
        // bh[][0] = 1st bit, bh[][1,2] = 2nd bit, bh[][3..6] = 3rd bit
        // bh[][0] is also a replacement priority, 0 = empty
        //  U8* get(U16 chk);  // Find element (0-6) matching checksum.
        // If not found, insert or replace lowest priority (not last).
    };
    U8 pad[B] ;
    __attribute__ ((noinline)) U8* get(U16 ch,int keep) {
        if (chk[last&15]==ch) return &bh[last&15][0];
        int b=0xffff, bi=0;

        for (int i=0; i<A; ++i) {
            if (chk[i]==ch) return last=last<<4|i, (U8*)&bh[i][0];
            int pri=bh[i][0];
            if (pri<b && (last&15)!=i && last>>4!=i) b=pri, bi=i;
        }
        return last=last<<4|bi|keep, chk[bi]=ch, (U8*)memset(&bh[bi][0], 0, 7);
    }
};

inline U32 getStateByteLocation(const int bpos, const int c0) {
    U32 pis = 0; //state byte position in slot
    const U32 smask = (U32(0x31031010) >> (bpos << 2)) & 0x0F;
    pis = smask + (c0 & smask);
    return pis;
}

#define MAXCXT 8
short st2_p0[4096];
short st2_p1[4096];
short rcpr[512]; //2-6 0-4
bool doCMprint=false;
template <const int A, const int B> // Warning: values 3, 7 for A are the only valid parameters
union  E1 {  // hash element, 64 bytes
    struct{ // this is bad uc
        U16 chk[A];  // byte context checksums
        U8 last;     // last 2 accesses (0-6) in low, high nibble
        U8 bh[A][7]; // byte context, 3-bit context -> bit history state
        // bh[][0] = 1st bit, bh[][1,2] = 2nd bit, bh[][3..6] = 3rd bit
        // bh[][0] is also a replacement priority, 0 = empty
        //  U8* get(U16 chk);  // Find element (0-6) matching checksum.
        // If not found, insert or replace lowest priority (not last).
    };
    U8 pad[B] ;
    __attribute__ ((noinline)) U8* get(U16 ch,int keep) {
        if (chk[last&15]==ch) return &bh[last&15][0];
        int b=0xffff, bi=0;

        for (int i=0; i<A; ++i) {
            if (chk[i]==ch) return last=last<<4|i, (U8*)&bh[i][0];
            int pri=bh[i][0];
            if (pri<b && (last&15)!=i && last>>4!=i) b=pri, bi=i;
        }
        return last=last<<4|bi|keep, chk[bi]=ch, (U8*)memset(&bh[bi][0], 0, 7);
    }
};

struct ContextMap3 {
    int C;  // max number of contexts
    U8* cp[MAXCXT];   // C pointers to current bit history
    U8* cp0[MAXCXT];  // First element of 7 element array containing cp[i]
    U32 cxt[MAXCXT];  // C whole byte contexts (hashes)
    U8* runp[MAXCXT]; // C [0..3] = count, value, unused, unused
    int cn;          // Next context to set by set()
    int result;
    short st1[4096];
    short *st2;
    short st32[256];
    int cms,cms3,cms4;
    int kep;
    const U8 *nn;
    E1<14,128> *ptr,*t;  // Double sized BH
    U32 tmask;
    int skip2;
    int skip3;
    U16 cxtMask;
    //state
    int cxtn[MAXCXT];    // Context of last prediction
    U32 *ts;             // cxt -> prediction in high 22 bits, count in low 10 bits
    int sti;
    inline U8  next(int i, int y) {
        return nn[ y + i*4];
    }

    inline int pre(const int state) {
        assert(state>=0 && state<256);
        U32 n0=next(state, 2)*3+1;
        U32 n1=next(state, 3)*3+1;
        return (n1<<12) / (n0+n1);
    }
    inline void update(const int i) {    
        U32 *p=&ts[cxtn[i]], p0=p[0];
        const int pr1=p0>>14;
        p[0]+=(x.y<<18)-pr1;
    }
    int set(const int c,int i) {  
        assert(c>=0 && c<256);
        const int  pr=ts[c]>>20; // predict from current state
        // look if new state is same as state before, if so skip this state
        // if first state, set it
        if (i==0) {
            cxtn[sti++]=c;
            return pr;
        }
        for (int j=0; j<sti; j++) {
            if (cxtn[j]==c) return pr; // skip if same
        }
        cxtn[sti++]=c;
        return pr;
    } 
    void upd() {
        for (int j=0; j<sti; j++) {
            update(j);
        }
    } 
    // Construct using m bytes of memory for c contexts(c+7)&-8
    void __attribute__ ((noinline)) Init(U32 m1, int c1, int s2, int s3,const U8 *nn1,int cs4,int k, const int u,short *st) {
        C=c1;
        int m=m1*2;
        tmask=((m>>7)-1); 
        cn=result=0;
        cxtMask=((1<C)-1)*2; // Inital zero contexts
        kep=k;
        alloc1(t,(m>>7)+64*2,ptr,128);  
        nn=nn1;        
        cms=s2;               // mix prediction mul value
        cms4=cs4;
        cms3=s3;
        skip2=u;
        assert(m>=64 && (m&m-1)==0);  // power of 2?
        assert(sizeof(E1<14,128>)==128);
        //state
        alloc(ts,256);
        for (int i=0; i<256; ++i){
            U32 n0=next(i, 2)*3+1;
            U32 n1=next(i, 3)*3+1;
            ts[i]=(((n1<<20) / (n0+n1)) << 12);
        }
        //cm
        for (int i=0; i<C; ++i) {
            cp0[i]=cp[i]=&t[0].bh[0][0];
            runp[i]=cp[i]+3;
        }

        st2=st;
        // precalc mix3 mixer inputs
        for (int i=0;i<4096;i++) {
            st1[i]=clp(sc(cms*stretch(i)));
        } 

        for (int s=0;s<256;s++) {
            int n0=-!next(s,2);
            int n1=-!next(s,3);
            int r=0;
            int sp0=0;
            if ((n1-n0)==1 ) sp0=0,r=1;
            if ((n1-n0)==-1 ) sp0=4095,r=1;
            if (r) {
                int st8 =clp(sc((cms4)*(pre(s)-sp0)));
                st32[s]=clp(sc((cms3)*stretch(pre(s))));
                if (s<8) st32[s]=st8;
                else st32[s]=(st8+st32[s])>>1;
            }else{
                st32[s]=0;
            }
        }
    }
    void reset() {
        memset((void*)ptr, 0, tmask*sizeof(E1<14,128>));
        for (int i=0; i<C; ++i) {
            cp0[i]=cp[i]=&t[0].bh[0][0];
            runp[i]=cp[i]+3;
        }
        cxtMask=((1<C)-1)*2;
    }
    void Free() {
        free(ts);
        free(ptr);
    }

    // Set the i'th context to cx
    inline void set(U32 cx) {
        int i=cn++;
        assert(i>=0 && i<C);
        cx=cx*987654323+i;  // permute (don't hash) cx to spread the distribution
        cx=cx<<16|cx>>16;
        cxt[i]=cx*123456791+i;
        cxtMask=cxtMask*2;
    }

    inline void sets() {
        int i=cn++;
        assert(i>=0 && i<C);
        cxtMask=cxtMask+1; cxtMask=cxtMask*2;
    }
    // Predict to mixer m from bit history state s, using sm to map s to
    // a probability.
    inline int mix3(const int s, int i) {
        if (s==0) {
            x.mxInputs1.add(0);
            if (skip2==1)x.mxInputs1.add(0);
            x.mxInputs1.add(0);
            return 0;
        } else {
            const int p1=set(s,i);
            x.mxInputs1.add(st1[p1]);
            if (skip2==1)x.mxInputs1.add(st2[p1]);
            x.mxInputs1.add(st32[s]);
            return 1;
        }
    }
    // Zero prediction
    inline void mix4() {
        x.mxInputs1.add(0);
        if (skip2==1) x.mxInputs1.add(0);
        x.mxInputs1.add(0);
        x.mxInputs1.add(0);
    }
    // Update the model with bit y1, and predict next bit to mixer m.
    int __attribute__ ((noinline)) mix() {
        // Update model with y
        result=0;
        upd(); // update statemap
        sti=0;
        for (int i=0; i<cn; ++i) {
            if ((cxtMask>>(cn-i))&1) {
                mix4();
            } else {
                if (cp[i]) {
                    assert(cp[i]>=&t[0].bh[0][0] && cp[i]<=&t[tmask].bh[14][6]);
                    assert(((long long)(cp[i])&127)>=29);
                    *cp[i]=next(*cp[i], x.y);
                }

                // Update context pointers
                int s=0;
                if (x.bpos>1 && runp[i][0]==0) {
                    cp[i]=0;
                } else {
                    const U16 chksum=(cxt[i]>>16)^i;
                    if ( x.bpos){     
                        if ( x.bpos==2 ||  x.bpos==5)cp0[i]=cp[i]=t[(cxt[i]+x.c0)&tmask].get(chksum,kep);
                        else cp[i]=cp0[i]+x.cmBitState;
                    } else {// default
                        cp0[i]=cp[i]=t[(cxt[i]+x.c0)&tmask].get(chksum,kep);
                        // Update pending bit histories for bits 2-7
                        if (cp0[i][3]==2) {
                            const int c=cp0[i][4]+256;
                            U8 *p=t[(cxt[i]+(c>>6))&tmask].get(chksum,kep);
                            p[0]=1+((c>>5)&1);
                            p[1+((c>>5)&1)]=1+((c>>4)&1);
                            p[3+((c>>4)&3)]=1+((c>>3)&1);
                            p=t[(cxt[i]+(c>>3))&tmask].get(chksum,kep);
                            p[0]=1+((c>>2)&1);
                            p[1+((c>>2)&1)]=1+((c>>1)&1);
                            p[3+((c>>1)&3)]=1+(c&1);
                            cp0[i][6]=0;
                        }
                        const U8 c1=x.c4;
                        // Update run count of previous context
                        if (runp[i][0]==0)  // new context
                        runp[i][0]=2, runp[i][1]=c1;
                        else if (runp[i][1]!=c1)  // different byte in context
                        runp[i][0]=1, runp[i][1]=c1;
                        else if (runp[i][0]<254)  // same byte in context
                        runp[i][0]+=2;
                        runp[i]=cp0[i]+3;
                    }
                    s = *cp[i];
                }
                // predict from bit context
                result=result+mix3(s, i);
                // predict from last byte in context
                int b=x.c0shift_bpos ^ (runp[i][1] >> x.bposshift);
                
                if (b<=1) {
                    b=b*256;   // predicted bit + for 1, - for 0
                    // count*2, +1 if 2 different bytes seen
                    x.mxInputs1.add(rcpr[runp[i][0]+b]);
                }
                else
                x.mxInputs1.add(0);
            }
        }
        if ( x.bpos==7) {
            assert(cn==0 || cn==C),cn=cxtMask=0;
        }
        return result;
    }
};

struct ContextMap4 {
    int C;  // max number of contexts
    U8* cp[MAXCXT];   // C pointers to current bit history
    U8* cp0[MAXCXT];  // First element of 7 element array containing cp[i]
    U32 cxt[MAXCXT];  // C whole byte contexts (hashes)
    U8* runp[MAXCXT]; // C [0..3] = count, value, unused, unused
    int cn;          // Next context to set by set()
    int result;
    short st32[256];
    short st9[256];
    int cms,cms3,cms4;
    int kep;
    const U8 *nn;
    E<3,32> *ptr,*t;  // Half sized BH
    U32 tmask;
    int skip2;
    U16 cxtMask;

    inline U8  next(int i, int y) {
        return nn[ y + i*4];
    }
    inline int pre(const int state) {
        assert(state>=0 && state<256);
        U32 n0=next(state, 2)*3+1;
        U32 n1=next(state, 3)*3+1;
        return (n1<<12) / (n0+n1);
    }

    // Construct using m bytes of memory for c contexts(c+7)&-8
    void __attribute__ ((noinline)) Init(U32 m, int c1, int s2, int s3,const U8 *nn1,int cs4,int k, const int u) {
        C=c1;
        tmask=((m>>6)-1); 
        cn=0;
        cxtMask=((1<C)-1)*2; // Inital zero contexts
        result=0;
        kep=k;
        alloc1(t,(m>>6)+64,ptr,64);  
        nn=nn1;        
        cms=s2;               // mix prediction mul value
        cms4=cs4;
        cms3=s3;
        skip2=u;
        assert(m>=64 && (m&m-1)==0);  // power of 2?
        assert(sizeof(E<3,32>)==32);

        for (int i=0; i<C; ++i) {
            cp0[i]=cp[i]=&t[0].bh[0][0];
            runp[i]=cp[i]+3;
        }

        // precalc mix3 mixer inputs
        for (int i=0;i<256;i++) {
            st9[i]=clp(sc(18*(pre(i) - 2048)));
        } 

        for (int s=0;s<256;s++) {
            int n0=-!next(s,2);
            int n1=-!next(s,3);
            int r=0;
            int sp0=0;
            if ((n1-n0)==1 ) sp0=0,r=1;
            if ((n1-n0)==-1 ) sp0=4095,r=1;
            if (r) {
                int st8 =clp(sc((cms4)*(pre(s)-sp0)));
                st32[s]=clp(sc((cms3)*stretch(pre(s))));
                if (s<8) st32[s]=st8;
                else st32[s]=(st8+st32[s])>>1;
            } else {
                st32[s]=0;
            }
        }
    }
    void reset() {
        memset((void*)ptr, 0, tmask*sizeof(E<3,32>));
        for (int i=0; i<C; ++i) {
            cp0[i]=cp[i]=&t[0].bh[0][0];
            runp[i]=cp[i]+3;
        }
        cxtMask=((1<C)-1)*2;
    }
    void Free() {
        free(ptr);
    }

    // Set the i'th context to cx
    inline void set(U32 cx) {
        int i=cn++;
        assert(i>=0 && i<C);
        cx=cx*987654323+i;  // permute (don't hash) cx to spread the distribution
        cx=cx<<16|cx>>16;
        cxt[i]=cx*123456791+i;
        cxtMask=cxtMask*2;
    }

    inline void sets() {
        int i=cn++;
        assert(i>=0 && i<C);
        cxtMask=cxtMask+1; cxtMask=cxtMask*2;
    }
    // Predict to mixer m from bit history state s, using sm to map s to
    // a probability.
    inline int mix3(const int s) {
        if (s==0){
            if (skip2==1) x.mxInputs1.add(0);
            x.mxInputs1.add(0);
            return 0;
        } else {
            if (skip2==1) x.mxInputs1.add(st9[s]);
            x.mxInputs1.add(st32[s]);
            return 1;
        }
    }

    // Zero prediction
    inline void mix4() {
        if (skip2==1) x.mxInputs1.add(0);
        x.mxInputs1.add(0);
        x.mxInputs1.add(0);
    }
    // Update the model with bit y1, and predict next bit to mixer m.
    int __attribute__ ((noinline)) mix() {
        // Update model with y
        result=0;
        for (int i=0; i<cn; ++i) {
            if ((cxtMask>>(cn-i))&1) {
                mix4();
            } else {
                if (cp[i]) {
                    assert(cp[i]>=&t[0].bh[0][0] && cp[i]<=&t[tmask].bh[3][6]);
                    assert(((long long)(cp[i])&31)>=7);
                    *cp[i]=next(*cp[i], x.y);
                }

                // Update context pointers
                int s = 0;
                if ( x.bpos>1 && runp[i][0]==0) {
                    cp[i]=0;
                } else {
                    const U16 chksum=(cxt[i]>>16)^i;
                    if ( x.bpos){     
                        if ( x.bpos==2 ||  x.bpos==5)cp0[i]=cp[i]=t[(cxt[i]+x.c0)&tmask].get(chksum,kep);
                        else cp[i]=cp0[i]+x.cmBitState;
                    } else {// default
                        cp0[i]=cp[i]=t[(cxt[i]+x.c0)&tmask].get(chksum,kep);
                        // Update pending bit histories for bits 2-7
                        if (cp0[i][3]==2) {
                            const int c=cp0[i][4]+256;
                            U8 *p=t[(cxt[i]+(c>>6))&tmask].get(chksum,kep);
                            p[0]=1+((c>>5)&1);
                            p[1+((c>>5)&1)]=1+((c>>4)&1);
                            p[3+((c>>4)&3)]=1+((c>>3)&1);
                            p=t[(cxt[i]+(c>>3))&tmask].get(chksum,kep);
                            p[0]=1+((c>>2)&1);
                            p[1+((c>>2)&1)]=1+((c>>1)&1);
                            p[3+((c>>1)&3)]=1+(c&1);
                            cp0[i][6]=0;
                        }
                        const U8 c1=x.c4;
                        // Update run count of previous context
                        if (runp[i][0]==0)  // new context
                        runp[i][0]=2, runp[i][1]=c1;
                        else if (runp[i][1]!=c1)  // different byte in context
                        runp[i][0]=1, runp[i][1]=c1;
                        else if (runp[i][0]<254)  // same byte in context
                        runp[i][0]+=2;
                        runp[i]=cp0[i]+3;
                    }
                    s = *cp[i];
                }
                // predict from bit context
                result=result+mix3(s);
                
                // predict from last byte in context
                int b=x.c0shift_bpos ^ (runp[i][1] >> x.bposshift);
                if (b<=1) {
                    b=b*256;   // predicted bit + for 1, - for 0
                    // count*2, +1 if 2 different bytes seen
                    x.mxInputs1.add(rcpr[runp[i][0]+b]);
                }
                else
                x.mxInputs1.add(0);
            }
        }
        if ( x.bpos==7) {
            assert(cn==0 || cn==C);
            cn=cxtMask=0;
        }
        return result;
    }
};


// APM maps a probability and a context into a new probability
// that bit y will next be 1.  After each guess it updates
// its state to improve future guesses.  Methods:
//
// APM a(N) creates with N contexts, uses 66*N bytes memory.
// a.p(pr, cx, rate=8) returned adjusted probability in context cx (0 to
//   N-1).  rate determines the learning rate (smaller = faster, default 8).
//   Probabilities are scaled 16 bits (0-65535).
template <const int B=8>
struct  APM {
    int index;     // last p, context
    static constexpr int S=(1<<B);
    static constexpr int mask=(1<<B)-1;
    U16 t[S*33];        // [N][33]:  p, context -> p
    
    int p(int pr=2048, int cxt=0, int rate=8, int y=0) {
        pr=stretch(pr);
        int g=(y<<16)+(y<<rate)-y*2;
        t[index]   += (g-t[index])   >> rate;
        t[index+1] += (g-t[index+1]) >> rate;
        const int w=pr&127;  // interpolation weight (33 points)
        index=((pr+2048)>>7)+(cxt&mask)*33;
        return clp1((t[index]*(128-w)+t[index+1]*w) >> 11);
    }

    // maps p, cxt -> p initially
    void Init() {
        index=0;
        for (int j=0; j<33; ++j) t[j]=squash((j-16)*128)*16;
        for (int i=33; i<S*33; ++i) t[i]=t[i-33];
    }
};


int buf(int i);
int bufr(int i);
int pos;

struct StationaryMap {
    U32 *Data;
    int Context, Mask, Stride, bCount, bTotal, B, N;
    U32 *cp;
    int Multiplier;
    void Init(int BitsOfContext, int InputBits=8, int mul=8,int Rate=0) {
        Multiplier=mul;
        N=((1ull<<BitsOfContext)*((1ull<<InputBits)-1));
        Context=0, Mask=(1<<BitsOfContext)-1, Stride=(1<<InputBits)-1, bCount=0, bTotal=InputBits, B=0; 
        assert(InputBits>0 && InputBits<=8);
        assert(BitsOfContext+InputBits<=24);
        alloc(Data,N);
        for (int i=0; i<N; ++i)
        Data[i]=(0x7FF<<20)|min(1023,Rate);
        cp=&Data[0];
    }
    void Free() {
        free(Data);
    }
    void set(U32 ctx) {
        Context = (ctx&Mask)*Stride;
        bCount=B=0;
    }
    void mix() {
        // update
        int Prediction;
        U32 p0=cp[0];
        int n=p0&1023, pr=p0>>13;  // count, prediction
        p0+=(n<1023);     
        p0+=(((x.y<<19)-pr))*dt[n]&0xfffffc00;  
        cp[0]=p0;
        // predict
        B+=(x.y && B>0);
        cp=&Data[Context+B];
        Prediction = (*cp)>>20;
        x.mxInputs1.add((stretch(Prediction)*Multiplier)/32);//      1/4    8/32
        x.mxInputs1.add(((Prediction-2048)*Multiplier)/(32*2));//    1/8    8/64
        bCount++; B+=B+1;
        if (bCount==bTotal)
        bCount=B=0;
    }
};

// Map 8 bit byte to 2 bit value (3 - upper 2 bits, adjusted)
static const U8 wrt_2b[256]={
2, 3, 1, 3, 3, 0, 1, 2, 3, 3, 0, 0, 1, 3, 3, 3, 
3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0, 3, 3, 3, 3, 
3, 2, 0, 2, 1, 3, 2, 1, 3, 3, 3, 3, 2, 3, 0, 2, // _!"#$%&'()*+,-./
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 2, 2, 3, 2, 2, // 0123456789:;<=>?
2, 2, 0, 0, 2, 3, 1, 2, 1, 2, 2, 2, 2, 2, 0, 0, // @ABCDEFGHIJKLMNO
2, 2, 2, 2, 2, 2, 2, 2, 3, 0, 2, 3, 2, 0, 2, 3, // PQRSTUVWXYZ[\]^_

1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, //  abcdefghijklmno
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, // pqrstuvwxyz{|}~
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

// Map 8 bit byte to 3 bit value (upper 3 bits, adjusted)
static const U8 wrt_3b[256]={
0, 0, 2, 0, 5, 6, 0, 6, 0, 2, 0, 4, 3, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
2, 4, 1, 4, 4, 7, 4, 7, 3, 7, 2, 2, 3, 5, 3, 1, // _!"#$%&'()*+,-./
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 5, 3, 3, 5, 5, // 0123456789:;<=>?
0, 5, 5, 7, 5, 0, 1, 5, 4, 5, 0, 0, 6, 0, 7, 1, // @ABCDEFGHIJKLMNO
3, 3, 7, 4, 5, 5, 7, 0, 2, 2, 5, 4, 4, 7, 4, 6, // PQRSTUVWXYZ[\]^_

5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7};

const U8 wrt_4b[256]={
 6, 0,12,15,12,15,14,14, 5, 3,14, 0,15,13, 8,13,
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
13, 5,15,11,10,12, 6,12, 0,11,14, 1, 1,10, 9, 8,
 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 9,11, 6, 1, 0, 4,
 9,10,10, 4, 5, 1, 4, 2,11, 8, 4, 1, 0,10,10, 5,
 4, 7,15, 4, 5,13, 0, 1, 4,12, 0, 1, 3, 3, 3,11,
 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 8, 0,11, 7,

 
 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2
 };

//wrt
#define COLON         'J' // :
#define SEMICOLON     'K' // ;
#define LESSTHAN      'L' // <
#define EQUALS        'M' // =
#define GREATERTHAN   'N' // >
#define QUESTION      'O' // ?
#define FIRSTUPPER     64 // @ - wrt first char in word is in upper case
#define SQUAREOPEN     91 // [
#define BACKSLASH      92 // '\'
#define SQUARECLOSE    93 // ]
#define CURLYOPENING  'P' // {
#define VERTICALBAR   'Q' // |
#define CURLYCLOSE    'R' // }
#define CHARSWAP


#define APOSTROPHE    39  // '
#define QUOTATION     34  // "
#define SPACE         32  // ' '
#define HTLINK        31  // http link
#define HTML          30  // 
#define LF            10  // 
#define ESCAPE        12  // 
#define UPPER         7   // Upper case word
#define TEXTDATA      96  // Any other char, probably text

// Vector for different contexts

template <typename T = int,const int S=256 >
struct vec {
    T cxt[S];
    static constexpr int capacity=S;
    int size=0;
};

template <typename T = int,const int S>
void vec_new(vec<T,S>* o) {
    o->size=0;
}
template <typename T = int,const int S>
int vec_size(vec<T,S> *o) {
    return o->size;
}
template <typename T  = int,const int S>
void vec_push( vec<T,S> *o, const T element){
    assert(o->size<S);
    o->cxt[o->size++]=element;
    o->size=o->size&(o->capacity-1); // roll over
}
template <typename T  = int,const int S>
int vec_at(vec<T,S> *o, const int index) {
    assert(index<S);
    return o->cxt[index];
}
template <typename T  = int,const int S>
T &vec_ref(vec<T,S> *o, const int index) {
    assert(index<S);
    return o->cxt[index];
}
template <typename T  = int,const int S>
void vec_i(vec<T,S> *o, const int index) {
    o->cxt[index]++;
}
template <typename T  = int,const int S>
void vec_pop(vec<T,S> *o) {
    assert(o->size<S);
    if (o->size>0) o->cxt[o->size]=0, o->size--; // no rollback
}
template <typename T  = int,const int S>
void vec_reset(vec<T,S> *o) {
    o->cxt[0]=0;
    o->size=0;
}
template <typename T  = int,const int S>
bool vec_empty(vec<T,S> *o) {
    return (o->size==0)?true:false;
}
template <typename T  = int,const int S>
int vec_prev(vec<T,S> *o) {
    assert(o->size>=0);
    return (o->size>1)?(o->cxt[o->size-2]):0; // no rollback
}
// This part is based on cmix BracketContext
template <typename T  = U8>
struct BracketContext {
    U32 context;           // bracket byte and distance
    vec<int,512> active;            // vector for brackets, max 512 elements
    vec<int,512> distance;          // vector for distance, max 512 elements
    const T *element;           
    int elementCount;
    bool doPop;            // set true for quotes
    int limit;
    T cxt,dst;

    void Init(const T *d,const int e,int pop=false,int l=(1<<(sizeof(T)*8))) {
        elementCount=e;
        element=d;
        context=cxt=dst=0;
        doPop=pop;
        limit=l;
        vec_new(&active);
        vec_new(&distance);
         Reset();
    }
    void __attribute__ ((noinline)) Reset() {
        vec_reset(&active);
        vec_reset(&distance);
        context=cxt=dst=0;
    }
    bool Find(int b){
        bool found=false;
        for (int i=0;i<elementCount;i=i+2) if (element[i]==b) {
            found=true;
            break;
        }
        return found;
    }
    bool FindEnd(int b,int c) {
        bool found=false;
        for (int i=0;i<elementCount;i=i+2) if (element[i]==b&&element[i+1]==c) found=true;
        return found;
    }
    int last() {
        return vec_prev(&active);
    }
    void __attribute__ ((noinline)) Update(int byte) {
        bool pop=false;
        if (!vec_empty(&active)) {
            if (FindEnd(vec_at(&active,vec_size(&active)-1) , byte) || vec_at(&distance,vec_size(&distance)-1) >= limit) {
                vec_pop(&active);
                vec_pop(&distance);
                pop=doPop;
            } else {
                vec_i(&distance,vec_size(&distance)-1);
            }
        }
        if (pop==false && Find(byte)) {
            vec_push( &active,byte);
            vec_push( &distance,0);
        }
        if (!vec_empty(&active)) {
            cxt=vec_at(&active,vec_size(&active)-1);
            dst=min(vec_at(&distance,vec_size(&distance)-1),(1<< (sizeof(T)*8))-1);
            context = (1<< ((sizeof(T)*8))) * cxt+dst;
        } else {
            context=cxt=dst=0;
        }
    }
};

bool isPre=false;
// Table/row & column context
struct Column {
    U32 linepos;
    U8 fc;
    vec<U8,1024*2> bytes; // max lenght 1024*2 chars
};
#define WIKIHEADER GREATERTHAN
#define WIKITABLE  '-'
struct ColumnContext {
    Column col[4];           // Content of last 3 + current row
    vec<U32,16*2> cell[4];   // Content of table row cell positions, max 4 rows, max 16*2 positions
    int rows;
    int cellCount,cells,abovecellpos,abovecellpos1;
    bool NL,isTemp;
    int limit;  // column lenght limit
    U8 nlChar;
    void Init( int l=31) {
        rows=abovecellpos=cellCount=abovecellpos1=cells=0;
        nlChar=LF;
        limit=l;
        for (int i=0;i<4;i++) vec_new(&col[i].bytes);
        for (int i=0;i<4;i++) vec_new(&cell[i]);
        NL=isTemp=false;
        resetCells();
        for (int i=0;i<4;i++) vec_reset(&col[i].bytes);
    }
   
    const U8 lastfc(int i=0) {
        return col[(rows-i)&3].fc;
    }
    void setlastfc(U8 f,int i=0) {
        col[(rows-i)&3].fc=f;
    }
    bool isNewLine() {
        return NL;
    }
    int  __attribute__ ((noinline)) collen(int i=0,int l=0) {
        return min((l?l:limit), vec_size(&col[(rows-i)&3].bytes)+1);
    }
    int nlpos(int i=0) {
        return col[(rows-i)&3].linepos;
    }
    U8  __attribute__ ((noinline)) colb(int i=1,int j=0,int l=0) {
        if (collen(0,l)<collen(i,l))
        return  vec_at(&col[(rows-i)&3].bytes,collen()-(1+j));
        else return 0;
    }
    void __attribute__ ((noinline)) Update(int byte,int b2=0) {
        // Start and end of table - this expects char { is swaped to {{
        if (b2==((CURLYOPENING<<16)+(CURLYOPENING<<8)+VERTICALBAR) ) nlChar=WIKITABLE;
        else if (b2==((VERTICALBAR<<16)+(CURLYCLOSE<<8)+CURLYCLOSE) ) nlChar=LF,resetCells();
        if (byte!=CURLYOPENING && (b2&0xff00)==(CURLYOPENING<<8) && (b2&0xff0000)!= (CURLYOPENING<<16)) 
        isTemp=true;
        else if ( isTemp==true && byte==CURLYCLOSE) isTemp=false;
        // Column
        NL=false;
        if (byte==LF) {
            vec_push( &col[rows].bytes,U8(byte));
            rows++;
            rows=rows&3;
            vec_reset(&col[rows].bytes); // reset new line.
            col[rows].fc=0;
            col[rows].linepos=x.blpos-1;
        }  else {
            vec_push(&col[rows].bytes,U8(byte)); // set new byte to line
            if (collen()==2) {
                col[rows].fc=min(byte,TEXTDATA);
                NL=true;
                if (col[rows].fc==GREATERTHAN && isPre==false) nlChar=WIKIHEADER;
                if (col[rows].fc==SQUAREOPEN && nlChar==WIKIHEADER) nlChar=LF;
            }
        }
        /*
        {|  Table start	It opens a table (and is required)
        |+  Table caption	It adds a caption
        |-  Table row	It adds a new row (but it is optional for the first row)
        !   Header cell	It adds a header cell, whose content can optionally be placed on a new line
        !!  Header cell (on the same line)	It adds a header cell on the same line
        |   Data cell	It adds a data cell, whose content can optionally be placed on a new line (see also the attribute separator)
        ||  Data cell (on the same line)	It adds a data cell on the same line
        |   Attribute separator	It separates a HTML attribute from cell or caption contents
        |}  Table end	It closes a table (and is required)
        */
        // Only  {| |- | || |} are implemented
        if (nlChar==WIKITABLE) {
            if ((b2&0xffff)==(WIKITABLE+VERTICALBAR*256)) {
                cells++;
                cells=cells&3;
                vec_reset(&cell[cells]); // reset new row.
                vec_push( &cell[cells],U32(x.blpos));
                cellCount=abovecellpos=abovecellpos1=0;
            }
            bool newcell=false;
            // Cells
            if ( (b2&0xffff)==(VERTICALBAR+VERTICALBAR*256) ||                // || 
             (b2&0xffff00)==((VERTICALBAR+LF*256)*256) ||                     // \n|x
            ((b2&0xffff00)==((VERTICALBAR+LF*256)*256) && byte!=VERTICALBAR)  // \n|yx  where y!=|
            ) vec_push( &cell[cells],U32(x.blpos)),cellCount++,newcell=true;
            // Advence above cell pos
            if (abovecellpos) {
                abovecellpos++;
                // When above cell is shorter reset
                if (abovecellpos>abovecellpos1) abovecellpos=abovecellpos1=0;
            }
            // If more then one cell get above cell based on current row cell
            if(newcell==true && cellsCount()>0) {
                // Get current above cell pos
                abovecellpos=cellPos(cellCount-1);
                abovecellpos1=cellPos(cellCount);
            }
        }
        // We have table of wikipeda article header containing time, etc ...
        if (nlChar==WIKIHEADER) {
            if ((b2&0xffff)==(WIKIHEADER+LF*256)) {
                //printf("%d %d\n",vec_size(&cell[cells]),vec_at(&cell[cells],vec_size(&cell[cells])-1));
                cells++;
                cells=cells&3;
                vec_reset(&cell[cells]); // reset new row.
                vec_push( &cell[cells],U32(x.blpos));
                cellCount=abovecellpos=abovecellpos1=0;
                //printf("\n");
            }else{
            
            bool newcell=false;
            // Cells
            if ((b2&0xff)==WIKIHEADER) vec_push(&cell[cells],U32(x.blpos)),cellCount++,newcell=true;
            // Advence above cell pos
            if (abovecellpos) {
                abovecellpos++;
                // When above cell is shorter reset
                if (abovecellpos>abovecellpos1) abovecellpos=abovecellpos1=0;
            }
            // If more then one cell get above cell based on current row cell
            if(newcell==true && cellsCount()>0) {
                //printf("%d  ",cellPos(cellCount-1));
                // Get current above cell pos
                abovecellpos=cellPos(cellCount-1);
                abovecellpos1=cellPos(cellCount);
            }
            }
        }
    }
    int cellsCount(int row=1) {
        return vec_size(&cell[(cells-row)&3]);
    }
    int cellPos(int cellID,int row=1) {
        int total=cellsCount(row)-1;
        total=min(total,cellID);
        return vec_at(&cell[(cells-row)&3],total);
    }
    void resetCells() {
        for (int i=0;i<4;i++) vec_reset(&cell[i]);
    }
};
// Keep track of main brackets
const U8 brackets[8]={'(',')', CURLYOPENING,CURLYCLOSE, '[',']', LESSTHAN,GREATERTHAN};
// Keep track of ' and " as quotes
const U8 quotes[4]={APOSTROPHE,APOSTROPHE,QUOTATION,QUOTATION};
// Keep track of first char including some brackets
const U8 fchar[20]={FIRSTUPPER,LF, TEXTDATA,LF, COLON,LF, LESSTHAN,GREATERTHAN,EQUALS,LF,SQUAREOPEN,SQUARECLOSE,CURLYOPENING,CURLYCLOSE,'*',LF,VERTICALBAR,LF,HTLINK,LF};
const U16 html[2]={'&'*256+'L','&'*256+'N'};

// Sentence & words context
struct WordsContext {
    vec<U16,64*4> sbytes; // List of bytes surrounded by a current word, max 64*4
    vec<U32,64*4> type;   // List of bytes surrounded by a current word, max 64*4
    vec<U32,64*4> stem;   // List of bytes surrounded by a current word, max 64*4
    vec<U8,64*4> capital;   // List of bytes surrounded by a current word, max 64*4
    vec<U32,64*4> codeword;   // List of bytes surrounded by a current word, max 64*4
    U32 fword,ftype;      // First word of a sentence
    U8 pbyte,tpbyte;             // Current byte before word
    int wordcount,upper;
    U32 codesum;
    int ref;
    bool paragraph;
    int maxc=0;
    int worInPar=0,worInLink=0;
    void Init() {
        vec_new(&sbytes);
        vec_new(&type);
        vec_new(&stem);
        Reset();
    }
    void Reset() {
        vec_reset(&sbytes);
        vec_reset(&type);
        vec_reset(&stem);
        vec_reset(&capital);
        vec_reset(&codeword);
        fword=pbyte=wordcount=upper=ftype=ref=codesum=tpbyte=0;
        paragraph=false;
        worInPar=0,worInLink=0;
    }
    void Set(U8 b,int a=0,U8 g=0) {
        pbyte=b;upper=a;tpbyte=g;
    }
    void  __attribute__ ((noinline)) Update(U32 w,U8 b, U32 t,U32 s,U32 cw=0,int worIn=0) {
        if (fword==0) fword=w;
        vec_push(&sbytes,U16(pbyte*256+b));  // Surrounding bytes
        vec_push(&type,t);
        vec_push(&stem,s);
        vec_push(&capital,U8(upper));
        vec_push(&codeword,cw);
        pbyte=tpbyte=0;wordcount++;
        codesum=codesum+cw;
        if (ftype==0 && t) ftype=t;
        if (worIn==1) worInLink++; 
        else  if (worIn==0) worInPar++;
    }
    void  __attribute__ ((noinline)) Remove() {
        const int num=vec_size(&stem);
        if (num) {
            vec_pop(&sbytes),vec_pop(&type),vec_pop(&stem),vec_pop(&capital),vec_pop(&codeword);
            if (wordcount)wordcount--;
        }
    }
    U32  __attribute__ ((noinline)) Word(int i=1) {
        const int num=vec_size(&stem);
        if (num>=i) return vec_at(&stem,num-(i));
        else return 0;
    }
    // Reverse order word
    U32  __attribute__ ((noinline)) WordR(int i=1) {
        int low=min(wordcount,i);
        return Word(wordcount-low);
    }
    U16  __attribute__ ((noinline)) sBytes(int i=1) {
        const int num=vec_size(&sbytes);
        if (num>=i) return vec_at(&sbytes,num-(i));
        else return 0;
    }
    U32  __attribute__ ((noinline)) Type(int i=1) {
        const int num=vec_size(&type);
        if (num>=i) return vec_at(&type,num-(i));
        else return 0;
    }
    U8  __attribute__ ((noinline)) Capital(int i=1) {
        const int num=vec_size(&capital);
        if (num>=i) return vec_at(&capital,num-(i));
        else return 0;
    }
    U32  __attribute__ ((noinline)) Code(int i=1) {
        const int num=vec_size(&codeword);
        if (num>=i) return vec_at(&codeword,num-(i));
        else return 0;
    }
    U32  __attribute__ ((noinline)) CodeR(int i,int j=1) {
        int low=min(wordcount,j);
        return Code(i+(wordcount-low));// i+, is it ok?
    }
    // Reverse order word of type t if any.
    U32  __attribute__ ((noinline)) LastR(int i,int j=1, U32 t=0) {
        int low=min(wordcount,j);
        return Last(i+(wordcount-low),t);// i+, is it ok?
    }
    // Return last word matching type, ... If not found return 0
    U32  __attribute__ ((noinline)) Last(int j=1, U32 t=0) {
        const int num=vec_size(&type);
        if (t==0) return Word(j);
        if (num>=j) {
            U32 typ=0;
            for (int i=j; i<num; i++) {
                typ=Type(i);
                if (typ&t) return Word(i);
            }
        }
        return Word(j);
    }
    U32  __attribute__ ((noinline)) LastIf(int j=1, U32 t=0) {
        const int num=vec_size(&type);
        if (t==0) return Word(j);
        if (num>=j) {
            U32 typ=0;
            for (int i=j; i<num; i++) {
                typ=Type(i);
                if (typ&t) return Word(i);
            }
        }
        return 0;
    }
    U32  __attribute__ ((noinline)) LastIdx(int j=1, U32 t=0) {
        const int num=vec_size(&type);
        if (t==0) return 0;
        if (num>=j) {
            
            U32 typ=0;
            for (int i=j; i<num; i++) {
                typ=Type(i);
                if (typ&t) return i;
            }
        }
        return 0;
    }
    void __attribute__ ((noinline)) removeWordsL(int len, U8 c,U8 d, const bool f=true) {
        if ((sBytes(1)&0xff)==d) {
            for (int i=1; i<len; i++) {
                if ((sBytes(i)>>8)==c) {
                    while( (sBytes(1)>>8)!=c) Remove();
                    if (f) Remove();
                    break;
                }
            }
        }
    }
    void __attribute__ ((noinline)) removeWordsR(int len, U8 c,U8 d, const bool f=true) {
        if ((sBytes(1)&0xff)==d) {
            for (int i=1; i<len; i++) {
                if ((sBytes(i)&0xff)==c) {
                    while( (sBytes(1)&0xff)!=c) Remove();
                    if (f) Remove();
                    break;
                }
            }
        }
    }
    U32  __attribute__ ((noinline)) Word0(int i=1) {
        const int num=vec_size(&stem);
        U8 lb=sBytes(i)&0xff;
        int idx=0;
        if (i==4) idx=37*47*53*83;
        else if (i==3) idx=47*53*83;
        else if (i==2) idx=53*83;
        else if (i==1) idx=83;
        if (lb==LF){
            if (i==3) idx=idx*37;
            else if (i==3) idx=idx*47;
            else if (i==2) idx=idx*53;
            else if (i==1) idx=idx*83;
        }
        if (num>=i) return vec_at(&stem,num-(i))*idx;
        else return 0;
    }
};

// Sentence Context hold upto 64 sentences.
// At init every sentence is empty. We can update, 
// reset all and request sentences.
// Similar sentences are calculated at runtime comparing
// words dictionary index sum differences. Asuming we have
// similar sentences with small differences, currently limited to 53.
#define SIMILARWORDS 64
struct SentenceContext {
    WordsContext sentence[SIMILARWORDS]; // List of sentences, max 64
    WordsContext empty;        // blank
    int sindex;
    U32 total;
    void Init() {
        for (int i=0; i<SIMILARWORDS; i++) sentence[i].Init();
        empty.Init();
        sindex=total=0;
    }
    void Reset() {
        for (int i=1; i<SIMILARWORDS; i++) sentence[i].Reset();
        sindex=total=0;
    }
    void  __attribute__ ((noinline)) Update(WordsContext *w) {
        if (w->wordcount) {
            memcpy(&sentence[sindex],w,sizeof(WordsContext));
            sindex=(sindex+1)&(SIMILARWORDS-1);
            total++;
        }
    }
    WordsContext *Sentence(int i) {
        return &sentence[(sindex-i)&(SIMILARWORDS-1)];
    }
    WordsContext __attribute__ ((noinline)) *SimilarSentence(WordsContext *wor,int wcount,U32 pres=53) {
        U32 isSimilar=0,isSimilarIdx=0;
        U32 codesum=0; // input
        if (wcount>=1) {
            for (int i=0; i<SIMILARWORDS; i++) {
                if (sentence[i].wordcount && sentence[i].wordcount<= wcount && sentence[i].wordcount>(wcount/2)) {
                    for (int k=0; k<wcount; k++) {
                        U32 curcode=wor->Code(wcount-k); //get word code of current sentance
                        // loop over and compare codewords
                        if (curcode) {
                            for (int j=0; j<sentence[i].wordcount; j++) {
                                U32 testcode=sentence[i].Code(sentence[i].wordcount-j);
                                if (testcode) {
                                    U32 diff=(curcode-testcode);
                                    if (diff==0) {
                                        codesum++;
                                        break;
                                    }
                                }
                            } 
                        }
                    }
                    if (codesum> isSimilar) isSimilar=codesum,isSimilarIdx=i+1;
                    codesum=0;
                }
            }
            if (isSimilar) {
                isSimilar=isSimilar*100/wcount;
                if (isSimilar<pres) isSimilar=isSimilarIdx=0;// if below 53% similar then skip
            }
        }
        if (isSimilarIdx) return &sentence[isSimilarIdx-1];
        else return &empty;
    }
};

inline U32 hash(U32 a, U32 b, U32 c=0xffffffff) {
    U32 h=a*110002499u+b*30005491u+c*50004239u; 
    return h^h>>9^a>>3^b>>3^c>>4;
}

inline int charSwap(int c){
    if (c>='{' && c<127) c+='P'-'{';
    else if (c>='P' && c<'T') c-='P'-'{';
    else if ( (c>=':' && c<='?') || (c>='J' && c<='O') ) c^=0x70;
    if (c=='X' || c=='`') c^='X'^'`';
    return c;
}

//                Stemming routines
// English affix stemmer, based on the Porter2 stemmer.
// This is mostly from paq8px with some modifications.

#define MAX_WORD_SIZE 64
class Word {
public:
  U8 Letters[MAX_WORD_SIZE];
  U8 Start, End,fEnd;
  U32 Hash, Type, Suffix, Preffix, NextW;
  Word(): Start(0), End(0), Hash(0), Type(0), Suffix(0), Preffix(0), NextW(0) {
    memset(&Letters[0], 0, sizeof(U8)*MAX_WORD_SIZE);
  }
  bool operator==(const char *s) const{
    size_t len=strlen(s);
    return ((size_t)(End-Start+(Letters[Start]!=0))==len && memcmp(&Letters[Start], s, len)==0);
  }
  bool operator!=(const char *s) const{
    return !operator==(s);
  }
  void operator+=(const char c){
    if (c>0 && End<MAX_WORD_SIZE-1){
      End+=(Letters[End]>0);
      Letters[End]=c;
    }
  }
  U8 operator[](U8 i) const{
    return (End-Start>=i)?Letters[Start+i]:0;
  }
  U8 operator()(U8 i) const{
    return (End-Start>=i)?Letters[End-i]:0;
  }
    U32 Length() const{
    if (Letters[Start]!=0)
      return End-Start+1;
    return 0;
  }

  bool ChangeSuffix(const char *OldSuffix, const char *NewSuffix){
    size_t len=strlen(OldSuffix);
    if (Length()>len && memcmp(&Letters[End-len+1], OldSuffix, len)==0){
      size_t n=strlen(NewSuffix);
      if (n>0){
        memcpy(&Letters[End-len+1], NewSuffix, min(MAX_WORD_SIZE-1,End+n)-End);
        End=min(MAX_WORD_SIZE-1, End-len+n);
      }
      else
        End-=len;
      return true;
    }
    return false;
  }
  bool MatchesAny(const char* a[], const int count) {
    int i=0;
    size_t len = (size_t)Length();
    for (; i<count && (len!=strlen(a[i]) || memcmp(&Letters[Start], a[i], len)!=0); i++);
    return i<count;
  }
  bool MatchesAnyP(const char* a[], const int count) {
    int i=0;
    NextW=0xff;
    size_t len = (size_t)Length();
    for (; i<count && (len!=strlen(a[i]) || memcmp(&Letters[Start], a[i], len)!=0); i++);
    if (i<count) NextW=i;
    return i<count;
  }
  bool EndsWith(const char *Suffix) const{
    size_t len=strlen(Suffix);
    return (Length()>len && memcmp(&Letters[End-len+1], Suffix, len)==0);
  }
  bool StartsWith(const char *Prefix) const{
    size_t len=strlen(Prefix);
    return (Length()>len && memcmp(&Letters[Start], Prefix, len)==0);
  }
};

enum EngWordTypeFlags {
  Verb                   = (1<<0),
  Noun                   = (1<<1),
  Adjective              = (1<<2),
  Plural                 = (1<<3),
  PastTense              = (1<<5)|Verb,
  PresentParticiple      = (1<<4)|Verb,
  AdjectiveSuperlative   = (1<<5)|Adjective,
  AdjectiveWithout       = (1<<6)|Adjective,
  AdjectiveFull          = (1<<7)|Adjective,
  AdverbOfManner         = (1<<8),
  Suffix                 = (1<<9),
  Prefix                 = (1<<10),
  Male                   = (1<<11),
  Female                 = (1<<13),
  Article                = (1<<14),
  Conjunction            = (1<<15),
  Adposition             = (1<<16),
  Number                 = (1<<17), // not used
  Preposition            = (1<<18), // not used
  ConjunctiveAdverb      = (1<<19),
  Pronoun                = (1<<20),
};
enum EngWordTypeFlagsNegation {
  Negation               = (1<<0),
  PrefixIrr              = (1<<1)|Negation,
  PrefixOver             = (1<<2),
  PrefixUnder            = (1<<3),
  PrefixUnn              = (1<<4)|Negation,
  PrefixNon              = (1<<5)|Negation,
  PrefixAnti             = (1<<6)|Negation,
  PrefixDis              = (1<<7)|Negation
};

enum EngWordTypeFlagsSuffix {
  SuffixNESS             = (1<<0),
  SuffixITY              = (1<<1)|Noun,
  SuffixCapable          = (1<<2),
  SuffixNCE              = (1<<3),
  SuffixNT               = (1<<4),
  SuffixION              = (1<<5),
  SuffixAL               = (1<<6)|Adjective,
  SuffixIC               = (1<<7)|Adjective,
  SuffixIVE              = (1<<8),
  SuffixOUS              = (1<<9)|Adjective,
};

// pronouns/determiner 
// https://en.wikipedia.org/wiki/English_pronouns
const char *Pronouns[14][5]={
{" "/*i*/,  "me",  "myself",     "mine",   "my"},
{"we",      "us",  "ourselves",  "ours",   "our"},
{"you",     "you", "yourself",   "yours",  "your"},
{"thou",    "thee","thyself",    "thine",  "thy"},
{"you",     "you", "yourselves", "yours",  "your"},
{"he",      "him", "himself",    "his",    "his"},
{"she",     "her", "herself",    "hers",   "her"},
{"it",      "it",  "itself",     " ",      "its"},
{"they",    "them","themself",   "theirs", "their"},
{"they",    "them","themselves", "theirs", "their"},
{"one",     "one", "oneself",    "one's",  "one's"},
{"who",     "whom", " ",         "whose",  "whose"},
{"what",    "what", " ",         " ",       " "},
{"which",   "which"," ",         " ",       " "},
};

#define NUM_VOWELS 6
const char Vowels[NUM_VOWELS]={'a','e','i','o','u','y'};
#define NUM_DOUBLES 9
const char Doubles[NUM_DOUBLES]={'b','d','f','g','m','n','p','r','t'};
#define NUM_LI_ENDINGS 10
const char LiEndings[NUM_LI_ENDINGS]={'c','d','e','g','h','k','m','n','r','t'};
#define NUM_NON_SHORT_CONSONANTS 3
const char NonShortConsonants[NUM_NON_SHORT_CONSONANTS]={'w','x','Y'}; 
#define NUM_VERB 23
const char *VerbWords1[NUM_VERB]={"has","had","have","was","were","may","might","must","shall","should","can","could","will","would","is","am","are","be","being","been","do","does","did"};
#define NUM_NUM 21
const char *Numbers[NUM_NUM]={"one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten","twenty","thirty","forty","fifty","sixty","seventy","eighty","ninety","hundred","thousand","million"};
#define  NUM_CONJ_WORDS  25
const char *ConjWords[NUM_CONJ_WORDS]={"for","and","nor","but","or","yet","so","than","as","that","if","when","because","while",
"where","after","though","whether","before","although","like","once","unless","now","except"};
#define  NUM_APO_WORDS  30  //Adposition
const char *ApoWords[NUM_APO_WORDS]={"in","during","at","on","since","until","above", "across", "against", "along", "among", "around",
   "behind", "below", "beneath", "beside", "between", "by", "down", "from",  "into", "near", "of", "off", "to", "toward", "under", "upon", "with", "within"};
#define  NUM_PREP_WORDS 5 //preposition
const char *PrepWords[NUM_PREP_WORDS]={"as","by","de","in","on"};
#define  NUM_CAVER_WORDS 2 //Conjunctive Adverb
const char *ConAdVerPrepWords[NUM_CAVER_WORDS]={"also","thus"};
#define  NUM_VERB_WORDS 12 
const char *VerbWords[NUM_VERB_WORDS]={"be","do","an","could","may","must","need","ought","shall","should","will","would"};
#define  NUM_MALE_WORDS  9
const char *MaleWords[NUM_MALE_WORDS]={"he","him","his","himself","man","men","boy","husband","actor"};
#define  NUM_FEMALE_WORDS  8
const char *FemaleWords[NUM_FEMALE_WORDS]={"she","her","herself","woman","women","girl","wife","actress"};
#define  NUM_ARTICLE_WORDS  3
const char *ArticleWords[NUM_ARTICLE_WORDS]={"a","an","the"};
#define NUM_SUFFIXES_STEP0 3
const char *SuffixesStep0[NUM_SUFFIXES_STEP0]={"'s'","'s","'"};
#define NUM_SUFFIXES_STEP1b 6
const char *SuffixesStep1b[NUM_SUFFIXES_STEP1b]={"eedly","eed","ed","edly","ing","ingly"};
const U32 TypesStep1b[NUM_SUFFIXES_STEP1b]={AdverbOfManner,0,PastTense,AdverbOfManner|PastTense,PresentParticiple,AdverbOfManner|PresentParticiple};
#define NUM_SUFFIXES_STEP2 22
const char *(SuffixesStep2[NUM_SUFFIXES_STEP2])[2]={
  {"ization", "ize"},
  {"ational", "ate"},
  {"ousness", "ous"},
  {"iveness", "ive"},
  {"fulness", "ful"},
  {"tional", "tion"},
  {"lessli", "less"},
  {"biliti", "ble"},
  {"entli", "ent"},
  {"ation", "ate"},
  {"alism", "al"},
  {"aliti", "al"},
  {"fulli", "ful"},
  {"ousli", "ous"},
  {"iviti", "ive"},
  {"enci", "ence"},
  {"anci", "ance"},
  {"abli", "able"},
  {"izer", "ize"},
  {"ator", "ate"},
  {"alli", "al"},
  {"bli", "ble"},
 
};
const U32 TypesStep2[NUM_SUFFIXES_STEP2]={
  Suffix,
  Suffix|Adjective,
  Suffix,
  Suffix,
  Suffix,
  Suffix|Adjective,
  AdverbOfManner,
  AdverbOfManner|Noun|Suffix,
  AdverbOfManner,
  Suffix,
  0,
  Noun|Suffix,
  AdverbOfManner,
  AdverbOfManner,
  Noun|Suffix,
  0,
  0,
  AdverbOfManner,
  0,
  0,
  AdverbOfManner,
  AdverbOfManner
};

const U32 TypesStep2Suffix[NUM_SUFFIXES_STEP2]={
  SuffixION,
  SuffixION|SuffixAL,
  SuffixNESS,
  SuffixNESS,
  SuffixNESS,
  SuffixION|SuffixAL,
  0,
  SuffixITY,
  0,
  SuffixION,
  0,
  SuffixITY,
  0,
  0,
  SuffixITY,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
};
#define NUM_SUFFIXES_STEP3 8
const char *(SuffixesStep3[NUM_SUFFIXES_STEP3])[2]={
  {"ational", "ate"},
  {"tional", "tion"},
  {"alize", "al"},
  {"icate", "ic"},
  {"iciti", "ic"},
  {"ical", "ic"},
  {"ful", ""},
  {"ness", ""},
  
};
const U32 TypesStep3[NUM_SUFFIXES_STEP3]={
Suffix|Adjective,
Suffix|Adjective,
0,
0,
Noun|Suffix,
Suffix|Adjective,
AdjectiveFull,
Suffix
};

const U32 TypesStep3Suffix[NUM_SUFFIXES_STEP3]={
SuffixION|SuffixAL,
SuffixION|SuffixAL,
0,
0,
SuffixITY,
SuffixAL,
0,
SuffixNESS
};
#define NUM_SUFFIXES_STEP4 20
const char *SuffixesStep4[NUM_SUFFIXES_STEP4]={"al","ance","ence","er","ic","able","ible","ant","ement","ment","ent","ou","ism","ate","iti","ous","ive","ize","sion","tion"};
const U32 TypesStep4[NUM_SUFFIXES_STEP4]={
  Suffix|Adjective,
  Suffix,
  Suffix,
  0,
  Suffix|Adjective,
  Suffix,
  Suffix,
  Suffix,
  0,
  0,
  Suffix,
  0,
  0,
  0,
  Suffix|Noun,
  Suffix|Adjective,
  Suffix,
  0,
  Suffix,
  Suffix
};

const U32 TypesStep4Suffix[NUM_SUFFIXES_STEP4]={
  SuffixAL,
  SuffixNCE,
  SuffixNCE,
  0,
  SuffixIC,
  SuffixCapable,
  SuffixCapable,
  SuffixNT,
  0,
  0,
  SuffixNT,
  0,
  0,
  0,
  SuffixITY,
  SuffixOUS,
  SuffixIVE,
  0,
  SuffixION,
  SuffixION
};
#define NUM_EXCEPTION_REGION1 3
const char *ExceptionsRegion1[NUM_EXCEPTION_REGION1]={"gener","arsen","commun"};
#define NUM_EXCEPTIONS1 19
const char *(Exceptions1[NUM_EXCEPTIONS1])[2]={
  {"skis", "ski"},
  {"skies", "sky"},
  {"dying", "die"},
  {"lying", "lie"},
  {"tying", "tie"},
  {"idly", "idle"},
  {"gently", "gentle"},
  {"ugly", "ugli"},
  {"early", "earli"},
  {"only", "onli"},
  {"singly", "singl"},
  {"sky", "sky"},
  {"news", "news"},
  {"howe", "howe"},
  {"atlas", "atlas"},
  {"cosmos", "cosmos"},
  {"bias", "bias"},
  {"andes", "andes"},
  {"texas", "texas"}
};
const U32 TypesExceptions1[NUM_EXCEPTIONS1]={
Noun|Plural,
Noun|Plural,
PresentParticiple,
PresentParticiple,
PresentParticiple,
AdverbOfManner,
AdverbOfManner,
Adjective,
Adjective|AdverbOfManner,
0,
AdverbOfManner,
Noun,
Noun,
0,
Noun,
Noun,
Noun,
Noun|Plural,
Noun};
#define NUM_EXCEPTIONS2 8
const char *Exceptions2[NUM_EXCEPTIONS2]={"inning","outing","canning","herring","earring","proceed","exceed","succeed"};
const U32 TypesExceptions2[NUM_EXCEPTIONS2]={Noun,Noun,Noun,Noun,Noun,Verb,Verb,Verb}; 

inline bool CharInArray(const char c, const char a[], const int len) {
  if (a==NULL)
    return false;
  int i=0;
  for (;i<len && c!=a[i];i++);
  return i<len;
}

class EnglishStemmer {
private:
  inline bool IsVowel(const char c) {
    return CharInArray(c, Vowels, NUM_VOWELS);
  }

  inline bool IsConsonant(const char c) {
    return !IsVowel(c);
  }

  inline bool IsShortConsonant(const char c) {
    return !CharInArray(c, NonShortConsonants, NUM_NON_SHORT_CONSONANTS);
  }

  inline bool IsDouble(const char c) {
    return CharInArray(c, Doubles, NUM_DOUBLES);
  }

  inline bool IsLiEnding(const char c) {
    return CharInArray(c, LiEndings, NUM_LI_ENDINGS);
  }
  inline void Hash(Word *W){
    (*W).Hash=0xb0a710ad;
    for (int i=(*W).Start;i<=(*W).End;i++)
      (*W).Hash=(*W).Hash*263*32+(*W).Letters[i];
  }
   
  U32 GetRegion(const Word *W, const U32 From) {
    bool hasVowel = false;
    for (int i=(*W).Start+From;i<=(*W).End;i++) {
      if (IsVowel((*W).Letters[i])) {
        hasVowel = true;
        continue;
      }
      else if (hasVowel)
        return i-(*W).Start+1;
    }
    return (*W).Length();
  }
  U32 GetRegion1(const Word *W){
    for (int i=0;i<NUM_EXCEPTION_REGION1;i++) {
      if ((*W).StartsWith(ExceptionsRegion1[i]))
        return strlen(ExceptionsRegion1[i]);
    }
    return GetRegion(W, 0);
  }
  bool SuffixInRn(const Word *W, const U32 Rn, const char *Suffix) {
    return ((*W).Start!=(*W).End && Rn<=(*W).Length()-strlen(Suffix));
  }
  bool EndsInShortSyllable(const Word *W) {
    if ((*W).End==(*W).Start)
      return false;
    else if ((*W).End==(*W).Start+1)
      return IsVowel((*W)(1)) && IsConsonant((*W)(0));
    else
      return (IsConsonant((*W)(2)) && IsVowel((*W)(1)) && IsConsonant((*W)(0)) && IsShortConsonant((*W)(0)));
  }
  bool IsShortWord(const Word *W){
    return (EndsInShortSyllable(W) && GetRegion1(W)==(*W).Length());
  }
  inline bool HasVowels(const Word *W) {
    for (int i=(*W).Start;i<=(*W).End;i++) {
      if (IsVowel((*W).Letters[i]))
        return true;
    }
    return false;
  }
  bool TrimStartingApostrophe(Word *W) {
    bool result=false;
    //trim all apostrophes from the beginning
    int cnt=0;
    while((*W).Start!=(*W).End && (*W)[0]==APOSTROPHE) {
      result=true;
      (*W).Start++;
      cnt++;
    }
    //trim the same number of apostrophes from the end (if there are)
    while((*W).Start!=(*W).End && (*W)(0)==APOSTROPHE) {
      if(cnt==0)break;
      (*W).End--;
      cnt--;
    }

    if ((*W)(0)=='-') {
      (*W).End--;
    }
    return result;
  }
  void MarkYsAsConsonants(Word *W) {
    if ((*W)[0]=='y')
      (*W).Letters[(*W).Start]='Y';
    for (int i=(*W).Start+1;i<=(*W).End;i++) {
      if (IsVowel((*W).Letters[i-1]) && (*W).Letters[i]=='y')
        (*W).Letters[i]='Y';
    }
  }
  bool ProcessPrefixes(Word *W) {
    if ((*W).StartsWith("irr") && (*W).Length()>5 && ((*W)[3]=='a' || (*W)[3]=='e'))
      (*W).Start+=2, (*W).Type|=Prefix, (*W).Preffix|=PrefixIrr;
    else if ((*W).StartsWith("over") && (*W).Length()>5)
      (*W).Start+=4, (*W).Type|=Prefix, (*W).Preffix|=PrefixOver;
    else if ((*W).StartsWith("under") && (*W).Length()>6)
      (*W).Start+=5, (*W).Type|=Prefix, (*W).Preffix|=PrefixUnder;
    else if ((*W).StartsWith("unn") && (*W).Length()>5)
      (*W).Start+=2, (*W).Type|=Prefix, (*W).Preffix|=PrefixUnn;
    else if ((*W).StartsWith("non") && (*W).Length()>(U32)(5+((*W)[3]=='-')))
      (*W).Start+=2+((*W)[3]=='-'), (*W).Type|=Prefix, (*W).Preffix|=PrefixNon;
    else if ((*W).StartsWith("anti") && (*W).Length()>6&& ((*W)[4]=='-'))
      (*W).Start+=4+((*W)[4]=='-'), (*W).Type|=Prefix, (*W).Preffix|=PrefixAnti;
    else if ((*W).StartsWith("dis") && (*W).Length()>5 && ((*W)[3]=='-'))
      (*W).Start+=2+((*W)[3]=='-'), (*W).Type|=Prefix, (*W).Preffix|=PrefixDis;
    else
      return false;
    return true;
  }
  bool ProcessSuperlatives(Word *W) {
    if ((*W).EndsWith("est") && (*W).Length()>4) {
      U8 i=(*W).End;
      (*W).End-=3;
      (*W).Type|=AdjectiveSuperlative;

      if ((*W)(0)==(*W)(1) && (*W)(0)!='r' && !((*W).Length()>=4 && memcmp("sugg",&(*W).Letters[(*W).End-3],4)==0)) {
        (*W).End-= ( ((*W)(0)!='f' && (*W)(0)!='l' && (*W)(0)!='s') ||
                   ((*W).Length()>4 && (*W)(1)=='l' && ((*W)(2)=='u' || (*W)(3)=='u' || (*W)(3)=='v'))) &&
                   (!((*W).Length()==3 && (*W)(1)=='d' && (*W)(2)=='o'));
        if ((*W).Length()==2 && ((*W)[0]!='i' || (*W)[1]!='n'))
          (*W).End = i, (*W).Type&=~AdjectiveSuperlative;
      }
      else{
        switch((*W)(0)) {
          case 'd': case 'k': case 'm': case 'y': break;
          case 'g': {
            if (!( (*W).Length()>3 && ((*W)(1)=='n' || (*W)(1)=='r') && memcmp("cong",&(*W).Letters[(*W).End-3],4)!=0 ))
              (*W).End = i, (*W).Type&=~AdjectiveSuperlative;
            else
              (*W).End+=((*W)(2)=='a');
            break;
          }
          case 'i': {(*W).Letters[(*W).End]='y'; break;}
          case 'l': {
            if ((*W).End==(*W).Start+1 || memcmp("mo",&(*W).Letters[(*W).End-2],2)==0)
              (*W).End = i, (*W).Type&=~AdjectiveSuperlative;
            else
              (*W).End+=IsConsonant((*W)(1));
            break;
          }
          case 'n': {
            if ((*W).Length()<3 || IsConsonant((*W)(1)) || IsConsonant((*W)(2)))
              (*W).End = i, (*W).Type&=~AdjectiveSuperlative;
            break;
          }
          case 'r': {
            if ((*W).Length()>3 && IsVowel((*W)(1)) && IsVowel((*W)(2)))
              (*W).End+=((*W)(2)=='u') && ((*W)(1)=='a' || (*W)(1)=='i');
            else
              (*W).End = i, (*W).Type&=~AdjectiveSuperlative;
            break;
          }
          case 's': {(*W).End++; break;}
          case 'w': {
            if (!((*W).Length()>2 && IsVowel((*W)(1))))
              (*W).End = i, (*W).Type&=~AdjectiveSuperlative;
            break;
          }
          case 'h': {
            if (!((*W).Length()>2 && IsConsonant((*W)(1))))
              (*W).End = i, (*W).Type&=~AdjectiveSuperlative;
            break;
          }
          default: {
            (*W).End+=3;
            (*W).Type&=~AdjectiveSuperlative;
          }
        }
      }
    }
    return ((*W).Type&AdjectiveSuperlative)>0;
  }
  bool Step0(Word *W) {
    for (int i=0;i<NUM_SUFFIXES_STEP0;i++) {
      if ((*W).EndsWith(SuffixesStep0[i])) {
        (*W).End-=strlen(SuffixesStep0[i]);
        (*W).Type|=Plural;
        return true;
      }
    }
    return false;
  }
  bool Step1a(Word *W) {
    if ((*W).EndsWith("sses")) {
      (*W).End-=2;
      (*W).Type|=Plural;
      return true;
    }
    if ((*W).EndsWith("ied") || (*W).EndsWith("ies")) {
      (*W).Type|=((*W)(0)=='d')?PastTense:Plural;
      (*W).End-=1+((*W).Length()>4);
      return true;
    }
    if ((*W).EndsWith("us") || (*W).EndsWith("ss"))
      return false;
    if ((*W)(0)=='s' && (*W).Length()>2){
      for (int i=(*W).Start;i<=(*W).End-2;i++) {
        if (IsVowel((*W).Letters[i])) {
          (*W).End--;
          (*W).Type|=Plural;
          return true;
        }
      }
    }
    if ((*W).EndsWith("n't") && (*W).Length()>4) {
      switch ((*W)(3)) {
        case 'a': {
          if ((*W)(4)=='c')
            (*W).End-=2;
          else
            (*W).ChangeSuffix("n't","ll");
          break;
        }
        case 'i': {(*W).ChangeSuffix("in't","m"); break;}
        case 'o': {
          if ((*W)(4)=='w')
            (*W).ChangeSuffix("on't","ill");
          else
            (*W).End-=3;
          break;
        }
        default: (*W).End-=3;
      }
      (*W).Type|=Prefix, (*W).Preffix|=Negation;// suffix as preffix
      return true;
    }
    if ((*W).EndsWith("hood") && (*W).Length()>7) {
      (*W).End-=4;
      return true;
    }
    return false;
  }
  bool Step1b(Word *W, const U32 R1) {
    for (int i=0;i<NUM_SUFFIXES_STEP1b;i++) {
      if ((*W).EndsWith(SuffixesStep1b[i])) {
        switch(i){
          case 0: case 1: {
            if (SuffixInRn(W, R1, SuffixesStep1b[i]))
              (*W).End-=1+i*2;
            break;
          }
          default: {
            U8 j=(*W).End;
            (*W).End-=strlen(SuffixesStep1b[i]);
            if (HasVowels(W)) {
              if ((*W).EndsWith("at") || (*W).EndsWith("bl") || (*W).EndsWith("iz") || IsShortWord(W))
                (*W)+='e';
              else if ((*W).Length()>2){
                if ((*W)(0)==(*W)(1) && IsDouble((*W)(0)))
                  (*W).End--;
                else if (i==2 || i==3){
                  switch((*W)(0)) {
                    case 'c': case 's': case 'v': {(*W).End+=!((*W).EndsWith("ss") || (*W).EndsWith("ias")); break;}
                    case 'd': {
                        static constexpr char nAllowed[4] = {'a', 'e', 'i', 'o'};
                        (*W).End+=IsVowel((*W)(1)) && (!CharInArray((*W)(2),nAllowed, 4)); 
                        break;
                    }
                    case 'k': {(*W).End+=(*W).EndsWith("uak"); break;}
                    case 'l': {
                      static constexpr char allowed1[10] = {'b', 'c', 'd', 'f', 'g', 'k', 'p', 't', 'y', 'z'};
                      static constexpr char allowed2[4] = {'a', 'i', 'o', 'u'};
                      (*W).End+= CharInArray((*W)(1),allowed1, 10) ||
                                (CharInArray((*W)(1),allowed2, 4) && IsConsonant((*W)(2)));
                      break;
                    }
                  }
                }
                else if (i>=4) {
                  switch((*W)(0)) {
                    case 'd': {
                      if (IsVowel((*W)(1)) && (*W)(2)!='a' && (*W)(2)!='e' && (*W)(2)!='o')
                        (*W)+='e';
                      break;
                    }
                    case 'g': {
                      static constexpr char allowed[7] = {'a', 'd', 'e', 'i', 'l', 'r', 'u'};  
                      if ( CharInArray((*W)(1),allowed, 7) || (
                         (*W)(1)=='n' && (
                          (*W)(2)=='e' ||
                          ((*W)(2)=='u' && (*W)(3)!='b' && (*W)(3)!='d') ||
                          ((*W)(2)=='a' && ((*W)(3)=='r' || ((*W)(3)=='h' && (*W)(4)=='c'))) ||
                          ((*W).EndsWith("ring") && ((*W)(4)=='c' || (*W)(4)=='f'))
                         )
                        ) 
                      )
                        (*W)+='e';
                      break;
                    }
                    case 'l': {
                      if (!((*W)(1)=='l' || (*W)(1)=='r' || (*W)(1)=='w' || (IsVowel((*W)(1)) && IsVowel((*W)(2)))))
                        (*W)+='e';
                      if ((*W).EndsWith("uell") && (*W).Length()>4 && (*W)(4)!='q')
                        (*W).End--;
                      break;
                    }
                    case 'r': {
                      if ((
                        ((*W)(1)=='i' && (*W)(2)!='a' && (*W)(2)!='e' && (*W)(2)!='o') ||
                        ((*W)(1)=='a' && (!((*W)(2)=='e' || (*W)(2)=='o' || ((*W)(2)=='l' && (*W)(3)=='l')))) ||
                        ((*W)(1)=='o' && (!((*W)(2)=='o' || ((*W)(2)=='t' && (*W)(3)!='s')))) ||
                        (*W)(1)=='c' || (*W)(1)=='t') && (!(*W).EndsWith("str"))
                      )
                        (*W)+='e';
                      break;
                    }
                    case 't': {
                      if ((*W)(1)=='o' && (*W)(2)!='g' && (*W)(2)!='l' && (*W)(2)!='i' && (*W)(2)!='o')
                        (*W)+='e';
                      break;
                    }
                    case 'u': {
                      if (!((*W).Length()>3 && IsVowel((*W)(1)) && IsVowel((*W)(2))))
                        (*W)+='e';
                      break;
                    }
                    case 'z': {
                      if ((*W).EndsWith("izz") && (*W).Length()>3 && ((*W)(3)=='h' || (*W)(3)=='u'))
                        (*W).End--;
                      else if ((*W)(1)!='t' && (*W)(1)!='z')
                        (*W)+='e';
                      break;
                    }
                    case 'k': {
                      if ((*W).EndsWith("uak"))
                        (*W)+='e';
                      break;
                    }
                    case 'b': case 'c': case 's': case 'v': {
                      if (!(
                        ((*W)(0)=='b' && ((*W)(1)=='m' || (*W)(1)=='r')) ||
                        (*W).EndsWith("ss") || (*W).EndsWith("ias") || (*W)=="zinc"
                      ))
                        (*W)+='e';
                      break;
                    }
                  }
                }
              }
            }
            else{
              (*W).End=j;
              return false;
            }
          }
        }
        (*W).Type|=TypesStep1b[i];
        return true;
      }
    }
    return false;
  }
  bool Step1c(Word *W) {
    if ((*W).Length()>2 && /*tolower*/((*W)(0))=='y' && IsConsonant((*W)(1))) {
      (*W).Letters[(*W).End]='i';
      return true;
    }
    return false;
  }
  bool Step2(Word *W, const U32 R1) {
    for (int i=0;i<NUM_SUFFIXES_STEP2;i++) {
      if ((*W).EndsWith(SuffixesStep2[i][0]) && SuffixInRn(W, R1, SuffixesStep2[i][0])) {
        (*W).ChangeSuffix(SuffixesStep2[i][0], SuffixesStep2[i][1]);
        (*W).Type|=TypesStep2[i],(*W).Suffix|=TypesStep2Suffix[i];
        return true;
      }
    }
    if ((*W).EndsWith("logi") && SuffixInRn(W, R1, "ogi")) {
      (*W).End--;
      return true;
    }
    else if ((*W).EndsWith("li")) {
      if (SuffixInRn(W, R1, "li") && IsLiEnding((*W)(2))) {
        (*W).End-=2;
        (*W).Type|=AdverbOfManner;
        return true;
      }
      else if ((*W).Length()>3) {
        switch((*W)(2)) {
            case 'b': {
              (*W).Letters[(*W).End]='e';
              (*W).Type|=AdverbOfManner;
              return true;              
            }
            case 'i': {
              if ((*W).Length()>4) {
                (*W).End-=2;
                (*W).Type|=AdverbOfManner;
                return true;
              }
              break;
            }
            case 'l': {
              if ((*W).Length()>5 && ((*W)(3)=='a' || (*W)(3)=='u')) {
                (*W).End-=2;
                (*W).Type|=AdverbOfManner;
                return true;
              }
              break;
            }
            case 's': {
              (*W).End-=2;
              (*W).Type|=AdverbOfManner;
              return true;
            }
            case 'e': case 'g': case 'm': case 'n': case 'r': case 'w': {
              if ((*W).Length()>(U32)(4+((*W)(2)=='r'))) {
                (*W).End-=2;
                (*W).Type|=AdverbOfManner;
                return true;
              }
            }
        }
      }
    }
    return false;
  }
  bool Step3(Word *W, const U32 R1, const U32 R2) {
    bool res=false;
    for (int i=0;i<NUM_SUFFIXES_STEP3;i++) {
      if ((*W).EndsWith(SuffixesStep3[i][0]) && SuffixInRn(W, R1, SuffixesStep3[i][0])) {
        (*W).ChangeSuffix(SuffixesStep3[i][0], SuffixesStep3[i][1]);
        (*W).Type|=TypesStep3[i],(*W).Suffix|=TypesStep3Suffix[i];
        res=true;
        break;
      }
    }
    if ((*W).EndsWith("ative") && SuffixInRn(W, R2, "ative")) {
      (*W).End-=5;
      (*W).Type|=Suffix,(*W).Suffix|=SuffixIVE;
      return true;
    }
    if ((*W).Length()>5 && (*W).EndsWith("less")) {
      (*W).End-=4;
      (*W).Type|=AdjectiveWithout;
      return true;
    }
    return res;
  }
  bool Step4(Word *W, const U32 R2){
    bool res=false;
    for (int i=0;i<NUM_SUFFIXES_STEP4;i++) {
      if ((*W).EndsWith(SuffixesStep4[i]) && SuffixInRn(W, R2, SuffixesStep4[i])) {
        (*W).End-=strlen(SuffixesStep4[i])-(i>17);
        if (i!=10 || (*W)(0)!='m')
          (*W).Type|=TypesStep4[i],(*W).Suffix|=TypesStep4Suffix[i];
        if (i==0 && (*W).EndsWith("nti")) {
          (*W).End--;
          res=true;
          continue;
        }
        return true;
      }
    }
    return res;
  }
  bool Step5(Word *W, const U32 R1, const U32 R2){
    if ((*W)(0)=='e' && (*W)!="here") {
      if (SuffixInRn(W, R2, "e"))
        (*W).End--;
      else if (SuffixInRn(W, R1, "e")) {
        (*W).End--;
        (*W).End+=EndsInShortSyllable(W);
      }
      else
        return false;
      return true;
    }
    else if ((*W).Length()>1 && (*W)(0)=='l' && SuffixInRn(W, R2, "l") && (*W)(1)=='l') {
      (*W).End--;
      return true;
    }
    return false;
  }
public:
  bool Stem(Word *W){
    /*if ((*W).Length()<2){
      Hash(W);
      return false;
    }*/
    bool res = TrimStartingApostrophe(W);

    if (ProcessPrefixes(W)) res = true;
    if (ProcessSuperlatives(W)) res = true;
    for (int i=0;i<NUM_EXCEPTIONS1;i++) {
      if ((*W)==Exceptions1[i][0]) {
        if (i<11){
          size_t len=strlen(Exceptions1[i][1]);
          memcpy(&(*W).Letters[(*W).Start], Exceptions1[i][1], len);
          (*W).End=(*W).Start+len-1;
        }
        Hash(W);
        (*W).Type|=TypesExceptions1[i];
        return (i<11);
      }
    }

    // Start of modified Porter2 Stemmer
    MarkYsAsConsonants(W);
    U32 R1=GetRegion1(W), R2=GetRegion(W,R1);
    if (Step0(W)) res = true;
    if (Step1a(W)) res = true;
    for (int i=0;i<NUM_EXCEPTIONS2;i++) {
      if ((*W)==Exceptions2[i]) {
        Hash(W);
        (*W).Type|=TypesExceptions2[i];
        return res;
      }
    }
    if (Step1b(W,R1)) res = true;
    if (Step1c(W)) res = true;
    if (Step2(W,R1)) res = true;
    if (Step3(W,R1,R2)) res = true;
    if (Step4(W,R2)) res = true;
    if (Step5(W,R1,R2)) res = true;

    for (U8 i=(*W).Start;i<=(*W).End;i++) {
      if ((*W).Letters[i]=='Y')
        (*W).Letters[i]='y';
    }
    if (!(*W).Type || (*W).Type==Plural) {
      if ((*W).MatchesAny(MaleWords, NUM_MALE_WORDS))
        res = true, (*W).Type|=Male;
      else if (W->MatchesAny(FemaleWords, NUM_FEMALE_WORDS))
        res = true, (*W).Type|=Female;
      else if (W->MatchesAny(ArticleWords, NUM_ARTICLE_WORDS))
        res = true, (*W).Type|=Article;
      else if (W->MatchesAny(ConjWords, NUM_CONJ_WORDS))
        res = true, (*W).Type|=Conjunction;
      else if (W->MatchesAny(ApoWords, NUM_APO_WORDS))
        res = true, (*W).Type|=Adposition;
      else if (W->MatchesAny(ConAdVerPrepWords, NUM_CAVER_WORDS))
        res = true, (*W).Type|=ConjunctiveAdverb;
      else if ((x.blpos<451531986) && W->MatchesAny(VerbWords1, NUM_VERB)) //77,06% disable 
        res = true, (*W).Type|=Verb;
      else if (W->MatchesAny(Numbers, NUM_NUM))
        res = true, (*W).Type|=Number;
      else {
          for (int i=0;i<14;i++) {
            if (W->MatchesAnyP(&Pronouns[i][0], 5)) {
                res = true, (*W).Type|=Pronoun;
                break;
            }
          }
      }
    }
    Hash(W);
    return res;
  }
};


// ===== XML model ======

#define CacheSize (1<<5)
struct XMLAttribute {
    U32 Name, Value, Length;
};

struct XMLContent {
    U32 Data, Length, Type;
};

struct XMLTag {
    U32 Name, Length;
    int Level;
    bool EndTag, Empty;
    XMLContent Content;
    struct XMLAttributes {
        XMLAttribute Items[4];
        U32 Index;
    } Attributes;    
};

struct XMLTagCache {
    XMLTag Tags[CacheSize];
    U32 Index;
};

enum ContentFlags {
    xText        = 0x001,
    xNumber      = 0x002,
    xDate        = 0x004,
    xTime        = 0x008,
    xURL         = 0x010,
    xLink        = 0x020,
    xCoordinates = 0x040,
    xTemperature = 0x080,
    xISBN        = 0x100
};

enum XMLState {
    xNone               = 0,
    xReadTagName        = 1,
    xReadTag            = 2,
    xReadAttributeName  = 3,
    xReadAttributeValue = 4,
    xReadContent        = 5,
    xReadCDATA          = 6,
    xReadComment        = 7
};

#define DetectContent() { \
  if ((x.c4&0xF0F0F0F0)==0x30303030){ \
    int i = 0, j = 0; \
    while ((i<4) && ( (j=(x.c4>>(8*i))&0xFF)>=0x30 && j<=0x39 )) \
      i++; \
\
    if (i==4 && ( ((c8&0xFDF0F0FD)==0x2D30302D && buf(9)>=0x30 && buf(9)<=0x39) || ((c8&0xF0FDF0FD)==0x302D302D) )) \
      (*Content).Type |= xDate; \
  } \
  else if (((c8&0xF0F0FDF0)==0x30302D30 || (c8&0xF0F0F0FD)==0x3030302D) && buf(9)>=0x30 && buf(9)<=0x39){ \
    int i = 2, j = 0; \
    while ((i<4) && ( (j=(c8>>(8*i))&0xFF)>=0x30 && j<=0x39 )) \
      i++; \
\
    if (i==4 && (x.c4&0xF0FDF0F0)==0x302D3030) \
      (*Content).Type |= xDate; \
  } \
\
  if ((x.c4&0xF0FFF0F0)==(0x30003030+COLON*256*256) && buf(5)>=0x30 && buf(5)<=0x39 && ((buf(6)<0x30 || buf(6)>0x39) || ((c8&0xF0F0FF00)==(0x30300000+COLON*256) && (buf(9)<0x30 || buf(9)>0x39)))) \
    (*Content).Type |= xTime; \
\
  if ((*Content).Length>=8 && (c8&0x80808080)==0 && (x.c4&0x80808080)==0) \
    (*Content).Type |= xText; \
\
  if ((c8&0xF0F0FF)==0x3030C2 && (x.c4&0xFFF0F0FF)==0xB0303027){ \
    int i = 2; \
    while ((i<7) && buf(i)>=0x30 && buf(i)<=0x39) \
      i+=(i&1)*2+1; \
\
    if (i==10) \
      (*Content).Type |= xCoordinates; \
  } \
\
  if ((x.c4&0xFFFFFA)==0xC2B042 && B!=0x47 && (((x.c4>>24)>=0x30 && (x.c4>>24)<=0x39) || ((x.c4>>24)==0x20 && (buf(5)>=0x30 && buf(5)<=0x39)))) \
    (*Content).Type |= xTemperature; \
\
  if (B>=0x30 && B<=0x39) \
    (*Content).Type |= xNumber; \
\
  if (lastCW==cwISBN && buf(4)==UPPER) \
    (*Content).Type |= xISBN; \
} 
U32 xlU1=0,xlU2=0,xlU3=0,xlU4=0;
bool isXML=false;
struct XMLModel1 {
    XMLTagCache Cache;
    XMLState State, pState;
    U32 c8, WhiteSpaceRun, pWSRun, IndentTab, IndentStep, LineEnding,lastState;
    U32 StateBH[8];

    void Init() {
        c8=0;
        Reset();
    }
    void Reset() {
        State=xNone, pState=xNone, 
        WhiteSpaceRun=0, pWSRun=0, IndentTab=0, IndentStep=2, LineEnding=2,lastState=0;
        memset(&Cache, 0, sizeof(XMLTagCache));
        memset(&StateBH, 0, sizeof(StateBH));  
    }
    int p() {
        xlU4=0;
        if (x.bpos==0) {
            U8 B=(U8)x.c4;
            XMLTag *pTag=&Cache.Tags[ (Cache.Index-1)&(CacheSize-1) ], *Tag = &Cache.Tags[ Cache.Index&(CacheSize-1)];
            XMLAttribute *Attribute=&((*Tag).Attributes.Items[ (*Tag).Attributes.Index&3]);
            XMLContent *Content=&(*Tag).Content;
            pState=State;
            c8=(c8<<8)|buf(5);
            if ((B==0x09 || B==0x20) && (B==(U8)(x.c4>>8) || !WhiteSpaceRun)) {
                WhiteSpaceRun++;
                IndentTab = (B==0x09);
            }
            else{
                if ((State==xNone || (State==xReadContent && (*Content).Length<=LineEnding+WhiteSpaceRun)) && WhiteSpaceRun>1+IndentTab && WhiteSpaceRun!=pWSRun) {
                    IndentStep=abs((int)(WhiteSpaceRun-pWSRun));
                    pWSRun=WhiteSpaceRun;
                }
                WhiteSpaceRun=0;
            }
            if (B==0x0A)
            LineEnding=1+((U8)(x.c4>>8)==0x0D);
            if(State!=xNone) lastState=x.blpos;
            switch (State) {
            case xNone: {
                    if (B==LESSTHAN){
                        State=xReadTagName;
                        memset(Tag, 0, sizeof(XMLTag));
                        (*Tag).Level=((*pTag).EndTag || (*pTag).Empty)?(*pTag).Level:(*pTag).Level+1;
                    }
                    if ((*Tag).Level>1)
                    DetectContent();
                    
                    xlU4=(hash(pState, State, ((*pTag).Level+1)*IndentStep - WhiteSpaceRun));
                    break;
                }
            case xReadTagName: {
                    if ((*Tag).Length>0 && (B==0x09 || B==0x0A || B==0x0D || B==0x20))
                    State = xReadTag;
                    else if ((B>127) || (B==COLON || B==FIRSTUPPER || B==UPPER || B==0x5F|| (B>='a' && B<='z')) || ((*Tag).Length>0 && (B==0x2D || B==0x2E || (B>='0' && B<='9')))) {
                        (*Tag).Length++;
                        (*Tag).Name=(*Tag).Name * 263 * 32 + (B&0xDF);
                    }
                    else if (B == GREATERTHAN) {
                        if ((*Tag).EndTag){
                            State=xNone;
                            Cache.Index++;
                        }
                        else
                        State = xReadContent;
                    }
                    else if (B!=0x21 && B!=0x2D && B!=0x2F && B!=0x5B) {
                        State=xNone;
                        Cache.Index++;
                    }
                    else if ((*Tag).Length==0) {
                        if (B==0x2F) {
                            (*Tag).EndTag=true;
                            (*Tag).Level=max(0,(*Tag).Level-1);
                        }
                        else if (x.c4==(LESSTHAN*256*256*256+0x212D2D)) { //LESSTHAN
                            State=xReadComment;
                            (*Tag).Level=max(0,(*Tag).Level-1);
                        }
                    }

                    if ((*Tag).Length==1 && (x.c4&0xFFFF00)==(LESSTHAN*256*256+0x2100)) {//LESSTHAN
                        memset(Tag, 0, sizeof(XMLTag));
                        State=xNone;
                    }
                    int i=1;
                    do {
                        pTag=&Cache.Tags[ (Cache.Index-i)&(CacheSize-1) ];
                        i+=1+((*pTag).EndTag && Cache.Tags[ (Cache.Index-i-1)&(CacheSize-1) ].Name==(*pTag).Name);
                    } while ( i<CacheSize && ((*pTag).EndTag || (*pTag).Empty) );

                    xlU4=(hash(pState*8+State, hash((*Tag).Name, (*Tag).Level),hash((*pTag).Name, (*pTag).Level!=(*Tag).Level) ) );
                    break;
                }
            case xReadTag: {
                    if (B==0x2F)
                    (*Tag).Empty=true;
                    else if (B==GREATERTHAN){
                        if ((*Tag).Empty){
                            State=xNone;
                            Cache.Index++;
                        }
                        else
                        State=xReadContent;
                    }
                    else if (B!=0x09 && B!=0x0A && B!=0x0D && B!=0x20) {
                        State = xReadAttributeName;
                        (*Attribute).Name=B;
                    }
                    xlU4=(hash(pState, State, hash((*Tag).Name, B, (*Tag).Attributes.Index)));
                    break;
                }
            case xReadAttributeName: {
                    if ((x.c4&0xFFF0)==(EQUALS*256+SPACE) && (B==0x22 || B==0x27)) {
                        State=xReadAttributeValue;
                        if ((c8&0xDFDF)==0x4852 && (x.c4&0xDFDF0000)==0x45460000)
                        (*Content).Type|=xLink;
                    }
                    else if (B!=0x22 && B!=0x27 && B!=EQUALS)
                    (*Attribute).Name=(*Attribute).Name*263*32+(B&0xDF);

                    xlU4=(hash(pState*8+State, (*Attribute).Name, hash((*Tag).Attributes.Index, (*Tag).Name, (*Content).Type )));
                    break;
                }
            case xReadAttributeValue: {
                    if (B==0x22 || B==0x27) {
                        (*Tag).Attributes.Index++;
                        State=xReadTag;
                    }
                    else{
                        (*Attribute).Value=(*Attribute).Value*263*32+B;
                        (*Attribute).Length++;
                        if (lastCW==cwHTTP && ((x.c4>>8)==(COLON*256*256+0x2F2F) )) // HTTP :// s://
                        (*Content).Type|=xURL;
                    }
                    xlU4=(hash(pState, State, hash((*Attribute).Name, (*Content).Type )));
                    break;
                }
            case xReadContent: {
                    if (B==LESSTHAN) {
                        State=xReadTagName;
                        Cache.Index++;
                        memset(&Cache.Tags[ Cache.Index&(CacheSize-1) ], 0, sizeof(XMLTag));
                        Cache.Tags[Cache.Index&(CacheSize-1) ].Level=(*Tag).Level+1;
                    }
                    else{
                        (*Content).Length++;
                        (*Content).Data=(*Content).Data*997*16+B;

                        DetectContent();
                    }
                    xlU4=(hash(pState, State,hash( (*Tag).Name, x.c4&0xC0FF )));
                    break;
                }
            case xReadComment: {
                    if ((x.c4&0xFFFFFF)==(0x2D2D00+GREATERTHAN)) { // -->
                        State=xNone;
                        Cache.Index++;
                    }
                    xlU4=(hash(pState, State));
                    break;
                }
            }
            StateBH[pState]=(StateBH[pState]<<8)|B;
            pTag=&Cache.Tags[ (Cache.Index-1)&(CacheSize-1) ];
            // set context if last state was less then 256 bytes ago
            isXML=(x.blpos-lastState)<64;
            xlU1=hash(State, (*Tag).Level, hash(pState*2+(*Tag).EndTag, (*Tag).Name));
            xlU2=hash((*pTag).Name, State*2+(*pTag).EndTag,hash( (*pTag).Content.Type, (*Tag).Content.Type));
            xlU1=hash(State*2+(*Tag).EndTag, (*Tag).Name,hash( (*Tag).Content.Type, x.c4&0xE0FF));
        }
        U8 s = ((StateBH[State]>>(28-x.bpos))&0x08) |
        ((StateBH[State]>>(21-x.bpos))&0x04) |
        ((StateBH[State]>>(14-x.bpos))&0x02) |
        ((StateBH[State]>>( 7-x.bpos))&0x01) |
        ((x.bpos)<<4);
        return (s<<3)|State;
    }
};
// Precalculated state to prediction for STA7
static const short pre1[256]={
    0,  -89,   88, -125,    0,    0,  124, -148,  -36,  -36,   35,   35,  147, -165,  -59,  -59,
    0,    0,   58,   58,  164, -178,  -76,  -76,  -23,  -23,   22,   22,   75,   75,  177, -189,
  -89,  -40,    0,   39,   88,  188, -198, -100,  -53,  -17,   16,   53,   99,  197, -207, -109,
  -64,  -30,    0,   30,   64,  109,  205, -214, -118,  -74,  -41,  -14,   13,   41,   73,  117,
  213, -220, -125,  -82,  -51,  -25,   24,   50,   81,  124,  219, -226, -131,  -89,  -59,  -34,
   33,   58,   88,  131,  225, -232, -137,  -96,  -66,  -42,   41,   66,   95,  136,  231, -237,
 -143, -102,  -73,  -50,   49,   72,  101,  142,  236, -241, -148, -107,  -79,  -56,   55,   78,
  106,  147,  240, -246, -152, -112,  -84,  -62,   61,   83,  111,  152,  244, -250, -157, -117,
  -89,  -67,   67,   88,  116,  156,  249, -253, -161, -121,  -94,  -72,   72,   93,  120,  160,
  252, -257, -165, -125,  -98,  -77,   76,   97,  124,  164,  255, -261, -168, -129, -102,  -81,
   80,  101,  128,  167,  259, -264, -172, -132, -106,  -85,   85,  105,  131,  171,  262, -267,
 -175, -136, -109,  -89,   88,  109,  135,  174,  265, -270, -178, -139, -113,  -93,   92,  112,
  138,  177,  268, -273, -181, -142, -116,  -96,   95,  115,  141,  180,  271, -275, -184, -145,
 -119,  -99,   99,  118,  144,  183,  273, -278, -186, -148, -122, -102,  102,  121,  147,  185,
  276, -280, -189, -151, -105,  105,  150,  188,  278, -283, -191, -153,  152,  190,  281, -194,
 -156,  155,  193, -196, -158,  157,  195, -160,  159, -163,  162, -165,  164, -167,  166, -169,
};

struct DirectStateMap {
    StateMap *sm;
    U32 *cxt;
    U32 mask;
    U8 *CxtState;
    int index;
    int count;
    const U8 *nn;
    Mix *mmm;  
    int pu;
    void __attribute__ ((noinline)) Init(int m,int c,const U8 *nn1) {
        nn=nn1;pu=0;
        mask=(1<<m)-1,index=0,count=c;
        assert(mask>=64 && ((mask+1)&mask)==0); 
        alloc(cxt,c);
        alloc(CxtState,(mask+1));
        alloc(sm,c);
        alloc(mmm,c);
        for (int i=0; i<count; i++) {
            sm[i].Init(256,nn1);
            mmm[i].Init(256);
        }
    }
    U8 next(int state, int y) {
        return nn[state*4+y];
    }
    void reset() {
        memset(CxtState, 0, mask);

    }
    void set(U32 cx,int y) {
        assert(cxt[index]>=0 && cxt[index]<=mask);
        assert(index<count);
        CxtState[cxt[index]]=next(CxtState[cxt[index]],y);       // update state
        cxt[index]=(cx)&mask;                                     // get new context
        const U8 state=CxtState[cxt[index]];
        sm[index].set(state);    // predict from new context
        if (index==0) pu=stretch(sm[index].pr);
        else mmm[index-1].update(x.y),  pu=clp(mmm[index-1].pp(pu,stretch(sm[index].pr),state));
        x.mxInputs1.add(pre1[state]);
        prediction_index--;
        index++;
    }
    void mix() {
        x.mxInputs1.add((pu)>>2);
        index=pu=0;
    }
};


// Brackets
static const U8 fcy[128]={
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 5, 0, 0, 0, 0, 6, 1, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 
2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};
// First char
static const U8 fcq[128]={
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 4, 5, 0, 0,
2, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0,
2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};

// Predictor
const U32 primes[14]={0, 257,251,241,239,233,229,227,223,211,199,197,193,191};
const U32 tri[4]={0,4,3,7}, trj[4]={0,6,6,12};

// Parameters
const U32 c_s[27]= {28, 32, 32, 32, 34, 31, 33, 33, 35, 35, 29, 32, 33, 34, 30, 36, 31, 32, 32, 32, 32, 32, 33, 32, 32, 32, 32};
const U32 c_s3[27]={43, 32, 32, 32, 34, 29, 32, 33, 37, 35, 33, 28, 31, 35, 28, 30, 33, 34, 32, 32, 32, 32, 32, 32, 32, 32, 32};
const U32 c_s4[27]={ 9,  8, 12, 12,  8, 12, 15,  8,  8, 12, 10,  7,  7,  8, 8, 13, 13, 14,  8,  8, 12, 12, 12, 12, 12, 12, 12};
 const int e_l[8]={1830, 1997, 1973, 1851, 1897, 1690, 1998, 1842};

const int MAXLEN=62; // longest allowed match + 1
U32 t[16]; 
int isParagraph;
U32 fails=0;
int c1,c2,c3;
U8 words,spaces,numbers;
U32 h,word0,wshift,x4,x5,isMatch,firstWord,linkword,senword;
U32 xword1=0,numberA=0;
U32 number0,number1,numlen0,numlen1,mybenum;
// First char context index, bracket/first char context index (max value 7)
U32 FcIdx=0,BrFcIdx;
U32 AH1=0,AH2=0x765BA55C; // APM hash
U32  failz=0, failcount=0;
int nl,nl1,col,fc,above,above1;
U32 t1[0x100];
U32 t2[0x10000];
int wp[0x10000];
int np[0x10000];
U16 ind3[0x2000000];
U32 context1_ind3=0,cxtind3=0;
U32 lastWT=0;
U32 cwSTR=0x10000,cwCOLON=0x10000;
// 3 bit stream 
U32 o3bState, n3bState, stream3bR, stream3b;
U32 stream3bMask=0,stream3bMask1=0,stream3bRMask1=0,stream3bRMask2=0;
// 2 bit stream 
U32 o2bState, n2bState, stream2bR, stream2b; 
U32 stream2bMask=0;
U32 o4bState, n4bState, stream4bR, stream4b;
int ordX,ordW; // Order x count, Order word count - max 6
U8  buffer[0x1000000];
enum {BMASK=0xffffff};
U8    cwbuf[0x1000];
enum {CBMASK=0xfff};
int cwpos=0;

Word StemWords[4];
Word *cWord, *pWord;
EnglishStemmer StemmerEN;
int StemIndex=0;
bool isMath=false;

int deccode=0;

bool isText=false;
bool skipM1=false;
U32 word00;
int utf8left=0;
U32 indirectBrByte=0,indirectByte=0,indirectWord0Pos=0, indirectWord=0,u8w=0,indirectNumberd0Pos=0;
U32 sVerb=0;
bool lastArt=false;
int dcw=0,dcwl=0;

bool inBR=false;
bool isHTTAG=false;
bool wasTag=false;
bool skipSeeExternal=false;
bool isCategory=false;

U32 lastPTOP=-1,pageParag=0,pageSent=0;
bool isLongTOP=false;
bool isPageStarted=false;
bool wasVerb=false;
bool wasNoun=false;
bool nestList=false;
U32 wasVerbH=0,wasNounH=0; 
U32 s5bByte[8]{};
U32 stream5b=0;
U32 s2Word0[8*16]{};
U32 stream6b=0;    
U8 indirect2[256]{};
U8 indirect3[256*256]{};

enum PageState {
    PNone     =0,
    PTemplate =1,
    PText     =2,
    PTopic    =4,
    PCategory =8,
};    
U32 PState=0,PStateH=0;

bool isNowiki=false;
int oldwt1=0;
U32 wt3b=0,wt3cxt=0,wt3cxtW=0,wt3cxtW1=0,wt4cxtW=0,wt4cxtW1=0;

int pr; // Our most important variable - final prediction

StateMap1 smA[3];
SmallStationaryContextMap scmA[3];   
Mixer1 mxA[18]; 
Mixer1 mxA1[6]; 
Mixer1 mxA2[4]; 
// Predictors are:
// medium state memory, per context max 7 uniqe contexts state sets
// for average sized contexts
ContextMap3 cmC[6]; 
ContextMap3 cmC44;
ContextMap3 cmcr[1+8];
ContextMap3 cmcr2[4];
// small state memory, per context max 3 uniqe contexts state sets
// for large amount of small contexts

// large state memory, per context max 14 uniqe contexts state sets
// for large amount of large contexts
ContextMap3 cmC2[21];
ContextMap4 cmC4[9];
ContextMap4 cmCR[9];

DirectStateMap dcsm;
DirectStateMap dcsm0;
DirectStateMap dcsm1;
DirectStateMap dcsm2;
DirectStateMap dcsmN;

StationaryMap maps1;
StationaryMap maps2;
APM<8>   apmA0;
APM<16>  apmA1;
APM<16>  apmA2;
APM<18>  apmA3;
APM<17>  apmA4;
APM<17>  apmA5;

RunContextMap rcmA[1];

BracketContext<U8> brcxt;
BracketContext<U8> qocxt;
BracketContext<U8> fccxt;
ColumnContext colcxt;
WordsContext worcxt; // decoded 
WordsContext worcxt1;
WordsContext worcxt2;
WordsContext worcxt3; // tag
WordsContext worcxt0; // undecoded


SentenceContext sencxt;
SentenceContext sencxtL; // lists '*'
SentenceContext sencxtT; // table
SentenceContext sencxtCL; // wikilinks

WordsContext *simiwor;
BracketContext<U16> htcxt;

XMLModel1 xml;

//int mxskip=0,mxskiptotal=0;

Mix mmmO[4];

void PredictorInit() { 
    n3bState=n2bState=0xffffffff;
    pr=2048;

    // Match
    smA[0].Init(1<<9,1023);
    smA[1].Init(1<<19,1023);
    smA[2].Init(1<<16,1023);
    
    for (int i=0;i<4;i++) mmmO[i].Init(256);

    scmA[0].Init(8);
    scmA[1].Init(9); 
    scmA[2].Init(8); 

    maps1.Init(16,8);
    maps2.Init(16,8);
    const int dccount=6;
    dcsm.Init(28,dccount-1, &STA7[0][0] );
    dcsm0.Init(28,dccount, &STA7[0][0] );
    dcsm1.Init(20,2, &STA7[0][0] );
    dcsm2.Init(26,3, &STA7[0][0] );
    dcsmN.Init(25,3, &STA7[0][0] );
    
    // Mixers      size,  shift, err, errmul 
    mxA[0].Init(  0x8000,  75,  8, 14);
    mxA[1].Init(   6*256,  75,  8, 14);
    mxA[2].Init(   6*256,  30,  1, 38);
    mxA[3].Init(   8*256,  31,  1, 34);
    mxA[4].Init(   6*256,  53,  1, 23);
    mxA[5].Init(  32*256,  79,  1, 24);
    mxA[6].Init(  0x4000,  75,  1, 20);
    mxA[7].Init(  0x4000,  55,  1, 24);
    mxA[8].Init( 0x20000,  55,  1, 24);
    mxA[9].Init( 0x20000,  55,  1, 24);
    mxA[10].Init(0x10000,  55,  1, 24);
    mxA[11].Init( 0x4000,  55,  1, 24);
    mxA[12].Init( 32*256,   6,  1,  4);
    mxA[13].Init( 32*256,   6,  1,  4);
    mxA[14].Init( 32*256,  55,  1, 24);
    mxA[15].Init( 16*256,  55,  1, 24);
    mxA[16].Init(   1024,   6,  1,  4);
    mxA[17].Init(   2048,   6,  1,  4);

    mxA1[0].Init(8*7*4,   6,  0,  4);
    mxA1[1].Init(    1,   6,  0,  4);
    mxA1[2].Init( 2048,  6,  1, 4);
    mxA1[3].Init( 2048,  6,  1, 4);
    mxA1[4].Init( 2048,  6,  1, 4);
    mxA1[5].Init( 2048,  6,  1, 4);
     
    mxA2[0].Init(  0x100,  30,  1, 14);
    mxA2[1].Init(  0x100,  30,  1, 14);
    mxA2[2].Init( 0x8000,  30,  1, 14);
    mxA2[3].Init(0x10000,  30,  1, 14);

    apmA0.Init();
    apmA1.Init();
    apmA2.Init();
    apmA3.Init();
    apmA4.Init();
    apmA5.Init();

    rcmA[0].Init(1*4096*4096,6);

    x.mxInputs1.ncount=544;
    x.mxInputs2.ncount=32;
    x.mxInputs4.ncount=16;

    // Provide inputs array info to mixers
    for (int i=0;i<18;i++)
        mxA[i].setTxWx(x.mxInputs1.ncount,&x.mxInputs1.n[0]);
    // Final mixer

    for (int i=0;i<5;i++)
        mxA1[i].setTxWx(x.mxInputs2.ncount,&x.mxInputs2.n[0]);
    
    mxA1[5].setTxWx(x.mxInputs4.ncount,&x.mxInputs4.n[0]);
    //
    for (int i=0;i<4;i++)
        mxA2[i].setTxWx(num_models,model_predictions1);
 
//                mem           c  rmul    prmul    smul   sta          s4      keep, st2 sta2pr
    cmC2[0].Init(  8*4096*4096, 3, c_s[0], c_s3[0],&STA6[0][0], c_s4[0],0xf0,1,&st2_p1[0]);
    cmC2[1].Init( 16*4096*4096, 1, c_s[1], c_s3[1],&STA6[0][0], c_s4[1],0xf0,0,&st2_p1[0]);
    cmC2[2].Init(  8*4096*4096, 1, c_s[2], c_s3[2],&STA6[0][0], c_s4[2],0xf0,0,&st2_p1[0]);
    cmC2[3].Init(  8*4096*4096, 1, c_s[3], c_s3[3],&STA6[0][0], c_s4[3],0xf0,0,&st2_p1[0]);
    cmC2[4].Init( 16*4096*4096, 2, c_s[4], c_s3[4],&STA6[0][0], c_s4[4],0xf0,1,&st2_p1[0]);
    cmC2[5].Init( 16*4096*4096, 6, c_s[5], c_s3[5],&STA6[0][0], c_s4[5],0xf0,1,&st2_p1[0]);
    cmC2[6].Init(      64*4096, 1, c_s[6], c_s3[6],&STA1[0][0], c_s4[6],0x00,1,&st2_p1[0]);
    cmC2[7].Init(  2*4096*4096, 1, c_s[7], c_s3[7],&STA5[0][0], c_s4[7],0xf0,1,&st2_p1[0]);
    cmC2[8].Init(  8*4096*4096, 4, c_s[8], c_s3[8],&STA4[0][0], c_s4[8],0x00,1,&st2_p1[0]);
    
    cmC4[0].Init(      32*4096, 2, c_s[9], c_s3[9],&STA6[0][0], c_s4[9],0x00,0);
    cmC4[1].Init(    8*32*4096, 3, c_s[10],c_s3[10],&STA7[0][0],c_s4[10],0x00,1);
    cmC4[2].Init(    8*32*4096, 4, c_s[11],c_s3[11],&STA2[0][0],c_s4[11],0x00,1);
    cmC4[4].Init(      32*4096, 6, c_s[12],c_s3[12],&STA7[0][0],c_s4[12],0x00,1);
   
    cmC[0].Init(    2* 16*4096, 7, c_s[13],c_s3[13],&STA2[0][0],c_s4[13],0x00,1,&st2_p1[0]);
    cmC44.Init(    1*4096*4096, 4, c_s[13],c_s3[13],&STA2[0][0],c_s4[13],0xf0,1,&st2_p1[0]);
     
    cmcr[0].Init(  1*4096*4096, 2, c_s[13],c_s3[13],&STA2[0][0],c_s4[13],0xf0,1,&st2_p1[0]);
    for (int i=0;i<8;i++)
        cmcr[i+1].Init(  2048*4096, 3, c_s[13],c_s3[13],&STA6[0][0],c_s4[13],0x00,0,&st2_p1[0]);
    
    for (int i=0;i<4;i++)
    cmcr2[i].Init( 1*4096*4096, 3, c_s[13],c_s3[13],&STA6[0][0],c_s4[13],0x00,0,&st2_p1[0]);
    
    cmC[1].Init(     64*2*4096, 3, c_s[14],c_s3[14],&STA5[0][0],c_s4[14],0xf0,0,&st2_p0[0]);
    cmC[2].Init(        2*4096, 2, c_s[15],c_s3[15],&STA2[0][0],c_s4[15],0xf0,0,&st2_p0[0]);
    
    cmC4[3].Init(     128*4096, 2, c_s[16],c_s3[16],&STA1[0][0],c_s4[16],0x00,0);
    
    cmC2[9].Init(  8*4096*4096, 4, c_s[17],c_s3[17],&STA6[0][0],c_s4[17],0xf0,1,&st2_p1[0]);
    cmC2[10].Init( 8*4096*4096, 6, c_s[18],c_s3[18],&STA5[0][0],c_s4[18],0xf0,1,&st2_p1[0]);
    cmC2[11].Init( 8*4096*4096, 5, c_s[19],c_s3[19],&STA5[0][0],c_s4[19],0xf0,1,&st2_p1[0]);
    cmC2[12].Init( 8*4096*4096, 2, c_s[20],c_s3[20],&STA6[0][0],c_s4[20],0xf0,1,&st2_p1[0]);
    cmC2[13].Init(16*4096*4096, 2, c_s[21],c_s3[21],&STA6[0][0],c_s4[21],0xf0,1,&st2_p1[0]);
    
    cmC[3].Init(       32*4096, 2, c_s[22],c_s3[22],&STA2[0][0],c_s4[22],0x00,1,&st2_p1[0]);
    
    cmC2[14].Init( 2*4096*4096, 1, c_s[23],c_s3[23],&STA6[0][0],c_s4[23],0xf0,1,&st2_p1[0]);
    cmC2[15].Init(   8*64*4096, 1, c_s[24],c_s3[24],&STA1[0][0],c_s4[24],0x00,0,&st2_p0[0]);
    
    cmC[4].Init(      512*4096, 1, c_s[25],c_s3[25],&STA1[0][0],c_s4[25],0xf0,1,&st2_p1[0]);
    cmC[5].Init(      512*4096, 1, c_s[26],c_s3[26],&STA1[0][0],c_s4[26],0xf0,1,&st2_p1[0]);
   
    cmC2[16].Init(   2048*4096, 1, c_s[17],c_s3[17],&STA6[0][0],c_s4[17],0xf0,1,&st2_p1[0]);
    cmC2[17].Init( 2*4096*4096, 2, c_s[17],c_s3[17],&STA6[0][0],c_s4[17],0xf0,1,&st2_p1[0]);
   
    cmC4[6].Init(    2*16*4096, 1, c_s[5], c_s3[5],&STA6[0][0], c_s4[5],0x00,0);
    
    cmC4[7].Init(      16*4096, 4, c_s[12],c_s3[12],&STA2[0][0],c_s4[12],0x00,1);
    cmC4[8].Init(      16*4096, 2, c_s[12],c_s3[12],&STA2[0][0],c_s4[12],0x00,1);
    
    cmC2[18].Init( 4*4096*4096, 5, c_s[5], c_s3[5],&STA6[0][0], c_s4[5],0xf0,1,&st2_p1[0]);
    cmC2[19].Init( 2*4096*4096, 1, c_s[17],c_s3[17],&STA6[0][0],c_s4[17],0xf0,1,&st2_p1[0]);
    cmC2[20].Init( 8*4096*4096, 3, c_s[5], c_s3[5],&STA6[0][0], c_s4[5],0xf0,1,&st2_p1[0]);
    
    cmCR[0].Init(    8*32*4096, 1, c_s[10],c_s3[10],&STA7[0][0],c_s4[10],0x00,1);
    cmCR[1].Init(      32*4096, 1, c_s[10],c_s3[10],&STA2[0][0],c_s4[10],0x00,1);
    cmCR[2].Init(   16*32*4096, 1, c_s[10],c_s3[10],&STA6[0][0],c_s4[10],0x00,0);
   
    brcxt.Init(&brackets[0],8);
    qocxt.Init(&quotes[0],4,true);
    fccxt.Init(&fchar[0],20);
    colcxt.Init();
    worcxt.Init();
    worcxt1.Init();
    worcxt2.Init();
    worcxt3.Init();
    worcxt0.Init();
    sencxt.Init();
    sencxtL.Init();
    sencxtT.Init();
    sencxtCL.Init();
    htcxt.Init(&html[0],2,false,0xfff);
    simiwor=sencxt.SimilarSentence(&worcxt,0);
    cWord=&StemWords[0], pWord=&StemWords[3];
    xml.Init();
}
  
void PredictorFree() {
    smA[0].Free(); 
    smA[1].Free(); 
    smA[2].Free(); 
    for (int i=0;i<11;i++) mxA[i].Free();
    for (int i=0;i<6;i++) cmC[i].Free();
    for (int i=0;i<8;i++) cmC4[i].Free();
    for (int i=0;i<18;i++) cmC2[i].Free();
    for (int i=0;i<44515;i++) free(dictW[i]);
    rcmA[0].Free();
}
  
int buf(int i){
    return buffer[(pos-i)&BMASK];
}
int bufr(int i){
    return buffer[i&BMASK];
}

// Match model 2
// based on paq8px v208
struct HashElementForMatchPositions { // sizeof(HashElementForMatchPositions) = 3*4 = 12
  #define mHashN   4
  U32 matchPositions[mHashN];
  void Add(int pos) {
    if (mHashN > 1) {
      memmove(&matchPositions[1], &matchPositions[0], (mHashN - 1) * sizeof(matchPositions[0]));
    }
    matchPositions[0] = pos;
  }
};

const int MINLEN_RM = 3; //minimum length in recovery mode before we "fully recover"
const int LEN1 = 5;      // order x
const int LEN2 = 7;      //
const int LEN3 = 9;
 
struct MatchInfo {
    U32 length;      // rebased length of match (length=1 represents the smallest accepted match length), or 0 if no match
    U32 index;       // points to next byte of match in buf, 0 when there is no match
    U32 lengthBak;   // allows match recovery after a 1-byte mismatch
    U32 indexBak;
    U8 expectedByte; // prediction is based on this byte (buf[index]), valid only when length>0
    bool delta;      // indicates that a match has just failed (delta mode)
    void Init() {
        length=0;
        index=0;
        lengthBak=0;
        indexBak=0;
        expectedByte=0;
        delta=false;
    }
    bool isInNoMatchMode() const {
      return length == 0 && !delta && lengthBak == 0;
    }

    bool isInPreRecoveryMode() const {
      return length == 0 && !delta && lengthBak != 0;
    }

    bool isInRecoveryMode() const {
      return length != 0 && lengthBak != 0;
    }

    U32 recoveryModePos() const {
      assert(isInRecoveryMode()); //must be in recovery mode
      return length - lengthBak;
    }

    U32 prio() {
      return
        (length != 0) << 31 |                          // normal mode (match)
        (delta) << 30 |                                // delta mode
        (delta ? (lengthBak>>1) : (length>>1)) << 24 | // the longer wins, halve
        (index&0x00ffffff);                            // the more recent wins
    }
    bool isBetterThan(MatchInfo* other) {
      return this->prio() > other->prio();
    }

    void update() {
      //printf("- pos %d %d  index %d  length %d  lengthBak %d  delta %d\n", x.blpos, x.bpos, index, length, lengthBak, delta ? 1 : 0);
      if (length != 0) {
        const int expectedBit = (expectedByte >> ((8 - x.bpos) & 7)) & 1;
        if (x.y != expectedBit) {
          if (isInRecoveryMode()) { // another mismatch in recovery mode -> give up
            lengthBak = 0;
            indexBak = 0;
          } else { //backup match information: maybe we can recover it just after this mismatch
            lengthBak = length;
            indexBak = index;
            delta = true; //enter into delta mode - for the remaining bits in this byte length will be 0; we will exit delta mode and enter into recovery mode on bpos==0
          }
          length = 0;
        }
      }

      if (x.bpos == 0) {
        // recover match after a 1-byte mismatch
        if (isInPreRecoveryMode()) { // just exited delta mode, so we have a backup
          //the match failed 2 bytes ago, we must increase indexBak by 2:
          indexBak++;
          if (lengthBak < MAXLEN) {
            lengthBak++;
          }
          if (bufr(indexBak) == c1) { //                     match continues -> recover 
            length = lengthBak;
            index = indexBak;
          } else { // still mismatch
            lengthBak = indexBak = 0; // purge backup (give up)
          }
        }
        // extend current match
        if (length != 0) {
          index++;
          if (length < MAXLEN) {
            length++;
          }
          if (isInRecoveryMode() && recoveryModePos() >= MINLEN_RM) { // recovery seems to be successful and stable -> exit recovery mode
            lengthBak = indexBak = 0; // purge backup
          }
        }
        delta = false;
      }
        //printf("  pos %d %d  index %d  length %d  lengthBak %d  delta %d\n", x.blpos, x.bpos, index, length, lengthBak, delta ? 1 : 0);
    }

    void registerMatch(const U32 pos, const U32 LEN) {
      assert(pos != 0);
      length = LEN - LEN1 + 1; // rebase
      index = pos;
      lengthBak = indexBak = 0;
      expectedByte = 0;
      delta = false;
    }
};

const int matchN=4; // maximum number of match candidates
MatchInfo matchCandidates[matchN];
U32 numberOfActiveCandidates=0;
HashElementForMatchPositions *mhashtable,*mhptr;
U32 mhashtablemask;
const int nST=3;
U32 ctx[nST];

bool isMMatch(const U32 pos, const int MINLEN) {
    for (int length = 1; length <= MINLEN; length++) {
      if (buf(length) != bufr(pos - length))
        return false;
    }
    return true;
}

void AddCandidates(HashElementForMatchPositions* matches, U32 LEN) {
    U32 i = 0;
    while (numberOfActiveCandidates < matchN && i < mHashN) {
      U32 matchpos = matches->matchPositions[i];
      if (matchpos == 0)
        break;
      if (isMMatch(matchpos, LEN)) {
        bool isSame = false;
        //is this position already registered?
        for (U32 j = 0; j < numberOfActiveCandidates; j++) {
          MatchInfo* oldcandidate = &matchCandidates[j];
          isSame = (oldcandidate->index == matchpos);
          if (isSame)
            break;
        }
        if (!isSame) { //don't register an already registered sequence
          matchCandidates[numberOfActiveCandidates].registerMatch(matchpos, LEN);
          numberOfActiveCandidates++;
        }
      }
      i++;
    }
  }
  
void MatchModel2update() {
  //update active candidates, remove dead candidates
  U32 n = max(numberOfActiveCandidates, 1);
  for (U32 i = 0; i < n; i++) {
    MatchInfo* matchInfo = &matchCandidates[i];
    matchInfo->update();
    if (numberOfActiveCandidates != 0 && matchInfo->isInNoMatchMode()) {
      numberOfActiveCandidates--;
      if (numberOfActiveCandidates == i)
        break;
      memmove(&matchCandidates[i], &matchCandidates[i + 1], (numberOfActiveCandidates - i) * sizeof(MatchInfo));
      i--;
    }
  }

  if( x.bpos == 0 ) {    
    U32 hash;
    HashElementForMatchPositions* matches;

    hash = t[LEN3];
    matches = &mhashtable[(hash& mhashtablemask)];
    if (numberOfActiveCandidates < matchN)
      AddCandidates(matches, LEN3); //longest
    matches->Add(pos);

    hash = t[LEN2];
    matches = &mhashtable[(hash& mhashtablemask)];
    if (numberOfActiveCandidates < matchN)
      AddCandidates(matches, LEN2); //middle
    matches->Add(pos);

    hash = t[LEN1];
    matches = &mhashtable[(hash& mhashtablemask)];
    if (numberOfActiveCandidates < matchN)
      AddCandidates(matches, LEN1); //shortest
    matches->Add(pos);
    
    hash =worcxt.Word(1);
    matches = &mhashtable[(hash& mhashtablemask)];
    if (numberOfActiveCandidates < matchN)
      AddCandidates(matches, LEN1); //shortest
    matches->Add(pos);

    for (U32 i = 0; i < numberOfActiveCandidates; i++) {
      matchCandidates[i].expectedByte = bufr(matchCandidates[i].index);
    }
  }
}

int MatchModel2mix() {
  MatchModel2update();

  for( int i = 0; i < nST; i++ ) { // reset contexts
    ctx[i] = 0;
  }
  
  int bestCandidateIdx = 0; //default item is the first candidate, let's see if any other candidate is better
  for (U32 i = 1; i < numberOfActiveCandidates; i++) {
    if (matchCandidates[i].isBetterThan(&matchCandidates[bestCandidateIdx]))
      bestCandidateIdx = i;
  }

  const U32 length = matchCandidates[bestCandidateIdx].length;
  const U8 expectedByte = matchCandidates[bestCandidateIdx].expectedByte;
  const bool isInDeltaMode = matchCandidates[bestCandidateIdx].delta;
  const int expectedBit = length != 0 ? (expectedByte >> (7 - x.bpos)) & 1 : 0;

  U32 denselength = 0; // 0..27
  if (length != 0) {
    if (length <= 16) {
      denselength = length - 1; // 0..15
    } else {
      denselength = 12 + ((length ) >> 2); // 16..27
    }
    ctx[0] = (denselength << 4) | (expectedBit << 3) | x.bpos; // 1..28*2*8
    ctx[1] = ((expectedByte << 11) | (x.bpos << 8) | c1) ;//+ 1;
    const int sign = 2 * expectedBit - 1;
    x.mxInputs1.add(sign * (length << 5));
  } else { // no match at all or delta mode
    x.mxInputs1.add(0);
  }

  if( isInDeltaMode ) { // delta mode: helps predicting the remaining bits of a character when a mismatch occurs
    ctx[2] = (expectedByte << 8) | x.c0;
  }

  for( int i = 0; i < nST; i++ ) {
    const U32 c = ctx[i];
    if( c != 0 ) {
         smA[i].set(c);
      const int p1 = smA[i].pr;
      const int st = stretch(p1);
      x.mxInputs1.add(st >> 2);
      x.mxInputs1.add((p1 - 2048) >> 3);
    } else {
      x.mxInputs1.add(0);
      x.mxInputs1.add(0);
    }
  }
  return length;
}

int buffer1(int i){
    return cwbuf[(cwpos-i)&CBMASK];
}


// Return value based on a word type. Must fit in 4 bits
int __attribute__ ((noinline)) getWT(U32 t) {
    if      (t& Verb) return 1;
    else if (t& Noun) return 2;
    else if (t& Adjective) return 3;
    else if (t& Male) return 4;
    else if (t& Female) return 5;
    else if (t& Article) return 6;
    else if (t& Conjunction) return 7;
    else if (t& Adposition) return 8;
    else if (t& ConjunctiveAdverb) return 9;
    else if (t& AdverbOfManner) return 11;
    else if (t& Suffix) return 12;
    else if (t& Prefix) return 13;
    else if (t& Plural) return 10;
    else if (t& Pronoun) return 2;
    else if (t) return 14;
    else return 15;
}

int __attribute__ ((noinline)) getWT3(U32 t) {
    if      (t&Verb) return 1;
    else if (t&Noun) return 2;
    else if (t&Adjective) return 3;
    else if (t&Plural) return 4;
    else if (t&PastTense) return 5;
    else if (t&PresentParticiple) return 6;
    else if (t&AdjectiveSuperlative) return 7;
    else if (t&AdjectiveWithout) return 8;
    else if (t&AdjectiveFull) return 9;
    else if (t&AdverbOfManner) return 10;
    else if (t&Suffix) return 11;
    else if (t&Prefix) return 12;
    else if (t&Male) return 13;
    else if (t&Female) return 14;
    else if (t&Article) return 15;
    else if (t&Conjunction) return 16;
    else if (t&Adposition) return 17;
    else if (t&Number) return 18;
    else if (t&Preposition) return 19;
    else if (t&ConjunctiveAdverb) return 20;
    else if (t&Pronoun) return 21;
    else return 0;
}

// Update string and stemm when string ends
void setbufstem(char c) {
    // Allow:
    // single ' at the beginning of word
    // - in middle or at end of word
    if ((c>='a' && c<='z') || (c==APOSTROPHE &&  c2!=APOSTROPHE) || (c=='-'&&  (*cWord).Length()>0) ) {
        (*cWord)+=c;
    }
    // Ignore wiki link and maintain current word: [dog]s
    // Disabled when in http link
    else if ((*cWord).Length()>0 && c==SQUARECLOSE && fccxt.cxt!=HTLINK && isParagraph) {
        inBR=true;
    } 
    else if ((*cWord).Length()>0) {
        bool res=StemmerEN.Stem(cWord);
        StemIndex=(StemIndex+1)&3;
        pWord=cWord;
        cWord=&StemWords[StemIndex];
        memset(cWord, 0, sizeof(Word));
        
        if ((*pWord).Type& Verb) sVerb=(*pWord).Hash;

        // This is not correct way. We assume that all words that have word 'the' before are nouns. This is not true.
        if (lastArt) {
            (*pWord).Type|=Noun;
        }
        if ( ((*pWord).Type==Article&&  buffer1(5)==SPACE&&  buffer1(4)=='t' &&  buffer1(3)=='h'&&buffer1(2)=='e')) {
            lastArt=true;
        } else  
            lastArt=false;

        U32 whash=(isMath?word0:(*pWord).Hash);
        lastWT=lastWT*16+getWT((*pWord).Type);

        if (isHTTAG) ;
        else if (worcxt.tpbyte==LESSTHAN && isHTTAG==false && qocxt.cxt==0) {
            isHTTAG=wasTag=true;
        }
        if (isHTTAG) {
            worcxt3.Update(word0,c1,0,whash,lastCW);
        }
        
        // Sentence, all words.
        worcxt.Update(word0,c1,(*pWord).Type,whash,lastCW,(brcxt.cxt==SQUAREOPEN || c1==SQUARECLOSE|| c2==SQUARECLOSE|| inBR)?1:0);
        inBR=false;
        
        // Paragraph, most words, exclude Conjunction etc.
        if (((*pWord).Type&(Conjunction+Article+Male+Female+Number+ConjunctiveAdverb))==0  && isHTTAG==false) worcxt1.Update(word0,c1,(*pWord).Type,whash);
        // Stream, words with type, exclude Conjunction etc.
        if (((*pWord).Type&(Conjunction+Article+Male+Female+Adposition+Number+AdverbOfManner+ConjunctiveAdverb ))==0  && (isHTTAG==false)) {  
            if ((*pWord).Type)  worcxt2.Update(word0,c1,(*pWord).Type,whash,lastCW);
        }
    }
    
    if (c==LF && wasTag) { 
        wasTag=isHTTAG=false;
    }
    if (c=='>') isHTTAG=false;
}

// Set decoded char to buffer and update stemmer string
void setbuf(char c){
    cwbuf[cwpos&CBMASK]=c;
    cwpos++; 
    setbufstem(c);
}

// Process partial or full codeword
// decode to text if found and update buffer
void __attribute__ ((noinline)) procWord() {
    if (dcwl>0) {
        if (dcwl==2) dcw=(dcw/256)+(dcw&255)*256;
        if (dcwl==3) dcw=((dcw/256)/256)+(dcw&0xff00)+(dcw&255)*256*256;
        decodeWord(dcw);
        cwSTR=lastCW;
        
        dcw=dcwl=0;
        const char *so1=&(*dictW[lastCW]);
        for (int i=0; i<dictWLen[lastCW]; ++i) {
            setbuf(so1[i]); 
        }
    }
}

void __attribute__ ((noinline)) updateSen(int i=0) {
    if (isParagraph) worcxt.paragraph=true;
    if (colcxt.lastfc(1)==SQUAREOPEN && c2==SQUARECLOSE) sencxt.Update(&worcxt),sencxtCL.Update(&worcxt);
    else if (colcxt.nlChar==WIKITABLE) sencxtT.Update(&worcxt);
    else if (colcxt.lastfc(1)!='*')sencxt.Update(&worcxt);
    else if (colcxt.lastfc(1)=='*' || colcxt.lastfc(1)=='#')sencxtL.Update(&worcxt);
    worcxt.Reset();
    wt3cxt=0;

    if (i==0 && colcxt.isTemp==false) worcxt1.Reset();
}

// Parse all contexts
void parseByte() {
    U32 j=0;
    const int c4=x.c4;
        // Skip mode based on a match lenght
        skipM1=isMatch>61;
        if (colcxt.isTemp || colcxt.nlChar==WIKITABLE) skipM1=isMatch>15;
        if (colcxt.lastfc()==SQUAREOPEN || worcxt.wordcount>1) skipM1=isMatch>15;
        if (colcxt.lastfc()=='*' ||colcxt.lastfc()=='#') skipM1=isMatch>15;
        if (colcxt.lastfc()==EQUALS) skipM1=isMatch>5;
        if (qocxt.context==APOSTROPHE) skipM1=isMatch>5;
        if (isParagraph && worcxt.wordcount>3) skipM1=isMatch>13;
        c3=c2;
        c2=c1;
        c1=c4&0xff;
        n2bState=wrt_2b[c1];
        n3bState=wrt_3b[c1];
        n4bState=wrt_4b[c1];

        stream2b=stream2b*4+n2bState;
        stream4b=stream4b*16+n4bState;
        buffer[pos&BMASK]=c1;
        pos++;
        // When 'text' tag ends force LF reset when line starts with text
        if ( c2==GREATERTHAN  && isText==true ) {isPageStarted=true;
            isText=false;
            PState=PText;
            PStateH=hash(PStateH,PState,0);
            if (c1!=LF){
                colcxt.Update(LF,0);
                updateSen();
                fc=isParagraph=firstWord=0;
                nl1=nl;
                nl=pos-2;
            }
        }
        if ((x.c4&0xffffff)==(((EQUALS*256)+EQUALS)*256+EQUALS)) isLongTOP=true;
        if (c1==CURLYCLOSE && c2==VERTICALBAR && colcxt.nlChar==WIKITABLE) {
            updateSen();
        }
        if (worcxt.wordcount<6 &&colcxt.isTemp && c1==CURLYCLOSE)worcxt.removeWordsL(8,CURLYOPENING,CURLYCLOSE);
        // Column context update
        colcxt.Update(c1,c4&0xffffff);
        // Bracket context update
        if (c1<'a' || c1==SQUAREOPEN) brcxt.Update(c1);            // advance bracket context only if no letters, so we do not get out of range
        if (c1==SPACE && c2==LESSTHAN) brcxt.Update(GREATERTHAN); // Probably math operator, ignore
        
        // Quote context update
        qocxt.Update(c1); 
        
        if (htcxt.cxt && c2=='L'&&(c1==SPACE || c1=='!' || c1<128 )) {
            htcxt.Update('&'*256+'N');// not an html tag
        } 
        htcxt.Update(c4&0xffff);

        // Look for byte stream x4 and order X context end.
        if (c1=='$' || c1==SQUARECLOSE || c1==VERTICALBAR || c1==')' || c1==SQUAREOPEN) {
            if (c1!=c2) {
                // update order X context hashes
                for (int i=13; i>0; --i)  
                    t[i]=t[i-1]*primes[i];
            }
            // Duplicate input byte, marks end
            x4=(x4<<8)+c2;
            stream2b=stream2b*4+n2bState;
             // Repeat, 2k, has problems with list, probably
            stream2bR=(stream2bR<<2)+n2bState;
            stream3bR=(stream3bR<<3)+n3bState;
        }
        // Update byte stream x4 and order X contexts
        x4=(x4<<8)+c1; 

        for (int i=13; i>0; --i)
            t[i]=t[i-1]*primes[i]+c1+i*256;
        
        words=words<<1;
        spaces=spaces<<1;
        numbers=numbers<<1;
        j=c1;
        if (((j-'a') <= ('z'-'a')) || (c1>127 && c2!=ESCAPE)) {
            if (word0==0 ) {
                if (isMath && c2=='/' && c3==LESSTHAN) isMath=false;
                U8 reChar=c2;
                // Skip first upper or upper word flag
                if (c2==FIRSTUPPER || c2==UPPER) {
                if (c3!=APOSTROPHE) reChar=c3;
                else if (buf(4)!=APOSTROPHE) reChar=buf(4);
                else if (buf(5)!=APOSTROPHE) reChar=buf(5);
                else if (buf(6)!=APOSTROPHE) reChar=buf(6);
                
                else reChar=c3;
            }
                else if (c2=='/' && c3==LESSTHAN) reChar=c3;  // </ to <
                 worcxt.Set(reChar,c2==FIRSTUPPER?1:0,c2==LESSTHAN?c2:((c3==LESSTHAN && c2=='/')?c3:0)); 
                 worcxt0.Set(reChar,c2==FIRSTUPPER?1:0); 
                 worcxt1.Set(reChar);
            }

            words=words|1;
            word0=word0*2104+j; //263*8
            word00=word0;
            h=word0*271;u8w=0;
            if (brcxt.cxt==SQUAREOPEN && fccxt.cxt!=HTLINK && fc!=HTML) linkword=linkword*2104+j;
            if (isParagraph && fccxt.cxt!=HTLINK && colcxt.isTemp==false) senword=senword*2104+j;
            // If ' or 0x27 (") is used for quotes, sometimes it has other meaning
            // remove quote context if any of the fallowing is true
            const int word3bit=(words&7);
            if (((word3bit==5) && (c2==APOSTROPHE)) ||                    // "x'x" where x is any letter in word
             ((word3bit==1) && (c3==SQUARECLOSE) && (c2==APOSTROPHE)) ||  // "]'x" where x is any letter in word
             ((word3bit==1) && (numbers&4)&&(c2==APOSTROPHE))             // "y'x" where y is number and x is any letter in word //(c3>='0' && c3<='9') 
             ) qocxt.Update(qocxt.cxt); 
            // Attempt to partialy decode word dictionary index 
            if (c1>127 && dcwl<3) {
                dcw=dcw*256+c1,dcwl++;
                if (x.blpos>6){
                    int dcw2=0; //dcw; // ignore first byte
                    if (dcwl==2)dcw2=(dcw/256)+(dcw&255)*256;
                    else if (dcwl==3)dcw2=((dcw/256)/256)+(dcw&0xff00)+(dcw&255)*256*256;
                    int i=decodeCodeWord(dcw2);
                    if (i>0) deccode=i; // Change only when found
                    assert(deccode>=0 && deccode<0x20000);
                }
            } else if (dcw) {
                // Codeword ended so decode it and put to the decoded text buffer
                procWord();
                // Set mixer(8) to last decoded dictionary index (befere main text ends)
                if (x.blpos<448131719)  deccode=lastCW; // for processed enwik9
                assert(deccode>=0 && deccode<0x20000);
            }
            // Set char to decoded text buffer
            if (c1==10 || c1==9 ||(c1>31 && c1<128)) {
                setbuf(charSwap(c1));
            } 
        } else if (((c1==ESCAPE && c2!=ESCAPE) || (c1>127 && c2==ESCAPE)) && (words&4)==4) {
            if (c1!=ESCAPE) {
                words=words|3;
                word0=word0*2104+j; //263*8
                word00=word0;
                h=word0*271;
            }
        } else {
            if (word0) {
                // Decode codeword if any and put to the decoded text buffer
                procWord();
                // Set mixer(8) to last decoded dictionary index (befere main text ends)
                if (x.blpos<448131719) deccode=lastCW; // for processed enwik9
                assert(deccode>=0 && deccode<0x20000);
            } else {
                // Set mixer(8) to non-codeword flag (bit 16) and last 8 bit2words
                deccode=0x10000+(stream2b&0xffff);
                assert(deccode>=0 && deccode<0x20000);
                lastCW=0;
            }
            // Set char to decoded text buffer
            if (c1==10 || c1==9 || (c1>31 && c1<128)) {
                setbuf(charSwap(c1));
            } 
            
            // Parse numbers: (number), (number.number) or (number,number) 
            if (c1>='0' && c1<='9') {
                numbers=numbers+1;np[numberA&0xffff]=pos;
                numberA=numberA*2104+j;
                if (mybenum && numlen1<=2) number0=number1,number1=0,numlen0=numlen1,numlen1=0;
                number0=number0*10+(c1&0x0f);
                numlen0=min(19,numlen0+1);mybenum=0;
            } else {
                if (numberA) {
                    worcxt0.Update(numberA,c1,Number,numberA) ,worcxt.Update(numberA,c1,Number,numberA,0,(brcxt.cxt==SQUAREOPEN||c1==SQUARECLOSE )?1:0);
                }
                numberA=0;
                if (numlen0 ||((numbers&0xf)==0)) {
                    number1=number0,numlen1=numlen0,number0=numlen0=0;
                }
                if (numlen1<=2 &&numlen1&&((numbers&5)==5) && numlen0==0 && c2=='.') mybenum=2;
                else if (numlen1<=2&&numlen1&&(numbers&2) && numlen0==0 && c1=='.') mybenum=1;
                else if (mybenum==1  && c1!='.') mybenum=0;
            }
            // If ' or 0x27 (") is used for quotes, sometimes it has other meaning
            // remove quote context if any of the fallowing is true
            const int word3bit=(words&7);
            if ( ((word3bit==4)&& (c1==SPACE) && (c2==APOSTROPHE)) ||     // "x' " where x is any letter in word
             ((c1==FIRSTUPPER) && (numbers&4) && (c2==APOSTROPHE)) ||     // "x'@" where x is number               // this is somehow semi good-bad //(c3>='0' && c3<='9')
             ((word3bit==4) && (c1==FIRSTUPPER) && (c2==APOSTROPHE)) ||   // "x'@" where x is any letter in word   // this is somehow semi good-bad
             ((word3bit==4) && (numbers&1)&&(c2==APOSTROPHE))             // "x'y" where y is number and x is any letter in word //(c3>='0' && c3<='9') 
              ) qocxt.Update(qocxt.cxt); 
            // Reset only when not '[word word ...'
            if (word00 && !(fccxt.cxt==SQUAREOPEN) ) word00=0;
            // Update word context
            if (word0){
                if (wt3cxt==0 || ((worcxt.Type(1)&Article)==Article)) wt3cxtW1=wt3cxtW, wt3cxtW=0;
                // Full sentance type hash
                wt3cxt=hash(wt3cxt,getWT3(worcxt.Type(1)),worcxt.wordcount);
                //Partial sentance codeword with type
                if ((worcxt.Type(1))!=0) wt4cxtW=hash(wt4cxtW,worcxt.Code(1),worcxt.Type(1));
                // Full sentance  codeword with type excluding Nouns
                if ((worcxt.Type(1)&Noun)==0 && (worcxt.Type(1))!=0) wt3cxtW=hash(wt3cxtW,worcxt.Code(1),worcxt.Type(1));
                // Skip some word tpyes in main word order
                if (((*pWord).Type&(ConjunctiveAdverb+Conjunction))==0) {
                    worcxt0.Update(word0,c1!=LF?c1:0,0,word0); 
                }
                
                if (worcxt.Type(1)==Number) {
                    stream3bR=(stream3bR<<7)+1;
                    stream3b=(stream3b<<7)+1;
                }
                if (firstWord==0 && fccxt.cxt!=SQUAREOPEN) {
                    firstWord=word0;
                }
                // Separate at Conjunction(ok)
                if (worcxt.Type()&Conjunction) {
                    stream3bR=stream3bR<<7;
                    stream3b=stream3b<<7;
                    if (isParagraph)senword=0;
                }
                
                if (worcxt.Type()&Article) {
                    stream3bR=(stream3bR<<7)+2;
                    stream3b=(stream3b<<7)+2;
                }
             
                // Add ok, PP problematic
                if (worcxt.Type()&Adposition || (isParagraph && (worcxt.Type()&PresentParticiple)) ) {
                    stream2bR=(stream2bR<<2)+(stream2bR&3);
                    stream2b=(stream2b<<2)+(stream2b&3);
                }
                // needs work
                if ((worcxt.Type()&AdverbOfManner)) {
                    if (isParagraph) {
                        worcxt.Remove();
                    }
                }
                // Article(good)
                if ((worcxt.Type()&Noun) && worcxt.Type(2)&Article) {
                    stream3bR=(stream3bR<<6)+1;
                    stream3b=(stream3b<<6)+1;
                    U16 sb=worcxt.sBytes(1);
                    U32 w=worcxt.Word(1);
                    U32 t=worcxt.Type(1);
                    U8 ca=worcxt.Capital(1);
                    U32 co=worcxt.Code(1);
                    worcxt.Remove();
                    worcxt.Remove();
                    worcxt.Set(sb>>8,ca);
                    worcxt.Update(w,c1,t,w,co,2);
                }
                // Reset all bit stream mask after a word
                stream3bRMask2=stream3bRMask1;
                stream3bMask1=stream3bMask;
                stream3bMask=stream2bMask=stream3bRMask1=0;
            }

            if (buffer1(6)==charSwap(LESSTHAN) && buffer1(5)== 't'&& isText==false && c1==SPACE  &&cwSTR==cwTEXT) isText=true,cwSTR=0x10000;
            
            if (buffer1(8)==charSwap(LESSTHAN) && isNowiki==false && cwSTR==cwNOWIKI) isNowiki=true;
            else if ((buffer1(9)=='/') && (c1==GREATERTHAN ) && isNowiki==true && cwSTR==cwNOWIKI) isNowiki=isPre=false,cwSTR=0x10000;
            
            U8 pChar=(worcxt.sBytes()>>8)&0xff;
            if (isMath && ( (c1==SPACE && colcxt.lastfc()!=COLON)|| c1==',') && c2==GREATERTHAN && cwSTR==cwMATH) isMath=false,cwSTR=0x10000;
            if (isMath && c1=='/' && c2==LESSTHAN && c3==GREATERTHAN && buffer1(4)=='h')isMath=false,cwSTR=0x10000;
                
            if (isNowiki==false && buffer1(6)==charSwap(LESSTHAN) && buffer1(5)=='m'&& isMath==false && c1!='.'&& buffer1(7)!='&'  && buffer1(8)!='&' && cwSTR==cwMATH) isMath=true;
            else if ((buffer1(6)=='/') && (c1==GREATERTHAN ||c1=='&') && isMath==true && cwSTR==cwMATH) isMath=false,cwSTR=0x10000;

            pChar=(worcxt.sBytes()>>8)&0xff;
            if ((buffer1(5)==charSwap(LESSTHAN)) && (c1==GREATERTHAN) && (buffer1(4)=='p') && isPre==false && cwSTR==cwPRE) isPre=true,cwSTR=0x10000;
            else if ((buffer1(5)=='/') && (c1==GREATERTHAN) &&(buffer1(4)=='p') && cwSTR==cwPRE) isPre=false,cwSTR=0x10000;
            // wikipedia article ended with page tag. Do reset on some contexts. 
            if ((buffer1(6)=='/') && (c1==GREATERTHAN) &&(buffer1(5)=='p') && cwSTR==cwPAGE) {
                isPre=isMath=isNowiki=false,colcxt.nlChar=LF;
                colcxt.resetCells();
                fccxt.Reset();
                brcxt.Reset();
                qocxt.Reset();
                htcxt.Reset();
                worcxt.Reset();wt3cxt=0;
                worcxt0.Reset();
                worcxt1.Reset();
                worcxt2.Reset();
                skipSeeExternal=isCategory=false;
                isPageStarted=false;
                PState=PStateH=0;
                 
                lastPTOP=lastPTOP*2;
                if (pageParag<2 && pageSent<5) lastPTOP++;
                
                cmC4[1].reset();
                cmC4[2].reset();
               
                cmC4[4].reset();
                cmC4[7].reset();
                cmC4[8].reset();
                cmC4[6].reset();
                cmC[0].reset();
                if ((lastPTOP&63)==63) cmC[1].reset(),lastPTOP++;
                pageParag=pageSent=0;
            }

            // Update word0 pos
            wp[word0&0xffff]=pos;

            word0=h=0; 
            if (linkword && c1==COLON) linkword=0;
            if (c1=='-'&& c2==SPACE) worcxt1.Reset(),sVerb=0;
            if (c1=='-'&& worcxt.wordcount==1 && worcxt.Type()==Number && colcxt.nlChar!=WIKITABLE) worcxt.Reset(),wt3cxt=0;

            // Paragraph or sentence related updates
            if (c1==SPACE) {
                spaces++;
            }
            else if (c1==LF) {
                if (wt4cxtW) wt4cxtW1=wt4cxtW,wt4cxtW=0;
                nl1=nl;
                nl=pos-1;
                stream3bR=(stream3bR<<7);
                stream2b=stream2b|0x3fc;
                words=0xfc;
                updateSen();
                stream2bR=stream2bR<<2;
                stream4b=stream4b|0xfff0; 
                if (colcxt.lastfc(1)==FIRSTUPPER) pageParag++;
                pageSent=pageSent+isParagraph;
                fc=isParagraph=firstWord=lastWT=0;
                if (c2==LF) isNowiki=false;
                if (!(colcxt.isTemp || colcxt.nlChar==WIKITABLE)) {
                    brcxt.Reset();
                    qocxt.Reset(); 
                }
                wt3cxt=0;
                wasVerb=wasNoun=nestList=false;
                wasVerbH=wasNounH=0;
            } else if (c1=='.' || c1==')' || c1==QUESTION) {
                lastWT=lastWT*16;
                stream3bR=stream3bR<<7;
                stream3b=stream3b<<7;
                words=words|0xfe;
                x5=(x5<<8)+(c4&0xff);
                stream2b=stream2b|204;
                stream4b=((stream4b&0xffff0)<<8)+(stream4b&0xf);
                stream2bR= stream2bR&0xffffffc0;
                if (c1=='.') {
                    if (fccxt.cxt!=SQUAREOPEN) wshift=1;
                    // Skip sentence reset if fallowing is true
                    if (!(fccxt.cxt==SQUAREOPEN || fccxt.cxt=='(' || colcxt.nlChar==WIKITABLE || colcxt.lastfc()=='*')) {
                        if (isParagraph) worcxt.paragraph=true;
                        sencxt.Update(&worcxt);
                        wasVerb=wasNoun=false;
                        wasVerbH=wasNounH=0;
                        worcxt.Reset();
                        wt3cxt=0;
                        if (wt4cxtW) wt4cxtW1=wt4cxtW,wt4cxtW=0;
                    }
                    senword=0; // Age words(stream) and reset word context
                }
                if ( c1==')') senword=0;
            }
            else if (c1==',') {
                wasVerb=wasNoun=false;
                wasVerbH=wasNounH=0;
                if (wt4cxtW) wt4cxtW1=wt4cxtW,wt4cxtW=0;
                words=words|0xfc;
                senword=0;
            }
            else if (c1=='(') {
                senword=0;
            }
            else if (c1==SEMICOLON && colcxt.nlChar!=WIKITABLE) {
                updateSen(1);
            }
            // Probably link, word list - this can probably be better
            else if (c1==COLON) {
                stream3b=(stream3b&0xfffffff8)+4;
                stream2b=stream2b|12;//      1100
                x5=(x5<<8)+(c4&0xff);
                 senword=0;
            }
            // Table or template - this can probably be better
            else if (c1==CURLYCLOSE || c1==CURLYOPENING) {
                words=words|0xfc;
                stream3bR=stream3bR&0xffffffc0;
                x5=(x5<<8)+(c4&0xff);
                stream3b= (stream3b&0xfffffff8)+3;
            }
            // Wiki link ended 
            else if (c1==SQUARECLOSE) {
                stream3b=(stream3b&0xfffffff8)+3;
                linkword=0;
            }
            // HTML - WIT
            else if (c1==LESSTHAN || c2=='&') {
                words=words|0xfc;
            }

            // List to paragraph
            else if ((c1=='-' && (colcxt.lastfc()=='*')) && brcxt.cxt!=SQUAREOPEN && isParagraph==0) {
                isParagraph=1;
                fc=FIRSTUPPER;
            }
            // Probably heading
            else if (c1==EQUALS) {
                stream3b=(stream3b&0xfffffff8)+4;
                c2='.'; // ok
                words=words*2;
            }
            // HTML - WIT
            if (c1=='!' && c2=='&') {// '&nbsp;' to '&!'  to ' ' WIT
                c1=SPACE;
                stream2b=(stream2b&0xfffffffc)+wrt_2b[SPACE];
                stream3b=(stream3b&0xfffffff8)+wrt_3b[SPACE];
            }
            else if ( (colcxt.lastfc()=='*') && (c1==',' || c1==SPACE) && c2==SQUARECLOSE && isParagraph==0) {
                isParagraph=1;
                fc=FIRSTUPPER;
            }
            if (nestList==false && colcxt.lastfc()=='*' && (c4&0xffff)==0x2a2a) nestList=true;
        }
        x5=(x5<<8)+(c4&0xff);
        // Update bit streams, serial and non-repeating
        // 2 bit stream, switch state if it is new
        if (o2bState!=n2bState) {
            stream2bR=(stream2bR<<2)+n2bState;
            o2bState=n2bState;
        }
        stream2bMask=(stream2bMask<<2)+3;
        
        // 3 bit stream, switch state if it is new
        if (o3bState!=n3bState) {
            stream3bR=(stream3bR<<3)+n3bState;
            stream3bRMask1=(stream3bRMask1<<3)+7;
            stream3bRMask2=(stream3bRMask2<<3)+7;
            o3bState=n3bState;
        }
        stream3b=(stream3b<<3)+n3bState;
        stream3bMask=(stream3bMask<<3)+7;
        stream3bMask1=(stream3bMask1<<3)+7;
        U8 brcontext=brcxt.cxt;
        
        // Map bracket or quote context to 1-7, 0 if not found. (3 bits)
        BrFcIdx=0;
        if (brcxt.context) BrFcIdx=fcy[brcontext];
        if (brcxt.context==0 && qocxt.context) BrFcIdx=fcy[qocxt.context>>8];
        
        // Column and first char 
        col=colcxt.collen();
        above=buffer[(nl1+col)&BMASK];
        above1=buffer[(nl1+col-1)&BMASK];
        // Filtered wiki. We ignore 10 as new line char, '>' marks new line.
        if (colcxt.nlChar==WIKIHEADER) {
           above=colcxt.colb(1,0);
           above1=colcxt.colb(1,1); 
        }
        if (colcxt.isNewLine()) {
                // Reset contexts when there are two empty lines
                if ((colcxt.nlpos(0)+2-colcxt.nlpos(1))<4) { 
                    fccxt.Reset();
                    htcxt.Reset(); 
                    worcxt0.Reset();
                    isLongTOP=false;
                }
                fc=colcxt.lastfc();
                // Reset first char context when there is >. For filtered wiki.
                if (fc==WIKIHEADER) fccxt.Reset();
                // Set paragraph
                if (fc==FIRSTUPPER) isParagraph=1;
                else isParagraph=0;
                if (fc!=SQUAREOPEN) isCategory=false;
                // Set new first char
                fccxt.Update(fc);

                if (fc==EQUALS)
                    PState=PTopic;
                else if (isParagraph)
                    PState=PText;
                else if (colcxt.isTemp)
                    PState=PTemplate;
                PStateH=hash(PStateH,PState,0);
        }

        if (col>2 && c1>FIRSTUPPER && isMath==false) {
            // Before updating first char context look
            if ( (fccxt.cxt==EQUALS || fccxt.cxt==VERTICALBAR || fccxt.cxt==COLON || fccxt.cxt==FIRSTUPPER || fccxt.cxt==CURLYOPENING) && (c1==CURLYCLOSE)) {
                while (fccxt.cxt==EQUALS || fccxt.cxt==VERTICALBAR || fccxt.cxt==COLON || fccxt.cxt==FIRSTUPPER) fccxt.Update(LF);
                fccxt.Update(c1);
            }
            //   If link or template ended then remove any vertical bars |.  [xx|xx] {xx|xx}
            if (fccxt.cxt==VERTICALBAR && (c1==SQUARECLOSE || c1==CURLYCLOSE)) {
                while (fccxt.cxt==VERTICALBAR) fccxt.Update(LF);
            }
            //   If html link ends with ]
            if ((fccxt.cxt==COLON || fccxt.cxt==HTLINK ) && c1==SQUARECLOSE) {
                while (fccxt.cxt==COLON || fccxt.cxt==HTLINK) fccxt.Update(LF);
            }
            if (c1<128)
                fccxt.Update(c1);
        }

        if (c1==COLON && (words&2)==2) cwCOLON=cwSTR;// copy word
        if ((c1==SPACE &&fccxt.cxt==COLON && colcxt.lastfc()!=COLON && colcxt.nlChar!=WIKITABLE) ) {
            // Remove if word before colon was not:
            if (cwCOLON!=cwIMAGE) {
                while (fccxt.cxt==COLON) fccxt.Update(LF);
            }
        }
       
        // wiki links with [catecory:...] or [wikipedia:...], remove word and reset first char colon context
        if (c1==COLON &&(cwCOLON==cwCATEGORY ||cwCOLON==cwUSER ||cwCOLON==cwWIKIPEDIA)) {
            if (cwCOLON==cwCATEGORY) {
                PState=PCategory,isCategory=true;
            }
            PStateH=hash(PStateH,PState,0);
            fccxt.Update(LF);
            worcxt.Remove();
            worcxt0.Remove();
        }

        if (c1==SPACE && c2==LESSTHAN) fccxt.Update(GREATERTHAN); // Probably math operator, ignore
        // Switch from possible category link to http link ( [word:// to [http:// )
        if (fccxt.cxt==COLON && c2=='/' && c1=='/') {
            fccxt.Update(LF);
            fccxt.Update(HTLINK);
        }
        // Wiki link in the beginning of line
        if (colcxt.lastfc(0)==SQUAREOPEN && c1==SPACE&& isParagraph==0) {
            if (c2==SQUARECLOSE || c3==SQUARECLOSE) {
                fc=FIRSTUPPER;
                isParagraph=1; // Paragraph
                // Line started with '[' and last chars were '] ', so probably rest of the line/paragraph fallows. 
                // Reset first char context and set new first char as paragraph start and continue.
                // Problem: inlink links also reset, like images that have description etc.
                fccxt.Reset();
                fccxt.Update(fc);
            }
        }
        // Fist char was space, look for another non-space char
        if (fc==SPACE  && c1!=SPACE) { 
            fc=min(c1,TEXTDATA);
            // Paragraph
            if (fc==FIRSTUPPER) isParagraph=1;
            else isParagraph=0;
            // Set new first char, we keep space from previous update
            fccxt.Update(fc);
        }
        const U8 fccontext =fccxt.cxt;
        if (BrFcIdx==0 &&fccxt.context) BrFcIdx=fcy[fccontext];
        FcIdx=fcq[fccontext];

        // List - needs fixme
        if (fc=='*' && c1!=SPACE) {
            fc=min(c1,TEXTDATA);
        } 
        if (fc=='&' && c1==LESSTHAN) fc=HTML;  // Inwiki html - WIT
        // Full paragraph probably starts after text tag, hint it.
        if (c2==GREATERTHAN && fc==LESSTHAN && c1==APOSTROPHE) fc=APOSTROPHE;
        // We have bold/italic first char/words. When words surrounded by ' end there is problably rest of paraghraph.
        // Reset first char context and set new first char as paragraph start and continue.
        // exeption: ignore in list
        if ((colcxt.lastfc(0)==APOSTROPHE||  (fc==APOSTROPHE &&(colcxt.lastfc(0)!='*')) ) && (c1==SPACE)) {
            if(c2==APOSTROPHE || c3==APOSTROPHE ) {
                fc=FIRSTUPPER;
                isParagraph=1;    // Paragraph
                fccxt.Reset();    // Not really needed
                fccxt.Update(fc); //
            }  
        }
        if ((fc!=FIRSTUPPER) && ((c4&0xffffff)==0x4a2f2f)) {//http link - fixme
            fc=HTLINK;
        }
        // Words surrounded by ()
        worcxt.removeWordsL(8,'(',')');
        worcxt1.removeWordsL(8,'(',')');
        // Words surrounded by [|  as wiki internal link: word [word word|word] word.
        worcxt0.removeWordsL(8,SQUAREOPEN,VERTICALBAR);
        worcxt.removeWordsL(8,SQUAREOPEN,VERTICALBAR);
        worcxt1.removeWordsL(8,SQUAREOPEN,VERTICALBAR);
        worcxt.removeWordsL(8,LESSTHAN,COLON);

        // Template - surrounded by = |
        if (colcxt.isTemp==true) worcxt.removeWordsR(10,EQUALS,VERTICALBAR);
        // Remove words between < >
        worcxt0.removeWordsL(8,LESSTHAN,GREATERTHAN);
        worcxt.removeWordsL(8,LESSTHAN,GREATERTHAN);
        worcxt1.removeWordsL(8,LESSTHAN,GREATERTHAN);
        

        // Indirect
        indirectWord=(c4>>8)&0xffff;
        t2[indirectWord]=(t2[indirectWord]<<8)|c1;
        indirectWord=c4&0xffff;
        indirectWord=indirectWord|(t2[indirectWord]<<16);
        indirectByte=(c4>>8)&0xff;
        t1[indirectByte]=(t1[indirectByte]<<8)|c1;
        indirectByte=c1|(t1[c1]<<8);
        
        t1[brcontext]=(t1[brcontext]<<2)|(stream2b&3); // this is wierd, also end is bad
        indirectBrByte=(stream3b&7)|(t1[brcontext]<<3);
        // 
        indirectWord0Pos=pos-wp[word0&0xffff];
        if (indirectWord0Pos>255) {
            indirectWord0Pos=256+(c1<<16);
        } else {
            indirectWord0Pos=indirectWord0Pos+(buf(indirectWord0Pos)<<8)+(c1<<16);
        }
        indirectNumberd0Pos=pos-np[numberA&0xffff];
        if (indirectNumberd0Pos>1024) {
            indirectNumberd0Pos=0;
        } else {
            indirectNumberd0Pos=(0*256<<16)+(buf(indirectNumberd0Pos)<<8)+(c1<<16);
        }
        if (indirectNumberd0Pos && numberA)indirectWord0Pos=indirectNumberd0Pos;
        //
        ind3[context1_ind3] = (cxtind3 * (1 << 5) + c1) & (0x2000000-1);//++
        context1_ind3 = (context1_ind3 * (1 << 5) + c1) & (0x2000000-1);
        cxtind3 = ind3[context1_ind3];
      
        // utf8
        if (c2==12) {
            if (utf8left==0) {
                if ((c1>>5)==6) utf8left=1,u8w=u8w*191+c1;
                else if ((c1>>4)==0xE) utf8left=2,u8w=u8w*191+c1;
                else if ((c1>>3)==0x1E) utf8left=3,u8w=u8w*191+c1;
                else utf8left=0; //ascii or utf8 error
            } else {
               utf8left--;
               if ((c1>>6)!=2) utf8left=0; 
           }
        }
        // looks like this is not used?
        if (isParagraph) {
            if (wasVerb==false && (worcxt.Type()&Verb)==Verb) wasVerb=true,wasVerbH=worcxt.Word();
            else if (wasNoun==false && (worcxt.Type()&Noun)==Noun) wasNoun=true,wasNounH=worcxt.Word();
        }

        h=h+numberA+c1;
        oldwt1=getWT3(worcxt2.Type(1));
        
        if (skipSeeExternal==false && colcxt.lastfc()==EQUALS) {
            if (worcxt.wordcount==2 && c1==EQUALS) {
                if (worcxt.Code(2)==cwEXTERNAL && worcxt.Code(1)==cwLINKS) {
                    skipSeeExternal=true;
                } else if (worcxt.Code(2)==cwSEE && worcxt.Code(1)==cwALSO) {
                    skipSeeExternal=true;
                }
            } else if (worcxt.wordcount==1 && c1==EQUALS) {
                if (worcxt.Code(1)==cwREFERENCES) {
                    skipSeeExternal=true;
                }else if (worcxt.Code(1)==cwBIBLIOGRAPHY) {
                    skipSeeExternal=true;
                } 
            }
        }
        if (isCategory==false && colcxt.lastfc()==SQUAREOPEN) {
            if (lastCW==cwCATEGORY && c1==COLON){
                isCategory=true;
            }
        }
}

int xmlS=0;
int mp0=0,mp1=0,mp2=0,mp3=0;
U8 mstate=0;

// Set contexts and predict
int modelPrediction() {
    int c=0;
    const int bpos=x.bpos;
    const int c4=x.c4;
    const int c0=x.c0;
    if (bpos==0) parseByte();
    xmlS=xml.p();
    if (bpos==0) {
        // Skip when line starts with spaces and it is space char
        if ((fc==SPACE && c1==SPACE) || (skipM1)) {
            cmC2[0].sets(); cmC2[0].sets(); cmC2[0].sets(); 
        } else {
            for (int i=3; i<6; ++i)
            cmC2[0].set(t[i]);
        }
        cmC2[1].set(t[6]);
        cmC2[2].set(t[8]);
        cmC2[3].set(t[13]);
        
        cmC[5].set((fccxt.context&0xff00)+c1+(stream2b&12)*256+((brcxt.cxt+ brcxt.last())<<24));
        
        if (qocxt.context)
            cmC[4].set((qocxt.context<<8)+c1);
        else
            cmC[4].set((brcxt.context<<8)+c1);
        // Contexts
        // Set run context with word(3), current byte and bit3word(1-5)
        rcmA[0].set(worcxt0.Word0(3)+c1+193 * (stream3b & 0xfff),c1);// needs work
        // Word stream cm(4-5)
        if (col<2 || fc==SPACE) {
            cmC2[4].sets(); 
            cmC2[4].sets();
            cmC2[17].sets();
        } else {
            cmC2[4].set(word00+(number0*191+numlen0)+u8w);// ~4
            // Disabled when:
            // & is line first char (HTML)
            // escaped UTF8
            if (colcxt.lastfc()=='&' || utf8left) 
                cmC2[4].sets();
            else 
                cmC2[4].set(h+worcxt0.Word0(1));
            
            if (brcxt.cxt==LESSTHAN || skipSeeExternal  || colcxt.nlChar==WIKITABLE) cmC2[17].sets(); 
            else cmC2[17].set(worcxt1.Word(1)*53+worcxt1.Word(2)*11+h+(lastWT&0xf));
        }
        if (c1==ESCAPE ||col<2 || utf8left || fc==SPACE) {
            cmC2[5].sets(); 
        } else {
            if (isMatch>61) cmC2[5].sets(); 
            else  cmC2[5].set(h+worcxt0.Word0(2)*71);
        }
        if (fc==SPACE || brcxt.cxt==LESSTHAN) {
            cmC2[5].sets(); cmC2[5].sets(); cmC2[5].sets(); cmC2[5].sets(); cmC2[5].sets();
            cmC2[20].sets(); cmC2[20].sets();cmC2[20].sets();
            cmC4[8].sets(); cmC4[8].sets();
        } else {
            // Last sentence word(4) that is not Adjective with last Adjectiv stream word in a line.
            cmC2[5].set(hash(worcxt.Word(4),worcxt1.Word(1)+h,(stream3b & 511)));
            // Last sentence word(4+) that is not Verb (when found 4+) with last Verb in a line.
            cmC2[5].set(hash(worcxt.Last(4,worcxt.Type(4)^Verb),sVerb+h,(stream3bR & 63)));
            cmC2[5].set(hash(worcxt.fword,worcxt1.Word(1)+h,(stream3b & 63)));
            cmC2[5].set(hash(worcxt2.Word(1),worcxt2.Word(2),word00+c1));
            // Look for last verb in paragraph if found set with word ()
            const U32 lastParVerb=worcxt2.LastIf(1,worcxt.Type(1)&Verb);
            if (lastParVerb) 
                cmC2[5].set(hash(lastParVerb,word00,c1));
            else 
                cmC2[5].sets();
            cmC4[8].set(hash(h,worcxt0.Word0(1)));
            cmC4[8].set(hash(worcxt1.Word(1),worcxt1.Word(2),h)); 

            cmC2[20].set(hash(wt3cxt,word0,(stream3b & 63)));             //hash of current sentance word types with current word and last 2 3bit words
            cmC2[20].set(hash(wt3cxt,worcxt2.Word(1),(stream3b & 511)));  //hash of current sentance word types with last worcxt2 word context and last 3 3bit words
            cmC2[20].set(hash(wt3cxtW,word0,(stream3b & 63)));            //hash of current sentance word codewords (no Nouns) with current word contex
        }
        // Sentence contexts
        WordsContext *lastwor,*lastwor1,*lastwor3;
        if (colcxt.lastfc()=='*') {
            lastwor=sencxt.Sentence(1);    // current word index in Sentence 1
            lastwor1=sencxtL.Sentence(2);  // current word index in Sentence 2
            lastwor3=sencxt.Sentence(1);   // current word index in Sentence 3
        } else if (colcxt.nlChar==WIKITABLE) {// table, FIXME for single multiline
            lastwor=sencxtT.Sentence(2); 
            lastwor1=sencxtT.Sentence(4); 
            lastwor3=sencxtT.Sentence(6); 
        } else {
            lastwor=sencxt.Sentence(1);  // current word index in Sentence 1
            lastwor1=sencxt.Sentence(2); // current word index in Sentence 2
            lastwor3=sencxt.Sentence(3); // current word index in Sentence 3
        }
        U32 lastWordMT=lastwor->LastR(1,worcxt.wordcount, Verb); // Last Verb in Sentence if found or Word(1+)
        U32 simNoun=0,simVerb=0;
        U32 xword=0;
        if (isParagraph==0) {
            if(colcxt.lastfc()=='*') {
                WordsContext *lastworL=sencxtL.Sentence(1);
                xword=lastworL->WordR(worcxt.wordcount);
            }
            else if (colcxt.lastfc()==SQUAREOPEN) {
                WordsContext *lastworCL=sencxtCL.SimilarSentence(&worcxt,worcxt.wordcount);
                xword=lastworCL->WordR(worcxt.wordcount);
            } else
                xword=lastwor1->WordR(worcxt.wordcount);
        }
        else xword=lastwor->LastR(1,worcxt.wordcount, Noun);

        xword1=lastwor3->WordR(worcxt.wordcount);

        cmC2[18].set(word00*1471+ lastwor->WordR(worcxt.wordcount)+lastWordMT +(stream3b & 511)*191+BrFcIdx); //40k
        cmC2[18].set(word00*1471+ xword+(stream3bR & 511)*191+BrFcIdx);//needs work, combined with above 17k ---57k

        if (worcxt.wordcount) {
            simiwor=sencxt.SimilarSentence(&worcxt,worcxt.wordcount);
            xword1=simiwor->WordR(worcxt.wordcount);
            simNoun=simiwor->LastR(1,worcxt.wordcount, Noun);
            simVerb=simiwor->LastR(1,worcxt.wordcount, Verb);
        }
        // List
        if (worcxt.wordcount && (colcxt.lastfc()=='*' || colcxt.lastfc(1)=='#')) {
            simiwor=sencxtL.SimilarSentence(&worcxt,worcxt.wordcount);
            xword1=simiwor->WordR(worcxt.wordcount);
        // Tables
        }else if (worcxt.wordcount && colcxt.nlChar==WIKITABLE) {
            simiwor=sencxtT.SimilarSentence(&worcxt,worcxt.wordcount);
            xword1=simiwor->WordR(worcxt.wordcount);
        // Wiki links
        }else if (worcxt.wordcount && colcxt.lastfc()==SQUAREOPEN) {
            simiwor=sencxtCL.SimilarSentence(&worcxt,worcxt.wordcount);
            xword1=simiwor->WordR(worcxt.wordcount);
        }
        
        cmC2[18].set(word0*83+worcxt1.Word(1)+xword1*191); // needs work 25k disable before census
        cmC2[18].set(h+worcxt.Word(2)*83+xword1*53+(stream3b & 511)+FcIdx+((indirectBrByte>>0)&0x7ff)*191); // seems ok 56k?
        cmC2[18].set(h+simVerb+simNoun+linkword+fccxt.cxt );//30k disable before census
        
        if (brcxt.cxt==LESSTHAN || (colcxt.lastfc()!=FIRSTUPPER && colcxt.lastfc()!='*' && colcxt.lastfc()!=APOSTROPHE)) cmC2[19].sets(); 
        else cmC2[19].set(h+worcxt0.Word0(3)+worcxt0.Word0(4));
        // end
        
        // current word and word(1) type upto preffix(not included), paragraph word(1) 
        cmC4[6].set(h+(worcxt.Type(1)&(0x1FF))+worcxt1.Word(1)); 
        
        if (isMatch>61) cmC2[6].sets(); else cmC2[6].set(((stream2b&15)<<16)+(t[2]&0xffff));  // o2 31 ka

        if (c1==ESCAPE || utf8left || fccxt.cxt==CURLYOPENING) 
            cmC2[7].set(0);
        else 
            cmC2[7].set(indirectBrByte);

        if (skipM1) {
            cmC2[8].sets();
            cmC2[8].sets();
            cmC2[8].sets();
        } else {
            cmC2[8].set(((indirectBrByte>>0)&0x7ff)*32 + ((stream4b & 0xfff0) << 16)+BrFcIdx);
            cmC2[8].set((stream3bR&0x3fffffff)*4+(stream2b&3));
            cmC2[8].set((fccxt.cxt*4) + ((stream3bR & 0x3ffff) << 9 )+BrFcIdx);
        }
        if (fccxt.cxt==HTLINK) 
            cmC2[8].sets();
        else
            cmC2[8].set((c4 & 0xffffff) + ((stream2b << 18) & 0xff000000));

        if (skipM1) {
            cmC4[0].sets();
            cmC4[0].sets();
        }else{
            cmC4[0].set(colcxt.lastfc(0) | (fccxt.cxt<< 15) | ((stream3b & 63) << 7)|(brcxt.cxt << 24) );
            cmC4[0].set((colcxt.lastfc(0) | ((c4 & 0xffffff) << 8)));
        }

        cmC4[1].set( (stream2b & 3) +word00*11);
        cmC4[1].set((c4 & 0xffff));
        cmC4[1].set(((fc<<11) | c1)+((stream2b&3)<< 18));

        cmC4[2].set((stream2b & 15)+((stream3b&7) << 6 ));
        cmC4[2].set(c1 | ((col * (c1==SPACE)) << 8)|((stream2b & 15) << 16));
 
        if (isCategory) cmC4[2].sets(); 
        else cmC4[2].set((wt4cxtW1*191)+word0);
        if (c1==ESCAPE || fc==SPACE || utf8left)  
            cmC4[2].sets();
        else 
            cmC4[2].set((91 * 83* worcxt.Word(1) + 89 * word0));
           
        if (fc==SPACE) 
            cmC4[4].sets();
        else
            cmC4[4].set((c1 + ((stream3b & 0xe38) << 6)) );

        cmC4[4].set(worcxt.fword*11+BrFcIdx);
        cmC4[4].set(c1+word0+number0*191);//00
        cmC4[4].set(((c4 & 0xffff) << 16) | (fccxt.cxt  << 8) |fc);
        cmC4[4].set(((stream3bR & 0xfff)<< 8)+((stream2b & 0xfc)));
        cmC4[4].set(h+PStateH);

        if (c1==ESCAPE) {
            cmC[0].sets();cmC[0].sets();cmC[0].sets();
            cmC[0].sets();cmC[0].sets();cmC[0].sets(); 
        } else {
            // Switch between word/paragraph or column mode
            if (isParagraph==1  ) {
                // Word
                if (fccxt.cxt==SQUAREOPEN) cmC[0].sets(); 
                else cmC[0].set(worcxt.fword*3191+(stream2b & 3));
                cmC[0].set(h+firstWord*89);
                cmC[0].set(word0*53+c1+BrFcIdx+PStateH);
            } else {
                // Column
                if (fccxt.cxt==SQUAREOPEN) cmC[0].sets(); 
                else cmC[0].set(above | ((stream3b & 0x3f)<<9) | (colcxt.collen()<<19)| ((stream2b & 3) << 16) );
                cmC[0].set(h+firstWord*89);
                if (fccxt.cxt==SQUAREOPEN) cmC[0].sets(); 
                else cmC[0].set(above | (c1 << 16)| ((col+numlen0+BrFcIdx)<<8)| (above1<<24)  );
            }
            if (colcxt.lastfc()=='*') {
                // List
                cmC[0].set((word0+( ( fccxt.cxt) << 8)) | ((BrFcIdx ) << 16));// or not add
                cmC[0].set(c1);
                cmC[0].set(word0+PStateH);//00
            
            } else {
                // Table
                cmC[0].set(wrt_2b[bufr(colcxt.abovecellpos)]|(fccxt.cxt<<8) | (BrFcIdx<<16));
                cmC[0].set(bufr(colcxt.abovecellpos) | (c1<<8) );
                cmC[0].set( word0+wrt_2b[bufr(colcxt.abovecellpos)] );
            }
        }
        // xml
        if (isXML==true) {
            cmC44.set(xlU1);// ?
            cmC44.set(xlU2);
            cmC44.set(xlU3);
            cmC44.set(xlU4);
        } else {
            cmC44.sets();
            cmC44.sets();
            cmC44.sets();
            cmC44.sets();
        }
               
        if (fc==SPACE|| skipM1 || skipSeeExternal || isPageStarted==false) {
            cmC[1].sets(); cmC[1].sets(); cmC[1].sets();
        } else {
            cmC[1].set((stream3b&0x7fff)*word0+BrFcIdx );
            cmC[1].set((x4&0xff0000ff) | ((stream3b&0xe07) << 8));
            cmC[1].set((indirectBrByte&0xffff) | ((stream3b&0x38) << 16));
        }
        // Indirect byte with sentence word(1) and current byte
        if (isMath) cmC[0].sets(); 
        else cmC[0].set((indirectByte& 0xff00)+257 * worcxt.Word(1)*53+(c1 ));// brcxt less?

        cmC[2].set((c1<<8) | (indirectByte>>2)| (fc<<16));  //
        cmC[2].set((c4 & 0xffff)+(c2==c3?1:0));
   
        cmC4[3].set((stream3b & stream3bMask)*256 | (stream2b&stream2bMask& 255) );
        if ( skipSeeExternal) cmC4[3].sets(); 
        else cmC4[3].set(x4);
        // Word stream. word(1) with first char context and last bit3word(1-x)
        cmC2[9].set(257*(*pWord).Hash+fccxt.cxt+193*(stream3b & stream3bMask));
        // First char, current byte and non repeating bit2word(1-6)
        cmC2[9].set(fc|((stream2bR & 0xfff) << 9) | ((c1  ) << 24));//end is good (lang)

        cmC2[16].set(  worcxt.fword*83+(stream2b & 15)*11+brcxt.cxt); // all category/language/image links (better as standalone)
        if (skipM1 || skipSeeExternal) cmC2[17].sets(); 
        else  cmC2[17].set( worcxt.Last(1,Verb)+worcxt.Word(1)*83+h);// brcxt less?

        cmC2[9].set((x4 & 0xffff00)+ brcxt.cxt+(fccxt.cxt<< 24));
        // Wikipedia has lot of links in form: [word word ...]. We collect context of whole link as singele word, no gaps.
        // Skip when html/xml tags
        if (linkword)        
            cmC2[9].set(linkword);
        else if (isMath) cmC2[9].sets();
        else if (senword)        
            cmC2[9].set(senword*1471+c1); // needs more work
        else {
            if (fc==HTML || brcxt.cxt==LESSTHAN) 
                cmC2[9].sets();
            else 
                cmC2[9].set(0);
        }
        
        cmC2[10].set(indirectByte);
        cmC2[10].set(((indirectByte& 0xffff00)>>4) | ((stream2b&stream2bMask & 0xf) )| ((stream3b & 0xfff) << 20));
        cmC2[10].set((x4 >>16) | ((stream2b & 255) << 24));
        if (c1>127) cmC2[10].set(( (((stream2b & 12)*256)+c1) << 11) | ((indirectWord & 0xffffff)>>16) );
        else cmC2[10].set((c1 << 11)| (BrFcIdx  << 8) | ((indirectWord & 0xffffff)>>16) );
        if (isMath) cmC2[10].sets(); 
        else cmC2[10].set((fccxt.cxt*4+BrFcIdx) | ((c4 & 0xffff)<< 9)| ((stream2b & 0xff) << 24)); 
        cmC2[10].set(((indirectWord >> 16) )| ((stream2b & 0x3c)<< 25 )| (((stream3b & 0x1ff))<< 16 ));

        cmC2[11].set(((words) )+((( spaces ))<< 8)+((stream2b&15)<< 16)+(((stream3bR>>3)&511)<< 21)+(isParagraph<<30));
        cmC2[11].set(c1 + ((stream3b<< 5) & 0x1fffff00));
        cmC2[11].set(stream2bR*16+BrFcIdx );
        // Indirect byte with bracket context and last non repeating bit2words(2-10)
        if (brcxt.cxt==LESSTHAN) cmC2[11].sets(); 
        else cmC2[11].set(((indirectByte& 0xffff)>>8) + ((64 * stream2bR) & 0x3ffff00)+(brcxt.cxt<< 25)); // end good
        // Byte from prvious word(0), pos if in range(255), indirect byte
        
        if ((fccxt.cxt==FIRSTUPPER && brcxt.cxt==SQUAREOPEN) || brcxt.cxt==LESSTHAN/* ||  (xmlS&7)!=xReadContent*/) {
            cmC2[11].sets();
        } else {
            cmC2[11].set((indirectWord0Pos )| ((indirectByte& 0xff00)<<16));// good 2,7k
        }
        // Byte stream of x4, msb of byte(4), 4 msb bits of byte(2,3) and full byte(1)
        if (brcxt.cxt==LESSTHAN||isCategory) cmC2[12].sets(); 
        else cmC2[12].set((x4&0x80f00000)+((x4&0x0000f0ff) << 12) );

        // Paragraph or column. 
        // In Paragraph: disabled when escaped utf8, html link, math
        // In Column: when col is max(31) use last two bytes only otherwise add above bytes
        if (brcxt.cxt==LESSTHAN )cmC2[12].set(h+worcxt3.Word(1) *53 *79+worcxt3.Word(2) *53*47 *71);
        else if (isParagraph==1) {
            // word
            if (c1==ESCAPE || fccxt.cxt==HTLINK || fccxt.cxt==CURLYOPENING || isMath|| isPre) {
                cmC2[12].sets();
            } else {
                cmC2[12].set(h+worcxt.Word(1) *53 *79+worcxt.Word(3) *53*47 *71);
            }
        } else {
            // Skip when html link, tag
            if (fccxt.cxt==HTLINK || brcxt.cxt==LESSTHAN || htcxt.cxt || isCategory){
                cmC2[12].sets();
            }
            //column
            else if (col==31) {
                cmC2[12].set(c4<<16);
            } else {
                cmC2[12].set(above | ((c4 &0xffff)<< 16)| (above1<< 8));
            }
        }
        // Word/centence. 
        if (c1==ESCAPE || utf8left || fccxt.cxt==CURLYOPENING || fccxt.cxt==HTLINK || fc==HTML || htcxt.cxt || 
         isCategory||fc==SPACE || isPre || c1=='&' || brcxt.cxt==LESSTHAN || isMath|| col<2 || (worcxt.sBytes(0)>>8)=='\\') {
            // Disabled when: 
            // escaped utf8, template (onliner) or table beginning, 
            // html link, html (tag), fist char space, &
            // start of possible html tag, col not started including first char
            cmC2[13].sets();
            cmC2[13].sets();
        } else {
            // Word/Centence with current word(0), word(1) and word(2)
            cmC2[13].set(worcxt.Word(1)*83*1471-(word0)*53+worcxt.Word(2));
            cmC2[13].set(h+worcxt.Word(2)*53*79+worcxt.Word(3)*53*47*71);
        }
       
        // Last byte in stream3bR and stream2b type with first char and bracket index
        cmC[3].set(((stream3bR&7)<< 10) + (stream2b&3)+fc*4+ (BrFcIdx<< 24));
        // Current word or number
        cmC[3].set( ((linkword?linkword:word0)*3301+ number0*3191));
        if (c1==ESCAPE || skipM1 || utf8left || fccxt.cxt==CURLYOPENING || fccxt.cxt==HTLINK || fc==SPACE || fc==HTML ||
         isCategory || brcxt.cxt==LESSTHAN|| col<2 || isMath || (worcxt.sBytes(0)>>8)=='\\') {
            // Disabled when: 
            // escaped utf8, template (onliner) or table beginning, 
            // html link, html (tag), fist char space,
            // start of possible html tag, col not started including first char
            cmC2[14].sets();//end - fc [ ?????
        } else {
            // Word/centence and non-repeating bit3words upto word(2) with bracket/firstchar index
            cmC2[14].set(BrFcIdx+ worcxt.Word(2) * (stream3bR&stream3bRMask2)+(worcxt.Type(1)&(0x1ff)));
        }
        // Local, small memory
        if (c1==ESCAPE|| utf8left || fc==SPACE || skipSeeExternal) {
            // Disabled when: escaped utf8, fist char space,
            cmC4[7].sets(); cmC4[7].sets(); cmC4[7].sets(); cmC4[7].sets();
        } else {
            cmC4[7].set(worcxt1.Word()+word00);
            cmC4[7].set(worcxt.Word(2)+word0*191+(stream3bR & 63));
            cmC4[7].set(word0*191+(stream3bR & 63));
            cmC4[7].set((indirectWord0Pos&0xffff)*191+ word0+(stream3bR & 63));
        }
        
        if (brcxt.cxt==SQUAREOPEN)
            cmCR[0].set((worcxt.Word(1)+word0));
        else if (fccxt.cxt==HTLINK)
            cmCR[0].sets();
        else 
            cmCR[0].set(0);

        cmCR[1].set((fccxt.cxt<< 15) | ((stream3b & 7) << 3)|(brcxt.cxt << 24));
        cmCR[2].set( hash(  worcxt.Word(), stream2b& 0xFC, c1));
 
        scmA[0].set(c1);
        scmA[1].set(stream3b&0x1ff);
        scmA[2].set(brcxt.cxt);
 
        if (wshift||c1==LF) {
            U16 sb=worcxt0.sBytes(1);
            U32 w=worcxt0.Word(1);
            worcxt0.Remove();
            worcxt0.Set(sb>>8,0);
            worcxt0.Update(w,LF,0,w);
            wshift=0;
            if (c1==LF)sVerb=0;
        }

        cmC2[15].set((BrFcIdx*256)+fc+((stream3bR&0xFFF)<< 16));

        U32 w4=(c4<<8)&0xff000000;
        // 3 bit stream of mapped chars, all chars above 127 are mapped to a-z range.
        stream5b=(stream5b<<3) | (c1>127?wrt_3b['a']:wrt_3b[c1]);
        const U8 buf2=(c4>>8)&0xff;
        const U8 buf3=(c4>>16)&0xff;
        //sparse indirect 1-byte contexts
        U8 g=indirect2[c1]; // context from indirect
        indirect2[buf2]=c1; // update indirect
        cmcr[0].set(hash(indirectBrByte, w4|g));
        //sparse indirect 2-byte contexts
        g=indirect3[buf2<<8 | c1]; // context from indirect
        indirect3[buf3<<8 | buf2]=c1; // update indirect
        cmcr[0].set(hash(indirectBrByte, w4|g, colcxt.lastfc()));

        //update s5bByte
        s5bByte[((stream5b>>3)&7)]=c4&0xffffff;
        for (int i=0; i<8; i++) {
            U32 k=s5bByte[i];
            cmcr[i+1].set(hash(BrFcIdx+(nestList==false?word0*191:0), k, deccode>>2));
            cmcr[i+1].set(hash(BrFcIdx, stream5b&0x7fff, k & 0xffff));
            cmcr[i+1].set(hash(BrFcIdx, ((stream5b&7)<<8) | (k & 0xff00ff),colcxt.lastfc()));
        }
        //update s2Word0
        s2Word0[BrFcIdx*16+((stream2b>>2)&15)]=word0;
        for (int i=0; i<4; i++) {//disable after 77.00%-89.53, 93.10-97.46% (enwik9)
            U32 k=s2Word0[BrFcIdx*16+((stream2b>>(2*i) )&15)];
            cmcr2[i].set(k);
            cmcr2[i].set(k+h);
            cmcr2[i].set(k+worcxt.Word());
        }
        // Some APM context
        AH1=hash((x5>>0)&255, (x5>>8)&255, (x5>>16)&0x80ff);
        AH2=hash(19,     x5&0x80ffff);
        
        maps1.set((word0*191));
        maps2.set(deccode>>2);
    }
    
    const int c0b=c0<<(8-bpos);
    dcsm.set((word0*191)*256+x.c0,x.y);
    dcsm.set((word0*191+worcxt.Word(1))*256+x.c0,x.y);
    dcsm.set((word0*191+worcxt.Word(2))*256+x.c0,x.y);
    dcsm.set((word0*191+worcxt.Word(3))*256+x.c0,x.y);
    dcsm.set((word0*191+worcxt.Word(4))*256+x.c0,x.y);
    
    dcsm2.set((word00*191+indirectBrByte)*256+x.c0,x.y);
    dcsm2.set((word00*191+worcxt0.Word(1)+indirectBrByte)*256+x.c0,x.y);
    dcsm2.set((word00*191+worcxt0.Word(2)+indirectBrByte)*256+x.c0,x.y);

    //indirectBrByte
    dcsm1.set( (indirectBrByte)*256+x.c0,x.y);
    dcsm1.set( (cxtind3)*191+x.c0,x.y);

    dcsm0.set((h+worcxt1.Word(1))*256+x.c0,x.y);
    dcsm0.set((h+worcxt1.Word(2))*256+x.c0,x.y);
    dcsm0.set((h+worcxt1.Word(3))*256+x.c0,x.y);
    dcsm0.set((h+worcxt1.Word(4))*256+x.c0,x.y);
    dcsm0.set((h+worcxt1.Word(5))*256+x.c0,x.y);
    dcsm0.set((h+worcxt1.Word(6))*256+x.c0,x.y);
    
    dcsmN.set((h+worcxt.Code(1)+fccxt.cxt )*256+x.c0,x.y);
    dcsmN.set((h+worcxt.Code(2)+fccxt.cxt)*256+x.c0,x.y);
    dcsmN.set((h+worcxt.Code(3)+fccxt.cxt)*256+x.c0,x.y);
    
    maps1.mix();
    maps2.mix();
    scmA[0].mix(sscmrate);
    scmA[1].mix(sscmrate);
    scmA[2].mix(sscmrate);
  
    isMatch=MatchModel2mix();
    
    // Order X
    ordX=0;
    if (cmC2[0].cxtMask) ordX=2;
    ordX=ordX+cmC2[0].mix();
    if (ordX==3) ordX=2; // low max 2
    ordX=ordX+cmC2[1].mix();
    ordX=ordX+cmC2[2].mix();
    ordX=ordX+cmC2[3].mix();
    ordW=cmC2[4].mix();
    ordW=ordW+cmC2[5].mix();
    if (ordW>3) ordW=3;
    cmC2[6].mix();
    cmC2[7].mix();
    cmC2[8].mix();
    cmC4[0].mix();
    cmC4[1].mix();
    cmC4[2].mix();
    cmC4[4].mix();
    int ordP=cmC[0].mix();
    cmC[1].mix();
    cmC[2].mix();
    cmC4[3].mix();
    cmC2[9].mix();
    cmC2[10].mix();
    cmC2[11].mix();
    cmC2[12].mix();
    cmC44.mix();

    // Order Word
    ordW=ordW+cmC2[13].mix();  
    cmC[3].mix();  
    ordW=ordW+cmC2[14].mix();   
    
    cmC2[15].mix();
    cmC[4].mix();
    cmC[5].mix();
    cmC2[16].mix();
    cmC2[17].mix();
    cmC2[19].mix();
    cmC4[6].mix();
    cmC4[7].mix();
    cmC4[8].mix();
    cmC2[18].mix();
    cmC2[20].mix();
    cmCR[0].mix();
    cmCR[1].mix();
    cmCR[2].mix();
    for (int i=0;i<9;i++)
        cmcr[i].mix();
    for (int i=0;i<4;i++)
        cmcr2[i].mix();
    
    rcmA[0].mix();
    dcsm.mix();
    dcsm0.mix();
    dcsm1.mix();
    dcsm2.mix();
    dcsmN.mix();
    
    mstate=STA2[mstate][x.y];
   
    AddPrediction((64));
    // Mixers

    mxA1[0].cxt=(ordX*8 + (BrFcIdx?1:0)*4 + (stream2b&3))*2+(words&1);
//    mxA1[1].cxt=0;
    mxA1[2].cxt=stream2b&0x3f;
    mxA1[3].cxt=(stream2b&3)*4+wrt_2b[c0b&255];
    mxA1[4].cxt=ordX*4+wrt_2b[c0b&255];
//    mxA1[5].cxt=0;
 
    mxA[17].cxt=PState*16+(stream2b&15);
    // mixer 0
    // at bpos=0   context is last 2 bit2word and 1 bit3word
    // at bpos=1-3 context is last 2 bit2word and first char/bracket index (max 7)
    // at bpos=4-7 context is last 1 bit2word, current bit2word from c0 and bit3word of current bracket or quote
    if (bpos==0)  mxA[0].cxt=(stream2bR&4095)*8 + (stream3b&7);
    else if (bpos>3) {
        c=wrt_2b[c0b&255];
        mxA[0].cxt=(((stream2b<<2)&4095)+c)*8+BrFcIdx;// no 4095
    } else    
        mxA[0].cxt=(stream2b&4095)*8 +BrFcIdx;
        
    // mixer 1
    // at bpos=0   context is was byte(3,4) a word and 2 bit3word
    // at bpos=1   context is bit 1xxxxxxx from c0, was byte(2) a word, bit pos max 5,last 1 bit3word and first char/bracket index (max 7)
    // at bpos=2   context is bit 11xxxxxx (bit pos 2) from c0, was byte(2) a word, bit pos max 5,last 1 bit3word and first char/bracket index (max 7)
    // at bpos=3   context is bit 111xxxxx (bit pos 3) from c0, was byte(2) a word, bit pos max 5,last 1 bit3word and first char/bracket index (max 7)
    // at bpos=4-7 context is bit current bit2word from c0, bit pos max 5,last 1 bit3word and first char/bracket index (max 7)
    if (bpos){
         c=c0b; 
         if (bpos==1) c=c+16 * (words*2& 4);
         else if (bpos>3)  c=wrt_2b[c0b&255]*64;
         c=(min(bpos,5))*256+(stream3bR&7)+FcIdx*8+(c&192); //FcIdx -> cm(12,1)
    }
    else c=(words&12)*16+(stream3bR&7)+BrFcIdx*8;
    mxA[1].cxt=c;
    
    // mixer 2
    // at bpos=0-7   context is was byte(3,4) a word, sum of context order(3-5,6,8) isState counts (max 5) and last 2 bit2word
    mxA[2].cxt=((4 * words) & 0xf0) + ordX*256 + (stream2b & 15);
    
    // mixer 6
    // at bpos=0   context is non-repeating 2 bit3word of byte(2,3), first char/bracket index (max 7) and last 1 bit2word
    // at bpos=1-3 context is non-repeating 2 bit3word of byte(2,3), was byte(1-3) a word and last 1 bit2word
    // at bpos=4-7 context is non-repeating 2 bit3word of byte(2,3), last 1 bit2word and last partialy decoded bit2word
    if (bpos>3) {
        c=wrt_2b[c0b&255];
        mxA[6].cxt=((stream3bR) & 0xff8)*4 + ((stream2b & 3)*4) + (c );
    }
    else if (bpos==0) 
        mxA[6].cxt=((stream3bR) & 0xff8)*4 + ((4 * FcIdx) ) + (stream2b & 3);
    else
        mxA[6].cxt=((stream3bR) & 0xff8)*4 + ((2 * words) & 0x1c) + (stream2b & 3);

    c=c0b;
    
    // mixer 3
    // at bpos=0   context is bit xxxxxxxx from c0, was byte(1-8) a word or space and bit pos
    // at bpos=1   context is bit 1xxxxxxx from c0, was byte(1-7) a word or space and bit pos
    // at bpos=2   context is bit 11xxxxxx from c0, was byte(1-6) a word or space and bit pos
    // at bpos=3   context is bit 111xxxxx from c0, was byte(1-5) a word or space and bit pos
    // at bpos=4   context is bit 1111xxxx from c0, was byte(1-4) a word or space and bit pos
    // at bpos=5   context is bit 11111xxx from c0, was byte(1-3) a word or space and bit pos
    // at bpos=6   context is bit 111111xx from c0, was byte(1-2) a word or space and bit pos
    // at bpos=7   context is bit 1111111x from c0, was byte(1)   a word or space and bit pos
    mxA[3].cxt=bpos*256 + (((( (numbers|words)<< bpos)&255)>> bpos) | (c&255));
    
    // mixer 8 - final mixer
    // at bpos=0-7   context is sum of context order(3-5,6,8) isState counts (max 5), 
    // bracket or quote state(0,1) 
    // last 1 bit2word
    // was last byte a word
    
    // mixer 4 
    // at bpos=0   context is bit xxxxxxxx from c0, first char type state(0,1) xxxx1xxx, 2 bit2word            1111xxxx
    // at bpos=1   context is bit 1xxxxxxx from c0, first char type state(0,1) xxxx1xxx, 1 bit3word            x111xxxx, bit pos xxxxx111
    // at bpos=2   context is bit 11xxxxxx from c0, first char type state(0,1) xxxx1xxx, 1 bit2word            xx11xxxx, bit pos xxxxx111
    // at bpos=3   context is bit 111xxxxx from c0, first char type state(0,1) xxxx1xxx, was byte(1) a word    xxx1xxxx, bit pos xxxxx111
    // at bpos=4   context is bit 1111xxxx from c0, first char type state(0,1) xxxx1xxx,                                 bit pos xxxxx111
    // at bpos=5   context is bit 11111xxx from c0, first char type state(0,1) xxxx1xxx (overflow, ok!)
    // at bpos=6   context is bit 111111xx from c0, first char type state(0,1) xxxx1xxx
    // at bpos=7   context is bit 1111111x from c0, first char type state(0,1) xxxx1xxx
    // at bpos=0-7 sum of context order(3-5,6,8) isState counts (max 5) and is match(0,1) 111 xxxxxxxx

    if (bpos) {
        if (bpos==1) {
            c=c + 16*(stream3b&7);
        }
        else if (bpos==2) {
            c=c + 16*(stream2b&3);
        }
        else if (bpos==3) {
            c=c + 16*(words&1);
        } else  {
            c=bpos + (c&0xf0);
        }
        if (bpos<5)
            c=bpos + (c&0xf0); 
    }else   c=16 * (stream2b&0xf);
    ordX=ordX-1;
    if (ordX<0)
       ordX=0;
    if (isMatch)
        ordX=ordX+1;
    mxA[4].cxt=c + ordX*256+ 8*isParagraph;

    // mixer 5
    // at bpos=0-7   context is sum of context words isState counts (max 6), first char index (0-7), 2 bit2word of byte(3,4) and 1 bit3word of byte(2)
    mxA[5].cxt=(ordW*256 + (stream2b&0xf0) + ((stream3b&0x38) >> 2))*4 + FcIdx;

    // mixer 7 
    // at bpos 0-2 bit3word(low 2 bits), first char/bracket index (max 7), was byte(1,2,3) a word, first char flag, is a match
    // at bpos 3-7 bit3word from c0, first char/bracket index (max 7), was byte(1,2,3) a word, first char flag, is a match
    if (bpos>2) 
        mxA[7].cxt=((stream3b&7)*8+wrt_3b[c0b&255])*256 +(BrFcIdx)*32 + (words&7)*4 + isParagraph+(isMatch?2:0);
    else
        mxA[7].cxt=((stream3b&63)*256 +(BrFcIdx)*16 + (words&7)*2 + isParagraph)|(isMatch?128:0);

    mxA[8].cxt=deccode;
    mxA[9].cxt=(stream3bR&511)*16*16+FcIdx*32+isParagraph*16+(lastWT&15);
    mxA[10].cxt=(x.c4&0xffff);
    mxA[11].cxt=(stream3b&0x3f)*256 +x.c0;  
    mxA[12].cxt=oldwt1+(stream2bR&255)*32;
    mxA[13].cxt=(stream3bR&511);
    mxA[14].cxt=(BrFcIdx*8+FcIdx)*8+wrt_3b[c0b&255];
    mxA[15].cxt=(numbers|words)*16+(stream2bR&15); 
    mxA[16].cxt=xmlS&1023; 
 
    x.mxInputs2.add(mxA[0].p1());
    x.mxInputs2.add(mxA[1].p1());
    x.mxInputs2.add(mxA[2].p1());
    x.mxInputs2.add(mxA[3].p1());
    x.mxInputs2.add(mxA[4].p1());
    x.mxInputs2.add(mxA[5].p1());
    x.mxInputs2.add(mxA[6].p1());
    x.mxInputs2.add(mxA[7].p1());
    x.mxInputs2.add(mxA[8].p1());
    x.mxInputs2.add(mxA[9].p1());
    x.mxInputs2.add(mxA[10].p1());
    x.mxInputs2.add(mxA[11].p1());
    x.mxInputs2.add(mxA[12].p1());
    x.mxInputs2.add(mxA[13].p1());
    x.mxInputs2.add(mxA[14].p1());
    x.mxInputs2.add(mxA[15].p1());
    if (isXML==true) x.mxInputs2.add(mxA[16].p1()); else  x.mxInputs2.add(0);
    x.mxInputs4.add(mxA[17].p1());
    //l2
    x.mxInputs4.add(mxA1[0].p1());
    x.mxInputs4.add(mxA1[1].p1());
    x.mxInputs4.add(mxA1[2].p1());
    x.mxInputs4.add(mxA1[3].p1());
    x.mxInputs4.add(mxA1[4].p1());
    return (mxA1[5].p());
}

int rate=6;
void update() {
    x.c0+=x.c0+x.y;

    if (x.c0>=256) {
        x.c4=(x.c4<<8)+(x.c0&0xff);	
        x.c0=1;
        ++x.blpos;
        // When last byte was predicted good/below error treshold then set new limits to mixer update
        // larger value means less updates and better speed.
        if ((fails&255)==0) {
            for (int i=0;i<18;i++) mxA[i].elim=max(256+63,mxA[i].elim+1);
        }else{ 
            for (int i=0;i<18;i++) mxA[i].elim=max(0,min(16,mxA[i].elim-1));
        }
        sscmrate=(x.blpos>14*256*1024);
        // APM update rate based on input file position
        rate=6 + (x.blpos>14*256*1024) + (x.blpos>28*512*1024);
    }

    x.bpos=(x.bpos+1)&7;
    x.bposshift=7-x.bpos;
    x.c0shift_bpos=(x.c0<<1)^(256>>(x.bposshift));
    x.cmBitState=getStateByteLocation(x.bpos,x.c0);
    for (int i=0;i<16;i++)    mxA[i].update(x.y);
    if (isXML) mxA[16].update(x.y);
        mxA[17].update(x.y);
    for (int i=0;i<6;i++) mxA1[i].update(x.y);
    for (int i=0;i<4;i++) mxA2[i].update(x.y);
    //printf("mixer 0 predictor count %d\n",x.mxInputs1.ncount);
    x.mxInputs1.ncount=0;
    //printf("mixer 1 predictor count %d\n",x.mxInputs2.ncount);
    x.mxInputs2.ncount=0;
    //printf("mixer 2 predictor count %d\n",x.mxInputs4.ncount);
    x.mxInputs4.ncount=0; 
    //printf("mixer a2 predictor count %d\n",num_models);

    // This part is from paq8hp12
    if (fails&0x00000080) --failcount;
    fails=fails*2;
    failz=failz*2;

    if (x.y) pr=4095-pr;
    if (pr>=e_l[x.bpos]) ++fails, ++failcount;
    if (pr>=848) ++failz;

    pr=modelPrediction();
    AddPrediction(stretch(pr));

    int pt, pu=(apmA0.p(pr, x.c0, 3,x.y)+7*pr+4)>>3, pv, pz=failcount+1;
   
    pz+=tri[(fails>>5)&3];
    pz+=trj[(fails>>3)&3];
    pz+=trj[(fails>>1)&3];
    if (fails&1) pz+=8;
    pz=pz/2;

    pu=apmA3.p(pu,   ((x.c0*2)^AH1), rate,x.y); 
    AddPrediction(stretch(pu));
    pv=apmA1.p(pr,   ((x.c0*8)^hash(29,failz&2047)), rate+1,x.y);
    AddPrediction(stretch(pv));
    // If fails use stream2b else non-repeating stream2b
    if (fails&255)
        pv=apmA4.p(pv, hash(x.c0,stream2b & 0xfffc,(stream3bR & 0x1ff)), rate,x.y);
    else
        pv=apmA4.p(pv, hash(x.c0,(stream2bR & 0xfffc)+0x10000,(stream3bR & 0x1ff)), rate,x.y);
    AddPrediction(stretch(pv));
    pt=apmA2.p(pr, ( (x.c0*32)^AH2), rate,x.y);
    AddPrediction(stretch(pt));
    pz=apmA5.p(pu,   ((x.c0*4)^hash(min(9,pz),x5&0x80ff)), rate,x.y);
    AddPrediction(stretch(pz));
    if (fails&255) pr=(pt*6+pu  +pv*11+pz*14 +31)>>5;
    else           pr=(pt*4+pu*7+pv*12+pz*9 +31)>>5;
    AddPrediction(stretch(pr));
    //l 4
    mxA2[0].cxt=c1;
    mxA2[1].cxt=c2;
    mxA2[2].cxt=(stream5b)&0x7fff;
    mxA2[3].cxt=stream2b&0xffff;
    mp0=mxA2[0].p1();
    mp1=mxA2[1].p1();
    mp2=mxA2[2].p1();
    mp3=mxA2[3].p1();
    
    pu=stretch(pr);
    
    mmmO[0].update(x.y),  pu=mmmO[0].pp(0,mp0,mstate);//clp this
    mmmO[1].update(x.y),  pu=mmmO[1].pp(pu,mp1,mstate);
    mmmO[2].update(x.y),  pu=mmmO[2].pp(pu,mp2,mstate);
    mmmO[3].update(x.y),  pu=mmmO[3].pp(pu,mp3,mstate);
    
    pr=(squash(clp(pu))+pr*3)>>2;
}

void put32(U32 x,FILE *f){
    fputc((x >> 24) & 255, f); 
    fputc((x >> 16) & 255, f); 
    fputc((x >> 8) & 255, f); 
    fputc(x & 255, f);
}
void put64(U64 x,FILE *out){
    putc((x >> 56) & 255,out);
    putc((x >> 48) & 255,out);
    putc((x >> 40) & 255,out);
    putc((x >> 32) & 255,out);
    putc((x >> 24) & 255,out); 
    putc((x >> 16) & 255,out); 
    putc((x >> 8) & 255,out); 
    putc(x & 255,out);
}

U64 get64(FILE *in) {
    return ((U64)getc(in) << 56) | ((U64)getc(in) << 48) | ((U64)getc(in) << 40) | ((U64)getc(in) << 32) | (getc(in) << 24) | (getc(in) << 16) | (getc(in) << 8) | (getc(in));
}
// Encoder
unsigned long long percent,input_bytes;
typedef enum {COMPRESS, DECOMPRESS} Mode;
clock_t start;
class Encoder {
private:
  const Mode mode;       // Compress or decompress?
  FILE* archive;         // Compressed data file
  U32 x1, x2;            // Range, initially [0, 1), scaled by 2^32
  U32 x0;                // Decompress mode: last 4 input bytes of archive
  FILE*alt;              // decompress() source in COMPRESS mode
  U32 scompsize;

  void code(int i=0) {
    int p=pr;
    p+=p==0;
    assert(p>0 && p<4096);
    U32 xmid=x1 + ((x2-x1)>>12)*p + (((x2-x1)&0xfff)*p>>12);
    assert(xmid>=x1 && xmid<x2);
    x.y=i;
    i ? (x2=xmid) : (x1=xmid+1);
    update();
    ResetPredictions();
    while (((x1^x2)&0xff000000)==0) {  // pass equal leading bytes of range
       fputc(x2>>24, archive);
      x1<<=8;
      x2=(x2<<8)+255;
    }
  }
  int decode() {
    int p=pr;
    p+=p==0;
    assert(p>0 && p<4096);
    U32 xmid=x1 + ((x2-x1)>>12)*p + (((x2-x1)&0xfff)*p>>12);
    assert(xmid>=x1 && xmid<x2);
    x0<=xmid ? (x2=xmid,x.y=1) : (x1=xmid+1,x.y=0);
    update();
    ResetPredictions();
    while (((x1^x2)&0xff000000)==0) {  // pass equal leading bytes of range
      x1<<=8;
      x2=(x2<<8)+255;
      x0=(x0<<8)+(fgetc(archive)&255);  // EOF is OK
    }
    return x.y;
  }

public:
  Encoder(Mode m, FILE* f);
  Mode getMode() const {return mode;}
  U32 size() const {return  ftell (archive);}  // length of archive so far archive->curpos()
  void flush();  // call this when compression is finished
  void setFile(FILE* f) {alt=f;}

  // Compress one byte
  void compress(int c) {
    assert(mode==COMPRESS);
#ifdef PLAINTEXT
        c=charSwap(c);
#endif
      for (int i=7; i>=0; --i)
        code((c>>i)&1);
        if (x.blpos % percent == 0) {
      double frac = 100.0 * pos / input_bytes;
      fprintf(stderr, "%.2f %d %d   %.2fh \n", frac, size(),size()-scompsize,(((double(clock()-start)/CLOCKS_PER_SEC)*100.0/frac)/60)/60);scompsize=size();
    }
  }

  // Decompress and return one byte
  int decompress() {
      int c=0;
      for (int i=0; i<8; ++i)
        c+=c+decode();
#ifdef PLAINTEXT
        c=charSwap(c);
#endif
      return c;
  }
  ~Encoder() { }
};

Encoder::Encoder(Mode m, FILE* f):
    mode(m), archive(f), x1(0), x2(0xffffffff), x0(0), alt(0) {      
    PredictorInit();
    scompsize=0;
    // x0 = first 4 bytes of archive
    if (mode==DECOMPRESS) {
        for (int i=0; i<4; ++i)
            x0=(x0<<8)+(fgetc (archive)&255);
    }
}

void Encoder::flush() {
    if (mode==COMPRESS)
        put32(x1,archive);  // Flush first unequal byte of range
}

//https://stackoverflow.com/questions/669438/how-to-get-memory-usage-at-runtime-using-c
/*#if defined(_WIN32)
#include <psapi.h>
#endif
size_t getPeakMemory() {
#if defined(_WIN32)
    PROCESS_MEMORY_COUNTERS info;
    GetProcessMemoryInfo( GetCurrentProcess( ), &info, sizeof(info) );
    return (size_t)info.PeakPagefileUsage; // recuested peak memory /PeakWorkingSetSize used memory
#elif defined(UNIX) 
    struct rusage rusage;
    getrusage( RUSAGE_SELF, &rusage );
    return (size_t)(rusage.ru_maxrss * 1024L); //not tested
#else
    return (size_t)0L;
#endif // -lpsapi
}
*/
size_t getPeakMemory(){return 0;}
  
int main(int argc, char** argv) {   
    // Check arguments
    if (argc!=4 || (argv[1][0]!='c' && argv[1][0]!='d')) {
        printf(
        "fxcm v%d file compressor (C) 2024, Kaido Orav.\n"
        "Licensed under GPL, http://www.gnu.org/copyleft/gpl.html\n"
        "\n"
        "To compress:   fxcm c input output \n"
        "To decompress: fxcm d input output \n",VERSION);
        return 1;
    }

    // Open input file
    FILE *in=fopen(argv[2], "rb");
    if (!in) exit(1);
    FILE *out=0;

    // Precalculate tabeles
    int o=2;
    for (int i=0; i<1024; ++i)
        dt[i]=4096/(o),o++;
    dt[1023]=1;

    // Stretch table
    for (int i=0; i<=4095; i++) {
        strt[i]=stretchc(i);
    }

    // Squash table
    for (int i=-2047; i<=2047; i++) {
        sqt[i+2047]=squashc(i);
    }
  
    InitIlog();
    x.Init();
    
    for (int i=0;i<4096;i++) {
        st2_p1[i]=clp(sc(13*(i - 2048)));
       // st2_p2[i]=clp(sc(14*(i - 2048)));
    } 
    // run pr
    // precalc int c=ilog(rc+1)<<(2+(~rc&1));
    for (int rc=0; rc<256; rc++) {
        int c=ilog[rc];
        c=c<<(2+(~rc&1));
        rcpr[rc+256]=clp(c);
        rcpr[rc]=clp(-c);
    }
    // Match model
    mhashtablemask=0x200000*2-1;
    alloc1(mhashtable,0x200000*2+32,mhptr,32);  
    // Generate state tables
    StateTable statetable;
    statetable.Init(28, 28, 31, 29, 23, 4, 17,&STA1[0][0]);
    statetable.Init(32, 28, 31, 28, 21, 5,  6,&STA2[0][0]);
    statetable.Init(31, 27, 30, 27, 24, 4, 27,&STA4[0][0]);
    statetable.Init(33, 31, 31, 24, 20, 4, 33,&STA5[0][0]);
    statetable.Init(28, 29, 30, 30, 23, 3, 22,&STA6[0][0]);
    statetable.Init(28, 29, 33, 23, 23, 6, 14,&STA7[0][0]);
    alloc1(model_predictions1,560+64,model_predictions1_ptr,64);

    // Load dictionary
    dosym();

    if (argv[1][0]=='c') {     
        // Compress
        // Encode header: signature, version, file size
        fseeko(in, 0, SEEK_END);
        U64 size=ftello(in);
        fseeko(in, 0, SEEK_SET);
        out=fopen(argv[3], "wb");
        if (!out) exit(1);
        fprintf(out, "xF%c", VERSION);
        put64(size,out);
        if ((size / 10000)==58645)
        percent = 1 + (size / 10000); 
        else percent = 1 + (size*4 / 1000); 
        input_bytes=size;
        Encoder e(COMPRESS, out); 
        
        int c;
        start=clock();
        fseeko(in, 0, SEEK_SET);
       
        while ((c=getc(in))!=EOF)
          e.compress(c);
        e.flush();
    } else {
        // Decompress
        // Check header signature, version, file size
        if (getc(in)!='x' || getc(in)!='F' || getc(in)!=VERSION)
            printf("Not a fx file."),exit(1);
        U64 size=get64(in);
        if (size==0) printf("Bad file size."),exit(1);
        out=fopen(argv[3], "wb");
        if (!out) exit(1);
        
        Encoder e(DECOMPRESS, in);
        while (size-->0)
            putc(e.decompress(), out);
    }
    PredictorFree();
    free(mhptr);
    printf("%ld -> %ld in %1.2f sec. %d MB memory.\n", 
    ftello(in), ftello(out), double(clock()-start)/CLOCKS_PER_SEC, U32(getPeakMemory()/1024)/1024);
    fclose(in);
    fclose(out);
    return 0;
}

