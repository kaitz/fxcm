    /* fxcm file compressor is based on paq8 model.

    Copyright (C) 2026 Kaido Orav

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


To compress:   fxcm_text c input output
To decompress: fxcm_text d input output

This compressor contains changed parts (mostly marked with comments) from fallowing open source programs:
paq8hp12
paq8px_v208
paq8pxv
cmix
fxcm

*/

#define VERSION 2

#include <cstdint>
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
#include <immintrin.h>
// 8, 16, 32 bit unsigned types (adjust as appropriate)
typedef unsigned char  U8;
typedef unsigned short U16;
typedef unsigned int   U32;
typedef unsigned long long int U64;
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

//////////////////////////// Array ////////////////////////////

template <class T> void alloc(T*&ptr, int c) {
  ptr=(T*)calloc(c, sizeof(T));
  if (!ptr) exit(1);//quit("Out of memory.\n");
}
 
// for aligned data
template <class T> void alloc1(T*&data, int c,T*&ptr,const int align=16) {
  ptr=(T*)calloc(c, sizeof(T));
  if (!ptr) exit(1);//quit("Out of memory.\n");
  data=(T*)(((uintptr_t)ptr+(align-1)) & ~(uintptr_t)(align-1));
  
}

struct BlockData {
    int y;        // Last bit, 0 or 1, set by encoder
    int c0;       // Last 0-7 bits of the partial byte with a leading 1 bit (1-255)
    U32 c4;       // Last 4 whole bytes, packed.
    U64 c8;
    int bpos;     // bits in c0 (0 to 7)
    int blpos;    // Relative position in block
    int bposshift;
    int c0shift_bpos;
    int cmBitState;
    struct Inputs{
        int ncount;     // mixer input count
        short *n,*ptr;
        void add(int p){ n[ncount++]=p; }
    };
    Inputs mxInputs[3]; // array of inputs, for two layers

    void Init(){
        y=0 ,c0=1, c4=0,bpos=0,blpos=0,bposshift=0,c0shift_bpos=0 ;
    }
};

BlockData x; //maintains current global data block

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

// Original table in paq8
//Statetable:
/*static const U8 STA0[256][4]={    
  {  1,  2, 0, 0},{  3,  5, 1, 0},{  4,  6, 0, 1},{  7, 10, 2, 0}, // 0-3
  {  8, 12, 1, 1},{  9, 13, 1, 1},{ 11, 14, 0, 2},{ 15, 19, 3, 0}, // 4-7
  { 16, 23, 2, 1},{ 17, 24, 2, 1},{ 18, 25, 2, 1},{ 20, 27, 1, 2}, // 8-11
  { 21, 28, 1, 2},{ 22, 29, 1, 2},{ 26, 30, 0, 3},{ 31, 33, 4, 0}, // 12-15
  { 32, 35, 3, 1},{ 32, 35, 3, 1},{ 32, 35, 3, 1},{ 32, 35, 3, 1}, // 16-19
  { 34, 37, 2, 2},{ 34, 37, 2, 2},{ 34, 37, 2, 2},{ 34, 37, 2, 2}, // 20-23
  { 34, 37, 2, 2},{ 34, 37, 2, 2},{ 36, 39, 1, 3},{ 36, 39, 1, 3}, // 24-27
  { 36, 39, 1, 3},{ 36, 39, 1, 3},{ 38, 40, 0, 4},{ 41, 43, 5, 0}, // 28-31
  { 42, 45, 4, 1},{ 42, 45, 4, 1},{ 44, 47, 3, 2},{ 44, 47, 3, 2}, // 32-35
  { 46, 49, 2, 3},{ 46, 49, 2, 3},{ 48, 51, 1, 4},{ 48, 51, 1, 4}, // 36-39
  { 50, 52, 0, 5},{ 53, 43, 6, 0},{ 54, 57, 5, 1},{ 54, 57, 5, 1}, // 40-43
  { 56, 59, 4, 2},{ 56, 59, 4, 2},{ 58, 61, 3, 3},{ 58, 61, 3, 3}, // 44-47
  { 60, 63, 2, 4},{ 60, 63, 2, 4},{ 62, 65, 1, 5},{ 62, 65, 1, 5}, // 48-51
  { 50, 66, 0, 6},{ 67, 55, 7, 0},{ 68, 57, 6, 1},{ 68, 57, 6, 1}, // 52-55
  { 70, 73, 5, 2},{ 70, 73, 5, 2},{ 72, 75, 4, 3},{ 72, 75, 4, 3}, // 56-59
  { 74, 77, 3, 4},{ 74, 77, 3, 4},{ 76, 79, 2, 5},{ 76, 79, 2, 5}, // 60-63
  { 62, 81, 1, 6},{ 62, 81, 1, 6},{ 64, 82, 0, 7},{ 83, 69, 8, 0}, // 64-67
  { 84, 71, 7, 1},{ 84, 71, 7, 1},{ 86, 73, 6, 2},{ 86, 73, 6, 2}, // 68-71
  { 44, 59, 5, 3},{ 44, 59, 5, 3},{ 58, 61, 4, 4},{ 58, 61, 4, 4}, // 72-75
  { 60, 49, 3, 5},{ 60, 49, 3, 5},{ 76, 89, 2, 6},{ 76, 89, 2, 6}, // 76-79
  { 78, 91, 1, 7},{ 78, 91, 1, 7},{ 80, 92, 0, 8},{ 93, 69, 9, 0}, // 80-83
  { 94, 87, 8, 1},{ 94, 87, 8, 1},{ 96, 45, 7, 2},{ 96, 45, 7, 2}, // 84-87
  { 48, 99, 2, 7},{ 48, 99, 2, 7},{ 88,101, 1, 8},{ 88,101, 1, 8}, // 88-91
  { 80,102, 0, 9},{103, 69,10, 0},{104, 87, 9, 1},{104, 87, 9, 1}, // 92-95
  {106, 57, 8, 2},{106, 57, 8, 2},{ 62,109, 2, 8},{ 62,109, 2, 8}, // 96-99
  { 88,111, 1, 9},{ 88,111, 1, 9},{ 80,112, 0,10},{113, 85,11, 0}, // 100-103
  {114, 87,10, 1},{114, 87,10, 1},{116, 57, 9, 2},{116, 57, 9, 2}, // 104-107
  { 62,119, 2, 9},{ 62,119, 2, 9},{ 88,121, 1,10},{ 88,121, 1,10}, // 108-111
  { 90,122, 0,11},{123, 85,12, 0},{124, 97,11, 1},{124, 97,11, 1}, // 112-115
  {126, 57,10, 2},{126, 57,10, 2},{ 62,129, 2,10},{ 62,129, 2,10}, // 116-119
  { 98,131, 1,11},{ 98,131, 1,11},{ 90,132, 0,12},{133, 85,13, 0}, // 120-123
  {134, 97,12, 1},{134, 97,12, 1},{136, 57,11, 2},{136, 57,11, 2}, // 124-127
  { 62,139, 2,11},{ 62,139, 2,11},{ 98,141, 1,12},{ 98,141, 1,12}, // 128-131
  { 90,142, 0,13},{143, 95,14, 0},{144, 97,13, 1},{144, 97,13, 1}, // 132-135
  { 68, 57,12, 2},{ 68, 57,12, 2},{ 62, 81, 2,12},{ 62, 81, 2,12}, // 136-139
  { 98,147, 1,13},{ 98,147, 1,13},{100,148, 0,14},{149, 95,15, 0}, // 140-143
  {150,107,14, 1},{150,107,14, 1},{108,151, 1,14},{108,151, 1,14}, // 144-147
  {100,152, 0,15},{153, 95,16, 0},{154,107,15, 1},{108,155, 1,15}, // 148-151
  {100,156, 0,16},{157, 95,17, 0},{158,107,16, 1},{108,159, 1,16}, // 152-155
  {100,160, 0,17},{161,105,18, 0},{162,107,17, 1},{108,163, 1,17}, // 156-159
  {110,164, 0,18},{165,105,19, 0},{166,117,18, 1},{118,167, 1,18}, // 160-163
  {110,168, 0,19},{169,105,20, 0},{170,117,19, 1},{118,171, 1,19}, // 164-167
  {110,172, 0,20},{173,105,21, 0},{174,117,20, 1},{118,175, 1,20}, // 168-171
  {110,176, 0,21},{177,105,22, 0},{178,117,21, 1},{118,179, 1,21}, // 172-175
  {110,180, 0,22},{181,115,23, 0},{182,117,22, 1},{118,183, 1,22}, // 176-179
  {120,184, 0,23},{185,115,24, 0},{186,127,23, 1},{128,187, 1,23}, // 180-183
  {120,188, 0,24},{189,115,25, 0},{190,127,24, 1},{128,191, 1,24}, // 184-187
  {120,192, 0,25},{193,115,26, 0},{194,127,25, 1},{128,195, 1,25}, // 188-191
  {120,196, 0,26},{197,115,27, 0},{198,127,26, 1},{128,199, 1,26}, // 192-195
  {120,200, 0,27},{201,115,28, 0},{202,127,27, 1},{128,203, 1,27}, // 196-199
  {120,204, 0,28},{205,115,29, 0},{206,127,28, 1},{128,207, 1,28}, // 200-203
  {120,208, 0,29},{209,125,30, 0},{210,127,29, 1},{128,211, 1,29}, // 204-207
  {130,212, 0,30},{213,125,31, 0},{214,137,30, 1},{138,215, 1,30}, // 208-211
  {130,216, 0,31},{217,125,32, 0},{218,137,31, 1},{138,219, 1,31}, // 212-215
  {130,220, 0,32},{221,125,33, 0},{222,137,32, 1},{138,223, 1,32}, // 216-219
  {130,224, 0,33},{225,125,34, 0},{226,137,33, 1},{138,227, 1,33}, // 220-223
  {130,228, 0,34},{229,125,35, 0},{230,137,34, 1},{138,231, 1,34}, // 224-227
  {130,232, 0,35},{233,125,36, 0},{234,137,35, 1},{138,235, 1,35}, // 228-231
  {130,236, 0,36},{237,125,37, 0},{238,137,36, 1},{138,239, 1,36}, // 232-235
  {130,240, 0,37},{241,125,38, 0},{242,137,37, 1},{138,243, 1,37}, // 236-239
  {130,244, 0,38},{245,135,39, 0},{246,137,38, 1},{138,247, 1,38}, // 240-243
  {140,248, 0,39},{249,135,40, 0},{250, 69,39, 1},{ 80,251, 1,39}, // 244-247
  {140,252, 0,40},{249,135,41, 0},{250, 69,40, 1},{ 80,251, 1,40}, // 248-251
  {140,252, 0,41}};  // 252, 253-255 are reserved
*/

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
        /*printf("{%3d,%3d,%2d,%2d},", ns[state*4], ns[state*4+1], ns[state*4+2], ns[state*4+3]);
        if (state%4==3) printf(" // %d-%d\n", state-3, state);*/
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
U8 STA[6][256][4];
enum {
    ST1=0,
    ST2=1,
    ST4=2,
    ST5=3,
    ST6=4,
    ST7=5
};
// Squash returns p = 1/(1 + exp(-d)), d scaled by 8 bits, p scaled by 12 bits
short sqt[4095];

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

int stretchc(int p) {
    assert(p >= 0 && p <= 4095);
    if (p == 0)p = 1;
    float f = p / 4096.0f;
    float d = log(f / (1.0f - f)) * 256.0f;
    int di = (int)round(d);
    if (di > 2047)di = 2047;
    if (di < -2047)di = -2047;
    return di;
}

inline short stretch(int p) {
    return strt[p];
}

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

#if !defined(__GNUC__)

#if (2 == _M_IX86_FP) // 2 if /arch:SSE2 was used.
# define __SSE2__
#elif (1 == _M_IX86_FP) // 1 if /arch:SSE was used.
# define __SSE__
#endif

#endif /* __GNUC__ */

#if defined(__AVX2__)
#include <immintrin.h>
#define OPTIMIZE "AVX2-"
#elif defined(__SSE4_1__)   
#include<smmintrin.h>
#elif   defined(__SSSE3__)
#include<tmmintrin.h>
#elif defined(__SSE2__) 
#include <emmintrin.h>
#define OPTIMIZE "SSE2-"

#elif defined(__SSE__)
#include <xmmintrin.h>
#define OPTIMIZE "SSE-"
#endif

// Vector product a*b of n signed words, returning signed integer scaled down by 8 bits.
// n is rounded up to a multiple of 8.

//static int dot_product (const short* const t, const short* const w, int n);

// Train n neural network weights w[n] on inputs t[n] and err.
// w[i] += ((t[i]*2*err)+(1<<16))>>17 bounded to +- 32K.
// n is rounded up to a multiple of 8.

//static void train (const short* const t, short* const w, int n, const int e);

#if defined(__MMX__)
typedef __m128i XMM;
#endif

struct Mixer1 { 
  int N, M;   // max inputs, max contexts, max context sets
  short*tx; // N inputs from add()  
  short* wx ; // N*M weights
  short *ptr;
  int cxt;  // S contexts
  int pr;   // last result (scaled 12 bits)
  int shift1; 
  int elim;
  int uperr;
  int err;
#if defined(__AVX2__)
    static inline int horizontal_sum_i32(YMM sum) {
        XMM lo = _mm256_castsi256_si128(sum);
        XMM hi = _mm256_extractf128_si256(sum, 1);
        XMM s = _mm_add_epi32(lo, hi);
        s = _mm_add_epi32(s, _mm_shuffle_epi32(s, _MM_SHUFFLE(1, 0, 3, 2)));
        s = _mm_add_epi32(s, _mm_shuffle_epi32(s, _MM_SHUFFLE(2, 3, 0, 1)));
        return _mm_cvtsi128_si32(s);
    }

    int dot_product (const short* const t, const short* const w, int n) {
        assert(n > 0);
        YMM sum0 = _mm256_setzero_si256();
        YMM sum1 = _mm256_setzero_si256();
        YMM sum2 = _mm256_setzero_si256();
        YMM sum3 = _mm256_setzero_si256();
        int i = 0;
        for (; i + 63 < n; i += 64) {
            const YMM t0 = _mm256_load_si256((const YMM*)&t[i]);
            const YMM w0 = _mm256_load_si256((const YMM*)&w[i]);
            const YMM t1 = _mm256_load_si256((const YMM*)&t[i + 16]);
            const YMM w1 = _mm256_load_si256((const YMM*)&w[i + 16]);
            const YMM t2 = _mm256_load_si256((const YMM*)&t[i + 32]);
            const YMM w2 = _mm256_load_si256((const YMM*)&w[i + 32]);
            const YMM t3 = _mm256_load_si256((const YMM*)&t[i + 48]);
            const YMM w3 = _mm256_load_si256((const YMM*)&w[i + 48]);
            sum0 = _mm256_add_epi32(sum0, _mm256_srai_epi32(_mm256_madd_epi16(t0, w0), 8));
            sum1 = _mm256_add_epi32(sum1, _mm256_srai_epi32(_mm256_madd_epi16(t1, w1), 8));
            sum2 = _mm256_add_epi32(sum2, _mm256_srai_epi32(_mm256_madd_epi16(t2, w2), 8));
            sum3 = _mm256_add_epi32(sum3, _mm256_srai_epi32(_mm256_madd_epi16(t3, w3), 8));
        }
        for (; i < n; i += 16) {
            const YMM tv = _mm256_loadu_si256((const YMM*)&t[i]);
            const YMM wv = _mm256_loadu_si256((const YMM*)&w[i]);
            sum0 = _mm256_add_epi32(sum0, _mm256_srai_epi32(_mm256_madd_epi16(tv, wv), 8));
        }
        sum0 = _mm256_add_epi32(sum0, sum1);
        sum2 = _mm256_add_epi32(sum2, sum3);
        return horizontal_sum_i32(_mm256_add_epi32(sum0, sum2));
    }

 
 void train (const short* const t, short* const w, int n, const int e) {
  assert(n == ((n + 15) & -16));
  if (e) {
    const YMM one = _mm256_set1_epi16 (1);
    const YMM err = _mm256_set1_epi16 (short(e));
       int i = 0;
            // increment exactly as in the original: (((t*2*err)>>16)+1)>>1
        #define MIX_INC(idx) \
            _mm256_srai_epi16( \
                _mm256_adds_epi16( \
                    _mm256_mulhi_epi16( \
                        _mm256_adds_epi16( \
                            _mm256_load_si256((const YMM*)&t[idx]), \
                            _mm256_load_si256((const YMM*)&t[idx])), \
                        err), \
                    one), \
                1)
            // One 64-B block (32 int16 lanes) = one cache line per check.
            for (; i + 31 < n; i += 32) {
                const YMM d0 = MIX_INC(i);
                const YMM d1 = MIX_INC(i + 16);
                const YMM any = _mm256_or_si256(d0, d1);
                if (_mm256_testz_si256(any, any)) continue;   // row stays clean
                _mm256_store_si256((YMM*)&w[i],    _mm256_adds_epi16(d0, _mm256_load_si256((const YMM*)&w[i])));
                _mm256_store_si256((YMM*)&w[i+16], _mm256_adds_epi16(d1, _mm256_load_si256((const YMM*)&w[i+16])));
            }
            for (; i < n; i += 16) {
                const YMM d0 = MIX_INC(i);
                if (!_mm256_testz_si256(d0, d0)) _mm256_store_si256((YMM*)&w[i], _mm256_adds_epi16(d0, _mm256_load_si256((const YMM*)&w[i])));//MIX_ADDST(i, d0);
            }
  }
}

#else

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
#endif 

  // Adjust weights to minimize coding cost of last prediction
  void update(int y) {
       err=((y<<12)-pr)*uperr/4;
      if (err>32767)
          err=32767;
      if (err<-32768)
          err=-32768;
      if(err>=-elim && err<=elim) err=0;
      train(&tx[0], &wx[cxt*N], N, err);
  }
 
  // predict next bit
  int p( ) {
    assert(cxt<M);
    int dp=dot_product(&tx[0], &wx[cxt*N], N)*shift1>>11;
    return pr=squash(dp);
  }
    int p1( ) {
    assert(cxt<M);
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
  void setTxWx(int n,short* mn){
    N=n;
    alloc1(wx,(N*M)+32,ptr,32);
    tx=mn; 
  }
  void Init(int m,  U32 s,U32 e,U32 ue){
    M=m,  cxt=0, shift1=s,elim=e,uperr=ue;err=0;
    pr=2048; //initial p=0.5
  }
  void Free(){
    free(ptr);
  }
 /* void Print(){
      // print N weights averaged over context
  printf("Mixer(%d,%d): ", N, M);
  for (int i=0; i<N; ++i) {
    int w=0;
    for (int j=0; j<M; ++j)
      w+=wx[j*N+i],printf("%d ",wx[j*N+i]);;
    printf("%d ", w/M);
  }
  printf("\n");
  }*/
};

// A StateMap maps a context to a probability.  Methods:

// Statemap sm(n) creates a StateMap with n contexts using 4*n bytes memory.
// sm.p(y, cx, limit) converts state cx (0..n-1) to a probability (0..4095).
//     that the next y=1, updating the previous prediction with y (0..1).
//     limit (1..1023, default 1023) is the maximum count for computing a
//     prediction.  Larger values are better for stationary sources.

static int dt[1024];  // i -> 16K/(i+i+3)

struct StateMap {
  int N;        // Number of contexts
  int cxt;      // Context of last prediction
  U32 *t;       // cxt -> prediction in high 22 bits, count in low 10 bits
  int pr;
  int mask;
  int limit; 
  const U8 *nn;
  int next(int i, int y){
      return nn[ y + i*4];
  }
  void Init(int n, int lim,const U8 *nn1){nn=nn1;
    N=n, cxt=0, pr=2048, mask=n-1,limit=lim;
    assert(ispowerof2(n));
    alloc(t,n);
    assert(limit>0 && limit<1024);
    if (N==256){
        for (int i=0; i<N; ++i){
            U32 n0=next(i, 2)*3+1;
            U32 n1=next(i, 3)*3+1;
            t[i]=(((n1<<20) / (n0+n1)) << 12);
        }
    }else{
        for (int i=0; i<N; ++i)
            t[i]=1<<31;
    }
  }
  void Free(){
    free(t);
  }
  inline void update(const int y) {    
    assert(y==0 || y==1);
    U32 *p=&t[cxt], p0=p[0];
    int n=p0&1023, pr1=p0>>12;  // count, prediction
    p0+=(n<limit);
    p0+=(((((y<<20)-pr1)))*dt[n]+512)&0xfffffc00;
    p[0]=p0;
  }
  // update bit y (0..1), predict next bit in context cx
  void set(const int c,const int y) {  
    assert(cxt>=0 && cxt<N);
    update(y);
    pr=t[cxt=(c&mask)]>>20;
  } 
  /*void print(){
      for (int i=0;i<N;i++){
          printf("%d\n",t[i]>>20);
      }
  }*/
}; 

inline short clp(int z){
    if (z<-2047){
        z=-2047;
    }else if (z>2047){
        z=2047;
    }
    return z;
}
inline short clp1(int z){
    if (z<0){
        z=0;
    }else if (z>4095){
        z=4095;
    }
    return z;
}
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
  void Free(){
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
  int mix(int m) {  // return run length
    x.mxInputs[m].add(p());
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

struct SmallStationaryContextMap {
  U16 *Data;
  int Context, Mask, Stride, bCount, bTotal, B,N;
  U16 *cp;

  void Init(int BitsOfContext,  int InputBits = 8)     {
    assert(InputBits>0 && InputBits<=8);
    Context=0, Mask=((1<<BitsOfContext)-1), 
    Stride=((1<<InputBits)-1), bCount=(0), bTotal=(InputBits), B=(0)  ;
    N=(1ull<<BitsOfContext)*((1ull<<InputBits)-1);
    alloc(Data,N);
    for (int i=0; i<N; ++i)
      Data[i]=0x7FFF;
    cp=&Data[0];
  }
  void Free(){
    free(Data);
  }
  void set(U32 ctx) {
    Context = (ctx&Mask)*Stride;
    bCount=B=0;
  }
  void __attribute__ ((noinline)) mix(int m ) {
  const int rate = 7; const int Multiplier = 1;const int Divisor = 4;
    *cp+=((x.y<<16)-(*cp)+(1<<(rate-1)))>>rate;
    B+=(x.y && B>0);
    cp = &Data[Context+B];
    int Prediction = (*cp)>>4;
    x.mxInputs[m].add((stretch(Prediction)*Multiplier)/Divisor);
    x.mxInputs[m].add(((Prediction-2048)*Multiplier)/(Divisor*2));
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
/*struct E {  // hash element, 64 bytes
    U16 chk[7];  // byte context checksums
    U8 last;     // last 2 accesses (0-6) in low, high nibble
    U8 bh[7][7]; // byte context, 3-bit context -> bit history state
      // bh[][0] = 1st bit, bh[][1,2] = 2nd bit, bh[][3..6] = 3rd bit
      // bh[][0] is also a replacement priority, 0 = empty
  //  U8* get(U16 chk);  // Find element (0-6) matching checksum.
      // If not found, insert or replace lowest priority (not last).
      inline U8* get(U16 ch,int keep) {
    
  if (chk[last&15]==ch) return &bh[last&15][0];
  int b=0xffff, bi=0;
 
  for (int i=0; i<7; ++i) {
    if (chk[i]==ch) return last=last<<4|i, (U8*)&bh[i][0];
    int pri=bh[i][0];
    if (pri<b && (last&15)!=i && last>>4!=i) b=pri, bi=i;
  }
  return last=last<<4|bi|keep, chk[bi]=ch, (U8*)memset(&bh[bi][0], 0, 7);
}
    
};*/

inline U32 getStateByteLocation(const int bpos, const int c0) {
  U32 pis = 0; //state byte position in slot
  const U32 smask = (U32(0x31031010) >> (bpos << 2)) & 0x0F;
  pis = smask + (c0 & smask);
  return pis;
}

#define MAXCXT 7
/*
struct ContextMap {
  int C;  // max number of contexts
  U8* cp[MAXCXT];   // C pointers to current bit history
  U8* cp0[MAXCXT];  // First element of 7 element array containing cp[i]
  U32 cxt[MAXCXT];  // C whole byte contexts (hashes)
  U8* runp[MAXCXT]; // C [0..3] = count, value, unused, unused
  StateMap *sm;    // C maps of state -> p
  int cn;          // Next context to set by set()
  int result;
  short rc1[512];
  short st1[4096];
  short st2[4096];
  short st32[256];
  short st8[256]; 
  int cms,cms2,cms3,cms4;
  int kep;
  const U8 *nn;
  E *ptr,*t;
  U32 tmask;
  
  const inline U8  next(int i, int y){
      return nn[ y + i*4];
  }

  int __attribute__ ((noinline)) mix(const int m) {return mix1(m,  x.c0,  x.bpos, (U8) x.c4,  x.y);}
  inline int pre(const int state) {
    assert(state>=0 && state<256);
    U32 n0=next(state, 2)*3+1;
    U32 n1=next(state, 3)*3+1;
    return (n1<<12) / (n0+n1);
  }

// Construct using m bytes of memory for c contexts(c+7)&-8
void __attribute__ ((noinline)) Init(U32 m, int c, int s3,const U8 *nn1,int cs4,int k){
    C=c&255;
    tmask=((m>>6)-1); 
    cn=0;
    result=0;
    kep=k;
    alloc1(t,(m>>6)+64,ptr,64);  
    nn=nn1;        
    int cmul=(c>>8)&255;          // run context mul value
    cms=(c>>16)&255;              // mix prediction mul value
    cms2=(U32(c)>>24)&255;
    cms4=cs4;
    cms3=s3;
    assert(m>=64 && (m&m-1)==0);  // power of 2?
    assert(sizeof(E)==64);
    alloc(sm,C);
    for (int i=0; i<C; i++) 
        sm[i].Init(256,1023,nn1);
    for (int i=0; i<C; ++i) {
        cp0[i]=cp[i]=&t[0].bh[0][0];
        runp[i]=cp[i]+3;
    }
    // precalc int c=ilog(rc+1)<<(2+(~rc&1));
    for (int rc=0;rc<256;rc++) {
        int c=ilog[rc];
        c=c<<(2+(~rc&1));
        if ((rc&1)==0) c=c*cmul/4;
        rc1[rc+256]=clp(c);
        rc1[rc]=clp(-c);
    }

    // precalc mix3 mixer inputs
    for (int i=0;i<4096;i++) {
        st1[i]=clp(sc(cms*stretch(i)));
        st2[i]=clp(sc(cms2*(i - 2048)));
    } 

    for (int s=0;s<256;s++) {
        int n0=-!next(s,2);
        int n1=-!next(s,3);
        int r=0;
        int sp0=0;
        if ((n1-n0)==1 ) sp0=0,r=1;
        if ((n1-n0)==-1 ) sp0=4095,r=1;
        if (r) {
            st8[s] =clp(sc((cms4)*(pre(s)-sp0)));
            st32[s]=clp(sc((cms3)*stretch(pre(s))));
            if (s<8) st32[s]=0;
        }else{
            st8[s] =0;
            st32[s]=0;
        }
    }
}

void Free() {
    for (int i=0; i<C; i++) {
        sm[i].Free();
    }
    free(sm);
    free(ptr);
}

// Set the i'th context to cx
inline void set(U32 cx) {
  int i=cn++;
  assert(i>=0 && i<C);
  cx=cx*987654323+i;  // permute (don't hash) cx to spread the distribution
  cx=cx<<16|cx>>16;
  cxt[i]=cx*123456791+i;
}

// Predict to mixer m from bit history state s, using sm to map s to
// a probability.
inline int mix3(const int y,const int m,const int s, StateMap& sm) {
  if (s==0){
    x.mxInputs[m].add(0); 
    x.mxInputs[m].add(0);
    x.mxInputs[m].add(0);
    x.mxInputs[m].add(0);
    x.mxInputs[m].add(32*2);
    return 0;
  }else{
    sm.set(s,y);
    const int p1=sm.pr;
    x.mxInputs[m].add(st1[p1]); // From StateMap
    x.mxInputs[m].add(st2[p1]);
    x.mxInputs[m].add(st8[s]);  // From state
    x.mxInputs[m].add(st32[s]);
    x.mxInputs[m].add(0);
    return 1;
  }
}

// Zero prediction
inline void mix4(int m) {
    x.mxInputs[m].add(0); 
    x.mxInputs[m].add(0);
    x.mxInputs[m].add(0);
    x.mxInputs[m].add(0);
    x.mxInputs[m].add(32*2);
    x.mxInputs[m].add(0);
}

// Update the model with bit y1, and predict next bit to mixer m.
// Context: cc=c0, bp=bpos, c1=buf(1), y1=y.
int mix1(const int m,const int cc,const int bp,const int c1,const int y1) {
  // Update model with y
   result=0;

  for (int i=0; i<cn; ++i) {
    if (cp[i]) {
      assert(cp[i]>=&t[0].bh[0][0] && cp[i]<=&t[tmask].bh[6][6]);
      assert(((long long)(cp[i])&63)>=15);
      *cp[i]=next(*cp[i], y1);
    }

    // Update context pointers
    int s = 0;
    if (bp>1 && runp[i][0]==0) {
     cp[i]=0;
    } else {
     U16 chksum=(cxt[i]>>16)^i;
     
     if (bp){     
       if (bp==2 || bp==5)cp0[i]=cp[i]=t[(cxt[i]+cc)&tmask].get(chksum,kep);
       else cp[i]=cp0[i]+getStateByteLocation(bp,cc);
    } else {// default
       cp0[i]=cp[i]=t[(cxt[i]+cc)&tmask].get(chksum,kep);
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

    result=result+mix3(y1,m, s, sm[i]);
    // predict from last byte in context
    int b=x.c0shift_bpos ^ (runp[i][1] >> x.bposshift);
    if (b<=1) {
       b=b*256;   // predicted bit + for 1, - for 0
       // count*2, +1 if 2 different bytes seen
	   x.mxInputs[m].add(rc1[runp[i][0]+b]);
    }
    else
      x.mxInputs[m].add(0);
  }
  if (bp==7) cn=0;
  return result;
}
};*/

// Context map for large contexts.  Most modeling uses this type of context
// map.  It includes a built in RunContextMap to predict the last byte seen
// in the same context, and also bit-level contexts that map to a bit
// history state.
//
// Bit histories are stored in a hash table.  The table is organized into
// 32/128-byte buckets alinged on cache page boundaries.  Each bucket contains
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

//#define MAXCXT 10
short st2_p0[4096];
short st2_p1[4096];

template <int A, int B>
static inline int FindChk(const uint16_t* chk, uint16_t ch) {
  const __m128i target = _mm_set1_epi16((short)ch);
  const __m128i c0 = _mm_loadu_si128((const __m128i*)chk);
  unsigned m0 = (unsigned)_mm_movemask_epi8(_mm_cmpeq_epi16(c0, target));
  if (m0) return (int)(__builtin_ctz(m0) >> 1);
  const __m128i c1 = _mm_loadu_si128((const __m128i*)(chk + 8));
  unsigned m1 = (unsigned)_mm_movemask_epi8(_mm_cmpeq_epi16(c1, target));
  m1 &= (1u << ((A - 8) * 2)) - 1u;
  if (m1) return 8 + (int)(__builtin_ctz(m1) >> 1);
  return -1;
}

template <const int A, const int B> // Warning: values 3, 7 for A are the only valid parameters
union  E1 {  // hash element, 128 bytes
  struct{
    U16 chk[A];  // byte context checksums
    U16 last4;   // last 4 accesses (0-14) in low, high nibble
    U8 bh[A][7]; // byte context, 3-bit context -> bit history state
      // bh[][0] = 1st bit, bh[][1,2] = 2nd bit, bh[][3..6] = 3rd bit
      // bh[][0] is also a replacement priority, 0 = empty
      // If not found, insert or replace lowest priority (not last).
    };
    U8 pad[B];
    // Find element matching checksum.
    __attribute__ ((noinline)) U8* get(U16 ch) {
        const int recent0=last4&15;
        if (chk[recent0]==ch) return &bh[recent0][0];
        const int recent1=(last4>>4)&15;
        const int recent2=(last4>>8)&15;
        const int recent3=(last4>>12)&15;
        const int hit=FindChk<A,B>(chk, ch);
        if (hit>=0) return last4=last4<<4|hit, (U8*)&bh[hit][0];
        int b=0xffff, bi=0;

        for (int i=0; i<A; ++i) {
            if (chk[i]==ch) return last4=last4<<4|i, (U8*)&bh[i][0];
            if (recent0!=i && recent1!=i && recent2!=i && recent3!=i) {
               int pri=bh[i][0];
               if (pri<b) b=pri, bi=i;
            }
        }
        last4=last4<<4|bi;
        chk[bi]=ch;
        return (U8*)memset(&bh[bi][0], 0, 7);
   }
};

short st1in[4096];
short g_rc1[5][512];
// Changes from cmix-obias
// Sink byte for branch-free dead-context handling in ContextMap3/4::mix().
// A dead cp[i] is redirected here so the bit-history
// read-modify-write and the state load run unconditionally; next() on any
// byte value has no side effects and the sink is never meaningfully read.

static U8 cmSink8;
// ContextMap3 outputs 3-4 predictions to mixer.
struct ContextMap3 {
  int C;            // max number of contexts
  U8* cp[MAXCXT];   // C pointers to current bit history
  U8* cp0[MAXCXT];  // First element of 7 element array containing cp[i]
  U32 cxt[MAXCXT];  // C whole byte contexts (hashes)
  U8* runp[MAXCXT]; // C [0..3] = count, value, unused, unused
  int cn;           // Next context to set by set()
  int result;
  short *rc1;
  short *st1;
  short *st2;
  short st32[256];
  const U8 *nn;
  E1<14,128> *ptr,*t;  // Use 128-byte bucket
  U32 tmask;
  int skip2;
  U16 cxtMask;
  //State to prediction
  int cxtn[MAXCXT];
  U32 *ts;             // State prediction array
  int sti;             // Unique states count
  int resetcount;
  U64 cmtotal,cmset;

inline U8 next(int i, int y) {
    return nn[y + i*4];
}

int pre(const int state) {
    assert(state>=0 && state<256);
    U32 n0=next(state, 2)*3+1;
    U32 n1=next(state, 3)*3+1;
    return (n1<<12) / (n0+n1);
}
// StateMap
void updates(int i) {
    assert(i<C);
    assert(cxtn[i]<256);
    U32 *p=&ts[cxtn[i]], p0=p[0];
    int pr1=p0>>14;
    p0+=(x.y<<18)-pr1;
    p[0]=p0;
}

void upd() {
    assert(sti<=C);
    for (int i=0;i<sti;++i) {
        updates(i);
    }
} 
// Construct using m bytes of memory for c contexts(c+7)&-8
void __attribute__ ((noinline)) Init(U32 m1, int c, int cmul, int cms3,const int nn1,int cms4,int k,int u,short *st){
    C=c;//&255;
    U32 m=m1*2;
    tmask=((m>>7)-1); 
    cn=0;
    cxtMask=((1<C)-1)*2; // Initial zero contexts
    result=0;
    resetcount=0;
    cmtotal=0,cmset=0;
    alloc1(t,(m>>7)+64*2,ptr,64);  
    nn=&STA[nn1][0][0];
    skip2=u;
    assert(m>=64 && (m&m-1)==0);  // power of 2?
    assert(sizeof(E1<14,128>)==128);
    // Alloc state memory and prefill values
    alloc(ts,256);
    for (int i=0; i<256; ++i){
        U32 n0=next(i, 2)*3+1;
        U32 n1=next(i, 3)*3+1;
        ts[i]=(((n1<<20) / (n0+n1)) << 12);
    }
    for (int i=0; i<MAXCXT; ++i)cxtn[i]=0;
    sti=0;
    // Init cm
    for (int i=0; i<C; ++i) {
        cp0[i]=cp[i]=&t[0].bh[0][0];
        runp[i]=cp[i]+3;
    }
    assert(cmul>1 && cmul<7);
    rc1=&g_rc1[cmul-2][0]; // Set c=ilog(rc+1)<<(2+(~rc&1)) for cm run
    st2=st;
    st1=&st1in[0];
    // Precalc outputs to the mixer
    for (int s=0; s<256; s++) {
        int n0=-!next(s,2);
        int n1=-!next(s,3);
        int r=0;
        int sp0=0;
        if ((n1-n0)==1) sp0=0,r=1;
        if ((n1-n0)==-1) sp0=4095,r=1;
        if (r) {
            int st8 =clp(sc((cms4)*(pre(s)-sp0)));
            st32[s]=clp(sc((cms3)*stretch(pre(s))));
            if (s<8) st32[s]=st8;
            else st32[s]=(clp(st8+st32[s]))>>1;
        } else {
            st32[s]=0;
        }
    }
}
// Reset cm hashtable
void __attribute__ ((noinline)) reset() {
    memset((void*)t, 0, (tmask+1)*sizeof(E1<14,128>));
    for (int i=0; i<C; ++i) {
        cp0[i]=cp[i]=&t[0].bh[0][0];
        runp[i]=cp[i]+3;
    }
    resetcount++;
}

void freem() {
    free(ptr);
    free(ts);
}

// Set the i'th context to cx
inline void set(U32 cx) {
  int i=cn++;
  assert(i>=0 && i<C);
  cx=cx*987654323+i;  // permute (don't hash) cx to spread the distribution
  cx=cx<<16|cx>>16;
  cxt[i]=cx*123456791+i;
  // The bit-0 bucket t[(cxt[i]+c0)&tmask] with c0==1 is dereferenced by
        // mix() well after all byte-boundary set() calls; prefetching both
        // 64B lines of the 128B bucket here overlaps the ~60 otherwise
        // latency-serialized DRAM accesses across all maps.
  const char* bucket=(const char*)&t[(cxt[i]+1)&tmask];
  _mm_prefetch(bucket, _MM_HINT_T0);
  _mm_prefetch(bucket+64, _MM_HINT_T0);
  cxtMask=cxtMask*2;
  cmtotal++,cmset++;
}

// Skip context
inline void sets() {
  int i=cn++;
  assert(i>=0 && i<C);
  cxtMask=cxtMask+1; cxtMask=cxtMask*2;cmtotal++;
}

// Update the model with bit x.y, and predict the next bit to mixer m.
int __attribute__ ((noinline)) mix(int m) {
    result=0;
    upd();
    sti=0;
    for (int i=0; i<cn; ++i) {
        const int live=-(int)(((cxtMask>>(cn-i))&1)==0);
        int s=0;
        if (live) {
            // Branch-free bit-history update: a dead cp[i] is redirected
            // to a sink byte (next() has no side effects).
            assert(!cp[i] || (cp[i]>=&t[0].bh[0][0] && cp[i]<=&t[tmask].bh[14][6]));
            //assert(!cp[i] || ((long long)(cp[i])&127)>=29);
            U8* const wp = cp[i] ? cp[i] : &cmSink8;
            *wp=next(*wp, x.y);

            // Update context pointers
            if (x.bpos) {
                const int dead=(int)(x.bpos>1) & (int)(runp[i][0]==0);
                if (x.bpos==2 || x.bpos==5) {
                    // get() has eviction/LRU side effects: never run it
                    // for dead contexts, so this fork stays a branch.
                    if (dead) cp[i]=0;
                    else {
                        const U16 chksum=(cxt[i]>>16)^i;
                        cp0[i]=cp[i]=t[(cxt[i]+x.c0)&tmask].get(chksum);
                    }
                } else {
                    // Pure pointer bump at bpos 1,3,4,6,7: compute it
                    // unconditionally and cmov-select against nullptr.
                    U8* const bumped=cp0[i]+x.cmBitState;
                    cp[i] = dead ? (U8*)0 : bumped;
                    if ( x.bpos==1 ||  x.bpos==4) {
                        // Mid-byte prefetch for the bpos 2/5 re-hash:
                        // its bucket t[(cxt+c0')&tmask] needs c0'=2*c0+y
                        // with y still unknown, so prefetch both successor
                        // candidates (adjacent 128B elements -> up to 4
                        // consecutive lines). Contexts with a zero run
                        // count will skip the re-hash; their prefetch is
                        // steered to the always-hot t[0..1] instead.
                        const U32 gate = (U32)0 - (U32)(runp[i][0]!=0);
                        const U32 p0i = (cxt[i]+2*x.c0)   & tmask & gate;
                        const U32 p1i = (cxt[i]+2*x.c0+1) & tmask & gate;
                        const char* const pf0=(const char*)&t[p0i];
                        const char* const pf1=(const char*)&t[p1i];
                        _mm_prefetch(pf0,    _MM_HINT_T0);
                        _mm_prefetch(pf0+64, _MM_HINT_T0);
                        _mm_prefetch(pf1,    _MM_HINT_T0);
                        _mm_prefetch(pf1+64, _MM_HINT_T0);
                    }
                }
                // Branch-free state load: dead cp[i] reads the sink byte
                // and masks the value to 0.
                U8* const rp = cp[i] ? cp[i] : &cmSink8;
                s=(int)*rp & -(int)(cp[i]!=0);
            } else { // default
                const U16 chksum=(cxt[i]>>16)^i;
                cp0[i]=cp[i]=t[(cxt[i]+x.c0)&tmask].get(chksum);
                // Update pending bit histories for bits 2-7
                if (cp0[i][3]==2) {
                    const int c=cp0[i][4]+256;
                    U8 *p=t[(cxt[i]+(c>>6))&tmask].get(chksum);
                    p[0]=1+((c>>5)&1);
                    p[1+((c>>5)&1)]=1+((c>>4)&1);
                    p[3+((c>>4)&3)]=1+((c>>3)&1);
                    p=t[(cxt[i]+(c>>3))&tmask].get(chksum);
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
                s = *cp[i];
            }
        }
        // Branch-free s==0/s!=0 unification: st32[0]==0 by construction and
        // ts[0] is never registered (only s!=0 ever enters cxtn) so its
        // ts[0]>>20 stays a constant valid index; masking st1/st2 with
        // -(s!=0) reproduces the old add(0) stream exactly.
        //const int live = -(int)(s!=0);
        const int p1=(int)(ts[s]>>20);
        // Dedup over cxtn[0..sti) in one compare; lanes >= sti masked off.
        // Registration side effect matches old set(): advance only when
        // fresh AND s!=0, store is dead otherwise.
        const YMM vc=_mm256_set1_epi32(s);
        const YMM arr=_mm256_loadu_si256((const YMM*)&cxtn[0]);
        const U32 dup=(U32)_mm256_movemask_ps(_mm256_castsi256_ps(_mm256_cmpeq_epi32(arr, vc))) & ((1u<<sti)-1u);
        const int fresh=(int)(dup==0) & (int)(s!=0);
        result+=s!=0;
        cxtn[sti]=s;
        sti+=fresh;
        x.mxInputs[m].add(st1[p1]&live);
        if (skip2==1) x.mxInputs[m].add(st2[p1]&live);
        x.mxInputs[m].add(st32[s]&live);
        // predict from last byte in context
        const int b=x.c0shift_bpos ^ (runp[i][1] >> x.bposshift);
        // Branchless: b is non-negative, (b&1)<<8 keeps the index in
        // rcpr's 512 entries, the -(b<=1) mask kills the b>1 case 
        x.mxInputs[m].add(rc1[runp[i][0] + ((b&1)<<8)] & -(int)(b<=1)&live);
    }
    if (x.bpos==7) assert(cn==0 || cn==C),cn=cxtMask=0;
    return result;
}
};

template <const int A, const int B> // Warning: values 3, 7 for A are the only valid parameters
union  E4 {  // hash element, 64 bytes
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
      __attribute__ ((noinline)) U8* get(U16 ch, const int keep) {
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

short g_st9[6][256];
// ContextMap4 outputs 2-3 predictions to the mixer.
// No StateMap
struct ContextMap4 {
    int C;            // max number of contexts
    U8* cp[MAXCXT];   // C pointers to current bit history
    U8* cp0[MAXCXT];  // First element of 7 element array containing cp[i]
    U32 cxt[MAXCXT];  // C whole byte contexts (hashes)
    U8* runp[MAXCXT]; // C [0..3] = count, value, unused, unused
    int cn;           // Next context to set by set()
    int result;
    short *rc1;
    short *st9;
    short st32[256];
    int kep;
    const U8 *nn;
    E4<3,32> *ptr,*t;  // Use 128-byte bucket
    U32 tmask;
    int skip2;
    U16 cxtMask;
    int resetcount;
    U64 cmset,cmtotal;
    inline U8  next(int i, int y){
        return nn[ y + i*4];
    }

     int pre(const int state) {
        assert(state>=0 && state<256);
        U32 n0=next(state, 2)*3+1;
        U32 n1=next(state, 3)*3+1;
        return (n1<<12) / (n0+n1);
    }

    // Construct using m bytes of memory for c contexts(c+7)&-8
    void __attribute__ ((noinline)) Init(U32 m, int c,int cmul, int cms3,const int nn1, int cms4, int u){
        C=c;//&255;
        tmask=((m>>6)-1); 
        cn=0;
        cxtMask=((1<C)-1)*2; // Inital zero contexts
        result=0;
        alloc1(t,(m>>6)+64,ptr,64);  
        nn=&STA[nn1][0][0];;        
        skip2=u;
        assert(m>=64 && (m&m-1)==0);  // power of 2?
        assert(sizeof(E4<3,32>)==32);
        resetcount=0;
        cmtotal=0,cmset=0;
        // Init cm
        for (int i=0; i<C; ++i) {
            cp0[i]=cp[i]=&t[0].bh[0][0];
            runp[i]=cp[i]+3;
        }
        assert(cmul>1 && cmul<7);
        rc1=&g_rc1[cmul-2][0]; // Precalc c=ilog(rc+1)<<(2+(~rc&1)) for cm run
        st9=&g_st9[nn1][0];
        for (int s=0; s<256; s++) {
            int n0=-!next(s,2);
            int n1=-!next(s,3);
            int r=0;
            int sp0=0;
            if ((n1-n0)==1) sp0=0,r=1;
            if ((n1-n0)==-1) sp0=4095,r=1;
            if (r) {
                int st8=clp(sc((cms4)*(pre(s)-sp0)));
                st32[s]=clp(sc((cms3)*stretch(pre(s))));
                if (s<8) st32[s]=st8;
                else st32[s]=clp(st8+st32[s])>>1;
            } else {
                st32[s]=0;
            }
        }
    }
    // Reset cm hashtable
    void __attribute__ ((noinline)) reset() {
        memset((void*)t, 0, (tmask+1)*sizeof(E4<3,32>));
        for (int i=0; i<C; ++i) {
            cp0[i]=cp[i]=&t[0].bh[0][0];
            runp[i]=cp[i]+3;
        }
        resetcount++;
    }
    // Set the i'th context to cx
    inline void set(U32 cx) {
        int i=cn++;
        assert(i>=0 && i<C);
        cx=cx*987654323+i;  // permute (don't hash) cx to spread the distribution
        cx=cx<<16|cx>>16;
        cxt[i]=cx*123456791+i;
        // The bit-0 bucket t[(cxt[i]+c0)&tmask] with c0==1 is dereferenced by
        // mix() well after all byte-boundary set() calls; prefetching both
        // 64B lines of the 128B bucket here overlaps the ~60 otherwise
        // latency-serialized DRAM accesses across all maps.
        const char* bucket=(const char*)&t[(cxt[i]+1)&tmask];
        _mm_prefetch(bucket, _MM_HINT_T0);
        _mm_prefetch(bucket+64, _MM_HINT_T0);
        cxtMask=cxtMask*2;
        cmset++; cmtotal++;
    }
    // Skip context
    inline void sets() {
        int i=cn++;
        assert(i>=0 && i<C);
        cxtMask=cxtMask+1; 
        cxtMask=cxtMask*2;
        cmtotal++;
    }

    // Update the model with bit y1, and predict next bit to the mixer m.
    // Context: cc=c0, bp=bpos, c1=buf(1), y1=y.
    int __attribute__ ((noinline))  mix(int m) {
        result=0;
        for (int i=0; i<cn; ++i) {
            const int live=-(int)(((cxtMask>>(cn-i))&1)==0);
            int s=0;
            if (live) {
                // Branch-free bit-history update: a dead cp[i] is redirected
                // to a sink byte (next() has no side effects).
                assert(!cp[i] || (cp[i]>=&t[0].bh[0][0] && cp[i]<=&t[tmask].bh[14][6]));
                assert(!cp[i] || ((long long)(cp[i])&31)>=7);
                U8* const wp=cp[i]?cp[i]:&cmSink8;
                *wp=next(*wp, x.y);
                // Update context pointers
                if (x.bpos) {
                    const int dead=(int)(x.bpos>1) & (int)(runp[i][0]==0);
                    if (x.bpos==2 || x.bpos==5) {
                        // get() has eviction/LRU side effects: never run it
                        // for dead contexts, so this fork stays a branch.
                        if (dead) cp[i]=0;
                        else {
                            const U16 chksum=(cxt[i]>>16)^i;
                            cp0[i]=cp[i]=t[(cxt[i]+x.c0)&tmask].get(chksum,kep);
                        }
                    } else {
                        // Pure pointer bump at bpos 1,3,4,6,7: compute it
                        // unconditionally and cmov-select against nullptr.
                        U8* const bumped=cp0[i]+x.cmBitState;
                        cp[i] = dead ? (U8*)0 : bumped;
                        if (x.bpos==1 ||  x.bpos==4) {
                            // Mid-byte prefetch for the bpos 2/5 re-hash:
                            // its bucket t[(cxt+c0')&tmask] needs c0'=2*c0+y
                            // with y still unknown, so prefetch both successor
                            // candidates (adjacent 128B elements -> up to 4
                            // consecutive lines). Contexts with a zero run
                            // count will skip the re-hash; their prefetch is
                            // steered to the always-hot t[0..1] instead.
                            const U32 gate = (U32)0 - (U32)(runp[i][0]!=0);
                            const U32 p0i = (cxt[i]+2*x.c0)   & tmask & gate;
                            const U32 p1i = (cxt[i]+2*x.c0+1) & tmask & gate;
                            const char* const pf0=(const char*)&t[p0i];
                            const char* const pf1=(const char*)&t[p1i];
                            _mm_prefetch(pf0,    _MM_HINT_T0);
                            _mm_prefetch(pf0+64, _MM_HINT_T0);
                            _mm_prefetch(pf1,    _MM_HINT_T0);
                            _mm_prefetch(pf1+64, _MM_HINT_T0);
                        }
                    }
                    // Branch-free state load: dead cp[i] reads the sink byte
                    // and masks the value to 0.
                    U8* const rp = cp[i] ? cp[i] : &cmSink8;
                        s = (int)*rp & -(int)(cp[i]!=0);
                } else { // default
                    const U16 chksum=(cxt[i]>>16)^i;
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
                    s = *cp[i];
                }
            }
            // predict from bit context
            if (skip2==1)x.mxInputs[m].add(st9[s]&live);
            x.mxInputs[m].add(st32[s]&live);
            result+=s!=0;
            // predict from last byte in context
            // Branchless: b is non-negative, (b&1)<<8 keeps the index in
            // rc1's 512 entries, the -(b<=1) mask kills the b>1 case 
            const int b=x.c0shift_bpos^(runp[i][1]>>x.bposshift);
            x.mxInputs[m].add(rc1[runp[i][0] + ((b&1)<<8)] & -(int)(b<=1)&live);
        }
        if (x.bpos==7) assert(cn==0 || cn==C),cn=cxtMask=0;
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

struct  APM {
    int index;     // last p, context
    U16 *t;        // [N][33]:  p, context -> p

    int p(int pr=2048, int cxt=0, int rate=8, int y=0) {
        pr=stretch(pr);
        int g=(y<<16)+(y<<rate)-y*2;
        t[index]   += (g-t[index])   >> rate;
        t[index+1] += (g-t[index+1]) >> rate;
        const int w=pr&127;  // interpolation weight (33 points)
        index=((pr+2048)>>7)+cxt*33;
        return (t[index]*(128-w)+t[index+1]*w) >> 11;
    }

    // maps p, cxt -> p initially
    void Init(int n){
        index=0;
        alloc(t,n*33);
        for (int j=0; j<33; ++j) t[j]=squash((j-16)*128)*16;
        for (int i=33; i<n*33; ++i) t[i]=t[i-33];
    }
    void Free(){
        free(t);
   }
};

static const U8 wrt_w[256]={
2, 3, 1, 3, 3, 0, 1, 2, 3, 3, 0, 0, 1, 3, 3, 3,
3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0, 3, 3, 3, 3,
3, 2, 0, 2, 1, 3, 2, 1, 3, 3, 3, 3, 2, 3, 0, 2, // _!"#$%&'()*+,-./
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 0, 0, // 0123456789:;<=>?
2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, // @ABCDEFGHIJKLMNO
2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 2, 0, 2, 3, // PQRSTUVWXYZ[\]^_
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, // `abcdefghijklmno
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 1, // pqrstuvwxyz{|}~
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
};

static const U8 wrt_t[256]={
0, 0, 2, 0, 5, 6, 0, 6, 0, 2, 0, 4, 3, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
2, 4, 1, 4, 4, 7, 4, 7, 3, 7, 2, 2, 3, 5, 3, 1, // _!"#$%&'()*+,-./
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 6, 0, 7, 1, // 0123456789:;<=>?
0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // @ABCDEFGHIJKLMNO
4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 7, 4, 6, // PQRSTUVWXYZ[\]^_
5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, // `abcdefghijklmno
5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 3, 3, 7, 4, 5, // pqrstuvwxyz{|}~
6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
};


#define COLON         58 // :
#define SEMICOLON     59 // ;
#define LESSTHAN      60 // <
#define EQUALS        61 // =
#define GREATERTHAN   62 // >
#define QUESTION      63 // ?
#define ATSIGN        64 // @
#define SQUAREOPEN    91 // [
#define BACKSLASH     92 // '\'
#define SQUARECLOSE   93 // ]
#define CURLYOPENING 123 // {
#define VERTICALBAR  124 // |
#define CURLYCLOSE   125 // }

#define APOSTROPHE    39  // '
#define QUOTATION     34  // "
#define SPACE         32  // ' '

// Vector for different contexts
static const int charSize = 32;
template <typename T = int>
struct vec {
    T* cxt;
    int capacity;
    int size;
    int maxsize;
    char mo;
};

template <typename T = int>
void vec_new(vec<T>* o, char a){
    o->cxt=static_cast<T*>( calloc(charSize, sizeof(T)));
    o->capacity=charSize;
    o->size=0;
}
template <typename T = int>
void vec_free(vec<T>* o){
    free(o->cxt);
}
template <typename T = int>
int vec_size(vec<T> *o){
    return o->size;
}
template <typename T  = int>
void vec_push(struct vec<T> *o, const T element){
    if(o->size>0 && o->size%o->capacity==0) {
        o->capacity=(o->size/charSize+1)*charSize;
        o->cxt=static_cast<T*>(realloc(o->cxt, o->capacity*sizeof(T)));
    }
    o->cxt[o->size++]=element;
}
template <typename T  = int>
int vec_at(vec<T> *o, const int index){
    return o->cxt[index];
}
template <typename T  = int>
void vec_i(vec<T> *o, const int index){
    o->cxt[index]++;
}
template <typename T  = int>
void vec_pop(vec<T> *o){
    o->cxt[o->size-1]=0;
    o->size--;
}
template <typename T  = int>
void vec_reset(vec<T> *o){
    o->cxt[0]=0;
    o->size=0;
}
template <typename T  = int>
bool vec_empty(vec<T> *o){
    return (o->size==0)?true:false;
}
template <typename T  = int>
int vec_prev(vec<T> *o){
    return (o->size>1)?(o->cxt[o->size-2]):0;
}
// This part is based on cmix BracketContext
template <typename T  = U8>
struct BracketContext {
    U32 context;           // bracket byte and distance
    vec<int> active;            // vector for brackets
    vec<int> distance;          // vector for distance
    const T *element;           
    int elementCount;
    bool doPop;            // set true for quotes
    bool stateChange;
    int limit;
    T cxt,dst;

    void Init(const U8*d,const int e,int pop=false,int l=(1<< (sizeof(T)*8))-1) {
        elementCount=e;
        element=d;
        context =cxt=dst= 0;
        doPop=pop;
        limit=l;
        vec_new(&active,'b');
        vec_new(&distance,'b');
        stateChange=false;
    }
    void Reset(){
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
    bool FindEnd(int b,int c){
        bool found=false;
        for (int i=0;i<elementCount;i=i+2) if (element[i]==b&&element[i+1]==c) found=true;
        return found;
    }
    int last(){
        return vec_prev(&active);
    }
    void __attribute__ ((noinline)) Update(int byte) {
        bool pop=false;
        stateChange=false;
        if (!vec_empty(&active)) {
            if (FindEnd(vec_at(&active,vec_size(&active)-1) , byte) || vec_at(&distance,vec_size(&distance)-1) >= limit) {
                vec_pop(&active);
                vec_pop(&distance);
                pop=doPop;
                stateChange=true;
            } else {
                vec_i(&distance,vec_size(&distance)-1);
            }
        }
        if (pop==false && Find(byte)) {
            vec_push( &active,byte);
            vec_push( &distance,0);
            stateChange=true;
        }
        if (!vec_empty(&active)) {
            cxt=vec_at(&active,vec_size(&active)-1);
            dst=min(vec_at(&distance,vec_size(&distance)-1),(1<< (sizeof(T)*8))-1);
            context = (1<< ((sizeof(T)*8))) * cxt+dst;
        } else {
            context =cxt=dst= 0;
        }
    }
    void Free(){
        vec_free(&active);
        vec_free(&distance);
    }
};

// Table/row & column context
struct Column {
    U32 linepos;
    U8 fc;
    vec<U8> bytes;
};

struct ColumnContext {
    Column col[4];         // Content of last 3 + current row
    vec<U32> cell[4];      // Content of table row cell positions, max 4 rows
    int rows;
    int cellCount,cells,abovecellpos,abovecellpos1;
    bool NL;
    int limit;  // column lenght limit
    U8 nlChar;
    void Init( int l=31) {
        rows=abovecellpos=cellCount=abovecellpos1=0;
        nlChar=10;
        limit=l;
        for (int i=0;i<4;i++) vec_new(&col[i].bytes,'c');
        for (int i=0;i<4;i++) vec_new(&cell[i],'c');
        NL=false;
    }
   
    U8 lastfc(int i=0){
        return col[(rows-i)&3].fc;
    }
    bool isNewLine(){
        return NL;
    }
    int collen(int i=0,int l=0){
        return min((l?l:limit), vec_size(&col[(rows-i)&3].bytes)+1);
    }
    int nlpos(int i=0){
        return col[(rows-i)&3].linepos;
    }
    U8 colb(int i=1,int j=0,int l=0){
        if (collen(0,l)<collen(i,l))
        return  vec_at(&col[(rows-i)&3].bytes,collen()-(1+j));
        else return 0;
    }
    void __attribute__ ((noinline)) Update(int byte,int b2=0) {
        // Not working, extra CURLYOPENING
        // Start and end of table
        if ( b2==((CURLYOPENING<<16)+ (CURLYOPENING<<8)+ VERTICALBAR) ) nlChar='-';
        if ( b2==((VERTICALBAR<<16)+ (CURLYCLOSE<<8)+CURLYCLOSE ) ) nlChar=10,resetCells();
        // Column
        NL=false;
        if (byte==10){
            vec_push( &col[rows].bytes,U8(byte));
            rows++;
            rows=rows&3;
            vec_reset(&col[rows].bytes); // reset new line.
            col[rows].fc=0;
            col[rows].linepos=x.blpos-1;
        }else{
            vec_push( &col[rows].bytes,U8(byte)); // set new byte to line
            if (collen()==2) {
                col[rows].fc=min(byte,96);
                NL=true;
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
        if (nlChar=='-'){
            if ((b2&0xffff)==('-'+VERTICALBAR*256)){
                cells++;
                cells=cells&3;
                vec_reset(&cell[cells]); // reset new row.
                vec_push( &cell[cells],U32(x.blpos));
                cellCount=abovecellpos=abovecellpos1=0;
            }
            bool newcell=false;
            // Cells
            if ( (b2&0xffff)==(VERTICALBAR+VERTICALBAR*256) ||                // || 
             (b2&0xffff00)==((VERTICALBAR+10*256)*256) ||                     // \n|x
            ((b2&0xffff00)==((VERTICALBAR+10*256)*256) && byte!=VERTICALBAR)  // \n|yx  where y!=|
            ) vec_push( &cell[cells],U32(x.blpos)),cellCount++,newcell=true;
            // Advence above cell pos
            if (abovecellpos ) {
                abovecellpos++;
                // When above cell is shorter reset
                if (abovecellpos>abovecellpos1) abovecellpos=abovecellpos1=0;
            }
            // If more then one cell get above cell based on current row cell
            if(newcell==true && cellsCount() >0){
                // Get current above cell pos
                abovecellpos=cellPos(cellCount-1);
                abovecellpos1=cellPos(cellCount);
            }
        }
        /*#ifndef TEXTMODE
        if (nlChar==GREATERTHAN){
            if ((b2&0xffff)==(GREATERTHAN+10*256)){
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
            if ( (b2&0xff)==(GREATERTHAN) ) vec_push( &cell[cells],U32(x.blpos)),cellCount++,newcell=true;
            // Advence above cell pos
            if (abovecellpos ) {
                abovecellpos++;
                // When above cell is shorter reset
                if (abovecellpos>abovecellpos1) abovecellpos=abovecellpos1=0;
            }
            // If more then one cell get above cell based on current row cell
            if(newcell==true && cellsCount() >0){
                //printf("%d  ",cellPos(cellCount-1));
                // Get current above cell pos
                abovecellpos=cellPos(cellCount-1);
                abovecellpos1=cellPos(cellCount);
            }
            }
        }
        #endif*/
    }
    int cellsCount(int row=1){
        return vec_size(&cell[(cells-row)&3]);
    }
    int cellPos(int cellID,int row=1){
        int total=cellsCount(row)-1;
        total=min(total,cellID);
        return vec_at(&cell[(cells-row)&3],total);
    }
    void resetCells(){
        for (int i=0;i<4;i++) vec_reset(&cell[i]);
    }
    void Free(){
        for (int i=0;i<4;i++) vec_free(&col[i].bytes);
        for (int i=0;i<4;i++) vec_reset(&cell[i]);
    }
};
// Keep track of main brackets
const U8 brackets[8]={'(',')', CURLYOPENING,CURLYCLOSE, '[',']', LESSTHAN,GREATERTHAN};
// Keep track of ' and " as quotes
const U8 quotes[4]={APOSTROPHE,APOSTROPHE,QUOTATION,QUOTATION};
// Keep track of first char including some brackets
const U8 fchar[22]={ATSIGN,10, 96,10, COLON,10, LESSTHAN,GREATERTHAN,EQUALS,10,SQUAREOPEN,SQUARECLOSE,CURLYOPENING,CURLYCLOSE,'*',10,VERTICALBAR,10,31,10, '&',';'};
//ATSIGN,10 dickens not good
// Sentence & words context
struct WordsContext {
    vec<U32> awords; // List of words
    vec<U16> sbytes; // List of bytes surrounded by a current word
    U32 fword;       // First word
    U8 pbyte;        // Current byte before word
    int wordcount;
    void Init() {
        vec_new(&awords,'w');
        vec_new(&sbytes,'w');
        wordcount=0;
    }
    void Reset(){
        vec_reset(&awords);
        vec_reset(&sbytes);
        fword=pbyte=0;
        wordcount=0;
    }
    void Set(U8 b){
        pbyte=b;
    }
    void Update(U32 w,U8 b) {
        if (fword==0) fword=w;
        vec_push(&awords,w);
        vec_push(&sbytes,U16(pbyte*256+b));  // Surrounding bytes
        pbyte=0;
        wordcount++;
    }
    void Remove(){
        int num=vec_size(&awords);
        if (num) vec_pop(&awords),vec_pop(&sbytes);
        if (wordcount)wordcount--;
    }
    U32 Word(int i=1){
        int num=vec_size(&awords);
        if (num>=i) return vec_at(&awords,num-(i));
        else return 0;
    }
    U16 sBytes(int i=1){
        int num=vec_size(&sbytes);
        if (num>=i) return vec_at(&sbytes,num-(i));
        else return 0;
    }
    void Free(){
        vec_free(&awords);
        vec_free(&sbytes);
    }
};

inline U32 hash(U32 a, U32 b, U32 c=0xffffffff) {
    U32 h=a*110002499u+b*30005491u+c*50004239u; 
    return h^h>>9^a>>3^b>>3^c>>4;
}

// Predictor

static const U32 primes[14]={0, 257,251,241,239,233,229,227,223,211,199,197,193,191};
static const U32 tri[4]={0,4,3,7}, trj[4]={0,6,6,12};

// Parameters
// These parameters were tuned befor version 1 and may be bad.
const U32 m_e[10]={8,8,8,1,1,1,1,1,1,0}; // mixer error
const U32 m_s[10]={194, 237, 204, 70, 54, 55,55, 70,55, 6};// mixer shift
const U32 m_m[10]={36,   69,  19, 34, 23, 24,24, 34,24,4};// mixer error mul

const U32 c_r[27]= { 3,  4,  6,  4,  6,  6,  2,  3,  3,  3,  6,  4,  3,  4,  5,  6,  2,  6,  4,  4,  4,  4,  4,  4,  4,  4,  4};  // contextmap run mul
const U32 c_s3[27]={43, 33, 34, 28, 34, 29, 32, 33, 37, 35, 33, 28, 31, 35, 28, 30, 33, 34, 32, 32, 32, 32, 32, 32, 32, 32, 32};
const U32 c_s4[27]={ 9,  8,  9,  5,  8, 12, 15,  8,  8, 12, 10,  7,  7,  8, 8, 13, 13, 14,  8,  8, 12, 12, 12, 12, 12, 12, 12};
 
int e_l[8]={1830, 1997, 1973, 1851, 1897, 1690, 1998, 1842};

const int MAXLEN=62; // longest allowed match + 1
U32 t[14],told[14]; 
U32 oldx4;

int c1,c2,c3;
U32 word0,word1,word2,word3,words,spaces,wshift,w4,w4r,w4br,x4,x5,ismatch,firstWord,number0,number1,numbers,numlen0,numlen1,mybenum;
U32 lastword0,lineSP;
U32 fqcxt=0;
U32 AH1=0,AH2=0x765BA55C;
int rate=6;
U32 fails=0, failz=0, failcount=0;
int nl,col,fc,fc1,nl1;
U32 t1[0x100];
U32 t2[0x10000];
int wp[0x10000];

U32 oState, wtype, ttype;
U32 nState;
U32 oStatew4,nStatew4; 

int ord,ord2;
U8 buffer[0x1000000];
enum {BMASK=0xffffff};
int pos;
int pr; 
bool isPage=false;
bool isPDFText=false;

StateMap smA[3];
SmallStationaryContextMap scmA[8];   
Mixer1 mxA[10]; 
Mixer1 mxAX[9]; 
ContextMap3 cmC3[16];
ContextMap3 cmCS[6];
ContextMap4 cmC4[5];
APM apmA[6];
RunContextMap rcmA[1];  
BracketContext<U8> brcxt;
BracketContext<U8> qocxt;
BracketContext<U8> fccxt;
ColumnContext colcxt;
WordsContext worcxt;

void PredictorInit() { 
    nState=nStatew4=0xffffffff;
    pr=2048;
    smA[0].Init(1<<9,1023,&STA[ST1][0][0]);//match
    smA[1].Init(1<<19,1023,&STA[ST1][0][0]);//match
    smA[2].Init(1<<16,1023,&STA[ST1][0][0]);//match

    for (int i=0;i<8;i++){ 
        scmA[i].Init(8); 
    }
    // Layer 0
    mxA[1].Init(2*256,m_s[1],m_e[1],m_m[1]);
    mxA[2].Init(6*256,m_s[2],m_e[2],m_m[2]);
    mxA[3].Init(6*256,m_s[3],m_e[3],m_m[3]);
    mxA[4].Init(8*256,m_s[4],m_e[4],m_m[4]);
    mxA[5].Init(6*256,m_s[5],m_e[5],m_m[5]);
    mxA[6].Init(7*256*4,m_s[6],m_e[6],m_m[6]);
    mxA[7].Init(8*256,m_s[7],m_e[7],m_m[7]);
    mxA[8].Init(8*256,m_s[8],m_e[8],m_m[8]);
    mxA[9].Init(1,m_s[9],m_e[9],m_m[9]);       // Layer 2 - final
    // Layer 1
    mxAX[1].Init(2*256,8,m_e[1],13);
    mxAX[2].Init(6*256,8,m_e[2],13);
    mxAX[3].Init(6*256,8,m_e[3],13);
    mxAX[4].Init(8*256,8,m_e[4],13);
    mxAX[5].Init(6*256,8,m_e[5],13);
    mxAX[6].Init(7*256*4,8,m_e[6],13);
    mxAX[7].Init(8*256,8,m_e[7],13);
    mxAX[8].Init(8*256,8,m_e[8],13);
    
    apmA[0].Init(256);
    apmA[1].Init(0x8000*2);
    apmA[2].Init(0x8000*2);
    apmA[3].Init(0x20000*2);
    apmA[4].Init(0x10000*2);
    apmA[5].Init(0x10000*2);
    rcmA[0].Init(1*4096*4096,6);

    x.mxInputs[0].ncount=254;
    x.mxInputs[1].ncount=9;
    x.mxInputs[2].ncount=9;
    
    for (int j=0;j<3;j++)  {
        x.mxInputs[j].ncount=(x.mxInputs[j].ncount+16)&-16;
        alloc1(x.mxInputs[j].n,x.mxInputs[j].ncount+32,x.mxInputs[j].ptr,32);
    }     
    // Provide inputs array info to mixers
    for (int i=1;i<9;i++)
         mxA[i].setTxWx(x.mxInputs[0].ncount,&x.mxInputs[0].n[0]);
    // Final mixer
    for (int i=1;i<9;i++)
         mxAX[i].setTxWx(x.mxInputs[1].ncount,&x.mxInputs[1].n[0]);
    mxA[9].setTxWx(x.mxInputs[2].ncount,&x.mxInputs[2].n[0]);
    cmC3[0].Init( 8*4096*4096,3,c_r[0],c_s3[0],ST6,c_s4[0],0xf0,1,&st2_p1[0]);
    cmC3[1].Init( 8*4096*4096,1,c_r[1],c_s3[1],ST6,c_s4[1],0xf0,0,&st2_p1[0]);
    cmC3[2].Init( 8*4096*4096,1,c_r[2],c_s3[2],ST6,c_s4[2],0xf0,0,&st2_p1[0]);
    cmC3[3].Init( 8*4096*4096,1,c_r[3],c_s3[3],ST6,c_s4[3],0xf0,0,&st2_p1[0]);
    cmC3[4].Init(8*4096*4096,2,c_r[4],c_s3[4],ST6,c_s4[4],0xf0,1,&st2_p1[0]);
    cmC3[5].Init(4*4096*4096,2,c_r[5],c_s3[5],ST6,c_s4[5],0xf0,0,&st2_p1[0]);
    cmC3[6].Init(     64*4096,1,c_r[6],c_s3[6],ST1,c_s4[6],0,1,&st2_p1[0]);
    cmC3[7].Init( 1*4096*4096,1,c_r[7],c_s3[7],ST5,c_s4[7],0xf0,1,&st2_p1[0]);
    cmC3[8].Init( 1*4096*4096/2,3,c_r[8],c_s3[8],ST4,c_s4[8],0,0,&st2_p1[0]);
    
    cmC4[0].Init(     32*4096,2,c_r[9],c_s3[9],ST6,c_s4[9],0); //
    cmC4[1].Init(   8*32*4096,4,c_r[10],c_s3[10],ST7,c_s4[10],1);//
    cmC4[2].Init(   8*32*4096,4,c_r[11],c_s3[11],ST2,c_s4[11],1);//
    cmC4[3].Init(   2*16*4096,6,c_r[12],c_s3[12],ST7,c_s4[12],1);//
    cmC4[4].Init(    128*4096,2,c_r[16],c_s3[16],ST1,c_s4[16],0);//
    
    cmCS[0].Init(    4*16*4096,7,c_r[13],c_s3[13],ST2,c_s4[13],0,1,&st2_p1[0]);
    cmCS[1].Init(    64*2*4096,3,c_r[14],c_s3[14],ST5,c_s4[14],0xf0,0,&st2_p0[0]);
    cmCS[2].Init(       2*4096,2,c_r[15],c_s3[15],ST2,c_s4[15],0xf0,0,&st2_p0[0]);
    
    cmCS[3].Init(       32*4096,2,c_r[22],c_s3[22],ST2,c_s4[22],0x00,1,&st2_p1[0]);
    cmCS[4].Init(      512*4096,1,c_r[25],c_s3[25],ST1,c_s4[25],0xf0,1,&st2_p1[0]);
    cmCS[5].Init(      512*4096,1,c_r[26],c_s3[26],ST1,c_s4[26],0xf0,1,&st2_p1[0]);
    
    cmC3[9].Init( 2*4096*4096,4,c_r[17],c_s3[17],ST6,c_s4[17],0xf0,0,&st2_p1[0]);
    cmC3[10].Init(2*4096*4096,6,c_r[18],c_s3[18],ST5,c_s4[18],0xf0,0,&st2_p1[0]);
    cmC3[11].Init(2*4096*4096,6,c_r[19],c_s3[19],ST5,c_s4[19],0xf0,0,&st2_p1[0]);
    cmC3[12].Init(1*4096*4096/2,2,c_r[20],c_s3[20],ST6,c_s4[20],0xf0,0,&st2_p1[0]);
    cmC3[13].Init(1*4096*4096/2,2,c_r[21],c_s3[21],ST6,c_s4[21],0xf0,1,&st2_p1[0]);
    cmC3[14].Init( 1*4096*4096,1,c_r[23],c_s3[23],ST6,c_s4[23],0xf0,1,&st2_p1[0]);
    cmC3[15].Init(   8*64*4096,1,c_r[24],c_s3[24],ST1,c_s4[24],0,0,&st2_p0[0]);
    
    brcxt.Init(&brackets[0],8,false,512);
    qocxt.Init(&quotes[0],4,true);
    fccxt.Init(&fchar[0],22,false,512*2);
    colcxt.Init();
    worcxt.Init();
}
  
void PredictorFree(){
    smA[0].Free(); 
    for (int i=0;i<8;i++) scmA[i].Free();
    for (int i=1;i<10;i++) mxA[i].Free();
    //for (int i=0;i<27;i++) cmC[i].Free();
    for (int i=0;i<6;i++) apmA[i].Free();
    rcmA[1].Free();
    free(x.mxInputs[0].ptr);
    free(x.mxInputs[1].ptr);    
    brcxt.Free();
    qocxt.Free();
    fccxt.Free();
    colcxt.Free();
    worcxt.Free();
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
  #define mHashN   3
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
    void Init(){
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
        (length != 0) << 31 | //normal mode (match)
        (delta) << 30 | //delta mode
        (delta ? (lengthBak>>1) : (length>>1)) << 24 | //the longer wins, halve
        (index&0x00ffffff); //the more recent wins
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
          }
          else { //backup match information: maybe we can recover it just after this mismatch
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
          }
          else { // still mismatch
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

bool isMatch(const U32 pos, const int MINLEN) {
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
      if (isMatch(matchpos, LEN)) {
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

    hash = x4;//t[LEN1];
    matches = &mhashtable[(hash& mhashtablemask)];
    if (numberOfActiveCandidates < matchN)
      AddCandidates(matches, LEN1); //shortest
    matches->Add(pos);

    for (U32 i = 0; i < numberOfActiveCandidates; i++) {
      matchCandidates[i].expectedByte = bufr(matchCandidates[i].index);
    }
  }
}

int MatchModel2mix(int m) {
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
    x.mxInputs[m].add(sign * (length << 5));
  } else { // no match at all or delta mode
    x.mxInputs[m].add(0);
  }

  if( isInDeltaMode ) { // delta mode: helps predicting the remaining bits of a character when a mismatch occurs
    ctx[2] = (expectedByte << 8) | x.c0;
  }

  for( int i = 0; i < nST; i++ ) {
    const U32 c = ctx[i];
    if( c != 0 ) {
         smA[i].set(c,x.y);
      const int p1 = smA[i].pr;
      const int st = stretch(p1);
      x.mxInputs[m].add(st >> 2);
      x.mxInputs[m].add((p1 - 2048) >> 3);
    } else {
      x.mxInputs[m].add(0);
      x.mxInputs[m].add(0);
    }
  }
  return length;
}

// Find bracket or quotes index.
U8 qy[6]={'(',CURLYOPENING, '[', LESSTHAN,QUOTATION,APOSTROPHE};
int FindQy(int b){
    int found=-1;
    for (int i=0;i<6;i++) {
        if (qy[i]==b) {
            found=i;
            break;
        }
    }
    return found;
}
// In first char context find
U8 qy1[10]= {96, COLON, LESSTHAN,EQUALS,'*',VERTICALBAR,CURLYOPENING,SQUAREOPEN,31,'&'};
int FindQy2(int b){
    int found=-1;
    for (int i=0;i<10;i++) {
        if (qy1[i]==b) {
            found=i;
            break;
        }
    }
    if (found==8) found=0; // 31 as 96
    if (found==9) found=0; //  & as 96
    return min(7,found+1);
}

int modelPrediction(int c0,int bpos,int c4) {
    int i,c;
    U32 h,j;
    
    if (bpos==0) {
        wshift=0;
        c3=c2;
        c2=c1;
        c1=c4&0xff;
        
        i=wrt_w[(c1)];
        nStatew4=i;
        w4=w4*4+i;
        buffer[pos&BMASK]=c1;
        pos++;
        j=c1;
        // PDF BT-ET block
        if (c4==0x0a42540a) {        // '\nBT\n'
            isPDFText=true;
        } else if (c4==0x0a45540a) { // '\nET\n'
            isPDFText=false;
            brcxt.Reset();
            fccxt.Reset();
            qocxt.Reset();
        }
        bool isUpper=false;
        bool isLetter=false;
        if ( ((j-'A') <= ('Z'-'A'))) j =j+32,isUpper=true;
        // Column content update
        colcxt.Update(c1,c4&0xffffff);
        // Bracket content update
        if (fccxt.cxt=='&' && !((j-'a') <= ('z'-'a'))) fccxt.Update(';'); //html entity
        if (c4==0x266C743B) brcxt.Update(LESSTHAN);             // '&lt;'
        else if (c4==0x2667743B) brcxt.Update(GREATERTHAN);     // '&gt;'
        else if (c1<128) brcxt.Update( c1 );

        cmCS[4].set((brcxt.context<<8)+c1);
        // Quote content update
        if ((x.c8&0xffffffffffff)==0x2671756F743B) qocxt.Update(QUOTATION); // '&quot;'
        else qocxt.Update(c1); 
        // Look for byte stream x4 and order X context end.
        if (isPDFText==false && (c1=='$' || c1==VERTICALBAR || c1==')')) {
                x4=(x4<<8)+c2;
        }
        
        if (isPDFText==false) {
            // Update byte bitstreams (enwik)
            if (c1!=SQUAREOPEN && c1!=SQUARECLOSE) x4=(x4<<8)+c1;
            if (c1!=SQUAREOPEN &&  (c1=='$' || c1==SQUARECLOSE || c1==VERTICALBAR ||
                c1==APOSTROPHE || c1==QUOTATION || c1==')' || c1=='(')) {
                w4=w4*4+nStatew4;
                w4r=(w4r<<2)+nStatew4;
               ttype=(ttype<<3)+wrt_t[(c1)];
               wtype=(wtype<<3)+wrt_t[(c1)];
            }
        } else { // PDF
            if (brcxt.cxt=='(' && c1!='(')
            x4=(x4<<8)+c1;
        }
        
        if (isPDFText && brcxt.cxt=='(' && c1!='(') {
            for (i=13; i>0; --i)
                told[i]=told[i-1]*primes[i]+c1+i*256;  
        } 
        if ((x.c4&0xffff)!=(32*256+32)) {
            for (i=13; i>0; --i)
                t[i]=t[i-1]*primes[i]+c1+i*256;
        }
        if (isPDFText && brcxt.cxt=='(' && c1!='(') {
            for (i=3; i<6; ++i)
                cmC3[0].set(told[i]);
            cmC3[1].set(told[6]);
            cmC3[2].set(told[8]);
        } else if (isPDFText && brcxt.cxt==0 ) {
            cmC3[0].sets(); cmC3[0].sets(); cmC3[0].sets();cmC3[1].sets();cmC3[2].sets();
        } else if ((fc==SPACE && c1==SPACE)) {
            cmC3[0].sets(); cmC3[0].sets(); cmC3[0].sets();cmC3[1].sets();cmC3[2].sets();
        } else {
            for (i=3; i<6; ++i)
                cmC3[0].set(t[i]);
            cmC3[1].set(t[6]);
            cmC3[2].set(t[8]);
        }
        cmC3[3].set(t[13]);

        nState=wrt_t[(c1)];
        words=words<<1;
        spaces=spaces<<1;
        numbers=numbers<<1;

        if (isPDFText && c1=='(' && number0>300) word0=lastword0=0,spaces|=1,words|=2;  // PDF - A number greater than 300 is probably a space between words.
        if (((j-'a') <= ('z'-'a')) || (c1>127 && c2!=12) || (isPDFText && brcxt.cxt=='(' && c1!='(')) {
            if (isPDFText==false && ((words&9)==8 && c2==c3 && c2==SQUARECLOSE)) { // wiki
                    word0=lastword0;
                    worcxt.Remove();
            }
            if (isPDFText && brcxt.cxt=='(' && lastword0) word0=lastword0,words|=2; // PDF
            
            lastword0=0;    
            if (word0==0){
                if (c2==ATSIGN)worcxt.Set(c3); else worcxt.Set(c2);
            }
            
            isLetter=true;
            words=words|1;
            word0=word0*2104+j; //263*8
            if (c4==0x70616765 && buf(5)=='/') isPage=true; // '/page'
            // ' or 0x27 is used for quotes, sometimes it has other meaning
            // remove quote content if any of the fallowing is true
            const int word3bit=(words&7);
            if (((word3bit==5) && (c2==APOSTROPHE))||                    // "x'x" where x is any letter in word
             ((word3bit==1) && (c3==SQUARECLOSE) && (c2==APOSTROPHE))||  // "]'x" where x is any letter in word
             ((word3bit==1) && (numbers&4)&&(c2==APOSTROPHE))            // "y'x" where y is number and x is any letter in word //(c3>='0' && c3<='9') 
             ) if (isPDFText==false) qocxt.Update(qocxt.cxt); 
        } else {
            // ' or 0x27 is used for quotes, sometimes it has other meaning
            // remove quote content if any of the fallowing is true
            const int word3bit=(words&7);
            if ( ((word3bit==4)&& (c1==SPACE) && (c2==APOSTROPHE))||            // "x' " where x is any letter in word
             ((c1==ATSIGN) && (numbers&4) && (c2==APOSTROPHE)) ||               // "x'@" where x is number               // this is somehow semi good-bad //(c3>='0' && c3<='9')
             ((word3bit==4) && (c1==ATSIGN) && (c2==APOSTROPHE))                // "x'@" where x is any letter in word   // this is somehow semi good-bad
              )qocxt.Update(qocxt.cxt); 
            if (word0){
                word3=word2*47;
                word2=word1*53;
                word1=word0*83;
                // Update word content
                worcxt.Update(word0,c1);
                if(firstWord==0) {
                    firstWord=word0;
                }
            }
            if (!(isPage && c2=='e' && c1==GREATERTHAN)) isPage=false;
            //enwik page reset
            if (isPage) {
                cmC4[1].reset();
                cmC4[2].reset();
                cmC4[4].reset();
                cmCS[0].reset();
                cmCS[2].reset();
                isPage=false;
            }
            //
            wp[word0&0xffff]=pos;
            if (c1=='+' || (c1==SQUARECLOSE && lastword0==0)) lastword0=word0;
            if (isPDFText && brcxt.stateChange && c1==')') lastword0=word0;
            word0=0;
            if (c1>='0' && c1<='9') {
                numbers=numbers+1;
                if(numbers&4 && c2==',') number0=number1,number1=0,numlen0=numlen1,numlen1=0;
                if (mybenum && numlen1<=2) number0=number1,number1=0,numlen0=numlen1,numlen1=0;
                number0=number0*10+(c1&0x0f);
                numlen0=min(19,numlen0+1);mybenum=0;
            } else {
                if (numlen0 ||((numbers&0xf)==0)) {
                    number1=number0,numlen1=numlen0,number0=numlen0=0;
                }
                if (numlen1<=2 &&numlen1&&((numbers&5)==5) && numlen0==0 && c2=='.') mybenum=2;
                else if (numlen1<=2&&numlen1&&(numbers&2) && numlen0==0 && c1=='.') mybenum=1;
                else if (mybenum==1  && c1!='.') mybenum=0;
            }
            
            if (c1==SPACE) {
                spaces++;
            }
            else if (c1==10 || (c1==CURLYOPENING  && isPDFText)) {
                worcxt.Reset();
                fc=fc1=firstWord=0;
                nl1=nl;
                nl=pos-1;
                wtype=(wtype<<3);
                w4=w4|0x3fc;
                words=0xfffffffc;
                w4r= w4r<<2;
            }
            else if ((c1=='.') || (isPDFText==false && c1==')') || c1==QUESTION) {
                wtype= wtype&0xffffffc0;
                ttype= ttype<<7;//&0xffffffc0;  ????
                words= words|0xfffffffe;
                x5=(x5<<8)+(c4&0xff);
                w4=w4|204;
                w4r= w4r&0xffffffc0;
                if (c1=='.') {
                    // Skip sentence reset if following is true
                    if (!(fccxt.cxt==SQUAREOPEN || fccxt.cxt=='(' || colcxt.lastfc()=='*')) {
                        if (worcxt.wordcount) {
                            worcxt.Reset(); 
                        }
                    }
                }
            }
            else if (c1==',') {
                words=words|0xfffffffc;
            }   
            else if (c1==COLON) {
                ttype= (ttype&0xfffffff8)+4;
                w4=w4|12;//      1100
                x5=(x5<<8)+(c4&0xff);
            }   
            else if (c1==CURLYCLOSE || c1==CURLYOPENING) {
                words = words | 0xfffffffc;
                wtype= wtype&0xffffffc0;
                x5=(x5<<8)+(c4&0xff);
                ttype= (ttype&0xfffffff8)+3;
            }   
            else if (c1==EQUALS) {
                ttype=(ttype&0xfffffff8)+4;
            }  
            if (x.c8==0x6D703B6E6273703B)  {//  'mp;nbsp;'
                spaces|=1;
            }
            if (c1=='.' ) wshift=1,worcxt.Reset();
            else if (c1==SEMICOLON)  worcxt.Reset();
        }
        if (isPDFText && brcxt.stateChange) ;
        else x5=(x5<<8)+(c4&0xff);
        // Switch state if it is new
        if (oStatew4!=nStatew4){
            w4r=(w4r<<2)+nStatew4;
            oStatew4=nStatew4;
        }
        // Switch state if it is new
        if (oState!=nState){
            wtype=(wtype<<3)+nState;
            oState=nState;
        }
        ttype=(ttype<<3)+nState;
     
        const U8 brcontext=(brcxt.cxt);
        
        rcmA[0].set(word3*53+c1+193 * (ttype & 0x7fff),c1);
        //Retrun index of bracket or quote in array, if found add 1 to result. value is in range 1-7, 0 if not found
        w4br=0;
        if(brcxt.context)w4br=FindQy(brcontext)+1;
        if(brcxt.context==0 &&qocxt.context)w4br=FindQy(qocxt.cxt)+1;
        // Column and first char 
        col=colcxt.collen();
        int above=buffer[(nl1+col)&BMASK];
        int above1=buffer[(nl1+col-1)&BMASK];
        if (colcxt.isNewLine()) {
                // Reset contexts when there are two empty lines
                if (((c4&0xffffff00)==0x0a0d0a00 && buf(5)==0x0d) || (c4&0xffff00)==0x0a0a00){ 
                    fccxt.Reset();
                    brcxt.Reset();
                    qocxt.Reset();
                    cmCS[2].reset();
                }
                U8 oldfc=colcxt.lastfc(1);
                fc=colcxt.lastfc();
                if (oldfc=='.' && isLetter==true) { // book2
                    fccxt.Reset();
                    brcxt.Reset();
                    qocxt.Reset();
                } 
                if (isUpper==true && isLetter==true) {
                    fc=ATSIGN;
                }
                // Set paragraph
                if (fc==ATSIGN) fc1=1;
                else fc1=0;
                // Set new first char
                fccxt.Update(fc);
        }

        if ((col>2 && isLetter==false) ) {
            // Before updating first char context look:
            // If link or template ended then remove any vertical bars |.  [xx|xx] {xx|xx}
            if ((fccxt.cxt)==VERTICALBAR && (c1==SQUARECLOSE || c1==CURLYCLOSE)) while( (fccxt.cxt)==VERTICALBAR) fccxt.Update(10);
            // If html link ends with space or ]
            if (((fccxt.cxt)==COLON || (fccxt.cxt)==31 ) && (c1==SPACE || c1==SQUARECLOSE)) while( (fccxt.cxt)==COLON || (fccxt.cxt)==31) fccxt.Update(10);
            fccxt.Update(c1);
        }
        // Switch from possible category link to http link
        if ((fccxt.cxt)==COLON && c2=='/' && c1=='/') fccxt.Update(10),fccxt.Update(31);
        // Wiki link in the beginning of line
        if (colcxt.lastfc(0)==SQUAREOPEN && c1==SPACE) { 
            if(c2==SQUARECLOSE || c3==SQUARECLOSE) {
                fc=ATSIGN;
                fc1=1; // Paragraph
                // Line started with '[' and last chars were '] ', so probably rest of the line/paragraph fallows. 
                // Reset first char context and set new first char as paragraph start and continue.
                // Problem: inlink links also reset, like images that have description etc.
                fccxt.Reset();
                fccxt.Update(fc);
            }
        }
        // Fist char was space, look for another non-space char
        if (fc==SPACE  && c1!=SPACE) { 
            fc=min(c1,96);
            if (isUpper==true && isLetter==true) fc=ATSIGN;
            // Paragraph
            if (fc==ATSIGN) fc1=1;
            else fc1=0;
            // Set new first char, we keep space from previous update
            fccxt.Update(fc);
        }
        const U8 fccontext =fccxt.cxt;
        if(w4br==0 &&fccxt.context)w4br=FindQy(fccontext)+1;
        fqcxt=FindQy2(fccontext);

        cmCS[5].set((fccxt.context&0xff00)+c1+(w4&12)*256+((brcontext+ brcxt.last())<<24));
        // We have bold/italic first char/words. When words surrounded by ' end there is problably rest of paraghraph.
        // Reset first char context and set new first char as paragraph start and continue.
        if (colcxt.lastfc(0)==APOSTROPHE && (c1==SPACE)) { 
            if(c2==APOSTROPHE || c3==APOSTROPHE ){
                fc=ATSIGN;
                fc1=1;            // Paragraph
                fccxt.Reset();    // Not really needed
                fccxt.Update(fc); //
            }  
        }
        if ((fc!=ATSIGN) && ((c4&0xffffff)==0x4a2f2f)) {//http link - fixme
            fc=31;
        }
        if (isPDFText==false) {
            // Words surrounded by ()
            if ((worcxt.sBytes(1)&0xff)==')') {
                bool isC=false;
                for (int i=1; i<8; i++) if ((worcxt.sBytes(i)>>8)=='(') isC=true;
                if (isC) {
                    while( (worcxt.sBytes(1)>>8)!='(') worcxt.Remove();
                    worcxt.Remove();
                }
            }
            // Words surrounded by [|  as wiki internal link: word [word word|word] word.
            if ((worcxt.sBytes(1)&0xff)==VERTICALBAR) {
                bool isC=false;
                for (int i=1; i<8; i++) if ((worcxt.sBytes(i)>>8)==SQUAREOPEN) isC=true;
                if (isC) {
                    while( (worcxt.sBytes(1)>>8)!=SQUAREOPEN) worcxt.Remove();
                    worcxt.Remove();
                }
            } 
        }
        h=word0*271+number0+c1;

        // Indirect
        U32   f=(c4>>8)&0xffff;
        t2[f]=(t2[f]<<8)|c1;
        f=c4&0xffff;
        f=f|(t2[f]<<16);
        U32   d=(c4>>8)&0xff;
        t1[d]=(t1[d]<<8)|c1;
        d=c1|(t1[c1]<<8);
        
        t1[brcontext]=(t1[brcontext]<<2)|(w4&3); // this is wierd, also end is bad
        U32   d4=(ttype&7)|(t1[brcontext]<<3);
        // Contexts
        // Word stream cm(4-5)
        cmC3[4].set(word0+(number0*191+numlen0));
        cmC3[4].set(h+word1);
    
        cmC3[5].set(h+ word2*71);
        cmC3[5].set(hash(worcxt.fword,word1+h,(ttype & 63)));

        cmC3[6].set(((ttype&0x3f)<<16)+(c4&0xffff));
   
        cmC3[8].set((c4 & 0xffffff) + ((w4 << 18) & 0xff000000));
        cmC3[8].set((wtype&0x3fffffff)*4+(w4&3));
        cmC3[8].set((fccontext*4) + ((wtype & 0x3ffff) << 9 )+w4br);
    
        cmC4[0].set(colcxt.lastfc(0) | (fccontext << 15) | ((ttype & 63) << 7)|(brcontext << 24) );
        cmC4[0].set((colcxt.lastfc(0) | ((c4 & 0xffffff) << 8)));
    
        cmC4[1].set( (w4 & 3) +word0*11);
        cmC4[1].set(x4 & 0xffff);
        cmC4[1].set(((fc << 11) | c1)+((w4 & 3)<< 18));
        cmC4[1].set(word0*11+(w4 & 3) );
    
        cmC4[2].set((w4 & 15)+((ttype & 7) << 6 ));
        cmC4[2].set(c1 | ((col * (c1 == SPACE)) << 8)|((w4 & 15) << 16));
 
        cmC4[2].set(fc1?firstWord:(fc<< 11));
        cmC4[2].set((91 * 83* worcxt.Word(1) + 89 * word0));
        
        cmC4[3].set((c1 + ((ttype & 0x38) << 6)));
        cmC4[3].set(worcxt.fword*11+w4br);
        cmC4[3].set(c1+word0+number0*191 );
        cmC4[3].set(((c4 & 0xffff) << 16) | (fccontext  << 8) |fc);
        cmC4[3].set(((wtype & 0xfff)<< 8)+((w4 & 0xfc)));
        cmC4[3].set(h);
        
        // Switch between word/paragraph or column mode
        if (fc1==1){
            // Word
            cmCS[0].set(worcxt.fword*3191+(w4 & 3));
            cmCS[0].set(h+firstWord*89);
            cmCS[0].set(word0*53+c1+w4br);
        }else{
            // Column
            cmCS[0].set(above | ((ttype & 0x3f) << 9) | (colcxt.collen() << 19)| ((w4 & 3) << 16) );
            cmCS[0].set(h+firstWord*89);
            cmCS[0].set(above | (c1 << 16)| ((col+numlen0+w4br) << 8)| (above1<< 24)  );
        }
        // List
        if (colcxt.lastfc()=='*' ){
            cmCS[0].set(( ( fccontext) << 8)  | ((w4br & 0xfff) << 16));
            cmCS[0].set(c1);
            cmCS[0].set(word0);
        }else{
        // Table
            cmCS[0].set(wrt_w[bufr(colcxt.abovecellpos)]|( ( fccontext) << 8)  | ((w4br & 0xff) << 16));
            cmCS[0].set(bufr(colcxt.abovecellpos)|( ( c1) << 8) );
            cmCS[0].set( word0+wrt_w[bufr(colcxt.abovecellpos)] );
        }
        
        cmCS[1].set(x4 & 0xff00ff);
        cmCS[1].set((x4 & 0xff0000ff) | ((ttype & 0xe07) << 8));

        
        if (fccontext==CURLYOPENING) cmC3[7].set(0);else 
        cmC3[7].set(d4);
        
        cmCS[1].set((d4& 0xffff) | ((ttype & 0x38) << 16));
        // Two indirect bytes with current byte
        cmCS[0].set((f& 0xffffff));

        cmCS[2].set((c1 << 8) | (d >> 2)| (fc << 16));  // fixme
        cmCS[2].set((c4 & 0xffff)+(c2==c3?1:0));
   
        cmC4[4].set((ttype & 0x3ffff) | ((w4 & 255) << 24));
        cmC4[4].set(x4);
        // Word stream. word(1) with first char context and last bit3word(1-5)
        cmC3[9].set(257 * word1+fccontext + 193 * (ttype & 0x7fff));
        // First char, current byte and non repeating bit2word(1-6)
        cmC3[9].set(fc|((w4r & 0xfff) << 9) | ((c1  ) << 24));//end is good (lang)
        cmC3[9].set( worcxt.fword*83+(w4 & 3)*11+brcontext); // all category/language/image links (better as standalone)
        cmC3[9].set((x4 & 0xffff00)+ brcontext+(fccontext<< 24));

        cmC3[10].set(d);
        cmC3[10].set(((d& 0xffff00)>>4) | ((w4 & 0xf) )| ((ttype & 0xfff) << 20));
        cmC3[10].set((x4 >>16) | ((w4 & 255) << 24));
        cmC3[10].set((c1 << 11) | ((f & 0xffffff)>>16) );
        cmC3[10].set((fccontext*4+w4br) | ((c4 & 0xffff)<< 9)| ((w4 & 0xff) << 24)); 
        cmC3[10].set(((f >> 16) )| ((w4 & 0x3c)<< 25 )| (((ttype & 0x1ff))<< 16 ));

        cmC3[11].set(((words& 0xff) )+((( spaces & 0xff))<< 8)+((w4&15)<< 16)+(((wtype>>3)&511)<< 21)+(fc1<<30));
        cmC3[11].set(c1 + ((ttype<< 5) & 0x1fffff00));
        cmC3[11].set(ttype*16+w4br);
        cmC3[11].set(hash(worcxt.Word(3),(ttype & 0x1ff)));
        // Indirect byte with bracket context and last non repeating bit2words(2-10)
        cmC3[11].set(((d& 0xffff)>>8) + ((64 * w4r) & 0x3ffff00)+(brcontext<< 25)); // end good
        // Byte from prvious word(0) pos if in range(255) with current byte and first char
        int dd=pos-wp[word0&0xffff];
        if (dd>255)
            dd=256 + (c1<<16);
        else dd=dd + (buf(dd)<< 8)+(c1 << 16);
        if (fccontext==64 && brcontext==SQUAREOPEN) cmC3[11].set(0);else
        cmC3[11].set((dd )| (fc << 24));

        // Byte stream of x4, msb of byte(4), 4 msb bits of byte(2,3) and full byte(1)
        cmC3[12].set((x4&0x80f00000)+((x4&0x0000f0ff) << 12) );
        // Paragraph or column. 
        // In Paragraph: disabled when escaped utf8 or html link
        // In Column: when col is max(31) use last two bytes only otherwise add above bytes
        if (fc1==1){
            // word
            if (fccontext==31 ) 
                cmC3[12].set(0);
            else
                cmC3[12].set(h+worcxt.Word(1) *53 *79+worcxt.Word(3) *53*47 *71);
        }else{
            //column
            if (col==31) 
                cmC3[12].set( c4 << 16);
            else
                cmC3[12].set(above | ((c4 &0xffff)<< 16)| (above1<< 8));
        }
        
        // Word/Centence with current word(0), word(1) and word(2)
        cmC3[13].set(worcxt.Word(1)*83*1471-word0*53+worcxt.Word(2));
        // Word/centence. Disabled when: escaped utf8, template (onliner) or table beginning, html link.
        if (fccontext==CURLYOPENING || fccontext==31) 
            cmC3[13].set(0);
        else
            cmC3[13].set(h+worcxt.Word(2) *53 *79+worcxt.Word(3) *53*47 *71);

        // Last byte in wtype and w4 type with first char and bracket index
        cmCS[3].set(((wtype&7)<< 10) + (w4&3)+fc*4+ (w4br<< 24));
        // Current word or number
        cmCS[3].set( (word0*3301+number0*3191 ));//3191
        // Word/centence and non-repeating bit3words(5) with bracket/firstchar index
        cmC3[14].set(w4br+ worcxt.Word(2) * (wtype & 0x7fff));

        scmA[0].set(c1);
        scmA[1].set(c2*(fc1));
        scmA[2].set((f&0xffffff)>>16);
        scmA[3].set(ttype&0x3f);
        scmA[4].set(w4&0xff);
        scmA[5].set(brcontext);
        scmA[6].set(fc1+ 2*((wtype&0x3f)) );
        scmA[7].set(fc);

        if (wshift||c1==10)  {
            word3=word3*47, word2=word2*53, word1=word1*83;
        }

        cmC3[15].set((w4br*256)+fc+(((wtype>>0)&0xFFF)<< 16));
        // Some APM context
        AH1=hash((x5>>0)&255, (x5>>8)&255, (x5>>16)&0x80ff);
        AH2=hash(19,     x5&0x80ffff);
    }
    const int c0b=c0<<(8-bpos);
    ismatch=MatchModel2mix(0);

    scmA[0].mix(0);
    scmA[1].mix(0);
    scmA[2].mix(0);
    scmA[3].mix(0);
    scmA[4].mix(0);
    scmA[5].mix(0);
    scmA[6].mix(0);
    if ((fccxt.cxt)==COLON && fccxt.last()==SQUAREOPEN) {
        x.mxInputs[0].add(0);
        x.mxInputs[0].add(0);
    }
    else
        scmA[7].mix(0);

    // order X
    ord=cmC3[0].mix(0);
    if (ord==3) ord=2; // low max 2
    ord=ord+cmC3[1].mix(0);
    ord=ord+cmC3[2].mix(0);
    ord=ord+cmC3[3].mix(0);
    ord2=cmC3[4].mix(0);
    ord2=ord2+cmC3[5].mix(0);
    if (ord2>3) ord2=3;
    cmC3[6].mix(0);
    cmC3[7].mix(0);
    cmC3[8].mix(0);
    cmC4[0].mix(0);
    cmC4[1].mix(0);
    cmC4[2].mix(0);
    cmC4[3].mix(0);
    cmCS[0].mix(0);
    cmCS[1].mix(0);
    cmCS[2].mix(0);
    cmC4[4].mix(0);
    cmC3[9].mix(0);
    cmC3[10].mix(0);
    cmC3[11].mix(0);
    cmC3[12].mix(0);
    // order Word
    if ((fccxt.cxt)==CURLYOPENING|| (fccxt.cxt)==31 ) 
        cmC3[13].mix(0);  // For speed, also improves compression ?
    else
        ord2=ord2+(cmC3[13].mix(0));  
    cmCS[3].mix(0);  
    ord2=ord2+cmC3[14].mix(0);   
    cmC3[15].mix(0);
    cmCS[4].mix(0);
    cmCS[5].mix(0);
    rcmA[0].mix(0);
    
    // Mixer
    
    // reference
    // if(bpos){	
    //   c=c0<<(8-bpos); if(bpos==1)c=c+c3/2;
    //   c=(min(bpos,5))*256+c1/32+8*(c2/32)+(c&192);
    // }
    // else c=c3/128+(c4>>31)*2+4*(c2/64)+(c1&240);

    // mixer 1
    // at bpos=0   context is last 2 bit2word and 1 bit3word
    // at bpos=1-3 context is last 2 bit2word and 1 bit3word of current bracket or quote
    // at bpos=4-7 context is last 1 bit2word, current bit2word from c0 and bit3word of current bracket or quote
    if (bpos==0)  mxA[1].cxt=(w4&63)*8 + (ttype&7);
    else if (bpos>3) {
        c=wrt_w[c0b&255];
        mxA[1].cxt=(((w4<<2)&63)+c)*8+w4br;
    } else    
        mxA[1].cxt=(w4&63)*8 +w4br;

    // mixer 2
    // at bpos=0   context is was byte(3,4) a word and 2 bit3word
    // at bpos=1   context is bit 1xxxxxxx from c0, was byte(2) a word, bit pos max 5,last 1 bit3word and 1 bit3word of current first char context
    // at bpos=2   context is bit 11xxxxxx (bit pos 2) from c0, was byte(2) a word, bit pos max 5,last 1 bit3word and 1 bit3word of current first char context
    // at bpos=3   context is bit 111xxxxx (bit pos 3) from c0, was byte(2) a word, bit pos max 5,last 1 bit3word and 1 bit3word of current first char context
    // at bpos=4-7 context is bit current bit2word from c0, bit pos max 5,last 1 bit3word and 1 bit3word of current first char context
    if (bpos){
         c=c0b; 
         if (bpos==1) c=c+16 * (words*2& 4);
         else if (bpos>3)  c=wrt_w[c0b&255]*64;
         c=(min(bpos,5))*256+(wtype&7)+fqcxt*8+(c&192); //fqcxt -> cm(12,1)
    }
    else c=((words>>0)&12)*16+(wtype&7)+w4br*8;
    mxA[2].cxt=c;
    
    // mixer 3
    // at bpos=0-7   context is was byte(3,4) a word, sum of context order(3-5,6,8) isState counts (max 5) and last 2 bit2word
    mxA[3].cxt=((4 * words) & 0xf0) + ord*256 + (w4 & 15);
    
    // mixer 7
    // at bpos=0-7   context is non-repeating 2 bit3word of byte(2,3), was byte(1-3) a word and last 1 bit2word
    //mxA[7].cxt=((wtype) & 0x1f8)*4 + ((2 * words) & 0x1c) + (w4 & 3);
    if (bpos>3) {
        c=wrt_w[c0b&255];
        mxA[7].cxt=((wtype) & 0x1f8)*4 + ((w4 & 3)*4) + (c );
    }
    else if (bpos==0) 
        mxA[7].cxt=((wtype) & 0x1f8)*4 + ((4 * w4br) ) + (w4 & 3);
    else
        mxA[7].cxt=((wtype) & 0x1f8)*4 + ((2 * words) & 0x1c) + (w4 & 3);
    c=c0b;
    // mixer 4
    // at bpos=0   context is bit xxxxxxxx from c0, was byte(1-8) a word or space and bit pos
    // at bpos=1   context is bit 1xxxxxxx from c0, was byte(1-7) a word or space and bit pos
    // at bpos=2   context is bit 11xxxxxx from c0, was byte(1-6) a word or space and bit pos
    // at bpos=3   context is bit 111xxxxx from c0, was byte(1-5) a word or space and bit pos
    // at bpos=4   context is bit 1111xxxx from c0, was byte(1-4) a word or space and bit pos
    // at bpos=5   context is bit 11111xxx from c0, was byte(1-3) a word or space and bit pos
    // at bpos=6   context is bit 111111xx from c0, was byte(1-2) a word or space and bit pos
    // at bpos=7   context is bit 1111111x from c0, was byte(1)   a word or space and bit pos
    mxA[4].cxt=bpos*256 + (((( (numbers|words)<< bpos)&255)>> bpos) | (c&255));
    // mixer 9 - final mixer
    // at bpos=0-7   context is sum of context order(3-5,6,8) isState counts (max 5), bracket or quote state(0,1) and last 1 bit2word
    mxA[9].cxt=0;//(ord*8 + (w4br?1:0)*4 + (w4&3));
    
    // mixer 5 
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
            c=c + 16*(ttype&7);
        }
        else if (bpos==2) {
            c=c + 16*(w4&3);
        }
        else if (bpos==3) {
            c=c + 16*(words&1);
        } else  {
            c=bpos + (c&0xf0);
        }
        if (bpos<5)
            c=bpos + (c&0xf0); 
    }else   c=16 * (w4&0xf);
    ord=ord-1;
    if (ord<0)
       ord=0;
    if (ismatch)
        ord=ord+1;
    mxA[5].cxt=c + ord*256+ 8*fc1;

    // mixer 6
    // at bpos=0-7   context is sum of context words isState counts (max 6), first char type state(0,1), 2 bit2word of byte(3,4) and 1 bit3word of byte(2)
    mxA[6].cxt=(ord2*256 + (w4&0xf0) + ((ttype&0x38) >> 2))*4 + fqcxt;

    // mixer 8 
    // at bpos 0-2 bit3word, 1 bit3word of current bracket or quote, was byte(1,2,3) a word, first char flag, is a match
    // at bpos 3-7 bit3word from c0, 1 bit3word of current bracket or quote, was byte(1,2,3) a word, first char flag, is a match
    if (bpos>2) 
        mxA[8].cxt= wrt_t[c0b&255]*256 +(w4br)*32 + (words&7)*4 + fc1+(ismatch?2:0);
    else
        mxA[8].cxt= (ttype&3)*256 +(w4br)*16 + (words&7)*2 + fc1+(ismatch?128:0);  // fixme
    // Layer 0
    x.mxInputs[1].add(mxA[1].p1());
    x.mxInputs[1].add(mxA[2].p1());
    x.mxInputs[1].add(mxA[3].p1());
    x.mxInputs[1].add(mxA[4].p1());
    x.mxInputs[1].add(mxA[5].p1());
    x.mxInputs[1].add(mxA[6].p1());
    x.mxInputs[1].add(mxA[7].p1());
    x.mxInputs[1].add(mxA[8].p1());
    // Layer 1
    for (int i=1;i<9;i++) mxAX[i].cxt=mxA[i].cxt;
    for (int i=1;i<9;i++) x.mxInputs[2].add(mxAX[i].p1());
    // Layer 2
    return mxA[9].p();
}

void update() {
    x.c0+=x.c0+x.y;
    if (x.c0>=256) {
        x.c8=(x.c8<<8)+(x.c0&0xff);
        x.c4=U32(x.c8);
        x.c0=1;
        ++x.blpos;
        // When last byte was predicted good/below error treshold then set new limits to mixer update
        // larger value means less updates and better speed.
        if ((fails&255)==0) {
            const int mx0Limit=256+(((fails>>8)&255)==0?32:0);
            //const int mx1Limit=16+(((fails>>8)&255)==0?16:0);
            for (int i=1;i<9;i++) mxA[i].elim=max(mx0Limit,mxA[i].elim+1);
            for (int i=1;i<9;i++) mxAX[i].elim=max(mx0Limit,mxA[i].elim+1);
        } else {
            for (int i=1;i<9;i++) mxA[i].elim=min(16,mxA[i].elim-1);
            for (int i=1;i<9;i++) mxAX[i].elim=min(16,mxA[i].elim-1);
        }
        // APM update rate based on input file position
        rate=6 + (x.blpos>14*256*1024) + (x.blpos>28*512*1024);
    }
    x.bpos=(x.bpos+1)&7;
    x.bposshift=7-x.bpos;
    x.c0shift_bpos=(x.c0<<1)^(256>>(x.bposshift));
    x.cmBitState=getStateByteLocation(x.bpos,x.c0);
    //mxA[0].update(x.y);
    for (int i=1;i<10;i++) mxA[i].update(x.y);
    for (int i=1;i<9;i++) mxAX[i].update(x.y);

    //if (x.blpos<2)printf("mixer 0 predictor count %d\n",x.mxInputs[0].ncount);
    x.mxInputs[0].ncount=0;
    //printf("mixer 1 predictor count %d\n",x.mxInputs[1].ncount);
    x.mxInputs[1].ncount=0;
    x.mxInputs[2].ncount=0;
    // This part is from paq8hp12
    if (fails&0x00000080) --failcount;
    fails=fails*2;
    failz=failz*2;

    if (x.y) pr=4095-pr;
    if (pr>=e_l[x.bpos]) ++fails, ++failcount;
    if (pr>=848) ++failz;
    
    pr=modelPrediction(x.c0,x.bpos,x.c4);
    int pt, pu=(apmA[0].p(pr, x.c0, 3,x.y)+7*pr+4)>>3, pv, pz=failcount+1;

    pz+=tri[(fails>>5)&3];
    pz+=trj[(fails>>3)&3];
    pz+=trj[(fails>>1)&3];
    if (fails&1) pz+=8;
    pz=pz/2;

    pu=apmA[3].p(pu,   ((x.c0*2)^AH1)&0x3ffff, rate,x.y);
    pv=apmA[1].p(pr,   ((x.c0*8)^hash(29,failz&2047))&0xffff, rate+1,x.y);
    // If fails use w4 else non-repeating w4
    if (fails&255)
        pv=apmA[4].p(pv, hash(x.c0,w4 & 0xfffc,(wtype & 0x1ff))&0x1ffff, rate,x.y);
    else
        pv=apmA[4].p(pv, hash(x.c0,(w4r & 0xfffc)+0x10000,(wtype & 0x1ff))&0x1ffff, rate,x.y);
    pt=apmA[2].p(pr, ( (x.c0*32)^AH2)&0xffff, rate,x.y);
    pz=apmA[5].p(pu,   ((x.c0*4)^hash(min(9,pz),x5&0x80ff))&0x1ffff, rate,x.y);
    
    if (fails&255) pr=(pt*6+pu  +pv*11+pz*14 +31)>>5;
    else           pr=(pt*4+pu*5+pv*12+pz*11 +31)>>5;

}

void Init() {
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
    
    for (int i=0; i<4096; i++) {
        st2_p1[i]=clp(sc(12*(i-2048)));
    }
    for (int i=0; i<4096; i++) {
        st1in[i]=clp(sc(32*stretch(i)));
    }
    // contextmap run mul tables for 2, 3, 4, 5, 6
    // c=ilog(rc+1)<<(2+(~rc&1)) 
    for (int i=2; i<7;i++) {
        for (int rc=0;rc<256;rc++) {
            int c=ilog[rc];
            c=c<<(2+(~rc&1));
            if ((rc&1)==0) c=c*i/4;
            g_rc1[i-2][rc+256]=clp(c);
            g_rc1[i-2][rc]=clp(-c);
        }
    }
    StateTable statetable;
    statetable.Init(28, 28, 31, 29, 23, 4, 17,&STA[0][0][0]); // 1
    statetable.Init(32, 28, 31, 28, 21, 5,  6,&STA[1][0][0]); // 2
    statetable.Init(31, 27, 30, 27, 24, 4, 27,&STA[2][0][0]); // 4
    statetable.Init(33, 31, 31, 24, 20, 4, 33,&STA[3][0][0]); // 5
    statetable.Init(28, 29, 30, 30, 23, 3, 22,&STA[4][0][0]); // 6
    statetable.Init(28, 29, 33, 23, 23, 6, 14,&STA[5][0][0]); // 7
    // CM4 st9
    for (int j=0; j<6; j++) {
        for (int i=0; i<256; i++) {
            U32 n0=STA[j][i][2]*3+1;
            U32 n1=STA[j][i][3]*3+1;
            int p=(n1<<12) / (n0+n1);
            g_st9[j][i]=clp(sc(18*(p - 2048)));
        }
    }
}

// Encoder
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

typedef enum {COMPRESS, DECOMPRESS} Mode;
class Encoder {
private:
  const Mode mode;       // Compress or decompress?
  FILE* archive;         // Compressed data file
  U32 x1, x2;            // Range, initially [0, 1), scaled by 2^32
  U32 x0;                // Decompress mode: last 4 input bytes of archive
  FILE*alt;              // decompress() source in COMPRESS mode
  U32 scompsize;
  // Compress bit y or return decompressed bit
  void code(int i=0) {
    int p=pr;
    p+=p==0;
    assert(p>0 && p<4096);
    U32 xmid=x1 + ((x2-x1)>>12)*p + (((x2-x1)&0xfff)*p>>12);
    assert(xmid>=x1 && xmid<x2);
    x.y=i;
    i ? (x2=xmid) : (x1=xmid+1);
    update();
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
    for (int i=7; i>=0; --i)
        code((c>>i)&1);
    if (!(x.blpos&0x3ffff))  printf("%d  %d\n",x.blpos,size()-scompsize), scompsize=size();
  }

  // Decompress and return one byte
  int decompress() {
      int c=0;
      for (int i=0; i<8; ++i)
        c+=c+decode();
      return c;
  }
  
  ~Encoder(){ }
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

#include <psapi.h>
size_t getPeakMemory(){
#if defined(_WIN32)
    PROCESS_MEMORY_COUNTERS info;
    GetProcessMemoryInfo( GetCurrentProcess( ), &info, sizeof(info) );
    return (size_t)info.PeakPagefileUsage; // recuested peak memory /PeakWorkingSetSize used memory
#elif defined(UNIX) 
    return (size_t)0L; //not tested
#else
    return (size_t)0L;
#endif // -lpsapi
}
  
int main(int argc, char** argv) {   
    // Check arguments
    if (argc!=4 || (argv[1][0]!='c' && argv[1][0]!='d')) {
        printf(
        "fxcm_text v%d file compressor (C)  2026, Kaido Orav.\n"
        "Licensed under GPL, http://www.gnu.org/copyleft/gpl.html\n"
        "\n"
        "To compress:   fxcm_text c input output \n"
        "To decompress: fxcm_text d input output \n",VERSION);
        return 1;
    }

    // Get start time
    clock_t start=clock();
    
    // Open input file
    FILE *in=fopen(argv[2], "rb");
    if (!in) exit(1);
    FILE *out=0;

    Init();

    // Match model
    mhashtablemask=0x200000*1-1;
    alloc1(mhashtable,0x200000*1+32,mhptr,32);  
    
    if (argv[1][0]=='c') {        
        // Compress
        // Encode header: signature, version, file size
        fseeko(in, 0, SEEK_END);
        U64 size=ftello(in);
        fseeko(in, 0, SEEK_SET);
        out=fopen(argv[3], "wb");
        if (!out) exit(1);
        fprintf(out, "tF%c", VERSION);
        put64(size,out);

        Encoder e(COMPRESS, out);
        int c;
        while ((c=getc(in))!=EOF)
            e.compress(c);
        e.flush();
    } else {
        // Decompress
        // Check header signature, version, file size
        if (getc(in)!='t' || getc(in)!='F' || getc(in)!=VERSION)
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



