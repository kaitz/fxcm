
#include <stdio.h>
#define NDEBUG
#include <assert.h>
#include <string.h>


class StateTable {
  int mdc; // maximum discount
  enum {B=5, N=64}; // sizes of b, t
  int b[6];  // x -> max y, y -> max x
  unsigned char ns[1024]; // state*4 -> next state if 0, if 1, n0, n1
  unsigned char t[N][N][2]={{{0}}};
  int num_states(int x, int y);  // compute t[x][y][1]
  void discount(int& x);  // set new value of x after 1 or y after 0
  void next_state(int& x, int& y, int b);  // new (x,y) after bit b
  void generate();  // compute t[x][y][1]
  int pre(int state) {  // initial probability of 1 * 2^23
    assert(state>=0 && state<256);
    return ((next(state,3)*2+1)<<14)/(next(state,2)+next(state,3)+1);
  }
public:
  int next(int state, int sel) {return ns[state*4+sel];}
  StateTable();
  Init(int s0,int s1,int s2,int s3,int s4,int s5,int s6);
} ;


int StateTable::num_states(int x, int y) {
  if (x<y) return num_states(y, x);
  if (x<0 || y<0 || x>=N || y>=N || y>=B || x>=b[y]) return 0;
  return 1+(y>0 && x+y<b[5]);
}

// New value of count x if the opposite bit is observed
void StateTable::discount(int& x) {
  int y=0;
  if (x>2){
    for (int i=1;i<mdc;i++) y+=x>=i;
    x=y;
  }
}

// compute next x,y (0 to N) given input b (0 or 1)
void StateTable::next_state(int& x, int& y, int b) {
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
void StateTable::generate() {
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
        printf("{%3d,%3d,%2d,%2d},", ns[state*4], ns[state*4+1],
          ns[state*4+2], ns[state*4+3]);
        if (state%4==3) printf(" // %d-%d\n", state-3, state);
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
      }
    }
  }
}

StateTable::StateTable() {
}

StateTable::Init(int s0,int s1,int s2,int s3,int s4,int s5,int s6) {
    b[0]=s0;b[1]=s1;b[2]=s2;b[3]=s3;b[4]=s4;b[5]=s5;mdc=s6;
    printf("\n Generating with state parameters p1 %d,p2 %d,p3 %d,p4 %d,p5 %d,p6 %d,p7 %d\n",s0,s1,s2,s3,s4,s5,s6);
    printf("Statetable:\n");
    generate(); 
    printf("\n");   
}
StateTable stable;


int main(int argc, char** argv) {
  if (argc!=8) {
    printf("Usage:\nstates p1 p2 p3 p4 p5 p6 p7\n");
    exit(1);
  }
  stable.Init(atoi(argv[1]) ,atoi(argv[2]) ,atoi(argv[3]) ,atoi(argv[4]) ,atoi(argv[5]), atoi(argv[6]),atoi(argv[7]) );
  
  return 0;
}
