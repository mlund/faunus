#define CPLUSPLUS

#include <iostream>
#include "slump.h"
#include "io.h"
#include "xtcio.h"
#include "container.h"

using namespace std;

class ioxtc : public iopart {
  private:
    vector<particle> load(string) {}
    rvec x[1000];
    int xd;
    float box[3][3], time, step;
  public:
    ioxtc(container &con, string file, string mode) : iopart(con) {
      time=step=0;
      xd=open_xtc("test.xtc", "w");
    }
    void save(vector<particle> &p) {
      for (unsigned short i=0; i<p.size(); i++) {
        x[i][0]=p[i].x/10;
        x[i][1]=p[i].y/10;
        x[i][2]=p[i].z/10;
      }
      write_xtc(xd,p.size(),step++,time++,box,x,1000.);
    }
    void close() { close_xtc(xd); }
};

int main() {
  slump s;
  species spc;
  vector<particle> p(10);
  ioxyz out(spc);
  rvec x[10];
  float box[3][3],len=30;
  int natoms=10, step=1, time=0, xd;
  xd=open_xtc("test.xtc", "w" );

  for (int l=0; l<100; l++) {
    for (int i=0; i<10; i++) {
      p[i].x=len*s.random_one();
      x[i][0]=p[i].x/10;
      p[i].y=len*s.random_one();
      x[i][1]=p[i].y/10;
      p[i].z=len*s.random_one();
      x[i][2]=p[i].z/10;
    }
    write_xtc(xd, natoms, step++, time++, box, x, 1000.);
  }
  out.save("test.xyz", p);
  close_xtc(xd);
}
