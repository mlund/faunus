#include <iostream>
#include "xtcio.h"

using namespace std;

int main() {
  float x[10];
  float box[3][3];
  int natoms=10;
  int magic=1995;
  int step=1;
  int time=1;
  int xd;
  xd=open_xtc("test.xtc", "w" );

  //write_xtc(XDR *xd, natoms, step,
  //    time,box,
}
