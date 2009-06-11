#include "faunus/iobabel.h"
using namespace std;
using namespace Faunus;
int main() {
  iobabel b;
  b.read("water.pdb");
  cout << "Size = " << b.p.size() << endl;
  for (int i=0; i<b.p.size(); i++)
    cout << b.p[i].x << " " << b.p[i].y << " " << b.p[i].z
      << " " << b.p[i].charge << " " << b.p[i].radius << endl;
}
