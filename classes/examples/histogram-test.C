#include "../histogram.h"

using namespace Faunus;

int main() {
  histogram h(0.1,-5,5);
  h.add(-4.9);
  h.add(4);
  h.add(-5.0);
  h.write("sletmig");
  float val=h.get(4); // val=0.3333
}
