// Converts an AAM file to PQR (for VMD)
#include "faunus/group.h"
#include "faunus/io.h"

using namespace Faunus;
using namespace std;

int main() {
  cell con(100);
  con.atom.load("../../../misc/faunatoms.dat");
  macromolecule protein;
  ioaam aam(con.atom);
  iopqr pqr(con.atom);
  protein.add( con, aam.load( "../twobody-hofmeister/lysozyme-ph4.7.aam" ));
  pqr.save("out.pqr", con.p);
}
