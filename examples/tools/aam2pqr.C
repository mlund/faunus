// Converts an AAM file to PQR (for VMD)
#include "group.h"
#include "io.h"

using namespace std;

int main() {
  cell con(100);
  macromolecule protein;
  ioaam aam(con);
  iopqr pqr(con);
  protein.add( con, aam.load( "../twobody-hofmeister/lysozyme-ph4.7.aam" ));
  pqr.save("out.pqr", con.p);
}
