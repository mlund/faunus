#define OPENMP_TEMPER

#include "faunus/faunus.h"
#include "faunus/potentials/pot_coulomb.h"
#include "faunus/potentials/pot_debyehuckelP3.h"
#include "faunus/bottles/npt_molecular.h"

/*!
 * \brief Replica exchange program - proof of concept.
 * \author Mikael Lund
 * \date Malmo 2010
 */

using namespace std;
using namespace Faunus;

typedef interaction<pot_debyehuckelP3> Tpot;
typedef npt_molecular<box,Tpot> Tbottle;

int main() {
  cout << faunus_splash() << slp.info() << endl;
  inputfile in("replica.conf");
  io fio;
  mcloop loop(in);
  Tbottle a("a"), b("b"), c("c");

  replicaexchange temper;
  bool temperBool=in.getboo("temper",true);

  a.prepare();
  b.prepare();
  c.prepare();

  cout << in.info();

  while (loop.macroCnt()) {
    while (loop.microCnt()) {
      #pragma omp parallel
      {
        #pragma omp sections
        {
          #pragma omp section
          a.microloop(100);
          #pragma omp section
          b.microloop(100);
          #pragma omp section
          c.microloop(100);
        }
      }
      if (temperBool) {
        switch (rand() % 2) {
          case 0:
            temper.swap(a,b);
            break;
          case 1:
            temper.swap(b,c);
            break;
        }
      }
    } // end of inner loop
    a.macroloop();
    b.macroloop();
    c.macroloop();

    a.save();
    b.save();
    c.save();

    cout << loop.timing();
  } // end of outer loop

  fio.writefile( "a.out", a.postinfo() );
  fio.writefile( "b.out", b.postinfo() );
  fio.writefile( "c.out", c.postinfo() );

  cout << loop.info() << in.info() << temper.info();
};
