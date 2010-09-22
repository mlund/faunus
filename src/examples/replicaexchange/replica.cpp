#define OPENMP_TEMPER

#include "faunus/faunus.h"
#include "faunus/potentials/pot_coulomb.h"
#include "faunus/potentials/pot_debyehuckelP3.h"
#include "faunus/bottles/npt_molecular.h"
#include "faunus/energy/coarsegrain.h"

/*!
 * \brief Replica exchange program - proof of concept.
 * \author Mikael Lund
 * \date Malmo 2010
 */

using namespace std;
using namespace Faunus;

typedef interaction_dipole<pot_debyehuckelP3> Tpot;
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
  
  fio.writefile( "a.out", a.preinfo() );
  fio.writefile( "b.out", b.preinfo() );
  fio.writefile( "c.out", c.preinfo() );

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
          //#pragma omp section
          //c.microloop(100);
        }
      }
      if (temperBool) {
        switch (rand() % 1) {
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

    cout << loop.timing();
  } // end of outer loop

  fio.writefile( "a.out", a.postinfo(), std::ios_base::app );
  fio.writefile( "b.out", b.postinfo(), std::ios_base::app );
  fio.writefile( "c.out", c.postinfo(), std::ios_base::app );

  cout << loop.info() << in.info() << temper.info();
};
