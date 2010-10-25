#define OPENMP_TEMPER

#include "faunus/faunus.h"
#include "faunus/potentials/pot_coulomb.h"
#include "faunus/potentials/pot_sqwell.h"
#include "faunus/bottles/nvt_sq.h"
#include "faunus/energy/coarsegrain.h"

/*!
 * \brief Replica exchange program - proof of concept.
 * \author Mikael Lund
 * \date Malmo 2010
 */

using namespace std;
using namespace Faunus;

typedef interaction<pot_sqwell> Tpot;
typedef nvt_sq<box,Tpot> Tbottle;

int main() {
  slump slp;
  cout << faunus_splash() << slp.info() << endl;
  inputfile in("replica_chains.conf");
  io fio;
  mcloop loop(in);
  Tbottle a("a"), b("b"), c("c"), d("d"), e("e"), f("f");
  vector<Tbottle*> bottles;
  replicaexchange temper;
  bool temperBool=in.getboo("temper",true);

  bottles.push_back(&a);
  bottles.push_back(&b);
  bottles.push_back(&c);
  bottles.push_back(&d);
  bottles.push_back(&e);
  bottles.push_back(&f);

  cout << in.info();

  while (loop.macroCnt()) {
    while (loop.microCnt()) {

    //  #pragma omp parallel for
      for (int i=0; i<bottles.size(); ++i)
        bottles[i]->microloop(50);

      if (temperBool) {
        int N=bottles.size()-1;
        int i=slp.rand() % N;
        int j=slp.rand() % N;
        int k=slp.rand() % N;
        temper.swap(*bottles[i], *bottles[i+1]);
        if (j!=i)
          temper.swap(*bottles[j], *bottles[j+1]);
        if (k!=j && k!=i)
          temper.swap(*bottles[k], *bottles[k+1]);
        //cout << i << " " << j << " " << k << endl;

     }

    } // end of inner loop

    for (int i=0; i<bottles.size(); ++i)
      bottles[i]->macroloop();

    cout << loop.timing();
  } // end of outer loop

  for (int i=0; i<bottles.size(); ++i)
    bottles[i]->finish();

  cout << loop.info() << in.info() << temper.info();
};
