#ifndef FAU_SWITCH_H
#define FAU_SWITCH_H

#include <faunus/moves/base.h>
#include <faunus/titrate.h>
#include "faunus/potentials/pot_debyehuckel.h"
#include "faunus/species.h"

namespace Faunus {

  // The real hamiltonian is computed relative to all particles being fType1
  // and a energy of 0

  class switch_type : public markovmove, public titrate_switch {
  private:
    struct data {
      average<double> charge;                                     //!< Average net-charge
    };
    average<double> frac_2;                                       //!< Fraction of fType2 particles in titrate_switch::sites
    unsigned int n;                                               //!< Particle number

  public:
    switch_type( ensemble&, container&, energybase&, inputfile&); //!< Initialize and find titratable sites
    string info();                                                //!< Info string
    double titrateall();                                          //!< Sum of gaussian energies for all sites
    void samplesites(container &);
  };

}//namespace
#endif
