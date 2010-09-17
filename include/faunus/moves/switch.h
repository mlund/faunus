#ifndef FAU_SWITCH_H
#define FAU_SWITCH_H

#include <faunus/moves/base.h>
#include <faunus/titrate.h>
#include "faunus/potentials/pot_debyehuckel.h"
#include "faunus/species.h"

namespace Faunus {

  class switch_type : public markovmove, public titrate_switch {
  private:
    struct data {
      average<double> charge;//!< Average net-charge
    };
    average<double> frac_2;
    unsigned int n;        //!< Particle number

  public:
    switch_type( ensemble&, container&, energybase&, inputfile&); //!< Initialize and find titratable sites
    string info();                //!< Info string
    double titrateall(); //!< Sum of gaussian energies for all sites
    void samplesites(container &);
  };

}//namespace
#endif
