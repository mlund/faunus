#ifndef FAU_BRANCHROTATION_H
#define FAU_BRANCHROTATION_H

#include "faunus/moves/base.h"
#include "faunus/moves/crankshaft.h"
#include "faunus/energy/springinteraction.h"

namespace Faunus {
    /*!
   * \brief Branch rotation for polymers
   * \author Chris Evers
   * \date Lund, 2011
   *
   * This class will rotate a branch of a polymer with an arbitrary number of monomers,
   * between a randomly choosen monomer (A) and an end terminus (T), around an axis between 
   * A and another randomly choosen monomer (B).
   */
  class branchRotation : public markovmove {
    private:
      unsigned short A;              //!< Begin of branch to be rotated and begin of rotation axis
      unsigned short B;              //!< End of rotation axis
      unsigned short T;              //!< End of branch to be rotated
      unsigned short len;            //!< Number of atoms to include in move
      vector<unsigned short> v;      //!< Particles between A and T to rotate;
      vectorRotation rot;            //!< Rotation class
      void setNumberOfMonomers(polymer &); //!< Randomly set number of monomers to rotate between 1 and the length of the polymer
      bool findEnds(polymer &);      //!< Randomly find A, B and T in polymer

    public:
      branchRotation(ensemble &, container &, energybase &, inputfile &);
      double move(polymer &);        //!< Do branch rotation
      double penaltymove(polymer &); //!< Do branch rotation, while taking into account masscenter position dependent potentials
      double move(polymer &, int);   //!< Do branch rotation n times
      string info();
  };
} // namespace
#endif
