#ifndef FAU_CRANKSHAFT_H
#define FAU_CRANKSHAFT_H

#include "faunus/moves/base.h"
#include "faunus/energy/springinteraction.h"

namespace Faunus {
  /*!
   * Vector rotation routines.
   * Boundary condition are respected.
   *
   * \author Mikael Lund
   * \date Canberra, 2009
   */
  class vectorRotation {
    private:
      point origin, u;
      double cosang, sinang, eb;
      double e1mcox, e1mcoy, e1mcoz;
      container* conPtr;

    public:
      vectorRotation(container &);        //!< Constructor - container needed for boundaries 
      void setAxis(point, point, double); //!< Set axis of rotation and angle
      point rotate(point);                //!< Rotate point around axis
  };

  /*!
   * \brief Crankshaft move for polymers
   * \author Mikael Lund
   * \date Canberra, 2009
   *
   * This class will do crank-shaft moves of an arbitrary number of monomers,
   * spanned by two randomly chosen end points. The reported mean-square displacement
   * represents the movement of the polymer mass center.
   */
  class crankShaft : public markovmove {
    private:
      unsigned short beg;            //!< First CS fix point
      unsigned short end;            //!< Last CS fix point
      unsigned short len;            //!< Number of atoms to include in crank-shaft move
      vector<unsigned short> v;      //!< Particles between beg and end to rotate;
      vectorRotation rot;            //!< Rotation class
      void setNumberOfMonomers();    //!< Randomly set number of monomers to rotate
      void setNumberOfMonomers(polymer &); //!< Randomly set number of monomers to rotate between 1 and the length of the polymer
      bool findEnds(polymer &);      //!< Find random end points to rotate around
      bool useMinMaxMonomers;        //!< Use minMonomers and maxMonomers from inputfile

    public:
      unsigned short minMonomers;    //!< Minimum number of monomers to move
      unsigned short maxMonomers;    //!< Maximum number of monomers to move

      crankShaft(ensemble &, container &, energybase &, inputfile &);
      double move(polymer &);        //!< Do crank-shaft move
      double penaltymove(polymer &); //!< Do crank-shaft move, while taking into account masscenter position dependent potentials
      double move(polymer &, int);   //!< Do crank-shaft move n times
      string info();
  };

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
