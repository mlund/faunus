#ifndef FAU_TRANSLATIONAL_H
#define FAU_TRANSLATIONAL_H

#include "faunus/moves/base.h"
#include "faunus/histogram.h"

namespace Faunus {

  //-----------Z-MOVE--------------------------------------
  /*! \brief Class to make moves of a single particle along z-axis
   *  
   *  This class will make a random walk along the z-axis for a
   *  macromolece
   *
   *  \todo Needs some testing and perhaps some optimization
   *  \note Replaced by dualmove(?)
   *  \author Bjoern Persson
   */
  class zmove : public markovmove {
    public:
      zmove( ensemble&, container&, energybase&);
      double move(macromolecule &);
      double z, zmax;
  };

  //----------------- DUAL MOVE --------------------------
  /*!
   * This move will symmetrically translate two macromolecules
   * along the line connecting their mass centers.
   *
   * \brief Symmetrically move two groups along z-axis
   * \author Mikael Lund
   */
  class dualmove : public markovmove {
    private:
      point v;
    public:
      histogram gofr;             //!< g(r) of the two group mass centers
      double r;                   //!< Current distance between group mass centers
      double rmax;                //!< Maximum allowed mass-center distance
      double rmin;                //!< Minimum allowed mass-center distance
      dualmove( ensemble&, container&, energybase&);
      void setup(inputfile &);
      void load(inputfile &, vector<macromolecule> &g, float=0);
      void direction(double, double, double);
      double move(macromolecule &, macromolecule &);
      string info();
  };

  //---------- TRANSLATE MOLECULE ----------------
  class translate : public markovmove {
    public: 
      point dpv;  //!< Displacement vector
      translate( ensemble&, container&, energybase&, inputfile&);
      double move(group &); 
  };

  //-----------MOVE----------------------------------------
  /*! \brief Random move of a macromolecule in a container 
   *  \todo generalize, it is cubix periodic boundry specific
   *  \author Bjoern Persson
   *
  class move : public markovmove {
    private:
      box *b;
    public:
      move( ensemble&, container&, energybase &, box &, float);
      bool mOve(macromolecule &);
  };*/

  //-------------- SALT MOVE ---------------------------------
  /*! \brief Move salt particles
   *  \author Mikael Lund
   */
  class saltmove : public markovmove {
    protected:
      void init();
      long double rsqr;  //!< Mean square displacement
    public:
      point dpv;                      //!< Displacement direction vector
      saltmove( ensemble &, container&, energybase& );
      saltmove( ensemble &, container&, energybase&, inputfile &);
      double move(group &, int);      //!< Move a single particle
      double move(group &);           //!< Loop over group particles (randomly)
      string info();
  };

  //-------------- TRANSLATE A MONOMER PARTICLE ----------------------------
  /*! \brief Translate monomer particle
   *  \author Mikael Lund
   */
  class monomermove : public saltmove {
    public:
      monomermove(ensemble &, container&, energybase&, inputfile&);
      double move(polymer &);      //!< Move a single monomer particle
  };
}//namespace
#endif
