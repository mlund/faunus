#ifndef FAU_SALTBATH_H
#define FAU_SALTBATH_H
#include "faunus/moves/base.h"

namespace Faunus {
  /*! \brief Class for Grand Canonical MC
   *  \author Bjorn Persson and Mikael Lund
   *  \date Lund, 2010
   *  \warning This is some general pice of shit
   *
   *  Class for handling an arbitrary salt solution of N components whereof M may 
   *  may have a fixed chemical potential rather then fixed number density.
   */
  class saltbath : public markovmove { 
    private:
      bool locked;              //!< Is class closed for external modification?
      grandcanonical *gcPtr;    //!< Pointer to GC class
      salt *sPtr;               //!< Pointer to Salt class (change to group)
      struct pair {
        int i, j;               //!< Indices for GC::*g  
        int ni,nj;              //!< Stochiometric ratio: ni:nj of a i-j salt
        average < double >  ai, aj;         //!< Number averages
        average < double >  inacc, outacc;  //!< Acceptances
      };
      vector<pair> pairs;       //!< All (electro neutral) ion pairs with reservoirs 
      pair *thispair;           //!< Choosen pair for a particular step in the Markov chain
      bool inserting;           //!< Is it a inserting step?
      particle I, J;            //!< Particles to insert
      vector<particle> trialin, tini, tinj;  //!< Particle vectors trailin=tini+tinj
      vector<unsigned int>      trialout, touti, toutj;  //!< Indices for attempted extraction
      void remove();
      void acc_ins();           //!< Accept insertion
      void acc_rem();           //!< Accept extraction
      void insert();
      void setChemPot(inputfile &); //!< Load and set atoms::list.chempot
      double metropolis();      //!< returns the full Boltzmann weight, taking mu's and ideal cont. in to account.
      int gcd(int, int);        //!< Find the greatest common divider of two integers (external library?)

    public:
      // Keep as much as possible private! Easier for the user of the class.
      vector<group> g;          //!< Vector of salt components (not known externaly, maybe the should?)
      saltbath( grandcanonical&, container&, energybase&, inputfile &, salt &);
      void searchsalt( salt&);
      double move();
      string info();
      void check(checkValue &); //!< Perform unittest
  };
}//namespace
#endif
