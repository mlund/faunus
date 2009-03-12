#ifndef FAU_ROSEMBLUTH_H
#define FAU_ROSEMBLUTH_H
#include <faunus/moves/base.h>
namespace Faunus {
  /*! \brief Rosenbluth polymer and salt insertion
   *  \author Christophe Labbez and Mikael Lund
   *  \date Dijon / Lund, 2009
   *  \warning ....
   *  \todo a lot...
   */
  class rosenbluth : public markovmove { 
    private:
      grandcanonical *gcPtr;    // Pointer to GC class
      struct data {
        group* gPtr;            // Pointer to group containing this species
        unsigned short valency; // Stoechiometry
        short charge;           // Particle charge
        vector<string> seq;     // Sequence of particle names
      };
      data polymer, counter;    // Store data for "polymer" and counter ions
      short index;              // Rosembluth index
      string bondtype;          // Insertion scheme or bond type
      unsigned int k;           // Rosembluth k-value
      double mu;                // Chemical potential
      void insert();
      void remove();
      vector<unsigned int> ins; // Particles inserted or deleted in last move.
    public:
      // Keep as much as possible private! Easier for the user of the class.
      rosenbluth( grandcanonical&, container&, energybase&, inputfile &, int);
      double move();
      string info();
  };
}//namespace
#endif
