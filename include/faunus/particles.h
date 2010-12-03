#ifndef FAU_particles_h
#define FAU_particles_h

#ifndef SWIG
#include "faunus/common.h"
#include "faunus/slump.h"
#endif

namespace Faunus {
  class point;
  class group;
  class particle;

  /*!
   * \brief Class for all particles.
   * \author Mikael Lund
   * \date Lund, 2004-2010
   *
   * This class contains all the particles of the system. Particles are
   * kept in two vectors, "p" and "trial". "trial" is used when performing
   * trial MC moves while "p" contains the current up-to-date set of
   * particles. Before and after any MC move the two vectors "p" and "trial"
   * MUST be in sync.
   *
   * This class can also take care of particle insertion in the system, for
   * example by grand canonical schemes. To make sure that all groups
   * (i.e. particle ranges) are kept in sync after an insertion, add
   * these groups to the "g" group pointer vector. The insertion/deletion methods will
   * then automatically increase or decrease their sizes.
   */

  class particles {
    private:
      slump slp;
    public:
      vector<particle> p;                   //!< The main particle vector
      vector<particle> trial;               //!< Trial particle vector. 
      vector<group*> g;                     //!< Pointers to all groups in the system.

      int push_back(const particle &);      //!< add particle to both "p" and "trial"
      bool insert(particle, unsigned int);  //!< Insert particle at pos n.
      bool remove(unsigned int);            //!< Remove particle n.
      double charge();                      //!< Sum all charges in particle vector
      double charge(const point &, double); //!< Sum all charges within a sphere region
      virtual bool clash(const particle &, const particle &); //!< Overlap between two particles
      virtual bool overlap(const particle &);       //!< Check for overlap w. particle
      bool overlap(const std::vector<particle> &);
      bool check_vector();                  //!< Check if p and trial are equal!
      int count(unsigned char, const point&,double);//!< Count particles of "type" within a sphere
      vector<unsigned char> list_of_species() const; //!< Get vector of particle types in the system
  };
}
#endif
