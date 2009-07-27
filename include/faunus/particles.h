#ifndef FAU_particles_h
#define FAU_particles_h

#include "faunus/slump.h"
#include "faunus/point.h"

namespace Faunus {
  /*!
   * \brief Class the contains the all particles including coordinates, trials etc.
   * \author Mikael Lund
   * \date Lund, 2004
   */
  class particles {
    private:
      slump slp;
    public:
      vector<particle> p;                   //!< The main particle vector
      vector<particle> trial;               //!< Trial particle vector. 
      vector<point> pmu;                    //!< Point dipole moment (optional)
      vector<point> trialmu;                //!< Point dipole moment trial (optional)

      int push_back(const particle &);      //!< add particle to both "p" and "trial"
      bool insert(particle, unsigned int);  //!< Insert particle at pos n.
      bool remove(unsigned int);            //!< Remove particle n.
      double charge();                      //!< Sum all charges in particle vector
      double charge(const point &, double); //!< Sum all charges within a sphere region
      virtual bool clash(const particle &, const particle &); //!< Overlap between two particles
      bool overlap(const particle &);       //!< Check for overlap w. particle
      bool overlap(const std::vector<particle> &);
      bool check_vector();                  //!< Check if p and trial are equal!
      int count(unsigned char, const point&,double);//!< Count particles of "type" within a sphere
      vector<unsigned char> list_of_species() const; //!< Get vector of particle types in the system
  };
}
#endif
