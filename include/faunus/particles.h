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

      int push_back(const particle &);      //!< add particle to both "p" and "trial"
      double charge();                      //!< Sum all charges in particle vector
      double charge(const point &, double); //!< Sum all charges within a sphere region
      bool overlap(const particle &);       //!< Check for overlap w. particle
      bool overlap(const std::vector<particle> &);
      bool check_vector();                  //!< Check if p and trial are equal!
      int count(particle::type, const point&,double);//!< Count particles of "type" within a sphere
      vector<particle::type> list_of_species() const; //!< Get vector of particle types in the system
  };
}
#endif
