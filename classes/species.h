#ifndef _species_h
#define _species_h
#include <string>
#include <vector>
#include "point.h"

/*!
 * \brief Specific data for particles
 * \author Mikael Lund
 */
class species {
  private:
    struct data {
      particle p;
      float pka;
      string name;
    };
    void set(particle::type,string,float,float,float,bool);

  public:
    species();
    vector<data> d; //!< Data is stored here.
    particle::type id(string);          //!< Transform name to particle
    particle get(particle::type);       //!< Get particle from id
    particle get(string);               //!< Get particle from string (name)
};
#endif
