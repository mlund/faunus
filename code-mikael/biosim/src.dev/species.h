#ifndef _species_h
#define _species_h
#include <string>
#include <vector>
#include <iostream>
#include "point.h"

/*!
 * \brief Specific data for particles
 * \author Mikael Lund
 */
class species {
 private:
  struct data {
    string name;
    float radius,charge,pka;
    bool hydrp;
  };
  void set(int,string,float,float,float,bool);

  public:
    species();
    vector<data> d;                     //!< Data is stored here.
    particle::type id(string);          //!< Transform species name to particle id
    double mw(int);
    double vol(double, double);         //!< Volume estimate
    double radius(double,double);       //!< Radius estimate
    particle get(string);
    string info(string);
};
#endif
