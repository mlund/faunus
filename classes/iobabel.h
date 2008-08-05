#ifndef _iobabel_h
#define _iobabel_h

#include <iostream>
#include <vector>
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include "point.h"

using namespace OpenBabel;

/*! \brief OpenBabel file support
 *  \author Mikael Lund
 *  \date Helsingor, 2008
 *
 *  This class offers support for OpenBabel - a set of programs and
 *  libraries for reading and writing a large number of molecular formats.
 *  The main purpose of this class is to bridge between the babel api and
 *  the faunus particle class approach.
 */
class iobabel {
  private:
    OBConversion obconv;
    OBMol mol;
    OBAtom *atom;
    particle a;
    vector3 v;    // OpenBabel vector
    double c[3];  // Temp. vector storage
  public:
    vector<particle> p;         //!< Placeholder for loaded data
    particle get(unsigned int); //!< Convert n'th babel atom to a particle
    void read(string,string);   //!< Read entire file
};
#endif
