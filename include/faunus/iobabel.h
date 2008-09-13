#ifndef FAU_IOBABEL_H
#define FAU_IOBABEL_H

#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/obmolecformat.h>
#include "faunus/point.h"

namespace Faunus {
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
      OpenBabel::OBConversion obconv;
      OpenBabel::OBMol mol;
      OpenBabel::OBAtom atom;
      OpenBabel::OBAtom *atomPtr;
      OpenBabel::vector3 v;    // OpenBabel vector
      particle a;   // Tmp particle
      double c[3];  // Temp. vector storage
      void p2atom(particle &); // Convert particle to ObenBabel atom
    public:
      vector<particle> p;         //!< Placeholder for loaded data
      particle get(unsigned int); //!< Convert n'th babel atom to a particle
      void read(string);          //!< Read entire file (autodetect format from extension)
      bool write(string,vector<particle> &);//!< Write coordinates (format from extension)
  };
}
#endif
