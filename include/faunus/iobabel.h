#ifndef FAU_IOBABEL_H
#define FAU_IOBABEL_H

#ifdef BABEL

#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/obmolecformat.h>
#include <openbabel/generic.h>
#include "faunus/point.h"
#include "faunus/species.h"

namespace Faunus {
  using OpenBabel::OBResidueIter;
  using OpenBabel::OBAtomAtomIter;
  /*! \brief OpenBabel file support
   *  \author Mikael Lund
   *  \date Helsingor, 2008
   *
   *  This class offers support for OpenBabel - a set of programs and
   *  libraries for reading and writing a large number of molecular formats.
   *  The main purpose of this class is to bridge between the babel api and
   *  the faunus particle class approach.
   *
   *  Particle types defined in "faunatoms.dat" are added to the OB
   *  element table so that external OB input routines automatically
   *  recognize these. Note that OB can have a maximum of 128 (one char)
   *  elements where approximately 116 of these are already occupied.
   *
   *  If the number of residues equals the number of particles, the
   *  residue name will be used to identify the particle according
   *  the list provided by Faunus::atoms
   *
   *  Recommended file format: Sybyl/mol2 but any should work.
   */
  class iobabel {
    private:
      vector<int> nb;
      OpenBabel::OBConversion obconv;
      OpenBabel::OBMol obmol;
      OpenBabel::OBAtom obatom;
      OpenBabel::OBAtom *obatomPtr;
      OpenBabel::vector3 v;  // OpenBabel vector
      particle a;   // Tmp particle
      double c[3];  // Temp. vector storage
      void p2atom(particle &); // Convert particle to ObenBabel atom
    public:
      iobabel();
      vector<particle> p;                         //!< Placeholder for loaded data
      particle get(int);                          //!< Convert i'th babel atom to a particle
      void read(string, bool=true);               //!< Read entire file (autodetect format from extension)
      bool write(string,const vector<particle> &);//!< Write coordinates (format from extension)
      vector<int> neighbors(int);                 //!< Get list of neighboring atoms
  };
}
#endif
#endif
