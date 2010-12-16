#ifndef FAU_IOBABEL_H
#define FAU_IOBABEL_H

#ifdef BABEL

#ifndef SWIG
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/obmolecformat.h>
#include <openbabel/generic.h>
#include "faunus/point.h"
#include "faunus/species.h"
#endif

namespace Faunus {
  using OpenBabel::OBResidueIter;
  using OpenBabel::OBAtomAtomIter;
  /*! \brief  OpenBabel file support
   *  \author Mikael Lund
   *  \date   Helsingor, 2008
   *
   *  This class offers support for OpenBabel - a set of programs and
   *  libraries for reading and writing a large number of molecular formats.
   *  The main purpose of this class is to bridge between the babel api and
   *  the faunus particle class approach.
   *
   *  The number of residues should equal the number of particles since
   *  the residue name is used to identify the particle (OB atom names
   *  are mapped to the periodic table elements).
   *
   *  Recommended file format: Sybyl/mol2 but any should work. The following shows a
   *  Sybyl/Mol2 file where "FauX" are the Faunus species used to identify the atoms.
   *  Hence "C" and "H" have no effect. The column just before "FauX" designates
   *  the residue numbers and should match the number of atoms. Partial charges are
   *  read from the last column.
   *
   *  \code
   *  @<TRIPOS>MOLECULE
   *  *****
   *  5 4 0 0 0
   *  SMALL
   *  GASTEIGER
   *  Energy = 0
   *
   *  @<TRIPOS>ATOM
   *  1 C          -0.4391    5.3926  -11.7502 C.3     1  FauC       -0.0776
   *  2 H           0.6309    5.3926  -11.7502 H       2  FauH        0.0194
   *  3 H          -0.7963    6.4012  -11.7502 H       3  FauH        0.0194
   *  4 H          -0.8096    4.8684  -10.8941 H       4  FauH        0.0194
   *  5 H          -0.7908    4.8949  -12.6297 H       5  FauH        0.0194
   *  @<TRIPOS>BOND
   *  1     1     2    1
   *  2     1     3    1
   *  3     1     4    1
   *  4     1     5    1
   *  \endcode
   */
  class iobabel {
    private:
      OpenBabel::OBConversion obconv;
      OpenBabel::OBMol obmol;
      OpenBabel::OBAtom obatom;
      OpenBabel::OBAtom *obatomPtr;
      OpenBabel::OBResidue *obresPtr;
      OpenBabel::vector3 v;         //!< OpenBabel vector
      void p2atom(particle &);      //!< Convert Faunus particle to ObenBabel atom
    public:
      iobabel();
      vector<particle> p;           //!< Placeholder for loaded data
      particle get(int);            //!< Convert i'th babel atom to a particle
      void read(string);            //!< Read entire file (autodetect format from extension)
      bool write(string,const vector<particle> &);//!< Write coordinates (format from extension)
      vector<int> neighbors(int);   //!< Get list of neighboring atoms
  };
}
#endif
#endif
