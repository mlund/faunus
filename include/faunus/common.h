#ifndef FAUNUS_COMMON_H
#define FAUNUS_COMMON_H

#ifndef SWIG
#include <string>
#include <vector>
#include <list>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cstdarg>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <complex>
#include <map>
#include <set>
#include <memory>
#include <numeric>
#endif


// Use explicit virtual override and final keywords (C++11)
#ifdef NO_EXPLICIT_OVERRIDE
  #define FOVERRIDE
#else
  #define FOVERRIDE override
#endif

namespace Faunus {
  //Keep this at a minimum, please!
  using std::string;
  using std::vector;
  using std::complex;
  using std::endl;
  using std::setw;
  using std::abs;    // *very* important if cstdlib.h or algorithm is included!
  using std::sqrt;
  using std::exp;
  using std::cout;
  using std::ostringstream;
  using std::shared_ptr;

  namespace Move {
    class Movebase;
  }

  namespace Geometry {
    class Geometrybase;
  }

  namespace Energy {
    class Energybase;
    class Hamiltonian;
    class EqEnergy;
    class HarmonicBond;
    class GouyChapman;
  }

  class Group;
  class GroupMolecular;
  class InputMap;
  class UnitTest;
  class Space;
  class Point;
  class PointParticle;
  class CigarParticle;

#ifdef CIGARPARTICLE
  typedef CigarParticle particle;
#else
  typedef PointParticle particle;
#endif
  typedef std::vector<particle> p_vec;

}//namespace

#endif
