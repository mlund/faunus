#ifndef FAUNUS_COMMON_H
#define FAUNUS_COMMON_H

#ifndef SWIG
#include <Eigen/StdVector>
#include <string>
#include <list>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <complex>
#include <map>
#include <random>
#include <set>
#include <memory>
#include <numeric>
#include <type_traits>
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
    class VectorRotate;
  }

  namespace Energy {
    class Bonded;
    class Energybase;
    class Hamiltonian;
    class EqEnergy;
    class HarmonicBond;
    class GouyChapman;
    class EnergyRest;
  }

  namespace Potential {
    class PairPotentialBase;
    class HardSphere;
  }

  class AtomData;
  class Group;
  class GroupMolecular;
  class InputMap;
  class UnitTest;
  class RandomBase;
  class Space;

  //class Point;
  class PointParticle;
  class DipoleParticle;
  class CigarParticle;

#if defined(CIGARPARTICLE)
  typedef CigarParticle particle;
#elif defined(DIPOLEPARTICLE)
  typedef DipoleParticle particle;
#else
  typedef PointParticle particle;
#endif

  typedef std::vector< particle, Eigen::aligned_allocator<particle> > p_vec;

  // FUNCTORS
  typedef std::function<double()> RandFunctor;

}//namespace

#endif
