#ifndef FAUNUS_COMMON_H
#define FAUNUS_COMMON_H

#ifndef SWIG
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cstdarg>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <complex>
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

  namespace Move {
    class movebase;
    class translate;
  }

  namespace Geometry {
    class geometrybase;
  }

  class inputfile;
  class space;
  class point;
  class pointparticle;
  class cigarparticle;
#ifdef CIGARPARTICLE
  typedef cigarparticle particle;
#else
  typedef pointparticle particle;
#endif

}//namespace
#endif
