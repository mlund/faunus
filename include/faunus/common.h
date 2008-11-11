#ifndef FAU_COMMON_H
#define FAU_COMMON_H
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cstdarg>
#include <cmath>
#include <algorithm>

namespace Faunus {
  //Keep this at a minimum, please!
  using std::string;
  using std::vector; // Perhaps try boost::vector?
  using std::endl;
  using std::setw;
  using std::abs;    // *very* important if cstdlib.h or algorithm is included!
  using std::sqrt;
  using std::exp;

  namespace phys {
    static double pi=acos(-1.),
                  e0=8.85419e-19,  //!< Permittivity of vacuum
                  kB=1.380658e-23, //!< Boltzmann's constant [J/K]
                  e=1.602177e-19,  //!< Electron unit charge [C]
                  Nav=6.022137e23; //!< Avogadro's number [1/mol]
  }
}

#endif
