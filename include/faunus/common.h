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
  using std::abs;    // *very* important if cstdlib.h or algorithm is included!
  using std::sqrt;
  using std::exp;
}
#endif
