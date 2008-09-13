#ifndef FAU_SCOPES_H
#define FAU_SCOPES_H
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cstdarg>
#include <cmath>
#include <algorithm>

#ifdef GROMACS
#ifndef CPLUSPLUS
#define CPLUSPLUS
#endif
#include "xtcio.h"
#endif

namespace Faunus {
  //Keep this at a minimum, please!
  using std::string;
  using std::vector; // Perhaps try boost::vector?
  using std::endl;
}
#endif
