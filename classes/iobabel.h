#ifndef _iobabel_h
#define _iobabel_h

#include <iostream>
#include <vector>
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include "point.h"

using namespace OpenBabel;

/*! \brief OpenBabel file import/export
 *  \author Mikael Lund
 *  \todo Unfinished!
 */
class iobabel {
  private:
    OBConversion obconv;
    OBMol mol;
    vector<particle> p;
  public:
    void read(string,string);
};

#endif
