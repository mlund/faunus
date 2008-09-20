#ifndef FAU_INPUTFILE_H
#define FAU_INPUTFILE_H

#include <faunus/common.h>

namespace Faunus {
  /*!
   * \brief Retrieve data from a formatted input file
   * \author Mikael Lund
   *
   * The input file is expected to have the following
   * format:
   *    KEYWORD VALUE (space separated)
   */
  class inputfile {
    public:
      inputfile(string);                     //!< Constructor
      string getstr(string, string="") const;//!< Get string value
      double getflt(string, double=0) const; //!< Get double value
      int getint(string, int=0) const;       //!< Get integer value
      bool getboo(string, bool=false) const; //!< Get boolean value
    private:
      struct dataformat { string name, val; };
      vector<dataformat> matrix;
      int findKey(string &) const;
  };
}
#endif

