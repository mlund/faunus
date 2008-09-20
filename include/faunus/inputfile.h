#ifndef FAU_INPUTFILE_H
#define FAU_INPUTFILE_H

#include <faunus/common.h>

namespace Faunus {
  /*!
   * \brief Retrieve parameters from a formatted input file
   * \author Mikael Lund
   *
   * The input file is expected to have the following
   * format:
   * \li  KEYWORD VALUE (space separated)
   *
   * Blank lines or lines containing [ or # are ignored.
   */
  class inputfile {
    private:
      struct dataformat { string name, val; };
      vector<dataformat> matrix;
      int findKey(string &) const;
    public:
      inputfile(string);                     //!< Constructor
      string getstr(string, string="") const;//!< Get string value
      double getflt(string, double=0) const; //!< Get double value
      int getint(string, int=0) const;       //!< Get integer value
      bool getboo(string, bool=false) const; //!< Get boolean value
      void add(string,string);               //!< Add an entry to the loaded list
      void add(string,double);               //!< Add an entry to the loaded list
  };
}
#endif

