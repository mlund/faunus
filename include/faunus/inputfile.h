#ifndef FAU_INPUTFILE_H
#define FAU_INPUTFILE_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <stdlib.h>

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
      inputfile(std::string);                    //!< Constructor
      std::string getstr(std::string, std::string="");//!< Get string value
      double getflt(std::string, double=0);      //!< Get double value
      int getint(std::string, int=0);            //!< Get integer value
      bool getboo(std::string, bool=false);      //!< Get boolean value
    private:
      struct dataformat { std::string name, val; };
      std::vector<dataformat> matrix;
      int findKey(std::string &);
  };
}
#endif

