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
      struct dataformat {
        string name;
        vector<string> val;
      };
      vector<dataformat> matrix;
      vector<string> calls;
      int findKey(string &);
      void record_call(string);
      string file;
    public:
      inputfile(string);                     //!< Constructor
      string getstr(string, string="");//!< Get string value
      double getflt(string, double=0); //!< Get double value
      int getint(string, int=0);       //!< Get integer value
      bool getboo(string, bool=false); //!< Get boolean value
      vector<string> getvec(string,string);  //!< Get vector of strings
      void add(string,string);               //!< Add an entry to the loaded list
      void add(string,double);               //!< Add an entry to the loaded list
      string info();                         //!< Show info
  };
}
#endif

