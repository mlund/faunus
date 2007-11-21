#ifndef _CONFIG_H
#define _CONFIG_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <stdlib.h>

using namespace std;

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
  inputfile(string);                    //!< Constructor
  string getstr(string, string="");     //!< Get string value
  double getflt(string, double=0);      //!< Get double value
  int getint(string, int=0);            //!< Get integer value
  bool getboo(string, bool=false);      //!< Get boolean value
private:
  struct dataformat { string name, val; };
  vector<dataformat> matrix;
  int findKey(string &);
};
#endif

