#ifndef _INPUTFILE_H
#define _INPUTFILE_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

using namespace std;

class inputfile {
public:
  inputfile(string);                    //constructor, opens file
  void setcatPot( string&);
  void updateInp( double, double, string); 
  string getStr(string);                //get string value
  double getDbl(string, double=0);      //get double value
  int getInt(string, int=0);            //get integer value
  bool getBool(string, bool=false);     //get boolean value
  inline string stringify( double x)
  {
    ostringstream o;
    o<<x;
    return o.str();
  }

private:
  struct dataformat {
     string name;
     string val;
  };
  vector<dataformat> matrix;
  int findKey(string &);
};

#endif
