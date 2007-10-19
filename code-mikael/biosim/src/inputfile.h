#ifndef _INPUTFILE_H
#define _INPUTFILE_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>

using namespace std;

class inputfile {
public:
  inputfile(string);                    //constructor, opens file
  string getStr(string);                //get string value
  double getDbl(string, double=0);      //get double value
  int getInt(string, int=0);            //get integer value
  bool getBool(string, bool=false);     //get boolean value

private:
  struct dataformat {
     string name;
     string val;
  };
  vector<dataformat> matrix;
  int findKey(string &);
};

#endif
