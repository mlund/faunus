#ifndef _pdbclass_h
#define _pdbclass_h

#include <iostream>
#include <string>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>
#include "point.h"

using namespace std;

class pdbclass {
public:
  struct format {
    string entry;
    int atomnumber;
    string atomname;
    string residuename;
    string chainid;
    string seqid;
    int residuenumber;
    double x,y,z;
    double unk1, unk2;
    string moleculename;
    int atomcount;
    float charge,radius;
  };
  int modelCnt, remarkCnt;
  ostringstream out, empty;

  pdbclass();
  void remark(double);
  void conect(int, int);
  void hetatm(format &);
  void atom(format &);
  void model(int=0);
  void endmdl();
  void save(string);
  void load_particles(vector<particle> &);
};
#endif

