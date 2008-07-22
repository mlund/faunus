#ifndef _peptide_h
#define _peptide_h

//Defines amino acids, peptides etc.
//M. Lund, 2004

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include "point.h"
#include "group.h"

using namespace std;

class datafmt {
  public:
    double radius,charge,pka;
    string name;
    bool hydrp;

    datafmt() {
      radius=charge=pka=0;
      hydrp=false;
    };
};

class aminoacid {
 public:
  enum id {FIRST=0,GLY,ALA,VAL,LEU,ILE,PHE,TRP,TYR,HIS,SER,THR,MET,CYS,
	   ASP,GLN,GLU,ASN,LYS,ARG,PRO,UNK,NTR,CTR,ION,CATION,A37,LAST};

  vector<datafmt> d;

  aminoacid();
  id getId(string);
  double vdW(id);
  double weightcalc(id);           //returns Mw of aminoacid
  double volume(double, double);   //...and a volume estimate
  double radius(double,double);    //...and a radius estimate
};

class peptide : private aminoacid {
 private:
  struct filefmt {
    double x,y,z,mw,charge;
    string atomname, aminoname;
    int atomnr, aminonr;
  };
  double charge(vector<particle> &);
  point center_of_mass(vector<particle> &);

 public:
  enum keys {AAM, ATOMIC, SPHERICAL, CHAIN};
  vector<particle> loadAAModel(string, keys=AAM);
  vector<particle> loadstructure(string, keys);
  bool saveAAModel(string, vector<particle> &, group &);
};

#endif
