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
#include "xydata.h"
#include "intrinsics.h"

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
  ASP,GLN,GLU,ASN,LYS,ARG,PRO,UNK,NTR,CTR,NA,K,CL,BR,I,ION,CATION,ANION,GHOST,
  RNH3,RNH4,RCOOH,RCOO,LAST};

  xydata<double> pairpot[LAST][LAST]; //pair-potentials
  vector<datafmt> d;
  double lB;

  aminoacid(double=7.1);
  id getId(string);
  double vdW(id);
  double weightcalc(id);           //returns Mw of aminoacid
  double volume(double, double);   //...and a volume estimate
  double radius(double,double);    //...and a radius estimate

  bool loadpmf(string);     //load pair potentials from disk
  void loadpmf(string,string);
  void showpmf();
  double energy(vector<particle> &);
  double energy(vector<particle> &, unsigned int);
  double energy(vector<particle> &, group &);
  double energy(vector<particle> &, group &, group &);

  inline double pow2(double x) { return x*x; };
  inline double pow6(double x) { return pow2(x*x*x); };
  inline double energy(particle &p1, particle &p2) {
    unsigned int i=p1.id,j=p2.id;
    if (i==GHOST || j==GHOST) return 0;
    if (i>j) swap(i,j);
    double r2=p1.sqdist(p2); 
    if (pairpot[i][j].xmax==0) {
      double qq=p1.charge*p2.charge;
      return pow6( pow2(p1.radius+p2.radius)/r2 )
        + (qq!=0) ? lB*qq/sqrt(r2) : 0;
    }
    r2=sqrt(r2);
    return (r2>pairpot[i][j].xmax) ?
      lB*p1.charge*p2.charge/r2 : // use Coulomb pot. outside data
      pairpot[i][j].x2y(r2);      // else use table.
  };
};

class peptide : public aminoacid {
 private:
  struct filefmt {
    double x,y,z,mw,charge;
    string atomname, aminoname;
    unsigned int atomnr, aminonr;
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
