#ifndef titrate_h
#define titrate_h

#include <vector>
#include <algorithm>
#include <iostream>
#include "point.h"
#include "peptide.h"
#include "slump.h"
#include "group.h"
#include "average.h"

using namespace std;

//! Class to perform proton titration of molecules
class titrate : private aminoacid, private slump {
 public:
  enum keywords {PROTONATED,DEPROTONATED,ANY,NOACID};
  vector<short int> sites, protons, neutrons;
  vector<average> q; //!< Stores the average charges of sites from titrate::sites
  vector<point> eqpos; //!< Stores the equilibrium (starting) positions of sites

  average nprot; //!< Average number of protons. Updated with titrate::samplesites
  
  struct action {
    keywords action;
    short int site;        
    short int proton;    
  };

  group sort(group &);
  titrate(double);
  titrate(vector<particle> &, group &, double);
  void init(vector<particle> &, group &);//!< Locate and initialize sites and protons
  action exchange(vector<particle> &);
  action exchange(vector<particle> &, action &);
  void update(vector<particle> &, group &); //!< Update proton vector, for GC etc...
  double check(action &, double, double);
  double idPref( action &, double);
  double energy(vector<particle> &, double, double, action &);
  double energy(vector<particle> &, double, action &);
  double sumsites();                     //!< Calculates total charge of titrateable sites
  void samplesites(vector<particle> &);  //!< Updates the average charge vector titrate::q
  void showsites(vector<particle> &);    //!< Print average charges of titrateable sites
  double applycharges(vector<particle> &); //!< Copy average charges to particles in the particle vector
  
  void info();

 private:
  double ph;  //!< System pH
  action takeFromBulk(vector<particle> &, short int, short int=-1);
  action moveToBulk(vector<particle> &, short int, short int=-1);
  keywords status( vector<particle> &, short int );
  short int random(vector<short int> &); //!< Pick a random item in a vector
};

#endif
