#ifndef titrate_h
#define titrate_h

#include <vector>
#include <algorithm>
#include <iostream>
#include "species.h"
#include "point.h"
#include "slump.h"
#include "group.h"
#include "average.h"

using namespace std;

//! Class to perform proton titration of molecules
class titrate : private slump {
  public:
    enum keywords {PROTONATED,DEPROTONATED,ANY,NOACID};
    vector<short int> sites, protons, neutrons;
    vector<average <double> >  q; //!< Stores the average charges of sites from titrate::sites
    vector<point> eqpos; //!< Stores the equilibrium (starting) positions of sites

    average<double> nprot; //!< Average number of protons. Updated with titrate::samplesites

    struct action {
      keywords action;
      short int site;        
      short int proton;    
    };

    group sort(group &);
    titrate(species &, double);
    titrate(species &, vector<particle> &, group &, double);
    void init(vector<particle> &, group &);//!< Locate and initialize sites and protons
    action exchange(vector<particle> &);
    action exchange(vector<particle> &, action &);
    double energy(vector<particle> &, double, action &);
    double sumsites();                     //!< Calculates total charge of titrateable sites
    void samplesites(vector<particle> &);  //!< Updates the average charge vector titrate::q
    void showsites(vector<particle> &);    //!< Print average charges of titrateable sites
    double applycharges(vector<particle> &);//!< Copy average charges to particles in the particle vector

    void info();

  private:
    species *spc;
    double ph;  //!< System pH
    action takeFromBulk(vector<particle> &, short int, short int=-1);
    action moveToBulk(vector<particle> &, short int, short int=-1);
    keywords status( vector<particle> &, short int );
    short int random(vector<short int> &); //!< Pick a random item in a vector
};

#endif
