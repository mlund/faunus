#ifndef titrate_h
#define titrate_h

#include <vector>
#include <algorithm>
#include <iostream>
#include "species.h"
#include "point.h"
#include "group.h"
#include "average.h"

using namespace std;

/*! \brief Class to perform proton titration of molecules
 *  \author Mikael Lund
 *  \todo Too much is public! Documentation is bad.
 */
class titrate {
  public:
    enum keywords {PROTONATED,DEPROTONATED,ANY,NOACID};
    vector<short int> sites, protons, neutrons;
    double ph;                          //!< System pH

    struct action {
      keywords action;
      short int site;        
      short int proton;    
    };

    titrate(species &, double);
    titrate(species &, vector<particle> &, group &, double);
    void init(vector<particle> &, group &);//!< Locate and initialize sites and protons
    action exchange(vector<particle> &);
    action exchange(vector<particle> &, action &);
    double sumsites();                     //!< Calculates total charge of titrateable sites
    double energy(vector<particle> &, double, action &);
    void samplesites(vector<particle> &);  //!< Updates the average charge vector titrate::q
    void showsites(vector<particle> &);    //!< Print average charges of titrateable sites
    double applycharges(vector<particle> &);//!< Copy average charges to particles in the particle vector
  private:
    void infos();
    species *spc;
    short int random(vector<short int> &); //!< Pick a random item in a vector
    average<float> nprot;               //!< Average number of protons. Updated with titrate::samplesites
    vector<average <float> >  q;        //!< Stores the average charges of sites from titrate::sites
    action takeFromBulk(vector<particle> &, short int, short int=-1);
    action moveToBulk(vector<particle> &, short int, short int=-1);
    keywords status( vector<particle> &, short int );
    group sort(group &);
};

/*! \brief Class to preform grand canonical titration where the 
 *         potential determining ion is coupeld with a bulk salt
 *  \author Bjorn Persson/Mikael Lund
 *  \todo Everything...
 */


class GCtitrate {
  public:
    double  CatPot, volume;         //!< Systems chemical potential of cation, system volume
    titrate tit;

    GCtitrate GCtitrate(species &, double, double, double);//!< Constructor
    GCtitrate GCtitrate(species &, vector<particle> &, group &, double, double, double); //!<Constroctor
    void updateVolume( double );            //!< Update the system volume, for isobaric ensemble
};

#endif
