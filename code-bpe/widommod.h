#ifndef _WIDOMMOD_H
#define _WIDOMMOD_H

#include "average.h"
#include "particle.h"
#include "group.h"
#include "simbox.h"
#include "interact.h"

using namespace std;

class Widommod : public Particle {
  private:
    vector<double> chel;        //< electrostatic
    vector<double> chhc;        //< hard collision
    vector<double> chex;        //< excess
    vector<double> chexw;       //< excess
    vector<double> chtot;       //< total
    vector< vector<double> > ewden;       //< charging denominator
    vector< vector<double> > ewnom;       //< charging nominator
    vector< vector<double> > chint;       //< charging integrand
    long long int cnt;                    //< count test insertions
    int ghostin;                //< ghost insertions
    int gspec;                  //< number of species to test

  public:
    Widommod(int, int);                 //< Constructor, number of species, number of test insertions
    
    vector<double> chid;                //< ideal term
//    vector<double> chelav;              
//    vector<double> chhcav;
//    vector<double> chexav;
//    vector<double> chtotav;
    vector<double> expuw;
    vector<int> ihc;
    vector<int> irej;

    void ghostInsert(vector<Particle> &, vector<Group> &, Simbox &, Interact &);                 // Ghost insertion
    void getWidomResults();
};

#endif
