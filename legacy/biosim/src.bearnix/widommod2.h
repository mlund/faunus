#ifndef _WIDOMMOD_H
#define _WIDOMMOD_H


#include "average.h"
#include "cell.h"
#include "point.h"
#include "group.h"
#include "interact.h"

using namespace std;

/****************
The class is bulit so that
each declared object returns 
the potential of one specis,
declared in the constructor
****************/

class Widommod : public particle {
  private:
    string name;        //< name of spec.
    double charge;      //< charge of spec.
    double radius;      //< radius of spec.
    double chel;        //< electrostatic
    double chhc;        //< hard collision
    double chex;        //< excess
    double chexw;       //< excess
    double chtot;       //< total
    vector<double>  ewden;       //< charging denominator
    vector<double>  ewnom;       //< charging nominator
    vector<double>  chint;       //< charging integrand
    long long int cnt;                    //< count test insertions
    int ghostin;                //< ghost insertions
    int gspec;                  //< number of species to test

  public:
    Widommod(int, double ,double ,string);                 //< Constructor, number of test insertions, radius and charge
    
    double chid;                //< ideal term
//    vector<double> chelav;              
//    vector<double> chhcav;
//    vector<double> chexav;
//    vector<double> chtotav;
    double expuw;
    int ihc;
    int irej;

    void ghostInsert(vector<particle> &, interact &, collision &, cell &);                 // Ghost insertion
    void getWidomResults();
    double exessChempot();
};

#endif
