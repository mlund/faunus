#ifndef HYPERSPHERE
#define HYPERSPHERE
#endif

/*
*/
#include "faunus/faunus.h"
#include <iostream>
#include <sstream>
#include <string>
#include <stdexcept>

using namespace Faunus;
using namespace std;

class hyperrdf : public FAUrdf {
  public:
    double R;
    hyperrdf(short species1, short species2, float res, float xmax) :
  FAUrdf(species1, species2, res, xmax) { };
    float get(float x) {
      float simsize=R*R*R*acos(-1.)*acos(-1.)*2.0,
            fact=2*std::acos(-1.)*R*R*R,
            volnew=((x+xres)/R-std::sin((x+xres)*2./R)/2.),
            volold=(x/R-std::sin(x*2./R)/2.);
      return histogram::get(x) * simsize / (volnew-volold) / fact;
    }
};

int main() {
  cout << faunus_splash();
  slump slump;                            // A random number generator
  inputfile in("dipolarfluid.conf");      // Read input file
  hypersphere con(in);                    // We want a hypersphere
  canonical nvt;                          // Use the canonical ensemble
  //interaction<pot_hypersphere> pot(in);   // 
  mcloop loop(in);                        // Keep track of time and MC loop

  hypergroup salt;                        // Group for dipoles
  hyperrdf saltrdf(atom["NA"].id, atom["NA"].id, .5, acos(-1.)*con.r);
  saltrdf.R=con.r;

  salt.add(con, atom["NA"].id, 1000);
  
  cout << con.info() << in.info();
  
  for (int macro=1; macro<=loop.macro; macro++) {//Markov chain 
    for (int micro=1; micro<=loop.micro; micro++) {
      switch (rand() % 2) {
        case 0:
          for (int i=0; i<con.p.size(); i++) {
            con.randompos(con.p[i]);
            con.trial[i]=con.p[i];
          }
          saltrdf.update(con);
          break;
        case 1:
          for (int i=0; i<con.p.size(); i++) {
            salt.displace(con, 0.1);
          }
          break;
      }
    } // End of inner loop
    cout << loop.timing();
  } // End of outer loop
  saltrdf.write("blah.dat");
}

