/*! 
 *
 * This example program calculates the probability of finding 
 * holes in a random point lattice
 * 
 *
 * \author Bjoern Persson
 * \date Malmoe, 2009
 * \include hole1.cpp
 */
#include "faunus/faunus.h"

using namespace std;
using namespace Faunus;                 // Access to Faunus classes

int main() {
  cout << faunus_splash();              // Show Faunus information
  inputfile in("hole.conf");            // Read input file
  mcloop loop(in);                      // Set loop lengths
  xyplane con( in );                      // Planare container
//  con.setlen( in.getflt("boxlen",-1) * sqrt(acos(-1)) );                      // Planare container
  
  // Variables

  vector<point> p(in.getint("points", 1000));

  // Analysis

  diskoverlap anl(p);

  cout << "#  Hole-seach program, finds the probability of inserting elipsoids"<<endl
       << "#  in a random point grid."<<endl
       << "#  ----------------------------------------------------------------"<<endl
       << "#  Point density (N/sigma-cube) = "<< p.size()/con.len/con.len<<endl
       << "#  Lengt of square side (sigma) = "<< con.len<<endl
       << "#  Number of points             = "<< p.size()<<endl<<"#"<<endl; 

  // Dummy integration
  while ( loop.macroCnt() ) {            
    while ( loop.microCnt() ) {
      con.randompos(p);
      anl.check(p);
    }                                   // END of micro loop
    anl.blockavg();                     // Collect block averaged data (based on the number of macrosteps)
    cout << loop.timing();              // Show progres
  }                                     // END of macro loop and simulation

  cout << loop.info()<<anl.info(); // Final information and results!
}

