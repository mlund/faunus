/*! \file Example how to simulate a protein at a surface
 *  \example wall.C
 */
#include <iostream>
#include "classes/io.h"
#include "classes/analysis.h"
#include "classes/container.h"
#include "classes/potentials.h"
#include "classes/countdown.h"
typedef pot_coulomb T_pairpot;
#include "classes/markovmove.C"
#include "classes/histogram.h"
#include "classes/io.h"

using namespace std;

int main() {

  double prot_dp=4;


  /*********************************
  SYSTEM
  *********************************/


  slump2 slp;
  cylinder cylinder(300.,50);		// Constuct the cylinder
  iopov povray(cylinder);		// Some pictures
  povray.cylinder(100.,50.);
  canonical nvt;			// We choose a canonical ensamble
  pot_setup cfg;			// Pair potential (change this)
  interaction<T_pairpot> pot(cfg);	// Functions of interactions
  countdown<int> clock(10);		// Estimate remaining simulation time

  macromolecule protein;		// Group for the protein
  ioaam aam(cylinder);			// File format of protein
  protein.add( cylinder, aam.load(
    "protein-example.aam"));		// Lock and load...
  protein.move(cylinder, -protein.cm);  // Place it in origo... 
  protein.accept(cylinder);
  protein.zmove(cylinder, cylinder.len/2.);      // ...and push it to the middle
  protein.accept(cylinder);                      // (accept start position)
  zmove zm(nvt, cylinder, pot, protein, 0.7);		

/*  group salt;				     // Group for mobile ions
  salt.add( cylinder, particle::NA, 11);     //+19
  salt.add( cylinder, particle::CL, 11);
  saltmove sm(nvt, cylinder, pot);*/	     // Class for salt movment?
//  chargereg tit(nvt,cylinder,pot,salt,7);    // Prepare to titrate, pH 7
  systemenergy sys(pot.energy(cylinder.p));  // System energy analysis

  cout << cylinder.info() /*<< tit.info()*/;     // Check out system
  cout << protein.masscenter(cylinder.p).z<<endl<<endl;
  /****************************************
  ANALYSIS                                
  ****************************************/

  histogram ldf(0.5,0,100);

  /****************************************
  START INTEGRATION
  ****************************************/


  for (int macro=1; macro<=10; macro++) {       // Markov chain
    for (int micro=1; micro<=3e4; micro++) {
//      sm.move(salt);                            // Displace salt particles
/*      if (tit.titrateall()) {                   // Titrate groups
        protein.charge(cylinder.p);                 // Re-calc. protein charge
        protein.dipole(cylinder.p);                 // Re-calc. dipole moment
      }*/
      zm.move(protein);                             // Translate protein
      ldf.add(protein.masscenter(cylinder.p).z); 
      sys+=/*sm.du+tit.du +*/zm.du;                   // Keep system energy updated
      cout << protein.masscenter(cylinder.p).z<<endl;
    }
    cout << "Macro step " << macro << " completed. ETA: " << clock.eta(macro);
    sys.update(pot.energy(cylinder.p));
  }

  /****************************************
  STOP CHAIN AND PRINT 
  ****************************************/
  protein.save("coord", cylinder);
  cout << sys.info() << zm.info() << /*sm.info() << tit.info()*/ // More information...
    /*<< salt.info() << */protein.info();


  io gofr;
//  gofr.writefile("wall.plot", ldf.show());
  
  ldf.show();

  povray.save("protein-example.pov", cylinder.p);   // Save POVRAY file
}

