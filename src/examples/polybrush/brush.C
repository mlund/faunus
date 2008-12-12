#include <iostream>
#include "analysis.h"
#include "container.h"
#include "potentials.h"
#include "mcloop.h"
typedef pot_hscoulomb2D T_pairpot;         // Specify pair potential
#include "markovmove.h"

using namespace std;

// Prelude

int main(int argc, char* argv[]) {
  string config = "brush.conf";
  if (argc==2) config = argv[1];
  slump slp;
  inputfile in(config);                // Read input file
  mcloop loop(in);                     // Set Markov chain loop lengths
  ratblock con(in.getflt("length_z"), in.getflt("length_xy"));      // Use a ratblock container
  canonical nvt;                       // Use the canonical ensemble
  pot_setup cfg(in);                   // Setup pair potential (default)
  interaction<T_pairpot> pot(cfg);     // Functions for interactions
  ioxyz xyz(con);
  ioaam aam(con);                      // Protein input file format is AAM
  iopov pov(con);

// Groups and particles
  brush brush(con, in.getint("Nchains"), in.getint("Mmonomers"), // Group for chains Note: fix the constructor
        in.getflt("rotdp"), in.getflt("kspring"), in.getflt("eqdist"),
        con.id(in.getstr("monomer")), con.len_xy, -con.len_z/2);
  macromolecule protein;               // Group for the protein
  protein.add( con,
      aam.load( in.getstr("protein")) );
  protein.center(con);                 // ..translate it to origo (0,0,0)
  protein.name=in.getstr("protein");     
  salt salt;                           // Group for salt and counter ions
  salt.add( con, in );                 //   Insert sodium ions


  particle tp;
  vector<particle> ts;

// Equlibrated?
  if (aam.load(con, "confout.aam")) {       // Load old config (if present)
    brush.loadgp("graftpoints.xyz");        // Load the graftpoints (not in particle vectors)
    protein.masscenter(con.p);
  } 

// Markovemoves
  brushmove bm(nvt, con, pot);         // Class for brush movements
  saltmove sm(nvt, con, pot);          // Class for salt movements
  sm.dp=in.getflt("saltdp");
  zmove z(nvt, con, pot, protein, 1);  // Class for protein z-movements (compare with free translation for efficiency?)
  z.dp=in.getflt("prot_zdp");
  z.min=in.getflt("z_min");
  z.max=in.getflt("z_max");
  chargereg tit(nvt,con,pot,salt,in.getflt("pH")); // Prepare titration. 
  macrorot mr(nvt, con, pot);          // Class for macromolecule rotation
  mr.dp=in.getflt("prot_rotdp");

// Tools of analysis
   histogram bdens(1,0,con.len_z);
   histogram enddens(1,0,con.len_z);
   histogram catdens(1,0,con.len_z);
   histogram andens(1,0,con.len_z);
   distributions dist(1,0,con.len_z);
   histogram protdens(1,0,con.len_z);

// The program
  double shint=0;
  pot_hs hs;
  for (int i=protein.beg; i<protein.end; i++)
    for (int j=i+1; j<=protein.end; j++) {
      double r2;
      r2=con.p[i].sqdist2D(con.p[j], con.len_xy);
      double& R2=r2;
      shint+=hs.hs(con.p[i], con.p[j], R2);
    }
  cout <<"shint = "<<shint<<endl;
  cout <<"brushenergy = "<< pot.brushenergy(con.p, brush)<<endl;
  cout <<"systemenergy = "<<pot.energy(con.p)<<endl;
  systemenergy sys(pot.energy(con.p) +pot.brushenergy(con.p, brush)-shint*pot.pair.f); // System energy analysis
  
  cout << con.info() << sys.info()     // Some information
       << pot.info();

  con.check_vector();

  for (int macro=1; macro<=loop.macro; macro++) {//Markov chain 
    for (int micro=1; micro<=loop.micro; micro++) {
      double Switch = slp.random_one(); {                 // Randomly chose move

        if(Switch<0.4)
          sys+=bm.move(brush);               // Crank and rattle chains
  
        if(Switch<0.5 &&Switch>=0.4) 
          sys+=sm.move(salt);                // Displace salt particles

        if(Switch<0.7 &&Switch>=0.5)         // Z-displacement of protein
          sys+=z.move(protein);

        if(Switch<0.8 &&Switch>=0.7) {         // Titrate proteins
          sys+=tit.titrateall();
          protein.charge(con.p);
          protein.dipole(con.p);
        }

        if(Switch<1 &&Switch>=0.8)
          sys+=mr.move(protein);
        
      }
      //Analysis
      if ( slp.random_one() <0.05) {
        for (int i = brush.beg; i<=brush.end; i++) {
          bdens.add(con.p[i].z+con.len_z/2);
        }
        for (int i = 0; i<brush.chains.size(); i++) {
          enddens.add(con.p[brush.chains[i].end].z+con.len_z/2);
        }
        for (int i =salt.beg;i<=salt.end; i++) {
          if(con.p[i].charge>0)
            catdens.add(con.p[i].z+con.len_z/2);
          if(con.p[i].charge<0)
            andens.add(con.p[i].z+con.len_z/2);
        }
        protdens.add(protein.cm.z+con.len_z/2);
        double sumq=0;
        for(int r=1;r<=con.len_z;r++) {
          for(int i=0;i<con.p.size();i++)
            if(con.p[i].z<r-con.len_z/2. && con.p[i].z>r-con.len_z/2.-1.)
              sumq+=con.p[i].charge;
          dist.add("Charge summation", r, sumq);
        }
        dist.add("Protein charge", protein.cm.z+con.len_z/2., protein.charge(con.p));    // Sample the protein charge
        dist.add("Protein dipole-z", protein.cm.z+con.len_z/2., protein.mu.z);           // - // - dipole z component
        dist.add("Protein dipole", protein.cm.z+con.len_z/2., protein.mu.len());         // - // - dipole 
        dist.add("Electrostatic energy", protein.cm.z+con.len_z/2.,                      // Sample total el.stat. energy
                  pot.energy(con.p)-shint*pot.pair.f-pot.brushenergy(con.p,brush));
        dist.add("Harmonic energy", protein.cm.z+con.len_z/2.,                            // Sample the bond energy (brush)
                  pot.brushenergy(con.p,brush));
      }

    }                                       // END of micro loop
    sys.update(pot.energy(con.p) + pot.brushenergy(con.p, brush)-shint*pot.pair.f);          // Update system energy
    aam.save("confout.aam", con.p);         // Save config. to disk
    xyz.save("confout.xyz", con.p);
    bdens.write("bdens.dat");
    cout << loop.timing(macro);             // Show progress
    con.check_vector();
    }                                         // END of macro loop

  cout << sys.info() << brush.info(con, con.len_xy*con.len_xy) <<bm.info()
       << salt.info(con) <<sm.info() <<protein.info(con) <<z.info() <<mr.info() <<tit.info();  // Print final results

// Termination, prepare restart, print analysis and pics
  aam.save("confout.aam", con.p);         // Save config. to disk
  xyz.save("confout.xyz", con.p);         // Save gp config. to disk
  ts.clear();
  for (int i=0; i<brush.chains.size(); i++) {
    tp=brush.chains[i].GP;
    ts.push_back(tp);
  }
  xyz.save("graftpoints.xyz", ts);
  for(int i=0;i<brush.chains.size();i++) {
    con.unboundary(brush.chains[i].GP, con.p[brush.chains[i].beg]);
    pov.connect(brush.chains[i].GP,con.p[brush.chains[i].beg],0.5);
    for(int j=brush.chains[i].beg;j<brush.chains[i].end;j++) {
      con.unboundary(con.p[j+1], con.p[j]);
      pov.connect(con.p[j],con.p[j+1],0.5); 
    }
  }
  pov.save("confout.pov", con.p);           
  bdens.write("bdens.dat");
  enddens.write("enddens.dat");
  catdens.write("catdens.dat");
  andens.write("andens.dat");
  dist.write("dist.dat");
  protdens.write("protedens.dat");
}



