/*! \page test_sqw_chains sqw_chains
 *
 * Simulate a number of flexible and titratable
 * polymers interacting via square wells
 * 
 *
 * \author Bjoern Persson
 * \date Lund, 2010
 * \include sqw_chains.cpp
 */
#include "faunus/faunus.h"
#include "faunus/energy/springinteraction.h"
#include "faunus/potentials/pot_sqwell.h"
#include "faunus/moves/reptation.h"
#include "faunus/moves/switch.h"

using namespace std;
using namespace Faunus;

int main() {
  // General setup
  cout << faunus_splash();
  slump slp;
  inputfile in("sqw_chains.conf");
  mcloop loop(in);
  canonical nvt;
  box con(in);
  springinteraction<pot_sqwell> pot(in);
  distributions dst;                     // Distance dep. averages

  // Handle polymers
  vector<polymer> polymers;              //the polymers
  vector<macromolecule> mpolymers;       //the polymers again but in a 
  point uv;                              //'shadow' vector 
#ifdef BABEL
  for(int i=0;i<in.getint("N_polymer"); i++){
    polymer pol;
    pol.babeladd( con, in );
    uv.ranunit(slp);
  //atom.reset_properties(con.p);
  //con.trial=con.p;
    pol.move(con, -pol.cm);              // Translate polymer to origo (0,0,0)
    pol.move(con, uv*con.len);           // .. and a random position
    pol.accept(con);                     // .. accept translation
    polymers.push_back(pol);
    mpolymers.push_back(pol);
  }
#endif
  crankShaft cs(nvt,con,pot,in);         // ...crankshaft
  macrorot mr(nvt, con, pot,in);         // ...rotate
  translate mt(nvt, con, pot, in);       // ...translate
  reptation rep(nvt, con, pot, in);      // ...reptate
  switch_type tit(nvt, con , pot, in);   // ...titrate
  clustertrans clt(nvt, con, pot, mpolymers);
  clt.dp=mt.dp;
  
  // File I/O
  io io;
  iopqr pqr;
  ioaam aam;
  ioxtc xtc(con.len);

  aggregation agg(con, polymers, in.getflt("cluster_def"));

  //Sloppy sync of radii and mw... blame it on openbabel
  for (int i=0; i<con.p.size(); i++){
    con.p[i].radius=atom[con.p[i].id].radius;
    con.p[i].mw=atom[con.p[i].id].mw;
    con.trial[i].radius=atom[con.p[i].id].radius;
    con.trial[i].mw=atom[con.p[i].id].mw;
  }

  // Load stored configuration?
  if(  aam.load(con, "confout.aam"))
    cout << "# Previus configuration loaded"<<endl;

  //Sloppy sync of radii and mw... blame it on openbabel
  for (int i=0; i<con.p.size(); i++){
    con.p[i].radius=atom[con.p[i].id].radius;
    con.p[i].mw=atom[con.p[i].id].mw;
    con.trial[i].radius=atom[con.p[i].id].radius;
    con.trial[i].mw=atom[con.p[i].id].mw;
  }
  // Calculate initial system energy
  systemenergy sys(
        pot.energy( con.p)); 

  cout << pot.energy(con.p)<<endl;
  pqr.save("confout.pqr",con.p);

  cout << con.info() << atom.info()
       << pot.info() 
       << in.info();

  int N=in.getint("N_polymer");

  //Parameters for input control of Markov chain
  double randy, picky, rot_f, crank_f, rep_f, trans_f, tit_f, clt_f, sum;
  rot_f=in.getflt("rot_f"), crank_f=in.getflt("crank_f"), rep_f=in.getflt("rep_f"), clt_f=in.getflt("clt_f");
  trans_f=in.getflt("trans_f"), tit_f=in.getflt("tit_f"),sum=crank_f+rep_f+rot_f+trans_f+tit_f;
  rot_f/=sum, crank_f/=sum, rep_f/=sum, trans_f/=sum, tit_f/=sum, clt_f/=sum;

  while ( loop.macroCnt() ) {
    while ( loop.microCnt() ) {
      randy=slp.random_one();
      picky=slp.random_one();
      if(trans_f>randy){
        sys+=mt.move(polymers[int(picky*N)]);         // translate polymers
      }
      if(trans_f+rot_f>randy && trans_f<randy){
        polymers[int(picky*N)].masscenter(con);
        sys+=mr.move(polymers[int(picky*N)]);         // rotate polymers
      }
      if(trans_f+rot_f+crank_f>randy && trans_f+rot_f<randy){
        sys+=cs.move(polymers[int(picky*N)]);         // crankshaft
      }
      if(trans_f+rot_f+crank_f+rep_f>randy && trans_f+rot_f+crank_f<randy){
        sys+=rep.move(polymers[int(picky*N)]);        // reptation
      }
      if(1.0-tit_f-clt_f<randy && 1.0-clt_f>randy){
        sys+=tit.titrateall();                        // titrate titratable sites
      }
      if(1.0-clt_f<randy){
        sys+=clt.move(mpolymers);                     // cluster translation of 'all' polymers
      }
    
      if (slp.random_one()>0.99)                      // On the fly analysis
        agg.count();
      //if (slp.random_one()>0.95) {
        //xtc.setbox(con.len,con.len,(*zhalfPtr)*2);
        //xtc.save( "traj.xtc", con.p, vg );
      //}
    } // END of micro loop

    sys.update(
        pot.energy( con.p));

    cout << loop.timing() << "Energy drift = " << sys.cur-sys.sum << " kT" << endl;     

    // Write files to disk
    std::ostringstream o;
    for (int i=0; i<polymers.size(); i++)
      o << polymers[i].getVMDBondScript() <<endl;
    io.writefile("vmdbonds.tcl", o.str());
    pqr.save("confout.pqr",con.p);
    aam.save("confout.aam", con.p);
    agg.write("aggregation.dat");

  } // END of macro loop and simulation

  cout << sys.info() << loop.info()
       << cs.info() << mr.info() << mt.info()
       << tit.info() << rep.info() <<clt.info();
}
