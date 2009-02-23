/*!\page test_spce SPCE/E
 *
 * \author Mikael Lund and Bjorn Persson
 * \include spce.cpp
 */
#include "faunus/faunus.h"
#include "faunus/potentials/pot_rfield.h"
#include "faunus/average_vec.h"

using namespace Faunus;
using namespace std;

class ionpmf {
  private:
    Faunus::slump slump;
    particle a,b;
  public:
    double dr,rmax;
    vector< average<double> > pmf,u1,u2;
    ionpmf() {
      rmax=11.;
      dr=0.1;
      a.charge=0.1;
      b.charge=-a.charge;
      pmf.resize( int(rmax/dr+.5) );
      u1.resize( int(rmax/dr+.5) );
      u2.resize( int(rmax/dr+.5) );
    }
    void insert(container &c, energybase &pot, double r=5.0) {
      double u=0,ut1,ut2;
      a.x=2.*r*slump.random_half()+0.000001;
      b.x=-a.x;
      r=abs(a.x-b.x);
      ut1=a.charge*pot.potential(c.p, a);
      ut2=b.charge*pot.potential(c.p, b);
      u=ut1+ut2;
      pmf[ int( r/dr+.5) ] +=exp(-u);
      u1[ int( r/dr+.5) ] +=exp(-ut1);
      u2[ int( r/dr+.5) ] +=exp(-ut2);
    }
};

class dgsolv {
  private:
    slump slp;
    particle a,b;
  public:
    int ref, cnt;
    vector< average_vec<double> > edgp;      
    vector< average_vec<double> > edgm;      
    vector< average_vec<double> > phi2;
    vector<int> itr;                  // particle vectro index
    vector<double> du;
    average_vec<double> Pref, Mref; //, Mrefstd, Prefstd;
    average<double> Prefaux, Mrefaux, Mrefstd, Prefstd;
    dgsolv(vector<int> pari ) {
      ref=pari[0];                    // use first as reference
      pari.erase( pari.begin() );
      itr=pari;
      edgp.resize( pari.size() );
      edgm.resize( pari.size() );
      phi2.resize( pari.size() );
      du.clear();
      a.charge=0.1;
      b.charge=-a.charge;
      cnt=0;
      std::cout << endl <<"# Maximum number of values in Bennet vector is "<<du.max_size() <<endl <<endl;

   //   Prefstd.maxcnt=1;
   //   Mrefstd.maxcnt=1;
    }
    void sample( container &c, sphericalimage<pot_test>  &pot) {
      cnt++;
      double uref, uitr, pref, pitr;
      a.charge=c.p[ref].charge;
      c.p[ref].charge=0.1;
      uref=pot.elenergy(c.p, ref);
      Pref+=exp(-uref);       //sign
      Mref+=exp(uref);
      Prefaux+=exp(-uref);    //sign
      Mrefaux+=exp(uref);
      c.p[ref].charge=a.charge;
      pot.updateimg(c.p[ref], ref);
      if ( du.size() < du.max_size() )
        du.push_back(uref);   //sign
      pref=pot.potential(c.p, ref);
      for (int i =0; i<itr.size(); i++) {
        a.charge=c.p[itr[i]].charge;
        c.p[itr[i]].charge=0.1;
        uitr=pot.elenergy(c.p, itr[i]);
       // std::cout << uitr<<endl;
        c.p[itr[i]].charge=a.charge;
        pot.updateimg(c.p[itr[i]], i);
        pitr=pot.potential(c.p ,i);
        edgp[i]+=exp(-uitr);
        edgm[i]+=exp(uitr);
        phi2[i]+=pref*pitr;
      }
    }
    void write() {
      // Calculate std
      Mrefstd+=-log(Mrefaux.avg());
      Prefstd+=-log(Prefaux.avg());
      Mrefaux.reset();
      Prefaux.reset();
      //
      std::cout << "#    Free energy of solvation analysis, ref to particle "<<ref<<endl
                << "#    "<<cnt<<" insertions at "<<itr.size()<<" points with +/-0.1e"<<endl
                << "#                         Plus          |   Minus"<<endl
                << "#     Ref. Part.   deltaG= "<<-log(Pref.avg()) <<" ("<< Prefstd.avg()<<" +/- "<<Prefstd.stdev()<<")"
                << " | "<<-log(Mref.avg()) <<" ("<<Mrefstd.avg()<<" +/- "<<Mrefstd.stdev()<<endl;
      for (int i =0; i<itr.size(); i++)
        std::cout << "#     Particle " <<itr[i] <<", deltaG = " <<-log(edgp[i].avg()) <<"      " //<< std::abs(log(edgp[i].stdev()))
                  <<" | "<<-log(edgm[i].avg())<<" "<<endl //<<std::abs(log(edgm[i].stdev()))<<endl
                  <<"#      Phiref dot Phi[i] = "<< phi2[i].avg() <<endl <<endl;
      if ( du.size() == du.max_size() )
        std::cout <<"# Bennet vector is full !!!"<<endl<<endl;
 
      ofstream dudist("widom.dat");
      if (dudist) {
        int s=du.size();
        for (int t=0;t<s;t++)
          dudist << t <<"  "<< du[t]<< endl;
        dudist.close();
      }

    }
};

class acidsolv {
  private:
    double H1,O1,O2,C1;
    int cnt;
  public:
    double uref, uitr, sc;
    vector< average<double> > dg;
    acidsolv(vector<group> g, double scale=0.1) {
      sc=scale;
      H1=-0.433;
      O1=-0.487;
      O2=-0.540;
      C1= 0.460;
      dg.resize(g.size());
      cnt=0;
    }
    void sample(container &c, energybase &pot, vector<group> &g) {
      cnt++;
      unprotonate(c,g[0]); //use first acid as ref
      uref=pot.energy(c.trial, g[0])-pot.energy(c.p, g[0]);
      protonate(c,g[0]);   //reset the charges
      for (int i=1; i<g.size(); i++) {
        unprotonate(c,g[i]);
        uitr=pot.energy(c.trial, g[i])-pot.energy(c.p, g[i]);
        protonate(c,g[i]); //reset the charges
        dg[i]+=exp(-uitr+uref);   //sample exp(uitr-uref);
      }
    }
    void unprotonate(container &con, group &g) {
      con.trial[g.beg].charge+=sc*C1;
      con.trial[g.beg+1].charge+=sc*O2;
      con.trial[g.beg+2].charge+=sc*O1;
      con.trial[g.beg+3].charge+=sc*H1;
    }
    void protonate(container &con, group &g) {
      con.trial[g.beg].charge-=sc*C1;
      con.trial[g.beg+1].charge-=sc*O2;
      con.trial[g.beg+2].charge-=sc*O1;
      con.trial[g.beg+3].charge-=sc*H1;
    }
    void write() {
      std::cout << "#     Free energy of solvation analysis, ref to first acid "<<endl
                << "#    "<<cnt<<" insertions "<<endl
                << "#                         "<<endl;
      for (int i =1; i<dg.size(); i++){
        std::cout << "#     Acid # " <<i+1 <<", deltaG = " <<-log(dg[i].avg()) <<" "<< std::abs(log(dg[i].stdev()))
                  <<endl;
      }
      std::cout <<endl;
    }

};

int main(int argc, char* argv[]) {
  slump slump;
  cout << faunus_splash();             // Faunus spam

  inputfile in("spce.conf");           // Read input file
#ifdef SPCE_MINIMAGE
  box con(in);                        // Use a spherical simulation container
  interaction<pot_testminim> pot(in);       // Specify pair potential
#else
  cell con(in);                        // Use a spherical simulation container
//  interaction<pot_test> pot(in);       // Specify pair potential
  sphericalimage<pot_test> pot(in);       // Specify pair potential
#endif
  io io;
  ioaam aam(con.atom);                 // Protein input file format is AAM
  iogro gro(con.atom,in);
  mcloop loop(in);                     // Set Markov chain loop lengths
  canonical nvt;                       // Use the canonical ensemble
  pot_rfield rfield(in);
  pot.pair.init(con.atom);
  ionpmf ip;

  macrorot mr(nvt, con, pot);          // Class for molecular rotation
  translate mt(nvt, con, pot);         // Class for molecular translation
  transrot mrt(nvt, con, pot);       // Class for combined molecular trans. and rot.

  // Load central protein
  macromolecule protein;               // Group for the protein
  protein.add( con,
      aam.load(in.getstr("protein"))); // Load protein structure
//  point skitp;
//  skitp=-protein.cm;
//  skitp.x+=0.001,
  protein.move(con, -protein.cm);      // ..translate it to origo (0,0,0)
  protein.accept(con);                 // ..accept translation

  // Add salt
  salt salt;                           // Group for salt and counter ions
  salt.add( con, in );                 //   Insert sodium ions
  saltmove sm(nvt, con, pot);          // Class for salt movements

  // Add water
  molecules sol(3);                    // We want a three point water model
  sol.name="SPC/E Solvent";
  vector<particle>
    water = aam.load("water.aam");     // Load bulk water from disk (typically from MD)
  con.atom.reset_properties(water);    // Set particle parameters according to Faunus
  sol.add(con,water,sol.numatom);      // Inject water into the cell - avoid salt and protein overlap
  water.clear();                       // Free the (large) bulk water reservoir

  std::cout <<"# Objects initiated"<< endl;


  // Measure electric potential at all protein atoms
/*  pointpotential ppot;
  for (int i=protein.beg; i<=protein.end; i++)
    ppot.add(con.p[i], con.atom[con.p[i].id].name );
*/
  // Distribution functions
  FAUrdf spccell(float(0.2), float(50.));
  FAUrdf nacell(float(0.2), float(50.));
  FAUrdf saltrdf(con.atom["NA"].id,con.atom["CL"].id,0.2,20.);
  FAUrdf catcat( con.atom["NA"].id,con.atom["NA"].id,0.2,20.);
  FAUrdf spcrdf( con.atom["OW"].id,con.atom["OW"].id,0.2,20.);
  FAUrdf acidwrdf( con.atom["C1"].id,con.atom["OW"].id,0.2,20.);

  distributions radorder(float(0.2), float(0.), float(in.getflt("distmax")));

  histogram kfuck(float(0.02), float(0.), float(40.));
  histogram dip(float(0.05), float(0.), float(500.));
  gfactor gf(con, sol);
  average_vec<double> Pexcess(1);
  average_vec<double> Mexcess(1);

  // Solvation analysis-----------------------
  vector<int> dgpoints;
  dgpoints.clear();
  for (int i=0; i<con.p.size(); i++)
    if (con.atom["GLU"].id==con.p[i].id)
      dgpoints.push_back(i);
  dgsolv solvation(dgpoints);

/*  vector<group> acids;
  group acid;
  acids.clear();
  for (int i=0; i<con.p.size(); i++)
    if (con.atom["C1"].id==con.p[i].id){
      acid.beg=i;
      acid.end=i+3;            //Every acid is assumed to consist of four atoms, carbonyl carbon, O=, O- and H-
      acids.push_back(acid);   //in the order above, see protonate in acidsolv i.g. Faunus atoms are C1,O2,O1 and H1
    }
  acidsolv acidsolv(acids);*/
  //-------------------------------------------

  mr.dp=in.getflt("mrdp",0.6);
  mt.dp=in.getflt("mtdp",0.3); //0.6
  sm.dp=0.4; //0.4
  mrt.dpr=in.getflt("mrtdpr",0.0);
  mrt.dpt=in.getflt("mrtdpt",0.3);

  aam.load(con, "confout.aam");        // Load old config (if present)

#ifdef SPCE_MINIMAGE
#else
  pot.updateimg(con.p);                // Initiate image charges
#endif

  systemenergy sys( 
      pot.internal(con.p, sol, sol.numatom) +
      //pot.internal(con.p, salt) +
      //pot.energy(con.p, protein, salt) +
      pot.energy(con.p, sol)
      );
  io.splash("spce.dat");

  cout << in.info() << salt.info()
    << con.info() << con.atom.info()
    << pot.info() << sol.info() << rfield.info();

  io.splash("../../../misc/faunatoms.dat");

  group head;                            // Group of first three particles [0:2]
  head.set(0,2);                         // (Used to swap w. SPC for parallization reasons)
  macromolecule m, n;
  double i_norm, scalar, drlen, dr1o, dr1n, dr2o, dr2n;
  point dr;
  int solsize, steptype;
  int eprint;
  eprint=int(0.001*in.getflt("microsteps"));
  point origo;
  vector<particle> check=con.p;
  point trace = con.p[sol.beg];
  average< double > imu;
  int cnt=0;

  // Switch paramters
  double trans, rots, trs, sums, randy;
  int switcher=-1;
  sums=in.getflt("t")+in.getflt("rot")+in.getflt("tr");
  trans=in.getflt("t")/sums;
  rots=in.getflt("rot")/sums;
  trs=in.getflt("tr")/sums;
  //

  for (int i=0; i<con.p.size(); i++)
    if (con.collision(con.p[i]))
      std::cout << "# Particle "<<i<<" is overlaping the cavity !"<<endl;

  while ( loop.macroCnt() ) {            // Markov chain 
    while ( loop.microCnt() ) {
      m=sol[ sol.random() ];
      m.cm=m.cm_trial=con.p[m.beg];
      randy=slump.random_one();
      if (randy<trans)
        switcher=0;
      if (randy>trans && randy <trans+rots)
        switcher=1;
      if (randy>trans+rots)
        switcher=2;
      switch (switcher) {              // Randomly choose move
        case 0:
          sys+=mt.move(m);             // Translate solvent
          cnt++;
          steptype=0;
          break;
        case 1:
          sys+=mr.move(m);             // Rotate solvent
          cnt++;
          steptype=1;
          break;
        case 2:
          sys+=mrt.move(m);
          cnt++;
          steptype=2;
          break;
        case 3:
         // sys+=sm.move(salt,1);          // Displace a salt particle
          break;
      }
#ifdef SPCE_MINIMAGE
#else
      pot.updateimg(con.p, m);
#endif

      if (slump.random_one()<0.05) {
//        ppot.sample(con,pot);
//        ip.insert(con,pot);
#ifdef SPCE_MINIMAGE
#else
        solvation.sample(con, pot);
#endif
//        acidsolv.sample(con, pot, acids); 
      }
      if (slump.random_one()<0.05) {
        saltrdf.update(con, salt);
        catcat.update(con, salt);
        nacell.update(con, origo, "NA");
        acidwrdf.update(con);
#ifdef SPCE_MINIMAGE
#else
        imu+=pot.image(con.p);
#endif
      }
      if (slump.random_one()<0.05) {  
        solsize=sol.size()/3;
//        for (int i=0; i<acids.size(); i++)
          for (int j=0; j<sol.size()/3; j++) {
            m=sol[j];
            m.cm=con.p[m.beg];
//            dr=con.vdist(m.cm,con.p[acids[i].beg]);
//            drlen=dr.len();
//            i_norm=1./(drlen*m.dipole(con.p));
//            scalar=dr.dot(m.mu)*i_norm;
//            radorder.add("Acid-Water, cos(theta)", drlen,scalar );
//            radorder.add("Acid-Water, cos(theta)^2", drlen, scalar*scalar);
            i_norm=1./(m.cm.len()*m.dipole(con.p));
            scalar=m.cm.dot(m.mu)*i_norm; 
            radorder.add("Scalar", m.cm.len(),scalar );
            radorder.add("Scalar squared", m.cm.len(), scalar*scalar); 
          }
        kfuck.add(gf.add(con.p, sol));
        m.beg=sol.beg;
        m.end=sol.end;
        dip.add(m.dipole(con.p));
      } 
      if (slump.random_one()<-0.01) {
        solsize=sol.size()/3;
        for (int j=0; j<1000; j++) {
        m=sol[ sol.random() ];
        do {n=sol[ sol.random() ]; }
          while (n.beg==m.beg);
        m.cm=con.p[m.beg];
        n.cm=con.p[n.beg];
        dr=con.vdist(m.cm,n.cm);
        drlen=dr.len();
        i_norm=1./(drlen*m.dipole(con.p));
        scalar=dr.dot(m.mu)*i_norm;
        radorder.add("Water-Water, cos(theta) radius vect with dip", drlen,scalar );
        radorder.add("Water-Water, cos(theta)^2 radius vect with dip", drlen, scalar*scalar);
        i_norm=1./(m.dipole(con.p)*n.dipole(con.p)); 
        scalar=n.mu.dot(m.mu)*i_norm;
        radorder.add("Water-Water, cos(theta) dip vect with dip", drlen,scalar );
        radorder.add("Water-Water, cos(theta)^2 dip vect with dip", drlen, scalar*scalar);
        }
      }
      if (slump.random_one()<0.001) {
        spccell.update(con, origo, "OW");
        spcrdf.update(con, sol);
      }
      if ( cnt%eprint==0 )
         sys.track();
    }//end of micro-loop 

    nacell.write("rdf-cell-NA.dat");
    spccell.write("rdf-cell-OW.dat");
    spcrdf.write("rdf-OW-OW.dat");
    catcat.write("rdf-Na-Na.dat");
    saltrdf.write("rdf-Na-Cl.dat");
    radorder.write("mu-OW.dat");
    acidwrdf.write("rdf-C1-OW.dat");
    dip.write("dipoledist.dat");
    sys.write();

    Pexcess+=-log(solvation.Pref.avg());
    Mexcess+=-log(solvation.Mref.avg());

    sys.update( 
        pot.internal(con.p, sol, sol.numatom) +
        //pot.internal(con.p, salt) +
        //pot.energy(con.p, protein, salt) +
        pot.energy(con.p, sol)
        );

    gro.save("confout.gro", con.p);      // Save config. to disk
    aam.save("confout.aam", con.p);      // Save config. to disk
    cout << loop.timing()                // Show progress
         << "#   Energy (kT): sum average drift = "
         << sys.sum << " " << sys.uavg.avg() << " "
         << std::abs(sys.cur-sys.sum) << endl;
//    cout << ppot.info();
    cout << "#   First water oxygen movement = " << con.dist(trace, con.p[sol.beg]) << endl;
/*    std::cout << "#   Reaction field energy  = "<<imu.avg()<<" , stdev = "<<imu.stdev() <<endl<<endl
              << "#   Excess chemical potential ref (p/n) = "<<Pexcess.avg() <<" / "<<Mexcess.avg()<<" , "
              << Pexcess.stdev()<<" , "<<Mexcess.stdev()<<endl<<endl;                  
*/

    solvation.write();
//    acidsolv.write();
    std::cout << gf.info();
    kfuck.write("kirkwood.dat");

    point slask;
    for (int i=0; i<con.p.size(); i++) {
      slask=con.vdist(con.p[i], check[i]);
      if (slask.len()<1e-14)
        std::cout <<"#  Particle "<<i<<" is static !!!"<<endl;
    }
    cout << endl<<endl;

  }//end of macro-loop
  cout << sys.info()
       << salt.info(con)
       << protein.info() << loop.info()  // Print final results
       << sm.info() << mr.info() << mt.info() <<mrt.info() <<pot.info();

#ifdef SPCE_MINIMAGE
#else
  io.writefile("img.dat", pot.printimg());
#endif

  // Print dipole intertion results
/*  particle a,b;
  a.charge=0.1;
  b.charge=-a.charge;
  for (float r=0; r<10.; r+=0.1) {
    a.x=.5*r;
    b.x=-a.x;
    cout
      << r << " "
      << -log(ip.pmf[ int(r/ip.dr+.5) ].avg()) << " "
      << -log(ip.u1[ int(r/ip.dr+.5) ].avg()) << " "
      << -log(ip.u2[ int(r/ip.dr+.5) ].avg()) << " ";
    double u = rfield.f * (
          rfield.selfenergy(a) + rfield.selfenergy(b)
          + rfield.pairpot(a,b) );
    cout << u << " " << u - rfield.f*a.charge*b.charge/a.dist(b) << endl;
  }
*/
}

