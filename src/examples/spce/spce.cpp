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

int main(int argc, char* argv[]) {
  slump slump;
  cout << faunus_splash();             // Faunus spam

  inputfile in("spce.conf");           // Read input file
#ifdef SPCE_MINIMAGE
  box con(in);                         // Use a spherical simulation container
  interaction<pot_testminim> pot(in);  // Specify pair potential
#else
  cell con(in);                        // Use a spherical simulation container
  sphericalimage<pot_test> pot(in);    // Specify pair potential
#endif
  io io;
  ioaam aam(con.atom);                 // Protein input file format is AAM
  iogro gro(con.atom,in);
  mcloop loop(in);                     // Set Markov chain loop lengths
  canonical nvt;                       // Use the canonical ensemble
  pot_rfield rfield(in);
  pot.pair.init(con.atom);

  macrorot mr(nvt, con, pot);          // Class for molecular rotation
  translate mt(nvt, con, pot);         // Class for molecular translation
  transrot mrt(nvt, con, pot);         // Class for combined molecular trans. and rot.
  multtr multtr(nvt, con, pot, in.getint("multinum"));
                                       // Class for mulitple combined molecualr trans. and rot.

  // Add water
  molecules sol(3);                    // We want a three point water model
  sol.name="SPC/E Solvent";
  vector<particle>
    water = aam.load("water.aam");     // Load bulk water from disk (typically from MD)
  con.atom.reset_properties(water);    // Set particle parameters according to Faunus
  sol.add(con,water,sol.numatom);      // Inject water into the cell - avoid salt and protein overlap
  water.clear();                       // Free the (large) bulk water reservoir

  std::cout <<"# Objects initiated"<< endl;


  // Distribution functions and analysis
  FAUrdf spccell(float(0.2), float(50.));
  FAUrdf nacell(float(0.2), float(50.));
  FAUrdf saltrdf(con.atom["NA"].id,con.atom["CL"].id,0.2,20.);
  FAUrdf catcat( con.atom["NA"].id,con.atom["NA"].id,0.2,20.);
  FAUrdf spcrdf1( con.atom["OW"].id,con.atom["OW"].id,0.2,20.);
  FAUrdf spcrdf2( con.atom["OW"].id,con.atom["HW1"].id,0.2,20.);
  FAUrdf spcrdf3( con.atom["HW1"].id,con.atom["HW1"].id,0.2,20.);
  FAUrdf spcinner1( float(0.2), float(50.), float(1.));
  FAUrdf spcinner2( float(0.2), float(50.), float(1.));
  FAUrdf spcinner3( float(0.2), float(50.), float(1.));

  distributions radorder(float(0.2), float(0.), float(in.getflt("distmax")));

  histogram kfuck(float(0.02), float(0.), float(40.));
  histogram dip(float(0.2), float(0.), float(500.));
  histogram dipc(float(0.2), float(-500.), float(500.));
  gfactor gf(con, sol);

  average<double> Phi11, Phi22, Phi12;
  //-------------------------------------------

  // Parametrize markov chain
  mr.dp=in.getflt("mrdp",0.6);
  mt.dp=in.getflt("mtdp",0.3); //0.6
  mrt.dpr=in.getflt("mrtdpr",0.0);
  mrt.dpt=in.getflt("mrtdpt",0.3);
  multtr.dpr=in.getflt("multrdp");
  multtr.dpt=in.getflt("multtdp");

  aam.load(con, "confout.aam");        // Load old config (if present)

#ifdef SPCE_MINIMAGE
#else
  pot.updateimg(con.p);                // Initiate image charges
#endif

  systemenergy sys( 
      pot.internal(con.p, sol, sol.numatom) +
      pot.energy(con.p, sol)
      );
  if (in.getboo("splash")==true)
    io.splash("README");

  cout << in.info() 
    << con.info() << con.atom.info()
    << pot.info() << sol.info() << rfield.info();
  if (in.getboo("splash")==true)
    io.splash("../../../misc/faunatoms.dat");

  // Help variables
  group head;                            // Group of first three particles [0:2]
  head.set(0,2);                         // (Used to swap w. SPC for parallization reasons)
  macromolecule m, n;
  double i_norm, scalar, drlen, dr1o, dr1n, dr2o, dr2n, mdip;
  point dr, p1, p2;
  int solsize, steptype;
  int eprint;
  eprint=int(0.001*in.getflt("microsteps"));
  point origo;
  vector<particle> check=con.p;
  point trace = con.p[sol.beg];
  average< double > imu;
  int cnt=0;
  vector<int> o;

  // Switch paramters
  double trans, rots, trs, multi, sums, randy;
  int switcher=-1;
  sums=in.getflt("t")+in.getflt("rot")+in.getflt("tr")+in.getflt("multi");
  trans=in.getflt("t")/sums;
  rots=in.getflt("rot")/sums;
  trs=in.getflt("tr")/sums;
  multi=in.getflt("multi")/sums;
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
      if (randy>trans+rots &&randy <trans+rots+trs)
        switcher=2;
      if (randy> trans+rots+trs)
        switcher=3;
      switch (switcher) {              // Randomly choose move
        case 0:
          sys+=mt.move(m);             // Translate solvent
          cnt++;
          steptype=0;
#ifdef SPCE_MINIMAGE
#else
          pot.updateimg(con.p, m);
#endif
          break;
        case 1:
          sys+=mr.move(m);             // Rotate solvent
          cnt++;
          steptype=1;
#ifdef SPCE_MINIMAGE
#else
          pot.updateimg(con.p, m);
#endif
          break;
        case 2:
          sys+=mrt.move(m);
          cnt++;
          steptype=2;
#ifdef SPCE_MINIMAGE
#else
          pot.updateimg(con.p, m);
#endif
          break;
        case 3:
          o=sol.pick(in.getint("multinum"));
          sys+=multtr.move(sol, o);
#ifdef SPCE_MINIMAGE
#else
          pot.updateimg(con.p, sol, o);
#endif
          break;
      }
      if (slump.random_one()<0.05) {
#ifdef SPCE_MINIMAGE
#else
        imu+=pot.image(con.p);
#endif
      }
      // Update hist. of dipolemoment, dipolar correlation of spce  and sample kirkwood factor
      if (slump.random_one()<0.005) {  
        solsize=sol.size()/3;
          for (int j=0; j<sol.size()/3; j++) {
            m=sol[j];
            m.cm=con.p[m.beg];
            i_norm=1./(m.cm.len()*m.dipole(con.p));
            scalar=m.cm.dot(m.mu)*i_norm; 
            radorder.add("Scalar", m.cm.len(),scalar );
            radorder.add("Scalar squared", m.cm.len(), scalar*scalar); 
          }
        kfuck.add(gf.add(con.p, sol));
        m.beg=sol.beg;
        m.end=sol.end;
        m.cm.x=0;
        m.cm.y=0;
        m.cm.z=0;
        mdip=m.dipole(con.p);
        dip.add(mdip);
        dipc.add(m.mu.x);
        dipc.add(m.mu.y);
        dipc.add(m.mu.z);
      } 
      // Update hist. spce dipole and radial vector correlation
      if (slump.random_one()<0.001) {
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
      // Sample potential correlations
#ifdef SPCE_MINIMAGE
#else
      if (slump.random_one()<0.001)
        for (int i=0; i<100; i++) {
          p1.ranunit(slump);
          p2.ranunit(slump);
          p1=p1*in.getflt("cellradius")*0.5;
          p2=p2*in.getflt("cellradius")*0.5;
          radorder.add("Potential correlation between two points in a central cavity of radius/2",
                       p1.dist(p2), pot.potential(con.p,p1)*pot.potential(con.p,p2));

        }
#endif
      // Update water rdf:s
      if (slump.random_one()<0.001) {
        spccell.update(con, origo, "OW");
        spcrdf1.update(con, sol);
        spcrdf2.update(con, sol);
        spcrdf3.update(con, sol);
        for (int i=0; i<10; i++) {
        m=sol[ sol.random() ];
          if (con.p[m.beg].dist(origo)<in.getflt("cellradius")*0.5) {
            spcinner1.update(con, con.p[m.beg], "OW");
            spcinner2.update(con, con.p[m.beg], "HW1");
          }
          if (con.p[m.beg+1].dist(origo)<in.getflt("cellradius")*0.5) { 
            spcinner2.update(con, con.p[m.beg+1], "OW");
            spcinner3.update(con, con.p[m.beg+1], "HW1");
            spcinner3.update(con, con.p[m.beg+1], "HW2");
          }
          if (con.p[m.beg+2].dist(origo)<in.getflt("cellradius")*0.5) 
            spcinner2.update(con, con.p[m.beg+2], "OW");
            spcinner3.update(con, con.p[m.beg+2], "HW1");
            spcinner3.update(con, con.p[m.beg+2], "HW2");
        }
      }
      // Trace system energy during simulation
      if ( cnt%eprint==0 )
         sys.track();
    }//end of micro-loop 

    spccell.write("rdf-cell-OW.dat");
    spcrdf1.write("rdf-OW-OW.dat");
    spcrdf2.write("rdf-OW-HW.dat");
    spcrdf3.write("rdf-HW-HW.dat");
    spcinner1.write("rdf-OW-OW-inner.dat");
    spcinner2.write("rdf-OW-HW-inner.dat");
    spcinner3.write("rdf-HW-HW-inner.dat");
    radorder.write("mu-OW.dat");
    dip.write("dipoledist.dat");
    dipc.write("dipolexdist.dat");
    kfuck.write("kirkwood.dat");
    sys.write();


    sys.update( 
        pot.internal(con.p, sol, sol.numatom) +
        pot.energy(con.p, sol)
        );

    gro.save("confout.gro", con.p);      // Save config. to disk
    aam.save("confout.aam", con.p);      // Save config. to disk
    cout << loop.timing()                // Show progress
         << "#   Energy (kT): sum average drift = "
         << sys.sum << " " << sys.uavg.avg() << " "
         << std::abs(sys.cur-sys.sum) << endl;
    std::cout << "#   Reaction field energy  = "<<imu.avg()<<" , stdev = "<<imu.stdev() <<endl;
    std::cout <<endl;
    //std::cout << gf.info();

    point slask;
    for (int i=0; i<con.p.size(); i++) {
      slask=con.vdist(con.p[i], check[i]);
    //  if (slask.len()<1e-14)
    //    std::cout <<"#  Particle "<<i<<" is static !!!"<<endl;
    }

  }//end of macro-loop
  cout << sys.info()
       << loop.info()  // Print final results
       << mr.info() << mt.info() <<mrt.info() <<multtr.info() <<pot.info() <<gf.info();

#ifdef SPCE_MINIMAGE
#else
//  io.writefile("img.dat", pot.printimg());
#endif
}

