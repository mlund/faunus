/*! \page test_wallprotein WallproteinDebyeHuckel
 *
 * Simulate a number of flexible and titratable
 * polymers in a salt solution in the presence
 * of a charged Gouy-Chapman surface with implicit
 * salt interactions given by Debye Huckel.
 *
 * \author Chris Evers
 * \date Lund, 2011
 * \include wpDH.cpp
 */
#include "faunus/faunus.h"
#include "faunus/energy/springinteraction.h"
#include "faunus/energy/externalpotential.h"
#include "faunus/potentials/pot_r12debyehuckel.h"
#include "faunus/energy/penaltyfunction.h"
#include "faunus/moves/membrane.h"

using namespace std;
using namespace Faunus;

int main() {
  // General
  cout << faunus_splash();
  inputfile in("wp.conf");
  slump slp;                                // load random number generator
  mcloop loop(in);                          // define simulation loops
  canonical nvt;                            // canonical ensemble

  // Simulation container
#ifdef NOSLIT
  cuboid con(in);
  #ifdef NOTAB                              // wpDH_noslit_notab
    typedef pot_r12debyehuckel Tpot;
  #else                                     // wpDH_noslit
    typedef pot_r12debyehuckel_hydrophobic_tab Tpot; 
  #endif
#else
  cuboidslit con(in);
  #ifdef NOTAB
    typedef pot_r12debyehuckelXY Tpot;      // wpDH_notab
  #else
    typedef pot_r12debyehuckelXY_hydrophobic_tab Tpot;
  #endif
  #ifdef HYDROPHIC                          // wpDH_hydrophobic
    typedef expot_hydrophobic Texpot;
  #else                                     // wpDH
//     typedef expot_gouychapman Texpot;
    typedef expot_gchydrophobic Texpot;
  #endif
#endif

  // Polymer
  polymer pol;
  pol.babeladd( con, in );          //  add from input
  con.trial=con.p;                  //  synchronize particle vector
  pol.masscenter(con);              //  update masscenter
  pol.move(con, -pol.cm+(con.slice_min+con.slice_max)*0.5);  // translate to the middle of the slice or the origo (0,0,0) if no slice defined
  pol.accept(con);                  //  accept translation
  cout << pol.info();

  // Potentials
#ifdef NOSLIT
  springinteraction<Tpot> pot(in); 
#else
  #ifdef PENALTY
    springinteraction_expot_penalty<Tpot, Texpot> pot(in, con, pol); 
    pot.pen.load("penalty.xy");  // Load penaltyfunction from disk 
    pot.pen.gofrload("gofr.xy"); // Load penaltyfunction from disk 
    pot.expot.update(con);
  #else
    springinteraction_expot<Tpot, Texpot> pot(in); 
    pot.expot.update(con);
  #endif
#endif
  // Distribution functions
  double *zhalfPtr = &con.len_half.z;             // half length of container in z-direction
  double iRee = pol.calcIdealRee( in.getflt( "harmonic_req", 5 ) ); // calculate ideal end-to-end distance
  if ( iRee > *zhalfPtr )
    cout << "! Warning estimated end-to-end distance " << iRee 
      << " AA is bigger than half cuboid_zlen " << *zhalfPtr << " AA" << endl;

  histogram2 gofr(0.1, 0, *zhalfPtr*2);  // radial distribution function
  histogram fQ(1 , -pol.nb.size(), pol.nb.size());  // polymer charge distribution function
  histogram fRg(   1, 0, 2*iRee    );   // radius of gyration distribution function
  histogram fRgz2( 1, 0, iRee*iRee );   // radius of gyration distribution function in z-direction
  histogram fRee(  1, 0, 4*iRee    );   // end-to-end distance distribution function
  histogram fReez( 1, -4*iRee, 4*iRee); // end-to-end distance distribution function in z-direction
  int polEnds = in.getint("pol_ends",0);// ends of polymer (for Ree calculations)

#ifdef NOSLIT
  histogram2 fzmax(.1,0, con.len.z);             // max(z) distribution function
  histogram2 internalGofz(1, -4*iRee, 4*iRee ); // internal g(z)
  histogram2 internalGofr(1, -4*iRee, 4*iRee ); // internal g(z)
#else
  distributions2 dst(.2, 0, *zhalfPtr*2);        // distance dependent averages
#endif
// old stuff
//   vector<histogram> fInternalGofx(10, histogram(.2, -con.len.z, con.len.z ));          // radius of gyration distribution function
//   vector<histogram> fInternalGofz(10, histogram(.2, -con.len.z, con.len.z ));           // radius of gyration distribution function
//   vector<histogram> fInternalGofr(10, histogram(.2, -*zhalfPtr, *zhalfPtr ));          // radius of gyration distribution function
//   xytable< int, average<double> > rij2(1,1,pol.nb.size()); 
//   histogram fBondr2(.01,0, 20 );                // bond length distribution function

  // Moves
  monomermove mm(nvt,con,pot,in);
  translate mt(nvt,con,pot,in);
  mt.dpv.x=mt.dpv.y=0; // no need to translate in xy direction
  macrorot mr(nvt,con,pot,in);
  crankShaft cs(nvt,con,pot,in);
  branchRotation br(nvt,con,pot,in);
  eqtitrate tit(nvt, con, pot, in, "eqtit_");

  // Particle input output
  io io;
  iopqr pqr;
  ioxtc xtc(con.len.z);
  ioqtraj qtraj;
  xtc.setbox(con.len.x,con.len.y,con.len.z);
  ioaam aam;

  if ( aam.load(con, "confout.aam") ) // load stored configuration
  {
    aam.save("confin.aam", con.p);
    pqr.save("confin.pqr", con.p);
  }
  pol.masscenter(con); // update masscenter
  point cm=pol.cm;
  cm.z=0;
  pol.move(con, -cm); // translate to the middle of the xy-plane
  if (con.slicecollision(pol.cm)==true) {
    pol.move(con, -pol.cm+(con.slice_min+con.slice_max)*0.5);  // translate to the middle of the slice or the origo (0,0,0) if no slice defined
    cout << "! Slice collision" << endl;
    cout << "# Translated polymer to middle of slice" << endl;
  }
  pol.accept(con); // accept translation

  if (in.getflt("eqtit_runfrac",0.5)<1e-3) {
    pol.loadCharges(in.getstr("pol_charges","q.in"), con.p); // load residue charges from file
    con.trial=con.p;
    cout << "! Loaded polymer charges from disk" << endl;
  }
  aam.save("confin.aam", con.p);

  // Initial system energy
  systemenergy sys(
    pot.energy( con.p, pol) +
    tit.intrinsicenergy(con.p) +
    pot.uself_polymer(con.p, pol) );

  cout << con.info() << atom.info()
    << pot.info() << tit.info() 
    << in.info();

  // Simulation loop
  while ( loop.macroCnt() ) {
    while ( loop.microCnt() ) {
      switch (rand() % 6) {     // random number between 0 and 5
        case 0:
          for (int i=0; i<pol.size(); i++)
            sys+=mm.move(pol);   // move monomers 
          pol.masscenter(con);
          gofr.add(*zhalfPtr-pol.cm.z);
#ifdef PENALTY
          sys+=pot.pen.update(*zhalfPtr-pol.cm.z);
#endif
          break;
        case 1:
          sys+=mt.move(pol);     // translate polymer
          gofr.add(*zhalfPtr-pol.cm.z);
#ifdef PENALTY
          sys+=pot.pen.update(*zhalfPtr-pol.cm.z);
#endif
          break;
        case 2:
          pol.masscenter(con);
          sys+=mr.move(pol);     // rotate polymer
          break;
        case 3:
#ifdef PENALTY
          sys+=cs.penaltymove(pol);     // crankshaft
          sys+=pot.pen.update(*zhalfPtr-pol.cm.z);
#else
          sys+=cs.move(pol);     // crankshaft
#endif
          gofr.add(*zhalfPtr-pol.cm.z);
          break;
        case 4:
#ifdef PENALTY
          sys+=br.penaltymove(pol);     // branch rotation
          sys+=pot.pen.update(*zhalfPtr-pol.cm.z);
#else
          sys+=br.move(pol);     // branch rotation
#endif
          gofr.add(*zhalfPtr-pol.cm.z);
          break;
        case 5:
          sys+=tit.move();       // titrate titratable sites
          pol.charge(con.p);
          break;
      }

      if (slp.random_one()>0.3) 
      {
        pol.masscenter(con);
        double z=*zhalfPtr-pol.cm.z;  // distance between masscenter and wall

        double charge=pol.charge(con.p);
        fQ.add(charge);

        point rg2Point = pol.sqmassgradius3D(con);
        double rg2 = rg2Point.len();
        double rg = sqrt(rg2);
        pol.rg2 += rg2;
        pol.rg += rg;
        fRg.add(rg);
        fRgz2.add(rg2Point.z);

        point reePoint = con.p[pol.beg+polEnds] - con.p[pol.end-polEnds];
        con.boundary(reePoint);
        double ree = reePoint.len();
        pol.ree2 += ree*ree;
        pol.ree += ree;
        fRee.add( ree );
        fReez.add( reePoint.z );

#ifdef NOSLIT // wpDH_noslit
        point rmax;
        for ( int i=pol.beg; i<=pol.end; i++ ) 
        {
          point t = con.p[i]-pol.cm; // vector to center of mass
          con.boundary(t); // periodic boundary (if any)
          if ( t.x > rmax.x )
            rmax.x=t.x;
          if ( t.y > rmax.y )
            rmax.y=t.y;
          if ( t.z > rmax.z )
            rmax.z=t.z;
          internalGofz.add( t.x );
          internalGofz.add( t.y );
          internalGofz.add( t.z );
          internalGofr.add( t.len() );
        }
        fzmax.add( rmax.x );
        fzmax.add( rmax.y );
        fzmax.add( rmax.z );

#else // wpDH
        dst.add("Q", z, charge );

        dst.add("Rg2",  z, rg2 );
        dst.add("Rg",   z, rg );
        dst.add("Rg2x", z, rg2Point.x );
        dst.add("Rg2y", z, rg2Point.y );
        dst.add("Rg2z", z, rg2Point.z );

        dst.add("Ree2", z, ree*ree );
        dst.add("Ree",  z, ree );
        dst.add("Ree2x", z, reePoint.x);
        dst.add("Ree2y", z, reePoint.y);
        dst.add("Ree2z", z, reePoint.z);

        for ( int i=pol.beg; i<=pol.end; i++ ) 
        {
          std::ostringstream s; s << "q" << i;
          dst.add( s.str(), *zhalfPtr-con.p[i].z, con.p[i].charge );
        }
#endif

// old stuff
//         // internal Gofr functions
//         size_t zint = round(z);
//         size_t zndx = (zint-2)/2; // what happens with negative vaules?
//         if ( zndx < fInternalGofr.size() && zint%2 == 0 && abs(z - zint) < .2 )
//         {
//           for (int i=pol.beg; i<pol.end; i++) 
//           {
//             point t = con.p[i] - pol.cm; // vector to center of mass
//             con.boundary(t); // periodic boundary (if any)
//             fInternalGofx[zndx].add( t.x );
//             fInternalGofz[zndx].add( t.z );
//             fInternalGofr[zndx].add( t.len() );
//           }
//         }

//         for (int i=pol.beg; i<=pol.end-1; i++) 
//         {
//           fBondr2.add(sqrt(con.sqdist(con.p[i],con.p[i+1]))  );
//           for (int j=i+1; j<=pol.end; j++) 
//           {
//             rij2(j-i) += con.sqdist( con.p[i], con.p[j] );
//           }
// 
//           point l = con.p[i] - con.p[i+1];
//           con.boundary(l);
//           lz.add( l.x * l.x );
//           lz.add( l.y * l.y );
//           lz.add( l.z * l.z );
//         }

      } // end of analysis

      if (slp.random_one()<in.getflt("traj_runfrac",0.05) ) 
      {
        xtc.save( "traj.xtc", con.p );
        qtraj.save( "q.traj", con.p );
      }
    } // END of micro loop

    sys.update(
      pot.energy( con.p, pol) +
      tit.intrinsicenergy(con.p) +
      pot.uself_polymer(con.p, pol) );

    int precision = cout.precision();
    cout << loop.timing() 
      << "#   Energy drift = " << sys.cur-sys.sum << " kT "
      << "(current energy = " << setprecision(3) << sys.cur << setprecision(precision) << " kT;"
      << " charge = " << con.charge() << ") " << endl;

    // Write files to disk
    aam.save("confout.aam", con.p);
    if ( in.getboo("write_files",true) ) 
    {
      io.writefile("vmdbonds.tcl", pol.getVMDBondScript());
      pqr.save("confout.pqr",con.p);

      gofr.write("gofr.out");
      gofr.dump("gofr.xy");
      fQ.dump("fluctQ.xy");
      fRg.dump("fluctRg.xy");
      fRgz2.dump("fluctRgz2.xy");
      fRee.dump("fluctRee.xy");
      fReez.dump("fluctReez.xy");

      tit.applycharges(con.trial);         // Set average charges on all titratable sites
      pol.saveCharges("q.out", con.trial); // Save average charges to disk
      con.trial=con.p;                     // Restore original charges

#ifdef NOSLIT
      fzmax.dump("fluctzmax.xy");
      internalGofr.dump("internalGofr.xy");
      internalGofz.dump("internalGofz.xy");
#else
      dst.write("dist.out");
      dst.cntwrite("cntdist.out");
      pot.expot.dump("expot.xy");
#endif
#ifdef PENALTY
      pot.pen.dump("penalty", loop.cnt_macro, "xy");
      pot.pen.gofrdump("gofr", loop.cnt_macro, "xy");
#endif
#ifndef NOTAB
      pot.pair.dump("pairpot.xy");
#endif

// old stuff
//       for ( size_t zndx=0; zndx<fInternalGofr.size(); zndx++)
//       {
//         int zint = 2+zndx*2;
//         ostringstream filename;
//         filename << "internalGofr" << zint << ".xy";
//         fInternalGofr[ zndx ].dump( filename.str() );
//         filename.str("");
//         filename << "internalGofx" << zint << ".xy";
//         fInternalGofx[ zndx ].dump( filename.str() );
//         filename.str("");
//         filename << "internalGofz" << zint << ".xy";
//         fInternalGofz[ zndx ].dump( filename.str() );
//       }
//       fBondr2.dump("bondr2.xy");
//       rij2.dumptodisk("rij2.xy");

    } // end of write_files
  } // END of macro loop and simulation

cout << sys.info() << loop.info() << mm.info()
  << cs.info() << br.info() << mr.info() 
  << mt.info() << tit.info(con.p) << pol.info(con) 
  << pot.info();
}
