/*
 * manybodyMPI.cpp
 *
 * This will simulate polymers in a rectangular slit container:
 * container. The default Hamiltonian is:
 * - Debye-Huckel electrostatics
 * - Lennard-Jones interactions
 * - Gouy-Chapman electrostatics (for slit containers)
 *
 * The following Monte Carlo Moves are implemented:
 * 1) Translation / rotation of rigid molecules
 * 2) Particle swap moves (proton titration etc.)
 * 3) Isobaric volume moves (NPT ensemble)
 * 4) Parallel Tempering using MPI
 */

#include <faunus/faunus.h>

using namespace Faunus;

#ifdef SLIT
typedef Geometry::Cuboidslit Tgeometry;
#else
//typedef Geometry::Sphere Tgeometry;
typedef Geometry::Cuboid Tgeometry;
#endif
typedef Potential::CombinedPairPotential<Potential::Harmonic, Potential::LennardJones> Tbondpot;
typedef Potential::DebyeHuckelLJ Tpairpot;
//typedef Potential::DebyeHuckelHS Tpairpot;
//typedef Potential::DebyeHuckelr12 Tpairpot;

int main(int argc, char** argv) {
  Faunus::MPI::MPIController mpi;
  mpi.cout << textio::splash();

  InputMap mcp(textio::prefix+"manybody.input");
  MCLoop loop(mcp);                    // class for handling mc loops
  FormatPQR pqr;                       // PQR structure file I/O
  FormatAAM aam;                       // AAM structure file I/O
  FormatTopology top;
  FormatXTC xtc(1000);                 // XTC gromacs trajectory format for only cuboid geomtries
  EnergyDrift sys;                     // class for tracking system energy drifts
  UnitTest test(mcp);

  Energy::Hamiltonian pot;
  auto nonbonded = pot.create( Energy::Nonbonded<Tpairpot,Tgeometry>(mcp) );
  auto bonded    = pot.create( Energy::Bonded() );
//  auto mfc       = pot.create( Energy::MeanFieldCorrection(mcp) );  //pot.create returns a smart pointer!!
  auto swpair    = pot.create( Energy::PairListID() );
  swpair->add( atom["HIS"].id, atom["Zn"].id, Potential::SquareWell(mcp) );  //add square-well poteintial between Zn+2 and HIS resiues
  swpair->add( atom["ZNHIS"].id, atom["HIS"].id, Potential::SquareWellShifted(mcp, "squarewell_ZNHIS") );  //add square-well poteintial between ZNHIS and HIS resiues
#ifdef SLIT
  auto gouy = pot.create( Energy::GouyChapman(mcp) );
  gouy->setPosition( nonbonded->geometry.len_half.z() ); // Place the surface xy plane at +z direction
  auto restricted = pot.create( Energy::RestrictedVolumeCM(mcp) );
#endif
  Space spc( pot.getGeometry() );

  // Add polymers
  vector<GroupMolecular> pol( mcp.get("polymer_N",0) );
  string polyfile = mcp.get<string>("polymer_file", "");
  atom["MM"].dp = 0.;
  int ii=1;
  mpi.cout << "Number of polymers: " << pol.size() << endl;
  for (auto &g : pol) {                              // load polymers
    aam.load(polyfile);
    Geometry::FindSpace f;
    f.find(*spc.geo, spc.p, aam.particles());        // find empty spot in particle vector
    g = spc.insert( aam.particles() );               // insert into space
    ostringstream o;
    o << "Polymer" << ii++;
    g.name=o.str();
    spc.enroll(g);
#ifdef SLIT
    //(*restricted).groups.push_back( &g );   // Restrict center of masses of all proteins
    restricted->groups.push_back( &g );       // Restrict center of masses of all proteins
#endif
    for (int i=g.front(); i<g.back(); i++)
      bonded->add(i, i+1, Tbondpot(mcp, "polymer_", "minuslj_")); // add bonds
  }
  Group allpol( pol.front().front(), pol.back().back() ); // make group w. all polymers
  
  // Add salts
  GroupAtomic salt(spc, mcp);
  salt.name="salt";
  mpi.cout << "Number of salts: " << salt.size() << endl;

  // atom["NTR"].dp = 10.;
  // atom["CTR"].dp = 10.;
  // atom["HIS"].dp = 10.;
  // atom["HNTR"].dp = 10.;
  // atom["HCTR"].dp = 10.;
  // atom["HHIS"].dp = 10.;

  //Move::Isobaric iso(mcp,pot,spc);
  Move::TranslateRotate gmv(mcp,pot,spc);
  Move::AtomicTranslation mv(mcp, pot, spc);
  Move::AtomicTranslation mv2(mcp, pot, spc, "ion");
  Move::SwapMove tit(mcp,pot,spc);
  Move::CrankShaft crank(mcp, pot, spc);
  Move::Pivot pivot(mcp, pot, spc);
  Move::ParallelTempering temper(mcp,pot,spc,mpi);
  Move::Reptation rep(mcp, pot, spc);

  Analysis::RadialDistribution<float,int> rdf(0.25);
#ifdef SLIT
  Analysis::LineDistribution<float,int> surfdist(0.25), surfmapall(0.25);
  //std::map< string , Analysis::LineDistribution<float,int> > surfmap, surfmapall;
  std::map< int , Analysis::LineDistribution<float,int> > surfmap_res;
  typedef Analysis::Table2D<double, Average<double> > Ttable;
  std::map<string , Ttable> dst_map;
#endif
  Analysis::ChargeMultipole mpol;
  Analysis::PolymerShape shape;
  spc.load(textio::prefix+"state");

  sys.init( Energy::systemEnergy(spc,pot,spc.p) );

  mpi.cout << atom.info() << spc.info() << pot.info() << tit.info()
    << textio::header("MC Simulation Begins!");

  while ( loop.macroCnt() ) {  // Markov chain 
    while ( loop.microCnt() ) {
      int k,i=rand() % 6;
      switch (i) {
        case 0: // translate and rotate molecules
          k=pol.size()+salt.size();
          while (k-->0) {
            int s=rand()%2;
            if (s==1) {
              gmv.setGroup( pol.at( rand() % pol.size() ) );
              sys+=gmv.move(); 
            }
            else {
              mv2.setGroup(salt); // translate salt particles 
              if (salt.size()>0)
                sys+=mv2.move( rand() % salt.size() ); 
            }
          }
#ifdef SLIT
          for (auto &g : pol) {
            //surfmap[g.name]( gouy->dist2surf(g.cm) )++;
            surfdist( gouy->dist2surf(g.cm) )++; // polymer mass center to GC surface histogram
	    for (int i=g.front(); i<=g.back(); i++) {
              surfmap_res[i]( gouy->dist2surf( spc.p[i] ) )++; // monomer to GC surface histogram
              surfmapall( gouy->dist2surf( spc.p.at(i) ) )++;
            }
          }
#endif
          for (auto i=pol.begin(); i!=pol.end()-1; i++)
            for (auto j=i+1; j!=pol.end(); j++)
              rdf( spc.geo->dist(i->cm,j->cm) )++; // sample polymer-polymer rdf
          break;
        case 1:
          mv.setGroup(allpol);
          sys+=mv.move( allpol.size() ); // translate monomers
          for (auto &g : pol) {
            g.setMassCenter(spc);
            shape.sample(g,spc);
          }
          break;
        case 2: // titration move
          sys+=tit.move();
          mpol.sample(pol,spc);
          break;
        case 3: // crankshaft move
          k=pol.size();
          while (k-->0) {
            crank.setGroup( pol[ rand() % pol.size() ] );
            sys+=crank.move();
          }
          for (auto &g : pol) {
            g.setMassCenter(spc);
            shape.sample(g,spc);
          }
          break;
        case 4: // pivot move
          k=pol.size();
          while (k-->0) {
            pivot.setGroup( pol[ rand() % pol.size() ] );
            sys+=pivot.move();
          }
          for (auto &g : pol) {
            g.setMassCenter(spc);
            shape.sample(g,spc);
          }
          break;
        case 5: // reptation move
          k=pol.size();
          while (k-->0) {
            rep.setGroup( pol[ rand() % pol.size() ] );
            sys+=rep.move();
          }
          for (auto &g : pol) {
            g.setMassCenter(spc);
            shape.sample(g,spc);
          }
          break;
        //case 5: // volume move
          //sys+=iso.move();
         // break;
        }

      if ( slump.runtest(0.00001) ) {
        xtc.setbox( nonbonded->geometry.len );
        xtc.save(textio::prefix+"traj.xtc", spc);
      }
      //if ( slump.runtest(0.1) )
        //(*mfc).sample(spc.p, )
#ifdef SLIT
      if ( slump.runtest(0.1) ) {
         for (auto &g : pol) {
           double d = gouy->dist2surf(g.cm);
           dst_map["Q"]( d )+=g.charge(spc.p);
           Point p=shape.vectorgyrationRadiusSquared(g, spc);
           dst_map["Rg2"]( d )+=p.x()+p.y()+p.z();
           dst_map["Rg2x"]( d )+=p.x();
           dst_map["Rg2z"]( d )+=p.z();
           dst_map["Ree2"]( d )+=spc.geo->sqdist( spc.p[g.front()], spc.p[g.back()] );
           //dst_map["<Energy>"]( d )+=sys.current();
           for (int i=g.front(); i<=g.back(); i++){
            ostringstream o;
            o << "qres" << i ;
            dst_map[o.str()]( d )+=spc.p[i].charge;
           }
         } 
      }
#endif

    } // end of micro loop

    temper.setCurrentEnergy( sys.current() );
    sys+=temper.move();

    sys.checkDrift( Energy::systemEnergy(spc,pot,spc.p) );

    mpi.cout << loop.timing();
    cout << loop.timing();

    rdf.save(textio::prefix+"rdf_p2p.dat");

#ifdef SLIT
    surfdist.save(textio::prefix+"surfdist.dat");
    surfmapall.save(textio::prefix+"surfall.dat");
    dst_map["Q"].save(textio::prefix+"netq.dat");
    //dst_map["<Energy>"].save(textio::prefix+"aveenergy.dat");

    for (auto &g : pol){   //Saving averages 
      //surfmap[g.name].save(textio::prefix+g.name+"surfdist.dat");
      std::ofstream f(textio::prefix+"res-surfdist.dat");
      std::ofstream f1(textio::prefix+"Rg-comp.dat");
      std::ofstream f2(textio::prefix+"qres.dat");
      f.precision(7);
      f1.precision(7);
      f2.precision(7);
      if (f && f1 && f2) {
        //f << "#g(z) of each residue from 1 to " << g.size() << endl;
        for (double d=0; d<=nonbonded->geometry.len.z(); d+=0.25){
          f1 << d; 
          if (dst_map["Rg2"](d).cnt > 0) f1 << "\t" << dst_map["Rg2"](d); //If counter is empty, then it prints shit lot of warning
          else f1 << "\t0";

          if (dst_map["Rg2x"](d).cnt > 0)f1 << "\t" << dst_map["Rg2x"](d); 
          else f1 << "\t0";

          if (dst_map["Rg2z"](d).cnt > 0)  f1 << "\t" << dst_map["Rg2z"](d); 
          else f1 << "\t0";

          if (dst_map["Ree2"](d).cnt > 0)  f1 << "\t" << dst_map["Ree2"](d); 
          else f1 << "\t0";
          f1 << endl;

          f << d; 
          f2 << d;
          for (int i=g.front(); i<=g.back(); i++){
            f << "\t" << surfmap_res[i](d);
            ostringstream o;
            o << "qres" << i ;
            if (dst_map[o.str()]( d ).cnt > 0 ) f2 << "\t" << dst_map[o.str()]( d );
            else f2 << "\t0"; 
          }
          f << endl;
          f2 << endl;
        }
      }
    }
#endif

    pqr.save(textio::prefix+"confout.pqr", spc.p);
    aam.save(textio::prefix+"confout.aam", spc.p);
    mcp.save(textio::prefix+"mdout.mdp");
    spc.save(textio::prefix+"state");


  } // end of macro loop

  mpi.cout << loop.info() << sys.info() << gmv.info() <<  mv.info() << mv2.info() 
    << crank.info() << pivot.info() << rep.info() << temper.info() 
    << tit.info() << mpol.info()  <<  shape.info() << spc.info();   

  tit.applycharges(spc.p);
  pqr.save(textio::prefix+"conf_aveq.pqr", spc.p);
  aam.save(textio::prefix+"conf_aveq.aam", spc.p);

}
