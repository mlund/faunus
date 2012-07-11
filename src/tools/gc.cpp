#include <faunus/faunus.h>

using namespace Faunus;

typedef Geometry::Cuboid Tgeometry;                // select simulation geometry
typedef Potential::CoulombSR<Tgeometry, Potential::Coulomb, Potential::LennardJones> Tpairpot;

/*! \brief Class for Grand Canonical MC
 *  \author Bjorn Persson and Mikael Lund
 *  \date Lund, 2010
 *  \warning This is some general pice of shit
 *
 *  Class for handling an arbitrary salt solution of N components whereof M may 
 *  may have a fixed chemical potential rather then fixed number density.
 */
class saltbath : public Move::Movebase { 
  private:
    bool locked;              //!< Is class closed for external modification?
    Space *gcPtr;    //!< Pointer to GC class
    Group *sPtr;               //!< Pointer to Salt class (change to group)
    struct pair {
      int i, j;               //!< Indices for GC::*g  
      int ni,nj;              //!< Stochiometric ratio: ni:nj of a i-j salt
      Average < double >  ai, aj;         //!< Number averages
      Average < double >  inacc, outacc;  //!< Acceptances
    };
    vector<pair> pairs;       //!< All (electro neutral) ion pairs with reservoirs 
    pair *thispair;           //!< Choosen pair for a particular step in the Markov chain
    bool inserting;           //!< Is it a inserting step?
    particle I, J;            //!< Particles to insert
    p_vec trialin, tini, tinj;  //!< Particle vectors trailin=tini+tinj
    vector<unsigned int>      trialout, touti, toutj;  //!< Indices for attempted extraction
    void remove();
    void acc_ins();           //!< Accept insertion
    void acc_rem();           //!< Accept extraction
    void insert();
    void setChemPot(InputMap &); //!< Load and set atoms::list.chempot
    double metropolis();      //!< returns the full Boltzmann weight, taking mu's and ideal cont. in to account.
    int gcd(int, int);        //!< Find the greatest common divider of two integers (external library?)

  public:
    // Keep as much as possible private! Easier for the user of the class.
    vector<Group> g;          //!< Vector of salt components (not known externaly, maybe the should?)
    saltbath(InputMap&, Energy::Hamiltonian&, Space&, string="saltbath");
    void searchsalt(Group&);
    string _info();
};
saltbath::saltbath(InputMap &in, Energy::Hamiltonian &e, Space &s, string pfx="saltbath") : Movebase(e,s,pfx) {
  name.assign("GRAND CANONICAL SALT RESERVOIR");
  cite.assign("n/a");
  runfraction=in.getflt("saltbath_runfrac",1.0);
  deltadp=0;
  dp=0;
  g.clear();
  gcPtr=&gc;
  gcPtr->gp.clear();                    //Perhaps not always desired but for now...
  sPtr=&s;
  // Find salt components and sync with gc
  searchsalt(s);
  // Fetch chemical potentials from input file and showel parameters in to atom
  // NOTE: setChemPot depends on initialization of gc::g
  setChemPot(in);

}

// Search a salt group and sort out all species within
// and adds pointers to these and main salt group to gc::gp
void saltbath::searchsalt( salt &salt) {
  particle ref = con->p[salt.beg];
  group w;
  w.beg=salt.beg;
  w.end=salt.beg-1;
  for (unsigned int i=salt.beg; i<=salt.end; i++)
    if (con->p[i].id==ref.id)
      w.end++;
    else {
      w.name=atom[ref.id].name;
      g.push_back(w);
      w.beg=w.end=i;
      ref.id=con->p[i].id;
    }
  w.name=atom[ref.id].name;
  g.push_back(w);
  for (int i=0; i<g.size(); i++)    
    gcPtr->gp.push_back(&g[i]);
  gcPtr->gp.push_back(sPtr);
}

// Find all electro neutral pair combinations in saltbath::g
// 
void saltbath::setChemPot(inputfile &in) {
  pairs.clear();
  vector<string> cps;
  assert(2>1);
  //cps = in.getvec("MUs", "empty");  // Get sets of species, val
  if (cps.size()%2!=0) {
    std::cerr << " Not one value for each name in ChemPot input list!!"<<endl;
  }
  for (int i=0; i<cps.size() ; i+=2)
    atom[cps[i]].chempot=atof(cps[i+1].c_str());
  // NOTE Atoms::find returns zero if name is not in the list.
  // Set pair-vector using gc::g
  pair t;
  string si, sj;
  for (int i=0; i<g.size()-1; i++)
    for (int j=i+1; j<g.size(); j++) {
      si=g[i].name, sj=g[j].name;
      if (atom[g[i].name].charge*atom[g[j].name].charge < 0 and 
          atom[si].chempot < 0 and atom[sj].chempot < 0) {
        t.i=i, t.j=j;
        int tni, tnj, gCd;
        tni=int( abs(atom[g[j].name].charge) ), tnj=int( abs(atom[ g[i].name ].charge));
        gCd=gcd(tni,tnj);
        t.ni=tni/gCd, t.nj=tnj/gCd;
        pairs.push_back(t);
      }
    }
}

double saltbath::metropolis() {
  double idfactor = 1;
  double deltau   = 0;
  if (inserting==true) {
    deltau-=thispair->ni*atom[g[thispair->i].name].chempot + thispair->nj*atom[g[thispair->j].name].chempot;
    for (int i=0; i<thispair->ni; i++)
      idfactor*=(g[thispair->i].size()+1+i)/con->getvolume();
    for (int i=0; i<thispair->nj; i++)
      idfactor*=(g[thispair->j].size()+1+i)/con->getvolume();
    deltau+=log(idfactor);

  } else {
    deltau+=thispair->ni*atom[g[thispair->i].name].chempot + thispair->nj*atom[g[thispair->j].name].chempot;
    for (int i=0; i<thispair->ni; i++)
      idfactor*=(g[thispair->i].size()-thispair->ni+1+i)/con->getvolume();
    for (int i=0; i<thispair->nj; i++)
      idfactor*=(g[thispair->j].size()-thispair->ni+1+i)/con->getvolume();
    deltau-=log(idfactor);

  }
  return deltau;
}

// Find salt particles to remove
void saltbath::remove() {
  bool uniqe=true;
  int ii=g[thispair->i].beg+g[thispair->i].size()*slp.random_one();
  touti.push_back(ii);
  trialout.push_back(ii);
  while (touti.size() != thispair->ni) {
    ii=g[thispair->i].beg+g[thispair->i].size()*slp.random_one();
    for (int i=0; i<trialout.size(); i++)
      if (ii==trialout[i])
        uniqe=false;
    if (uniqe==true) {
      touti.push_back(ii);
      trialout.push_back(ii);
    }
    uniqe=true;
  }
  int jj=g[thispair->j].beg+g[thispair->j].size()*slp.random_one();
  toutj.push_back(jj);
  trialout.push_back(jj);
  while (toutj.size() != thispair->nj) {
    jj=g[thispair->j].beg+g[thispair->j].size()*slp.random_one();
    for (int i=0; i<toutj.size(); i++)
      if (jj==trialout[i])
        uniqe=false;
    if (uniqe==true) {
      toutj.push_back(jj);
      trialout.push_back(jj);
    }
    uniqe=true;
  }
}

void saltbath::acc_ins() {
  for (int i=0; i<trialin.size(); i ++)
    if (gcPtr->insert(*con, trialin[i]) == false)
      std::cerr << "GC insertion failed!!!"<<endl;
}

void saltbath::acc_rem() {
  std::sort(trialout.begin(), trialout.end() );
  //    for (int i=0; i<trialout.size(); i++) {
  //      std::cout << trialout[i]<<endl;
  //    }
  //    std::cout << "finished sorting"<<endl;
  for (int i=0; i<trialout.size(); i++) {
    if (gcPtr->erase(*con, (trialout[i]) ) == false) 
      std::cerr << "GC erase failed!!!"<<endl;
    for (int j=i+1; j<trialout.size(); j++)
      trialout[j]--;
  }
}

/*!
 * \todo Cell overlap test missing
 */
double saltbath::move() {  //Two sorted vectors are used to simplify insertions
  trialin.clear();         //while the combined is used to simplify energy evalutations
  tini.clear();
  tinj.clear();
  trialout.clear();
  touti.clear();
  toutj.clear();
  thispair=&pairs[ int(slp.random_one()*pairs.size()) ];
  du=uold=unew=0;
  cnt++;
  // 1. Randomly insert or remove.
  if (slp.random_one()<0.5) { // Insert
    inserting=true;
    I=atom(g[thispair->i].name);
    J=atom(g[thispair->j].name);
    for (int i=0; i<thispair->ni; i++) {
      con->randompos(I);
      tini.push_back(I);
      trialin.push_back(I);
    }
    for (int i=0; i<thispair->nj; i++) {
      con->randompos(J);
      tinj.push_back(J);
      trialin.push_back(J);
    }
    for (int i=0; i<trialin.size(); i++) {
      unew+=pot->energy(con->p, trialin[i]);
      for (int j=i+1; j<trialin.size(); j++)
        unew+=pot->energy(trialin[i],trialin[j]);
    }
  } else {                  // Remove
    if (g[thispair->i].size()<thispair->ni || g[thispair->j].size()<thispair->nj) {
      du=0;
      thispair->ai += double(g[thispair->i].size() );
      thispair->aj += double(g[thispair->j].size() );
      return du;
    }
    inserting=false;
    remove();
    for (int i=0; i<trialout.size(); i++) {
      uold+=pot->energy(con->p,trialout[i]);
      for (int j=i+1; j<trialout.size(); j++) {
        uold-=pot->energy(con->p[trialout[i]], con->p[trialout[j]]);
      }
    }
  }
  du = unew-uold;

  if( ens->metropolis(metropolis()+du)==true) {
    if (inserting==true) {
      acc_ins();
      thispair->inacc+=100.;
    }
    else {
      acc_rem();
      thispair->outacc+=100.;
    }
    rc=OK;
    utot+=du;
    naccept++;
    //accept move
    thispair->ai += double(g[thispair->i].size() );
    thispair->aj += double(g[thispair->j].size() );
    return du;
  } else rc=ENERGY;
  //reject move
  du=0;
  thispair->ai += double(g[thispair->i].size() );
  thispair->aj += double(g[thispair->j].size() ); 
  if (inserting==true) {
    thispair->inacc +=0.;
  } else {
    thispair->outacc+=0.; 
  }
  return du;
}

string saltbath::info() {
  std::ostringstream o;
  if (runfraction>0) {
    o << markovmove::info()
      << "#   Salt pairs -- " << pairs.size() << " combination(s):" << endl
                                                                       << "#     "
                                                                       << setw(5) << std::left << "i" << setw(5) << "j" << setw(7)
                                                                       << "ratio" << setw(9) << "ci/M" << setw(9) << "cj/M"
                                                                       << setw(6) << "ni" << setw(6) << "nj" << setw(5) << "flux"
                                                                       << setw(8) << "mui/kT" << setw(8) << "muj/kT"
                                                                       << endl;
    for (int i=0; i<pairs.size(); i++) {
      o << "#     "
        << setw(5) << g[pairs[i].i].name << setw(5) << g[pairs[i].j].name
        << setw(3) << std::right << pairs[i].ni << ":" << setw(3) << std::left << pairs[i].nj;
      if (cnt!=0) {
        double ci=pairs[i].ai.avg()/con->getvolume()*1e27/pyc.Nav;
        double cj=pairs[i].aj.avg()/con->getvolume()*1e27/pyc.Nav;
        o << setiosflags(std::ios::fixed);
        o.precision(5);
        o << setw(9) << ci << setw(9) << cj;
      }
      //double muidi=log( pairs[i].ai.avg()/con->getvolume() ),
      //       muidj=log( pairs[i].aj.avg()/con->getvolume() );
      //       muexi=atom[ g[ pairs[i].i ].name ].chempot - muidi,
      //       muexj=atom[ g[ pairs[i].j ].name ].chempot - muidj;
      o << setw(6) << g[pairs[i].i].size() << setw(6) << g[pairs[i].j].size();
      o.precision(1);
      o << setw(5) << pairs[i].inacc.avg()-pairs[i].outacc.avg();
      o.precision(3);
      o << setw(8) << atom[g[pairs[i].i].name].chempot << setw(8) << atom[g[pairs[i].j].name].chempot;
      o << endl;
    }
  }
  return o.str();
}

void saltbath::check(checkValue &test) {
  markovmove::check(test);
  test.check("SALTBATH_numpairs", pairs.size(), 1e-5);
  for (int i=0; i<pairs.size(); i++) {
    std::ostringstream o;
    o << "SALTBATH" << i+1; 
    test.check(o.str() + "_ni", pairs[i].ai.avg() );
    test.check(o.str() + "_nj", pairs[i].aj.avg() );
    double flux=pairs[i].inacc.avg() - pairs[i].outacc.avg();
    test.check(o.str() + "_flux", flux);
    test.smallerThan(o.str() + "_flux", flux, 1.);
  }
}

int saltbath::gcd(int a, int b) {
  while(1) {
    a = a%b;
    if (a==0)
      return b;
    b = b%a;
    if (b==0)
      return a;
  }
}

int main() {
  cout << textio::splash();
  atom.includefile("atomlist.inp");    // load atom properties
  InputMap mcp("gc.inp");
  MCLoop loop(mcp);                    // class for handling mc loops
  FormatPQR pqr;                           // PQR structure file I/O
  EnergyDrift sys;                     // class for tracking system energy drifts

  Energy::Hamiltonian pot;
  auto nonbonded = pot.create( Energy::Nonbonded<Tpairpot>(mcp) );
  auto bonded    = pot.create( Energy::Bonded() );
  Space spc( pot.getGeometry() );

  // Handle particles
  GroupAtomic salt(spc, mcp);
  salt.name="Salt";
  Move::AtomicTranslation mv(mcp, pot, spc);  // Particle move class
  mv.setGroup(salt);

  GroupMolecular test;
  test.beg=spc.p.size();
  spc.insert("Na", 2);
  test.last=spc.p.size()-1;
  test.name="Test";
  spc.p[test.beg].x=0;
  spc.p[test.beg].y=0;
  spc.p[test.beg].z=0;
  spc.p[test.last].x=2;
  spc.p[test.last].y=2;
  spc.p[test.last].z=2;
  spc.trial=spc.p;
  spc.enroll(test);

  spc.load("space.state");
  salt.setMassCenter(spc);
  test.setMassCenter(spc);

  //bonded->bonds.add(0,1, Potential::Harmonic(0.2, 10.0));
  //bonded->bonds.add(1,2, Potential::Harmonic(0.3,  5.0));

  // Particle titration
  Move::SwapMove tit(mcp,pot,spc);
  Move::Isobaric iso(mcp,pot,spc);
  Move::TranslateRotate gmv(mcp,pot,spc);

  // Widom particle insertion
  Analysis::Widom widom(spc, pot);
  widom.addGhost(spc);

  //FastaSequence fasta;
  //Group protein1 = fasta.insert( "AAK", spc, bonded->bonds );

#define UTOTAL \
  pot.g_external(spc.p, test)\
  + pot.g_internal(spc.p, salt)  + pot.g_external(spc.p, salt)\
  + pot.g2g(spc.p,salt,test) + pot.external()

  //#define UTOTAL pot.g_internal(spc.p, salt)  + pot.g_external(spc.p, salt) + pot.external()

  sys.init( UTOTAL );

  cout << atom.info() << spc.info() << pot.info() << mv.info()
    << textio::header("MC Simulation Begins!");

  while ( loop.macroCnt() ) {  // Markov chain 
    while ( loop.microCnt() ) {
      sys+=mv.move();
      if (slp_global.randOne()>0.9)
        widom.sample();
      sys+=tit.move();
      gmv.setGroup(test);
      sys+=gmv.move();
      sys+=iso.move();
    }
    sys.checkDrift( UTOTAL );
    cout << loop.timing();
  }

  pqr.save("confout.pqr", spc.p);
  spc.save("space.state");

  cout << mv.info() << sys.info() << loop.info() << widom.info()
    << tit.info() << gmv.info() << iso.info();
}
