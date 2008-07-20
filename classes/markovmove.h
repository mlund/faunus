#ifndef _markovmove_h
#define _markovmove_h

#include <cstdarg>
#include "container.h"
#include "potentials.h"
#include "ensemble.h"
#include "titrate.h"
#include "slump.h"
#include "io.h"
#include "histogram.h"
//typedef pot_coulomb T_pairpot;

/*! \brief Base class for MC moves
 *  \author Mikael Lund
 *  \todo Perhaps the T_pairpot could be made more elegant?
 *
 *  This class will keep track of the returned energy, if the
 *  move was successful of not and can provide statistics about the
 *  accepted and rejected moves.
 *
 *  Each derived class should provide a simple move() function that will
 *  perform a trial move, calculate the energy and either accept or reject
 *  the move. Unsuccessful moves are automatically undone.
 *
 *  The used pair-potential is identified using the
 *  type \verb T_pairpot \endverbatim that must be
 *  defined before processing the source code. For example,
 *
 *  #include "potentials.h"
 *  typedef pot_coulomb T_pairpot
 *  #include "markovmove.C"
 *  \endcode
 */
class markovmove {
  protected:
    slump slp;
    float runfraction;
    double uold, unew, deltadp;
    unsigned long long int cnt, naccept;
    string name;                        //!< Arbitrary name for the move
    string cite;                        //!< Reference to additional info. (article, url etc.)
    container *con;
    ensemble *ens;
    interaction<T_pairpot> *pot;
    average<float> dpsqr;               //!< Average displacement squared
  public:
    enum keys {OK, ENERGY, HC};
    keys rc;                            //!< Return code from move() functions
    double dp,                          //!< Displacement parameter
           du,                          //!< Energy change of last move
           utot;                        //!< Sum of energy changes for all moves
    float accepted();                   //!< Return fraction of accepted moves
    bool run(float);                    //!< Probability
    void adjust_dp(float=30, float=40); //!< Adjust displacement parameter
    virtual string info();              //!< Show info about group 
    markovmove(ensemble &e, container &c, interaction<T_pairpot> &inter) {
      du=utot=dp=deltadp=0;
      cnt=naccept=0;
      ens=&e;
      con=&c;
      pot=&inter;
      runfraction=1;
    }
};
string markovmove::info() {
  ostringstream o;
  o << endl
    << "# " << name << ":" << endl;
  if (cnt>0) {
    o << "#   Acceptance          = " << accepted()*100 << endl
      << "#   Number of trials    = " << cnt << endl
      << "#   Pct. of Markov steps= " << runfraction*100 << endl
      << "#   Total energy change = " << utot << endl;
    if (dp!=0) {
      o << "#   Displacement param. = " << dp << endl
        << "#   Average displacement= " << sqrt(dpsqr.avg()) << endl;
    }
  }
  if (cite.empty()==false)
    o << "#   More information:     " << cite << endl;
  return o.str();
}

/*! \brief Move salt particles
 *  \author Mikael Lund
 */
class saltmove : public markovmove {
  public:
    saltmove( ensemble &, container&, interaction<T_pairpot>& );
    double move(group &, int);  //!< Move a single particle
    double move(group &);         //!< Loop over group particles (randomly)
};

/*! \brief Random move of a macromolecule in a container 
 *  \todo generalize, it is cubix periodic boundry specific
 *  \author Bjoern Persson
 */
class move : public markovmove {
  private:
    box *b;
  public:
    move( ensemble&, container&, interaction<T_pairpot>&, box &, float);
    bool mOve(macromolecule &);
};
/*! \brief Move one group parallel z-axis
 *  \author Bjoern Persson
 */
class zmove : public markovmove {
  public:
    zmove( ensemble&, container&, interaction<T_pairpot>&, macromolecule&, float);
    bool move(macromolecule &);
    float z;
};

/*!
 * This move will symmetrically translate two macromolecules
 * along the line connecting their mass centers.
 *
 * \brief Symmetrically move two groups along z-axis
 * \author Mikael Lund
 */
class dualmove : public markovmove {
  private:
    point v;
  public:
    histogram gofr;             //!< g(r) of the two group mass centers
    double r;                   //!< Current distance between group mass centers
    double rmax;                //!< Maximum allowed mass-center distance
    double rmin;                //!< Minimum allowed mass-center distance
    dualmove( ensemble&, container&, interaction<T_pairpot>&);
    void setup(inputfile &);
    void load(inputfile &, vector<macromolecule> &g, float=0);
    void direction(double, double, double);
    double move(macromolecule &, macromolecule &);
    string info();
};

class translate : public markovmove {
  public: 
    translate( ensemble&, container&, interaction<T_pairpot>&);
    double move(macromolecule &); 
};

/*!
 * \brief Fluctuate the volume against an external pressure 
 * \author Bjoern Persson and Mikael Lund
 * \todo Should take an arbritray vector of groups and implement a more efficent energycalculation.
 * Salt scaling not yet implemented
 * \note Tested and found valid using debyehuckel_P3 and point charges.
 * \date Lund/Canberra, Jan-Feb 2008
 *
 * This move is intended to maintain constant pressure in a canonical
 * Monte Carlo simulation. It will perform volume fluctuations and scale
 * particle positions according to the scaling algorithm specified for the container.
 */
class isobaric : public markovmove {
  public:
    isobaric( ensemble&, container&, interaction<T_pairpot>&, double, double, int);
    string info();
    void move(group &, unsigned short);  //!< Pressure scale group
    void move(unsigned short, group &,...);
    double move(vector<macromolecule> &);
    average<float> vol, vol2, len, len2;
    double pen;
    vector<double> penalty;
    void loadpenaltyfunction(string);
    void printpenalty(string, histogram &);
    void printpenalty(string);
  private:
    unsigned short N;
    interaction<T_pairpot> trialpot; // Copy of potential class for volume changes
    double P, dV, dh; // Pressure, volume difference and hamiltonian difference
    double newV;      // New volume
    void newvolume();
    void accept();
    vector<double> newpenalty;
};

/*! \brief Rotate group around its mass-center.
 *  \author Mikael Lund
 *  \date Prague 2007
*/
class macrorot : public markovmove { 
  public:
    macrorot( ensemble&, container&, interaction<T_pairpot>&);
    double move(macromolecule &);
};

/*! \brief Titrate all titrateable sites
 *  \author Mikael Lund
 *  \note "titrate" used to be private. Changed because iopqr::save()
 */
class chargereg : public markovmove, public titrate {
  public:
    chargereg( ensemble&, container&, interaction<T_pairpot>&, group&, float);
    double titrateall();
    string info();
};

/*!
 * \brief Grand Canonical titation of all sites
 * \author Bjoern Persson
 * \todo Untested, in principle it must be supplemented with grand canonical salt
 */
class HAchargereg : public chargereg {
  public: 
    HAchargereg( ensemble&, container&, interaction<T_pairpot>&, group&, float, float);
    string info();
  private:
    double energy(vector<particle> &, double, action &); //!< New titrate energy function
    double CatPot;  //!< Chemical potential of coupled cation
};

/*!
 * This type of move will attempt to move collective sets of macromolecules that
 * obeys some criteria (here a hardcore overlap) with a symmetric transition
 * matrix (no flow through the clusters).
 * \author Bjoern Persson
 * \brief Quick and dirty translational cluster move
 * \warning This rutine is only compatible for systems containing a set of macromolecules!
 *
 */
class clustertranslate : public markovmove {
  public: 
    clustertranslate( ensemble&, container&, interaction<T_pairpot>&, double);
    double move(vector<macromolecule> &); 
    vector<int> cluster;    //! Index vector for macromolecule number in cluster
    vector<int> free;       //! Index vector for 'free' macromolecules
    vector<int> nacc;
    vector<int> ntrial ;
    int np;                 //! Number of particles in config.
    bool decreasing(int, int);
    void sort(vector<macromolecule> &);  //function to sort macromolecules in to two classes, cluster and free
    void flowcheck(vector<macromolecule> &); //ensure detailed balance
    hardsphere coll;
    double sep;             // separation parameter
    int FLOW;               // control variable
    string info();
};

/*!
 * This type of move will attempt to move collective sets of macromolecules that
 * obeys some criteria (here a hardcore overlap) with a symmetric transition
 * matrix (no flow through the clusters).
 * \author Bjoern Persson
 * \brief Quick and dirty rotational cluster move
 * \warning This rutine is only compatible for systems containing a set of macromolecules!
 *
 */
class clusterrotate : public markovmove {
  public: 
    clusterrotate( ensemble&, container&, interaction<T_pairpot>&);
    double move(vector<macromolecule> &); 
    vector<int> cluster;    //! Index vector for macromolecule number in cluster
    vector<int> free;       //! Index vector for 'free' macromolecules
    vector<int> nacc;
    vector<int> ntrial ;
    int np;                 //! Number of particles in config.
    bool decreasing(int, int);
    void sort(vector<macromolecule> &);  //function to sort macromolecules in to two classes, cluster and free
    void flowcheck(vector<macromolecule> &); //ensure detailed balance
    hardsphere coll;
    double sep, angle;             // separation parameter
    int FLOW;               // control variable
    string info();
};


class DHchargereg : public chargereg {

};
#endif

//--------------- MARKOV MOVE ---------------------
bool markovmove::run(float p) { return (p>slp.random_one())?true:false; }
float markovmove::accepted() { return naccept/float(cnt); }

/*!
 * This function will adjust the displacement parameter in a way
 * that the acceptance ration lies within a certain tange. Useful
 * for equilibration runs -- do not use it in production runs!
 * \param max Maximum percentage of accepted moves
 * \warning This violates the detailed balance criteria!
 * \param min Minimum percentage of accepted moves
 * \author Mikael Lund
 * \todo Specify a maxmimum dp
 */
void markovmove::adjust_dp(float min, float max) {
  float a=accepted()*100.;
  if (a>max) dp+=deltadp;
  if (a<min) dp-=deltadp;
  if (dp<=0) dp=deltadp;
}

//-------------- SALT MOVE ---------------------------------
saltmove::saltmove(
    ensemble &e, container &c, interaction<T_pairpot> &i ) : markovmove(e,c,i) {
  dp=30;
  deltadp=2;
  name="SALT DISPLACEMENTS";
  runfraction=0.5;
}

/*! \param group Group containing mobile ions
 */
double saltmove::move(group &g) {
  du=0;
  if (slp.runtest(runfraction)==false)
    return du;
  double sum=0;
  for (unsigned short i=0; i<g.size(); i++) {
    move(g, 1);
    sum+=du;
  }
  du=sum;
  return du;
}
double saltmove::move(group &g, int n) {
  du=0;
  cnt++;
  n=g.displace(*con, dp); 
  if (con->collision( con->trial[n] )==true)
    rc=HC;
  else {
    //#pragma omp parallel
    {
      //#pragma omp sections
      {
        //#pragma omp section
        { uold = pot->energy(con->p, n);   }
        //#pragma omp section
        { unew = pot->energy(con->trial,n);   }
      }
    }
    du = unew - uold;
    if (ens->metropolis(du)==true) {
      rc=OK;
      utot+=du;                               // track energy changes
      dpsqr+=con->sqdist(con->p[n],con->trial[n]); // track avg. displacement
      naccept++;                              // accept counter
      con->p[n] = con->trial[n];              // Accept move
      return du;
    } else rc=ENERGY;
  }
  du=0;
  dpsqr+=0;
  con->trial[n] = con->p[n];
  return du;
}

//---------- TRANSLATE GROUP ----------------
translate::translate( ensemble &e,
    container &c, interaction<T_pairpot> &i ) : markovmove(e,c,i) {
  name = "MACROMOLECULAR TRANSLATION";
  runfraction=1.0;
  deltadp=1.;
  dp=10.;
};
double translate::move(macromolecule &g) {
  du=0;
  cnt++;
  point p;
  p.x=dp*slp.random_half();
  p.y=dp*slp.random_half();
  p.z=dp*slp.random_half();
  g.move(*con, p); 
  if (con->collision(g.cm_trial)==true) {
    rc=HC;
    dpsqr+=0;
    g.undo(*con);
    return du;
  }
  //#pragma omp parallel
  {
    //#pragma omp sections
    {
      //#pragma omp section
      { uold = pot->energy(con->p, g);   }
      //#pragma omp section
      { unew = pot->energy(con->trial, g);   }
    }
  }
  du = unew-uold;
  if (ens->metropolis(du)==true) {
    rc=OK;
    utot+=du;
    dpsqr+=con->sqdist( g.cm, g.cm_trial );
    naccept++;
    g.accept(*con);
    return du;
  } else rc=ENERGY;
  du=0;
  dpsqr+=0;
  g.undo(*con);
  return du;
}

//---------- ROTATE GROUP AROUND CM ---------
macrorot::macrorot( ensemble &e,
    container &c, interaction<T_pairpot> &i ) : markovmove(e,c,i)
{
  name = "MACROMOLECULAR ROTATION";
  runfraction=1.0;
  deltadp=0.1;
  dp=1.0;
};

/*!
 * \todo Cell overlap test missing
 */
double macrorot::move(macromolecule &g) {
  du=0;
  cnt++;
  g.rotate(*con, dp); 
  //insert cell overlap test
  //#pragma omp parallel
  {
    //#pragma omp sections
    {
      //#pragma omp section
      { uold = pot->energy(con->p, g);   }
      //#pragma omp section
      { unew = pot->energy(con->trial, g);   }
    }
  }
  du = unew-uold;
  if (ens->metropolis(du)==true) {
    rc=OK;
    utot+=du;
    naccept++;
    g.accept(*con);
    return du;
  } else rc=ENERGY;
  du=0;
  g.undo(*con);
  return du;
}
//---------- CHARGE REG ---------------------
string chargereg::info() {
  ostringstream o;
  o <<  markovmove::info()
    << "#   pH (concentration)  = " << ph << endl
    << "#   Titrateable sites   = " << sites.size() << endl
    << "#   Number of protons   = " << protons.size() << endl;
  return o.str();
}
chargereg::chargereg(ensemble &e,
    container &c,
    interaction<T_pairpot> &i,
    group &g,
    float ph ) : markovmove(e,c,i), titrate(c,c.p,g,ph)
{
  name="PROTON TITRATION";
  cite="Biochem. 2005, 44, 5722-5727.";
  runfraction=0.2;
  con->trial = con->p;
}

/*! \brief Exchange protons between bulk and titrateable sites.
 *
 *  This move will randomly go through the titrateable sites and
 *  try to exchange protons with the bulk. The trial energy is:
 */
double chargereg::titrateall() {
  du=0;
  if (slp.runtest(runfraction)==false)
    return du;
  action t;
  double sum=0;
  for (unsigned short i=0; i<sites.size(); i++) {
    cnt++;
    t=exchange(con->trial);
    //#pragma omp parallel
    {
      //#pragma omp sections
      {
        //#pragma omp section
        { uold = pot->potential( con->p, t.site ) * con->p[t.site].charge
          + pot->potential( con->p, t.proton ) * con->p[t.proton].charge
            - con->p[t.site].potential(con->p[t.proton] )*con->p[t.proton].charge;
        }
        //#pragma omp section
        { unew = pot->potential(con->trial,t.site)*con->trial[t.site].charge 
          + pot->potential(con->trial,t.proton)*con->trial[t.proton].charge
            - con->trial[t.site].potential(con->trial[t.proton] )*con->trial[t.proton].charge;
        }
      }
    }
    du = (unew-uold)*pot->pair.f;

    //uold = pot->energy(con->p, t.site)
    //  +    pot->energy(con->p, t.proton)
    //  -    pot->pair.pairpot(con->p[t.site], con->p[t.proton] ) * pot->pair.f;
    //unew = pot->energy(con->trial, t.site)
    //  +    pot->energy(con->trial, t.proton)
    //  -    pot->pair.pairpot(con->trial[t.site], con->trial[t.proton] ) * pot->pair.f;
    //du = (unew-uold);

    if (ens->metropolis( energy(con->trial, du, t) )==true) {
      rc=OK;
      utot+=du;
      naccept++;
      con->p[t.site].charge   = con->trial[t.site].charge;
      con->p[t.proton].charge = con->trial[t.proton].charge;
      sum+=du;
    } else {
      rc=ENERGY;
      exchange(con->trial, t);
    }
    samplesites(con->p);  // Average charges on all sites
  }
  return sum;
}

//-----------MOVE----------------------------------------
/*! \brief Class to prefom a random walk of a macromolecule in space
 *  \note Replaced by translate(?)
 */
move::move(
   ensemble &e, container &c, interaction<T_pairpot> &i 
   , box &BOX,float rf) : markovmove(e, c, i) {
  runfraction=rf;
  dp=100;
  deltadp=1;
  name="MACROMOLECULE MOVE";
  b=&BOX;
}

bool move::mOve(macromolecule &g) {
  du=0;
  if (slp.runtest(runfraction)==false)
    return false;
  cnt++;
  point rand;
  rand.x = dp*slp.random_half();
  rand.y = dp*slp.random_half();
  rand.z = dp*slp.random_half();
  g.move(*con, rand);
  for (int i=g.beg; i<(g.size()+g.beg); i++) { 
    b->boundary(con->trial[i]);
    if (con->collision( con->trial[i] )==true) 
      rc=HC;
    }
  if (rc==HC) {
    for (int i=g.beg; i<(g.size()+g.beg); i++) 
      con->trial[i] = con->p[i];
    cout << "rejected"<<endl;
    return false; }
  else {
    //#pragma omp parallel
    {
      //#pragma omp sections
      {
        //#pragma omp section
        { uold = pot->energy(con->p, g);   }
        //#pragma omp section
        { unew = pot->energy(con->trial,g);   }
      }
    }
    du += unew - uold; }
    if (ens->metropolis(du)==true) {
      rc=OK;
      utot+=du;
      naccept++;
      for (int i=g.beg; i<(g.beg+g.size()); i++)
        con->p[i] = con->trial[i];
      return true;
    } else rc=ENERGY;
    
  du=0;
  for (int i=g.beg; i<(g.size()+g.beg); i++) 
    con->trial[i] = con->p[i];
  return false;
}

//-----------Z-MOVE--------------------------------------
/*! \brief Class to make moves of a single particle along z-axis
 *  
 *  This class will make a random walk along the z-axis for a
 *  macromolece
 *
 *  \todo Needs some testing and perhaps some optimization
 *  \note Replaced by dualmove(?)
 */

zmove::zmove(
    ensemble &e, container &c, interaction<T_pairpot> &i, macromolecule &g 
    ,float rf) : markovmove(e,c,i) {
  runfraction=rf;
  dp=8;
  deltadp=1;
  name="MACROMOLECULE Z-DISPLACEMENTS";
}     

bool zmove::move(macromolecule &g) {
  du=0;
  if (slp.runtest(runfraction)==false)
    return false;
  cnt++;
  z=2*dp*slp.random_half();
  g.zmove(*con, z);
  for (int i=g.beg; i<(g.size()+g.beg); i++) { 
    if (con->collision( con->trial[i] )==true) 
      rc=HC;
    }
  if (rc==HC) {
    for (int i=g.beg; i<(g.size()+g.beg); i++) 
      con->trial[i] = con->p[i];
    return false; }
  else {
    //#pragma omp parallel
    {
      //#pragma omp sections
      {
        //#pragma omp section
        { uold = pot->energy(con->p, g);   }
        //#pragma omp section
        { unew = pot->energy(con->trial,g);   }
      }
    }
    du += unew - uold; }
    if (ens->metropolis(du)==true) {
      rc=OK;
      utot+=du;
      naccept++;
      for (int i=g.beg; i<(g.beg+g.size()); i++)
        con->p[i] = con->trial[i];
      return true;
    } else rc=ENERGY;
    
  du=0;
  for (int i=g.beg; i<(g.size()+g.beg); i++) 
    con->trial[i] = con->p[i];
  return false;
}

//------- DUAL MOVE --------------
dualmove::dualmove( ensemble &e,
    container &c, interaction<T_pairpot> &i ) : markovmove(e,c,i), gofr(0.1,0.,100.)
{
  name = "SYMMETRIC 1D GROUP TRANSLATION";
  cite = "Biophys J. 2003, 85, 2940";
  runfraction=1.0;
  deltadp=1.;
  dp=1.;
  direction(0,0,1); // move only in z direction
  rmin=0;
  rmax=pow(double(c.getvolume()), 1./3)/2.; // rough estimate from volume
}

/*! Load the following parameters from an input
 * object: displacement (dm_dp), minimum separation
 * (dm_minsep), maximum separation (dm_maxsep)
 */
void dualmove::setup( inputfile &in ) {
  dp = in.getflt("dm_dp",dp);
  rmin = in.getflt("dm_minsep",rmin);
  rmax = in.getflt("dm_maxsep",rmax);
}

/*! Specify unit vector that determines which coordinates
 *  that will be moved. The default is (0,0,1) meaning that
 *  the macromolecules will be moved in the z direction only.
 */
void dualmove::direction(double x, double y, double z) {
  v.x=x;
  v.y=y;
  v.z=z;
}

/*!
 * Loads two macromolecules from disk and place them symmetrically
 * around the cell origin, separated by a specified distance.
 *
 * \param in Inputfile object that contains the protein filenames
 * \param g Macromolecular vector (will be erased!)
 * \param dist Initial distance between the protein mass-centers
 * \note The macromolecule vector is erased before any proteins are loaded!
 */
void dualmove::load(inputfile &in, vector<macromolecule> &g, float dist) {
  if (dist==0)
    dist=rmax;
  g.clear();
  point a = v*(dist/2.);
  ioaam aam(*con);
  aam.load(*con, in, g); 
  g[0].move(*con, -( g[0].cm + a ));    // ...around the cell origo
  g[1].move(*con, -( g[1].cm - a ));    // ...along the z-axis
  g[0].accept(*con);
  g[1].accept(*con);
}

string dualmove::info() {
  ostringstream o;
  o <<  markovmove::info()
    << "#   Min/max separation  = " << rmin << " " << rmax << endl;
  return o.str();
} 

double dualmove::move(macromolecule &g1, macromolecule &g2) {
  du=0;
  if (dp==0)
    return du;
  cnt++;
  r=con->dist(g1.cm, g2.cm);
  point p;
  group g12=g1+g2;
  p.x=v.x*dp*slp.random_half();
  p.y=v.y*dp*slp.random_half();
  p.z=v.z*dp*slp.random_half();
  g1.move(*con, p);
  g2.move(*con,-p);
  double rtrial = con->dist(g1.cm_trial,g2.cm_trial);
  if (con->collision(g1.cm_trial)==true ||
      con->collision(g2.cm_trial)==true || 
      rtrial > rmax || rtrial < rmin) {
    rc=HC;
    dpsqr+=0;
    g1.undo(*con);
    g2.undo(*con);
    gofr.add(r);
    return du;
  }
  //#pragma omp parallel
  {
    //#pragma omp sections
    {
      //#pragma omp section
      { uold = pot->energy(con->p,g12) + pot->energy(con->p,g1,g2);   }
      //#pragma omp section
      { unew = pot->energy(con->trial,g12) + pot->energy(con->trial,g1,g2);   }
    }
  }
  du = unew-uold;
  if (ens->metropolis(du)==true) {
    rc=OK;
    utot+=du;
    r=con->dist(g1.cm, g2.cm);
    gofr.add(r);
    dpsqr+=4.*con->sqdist(g1.cm,g1.cm_trial);
    naccept++;
    g1.accept(*con);
    g2.accept(*con);
    return du;
  } else
    rc=ENERGY;
  du=0;
  dpsqr+=0;
  g1.undo(*con);
  g2.undo(*con);
  gofr.add(r);
  return du;
}

//---------- INHERITED GCTITRATE ---------
/*!
 * \param ph pH value
 * \param mu Proton excess chemical potential
 */
HAchargereg::HAchargereg(ensemble &e,
    container &c,
    interaction<T_pairpot> &i,
    group &g, float ph, float mu ) : chargereg(e,c,i,g,ph)
{
  name="GC PROTON TITRATION...";
  cite="Labbez+Jonsson....";
  CatPot=mu;
}

string HAchargereg::info() {
  ostringstream o;
  o << chargereg::info();
  o << "#   Excess chem. pot.   = " << CatPot << endl;
  return o.str();
}

double HAchargereg::energy( vector<particle> &p, double du, titrate::action &a ) {
  int i=p[a.site].id;
  if (a.action==PROTONATED)
    return du+( log(10.)*( ph - con->d[i].pka ) )+CatPot-log(protons.size()/con->getvolume()) ;
  else
    return du-( log(10.)*( ph - con->d[i].pka ) )-CatPot+log((protons.size()+1)/con->getvolume());
}

//---------- ISOBARIC ----------
/*!
 * \note If your pair-potential depends on the geometry of the 
 *       container (minimum image, for example) make sure that
 *       the pair-potential has a correct implementation of the
 *       setvolume() function. To allow for parallization this 
 *       class will instantiate a new (private) interaction class
 *       and in this way allow for two different container sizes.
 *
 * \param pressure bulk pressure in A^-3
 * \param i Potential class. Make sure that the T_pairpot::setvolume()
 *        function is appropriate.
 */
isobaric::isobaric (
  ensemble &e,
  container &c,
  interaction<T_pairpot> &i, double pressure, double PEN,int maxsize ) : markovmove(e,c,i), trialpot(i) {
    name="ISOBARIC";
    cite="none so far";
    P=pressure;
    runfraction=0.10;
    dp=100; 
    N=0;
    penalty.clear();
    penalty.resize(maxsize,0);
    pen=PEN;
    newpenalty.clear();
    newpenalty.resize(maxsize,0);
}
void isobaric::printpenalty(string file, histogram &ld) {
  ofstream f(file.c_str());
  if (f) {
    int j=penalty.size();
    f.precision(12);
    for (int i=0; i<j; i++)
        f << i <<" "<< penalty[i]+newpenalty[i]/cnt+log(ld.get(i))<<endl;
    f.close();
  }
}
void isobaric::printpenalty(string file) {
  ofstream f(file.c_str());
  if (f) {
    int j=penalty.size();
    f.precision(12);
    for (int i=0; i<j; i++)
        f << i <<" "<< penalty[i]<<endl;
    f.close();
  }
}
void isobaric::loadpenaltyfunction(string file) {
  vector<string> v;
  string s;
  v.clear();
  ifstream f(file.c_str() );
  if (f) {
    while (getline(f,s))
      v.push_back(s);
    f.close();
    cout << "# Penalty function loaded!! Let's go biased" << endl;
  }
  else cout << "# WARNING! Penalty function " << file << " NOT READ!\n";
  short c=v.size();
  double val;
  int bin;
  for (short i=0;i<c;i++) {
    val=0;
    bin=0;
    stringstream o;
    o << v[i];
    o >> bin >> val;
    penalty[bin]=val;
  }
}
string isobaric::info() {
  ostringstream o;
  o << markovmove::info();
  o << "#   External pressure   = "<< P <<"(A^-3) = "<<P*1660<<" (M)"<< endl
    << "#   Average volume      = "<< vol.avg() << " A^3" << endl
    << "#     sqrt(<V^2>-<V>^2) = "<< sqrt(vol2.avg()-pow(vol.avg(),2)) << endl
    << "#   Average box length  = "<< len.avg() << " A^3" << endl
    << "#     sqrt(<l^2>-<l>^2) = "<< sqrt(len2.avg()-pow(len.avg(),2)) << endl;
   return o.str();
}
void isobaric::newvolume() {
  newV = exp(log(con->getvolume())  // Random move in ln(V)
      + slp.random_half()*dp);
  trialpot.pair.setvolume(newV);    // Set volume for trial-pair-potential
  dV=newV - con->getvolume();
}
void isobaric::accept() {
  naccept++;
  rc=OK;
  utot+=du;
  dpsqr+=pow(newV,1./3.);
  con->setvolume(newV);         // save the trial volume in the container
  pot->pair.setvolume(newV);    // ...AND in the pair-potential function
  vol += con->getvolume();
  vol2+= pow(con->getvolume(),2);
  len += pow(con->getvolume(),1./3);
  len2+= pow(con->getvolume(),2./3);
}
/*!\param n Number of groups
 * \param g List of groups
 * \note Unfinished
 */
void isobaric::move(unsigned short n, group &g,...) {
  /*
  group first=g; // check if first group is included!!
  double newlen=pow(newV,1./3);
  group *gPtr;
  va_list ap;
  va_start(ap, g);
  while (n>0) { 
    //gPtr=&va_arg(ap, group);
    //gPtr->isobaricmove(*con, newlen);
    //cout << gPtr->name << endl;
    //N+=gPtr->nummolecules();
    n--;
  }
  */
}
double isobaric::move(vector<macromolecule> &g) {
  cnt++;
  du=0;
  newvolume();                                // Generate new trial volume - prep. trial potential
  unsigned short i=g.size();
  for (unsigned short n=0; n<i; n++)          // Loop over macromolecules
    g[n].isobaricmove(*con, pow(newV,1./3));  // ..and scale their mass-centers
  uold = pot->energy(con->p);                 // calc. old energy with original potential class
  unew = trialpot.energy(con->trial);         // calc. new energy using a COPY of the potential class
  du = unew-uold;
  dh = du + P*dV-(i+1)*log( newV/con->getvolume() );
  if (pen!=0) 
    dh += penalty[int(pow(newV,1./3))]-penalty[int(pow(con->getvolume(),1./3))];
  if (ens->metropolis(dh)==true) {
    accept();
    for (unsigned short n=0; n<i; n++) 
      g[n].accept(*con);
    if (pen!=0)
      newpenalty[int(pow(con->getvolume(),1./3))]+=pen;
    return du;
  } else rc=ENERGY;
  du=0;
  dpsqr+=0;
  for (unsigned short n=0; n<i; n++)
    g[n].undo(*con);
  vol += con->getvolume();
  vol2+= pow(con->getvolume(),2);
  if (pen!=0)
    newpenalty[int(pow(con->getvolume(),1./3))]+=pen;
  len += pow(con->getvolume(),1./3);
  len2+= pow(con->getvolume(),2./3);
  return du;
} 

//---------- TRANSLATE CLUSTER OF GROUP ----------------
clustertranslate::clustertranslate( ensemble &e,
    container &c, interaction<T_pairpot> &i, double S ) : markovmove(e,c,i) {
  name = "MACROMOLECULAR CLUSTER TRANSLATION";
  runfraction=0.5;
  deltadp=1.;
  dp=10.;
  sep=S;
  nacc.assign(30,0);
  ntrial.assign(30,0);
};

double clustertranslate::move(vector<macromolecule> &g) {
  du=0;
  if (slp.runtest(runfraction)==false)
    return du;
  unew=0;
  uold=0;
  FLOW=0;
  cnt++;
  point p;
  p.x=dp*slp.random_half();
  p.y=dp*slp.random_half();
  p.z=dp*slp.random_half();
  sort(g);


  short int iend=cluster.size();
  short int jend=free.size();

  ntrial[iend]++;

  for (short i=0; i<iend; i++)
    g[cluster[i]].move(*con, p);
  flowcheck(g);
  for (short i=0; i<iend; i++) {
    if (con->collision(g[i].cm_trial)==true || FLOW!=0) {
      rc=HC;
      dpsqr+=0;
      for (short i=0; i<iend; i++)
        g[i].undo(*con);
      return du;
    }
  }
  double Uold;
  double Unew;
  for (short i=0; i<iend; i++)
    for (short j=0; j<jend; j++) {
       Uold+= pot->energy( con->p, g[cluster[i]], g[free[j]] );
       Unew+= pot->energy( con->trial, g[cluster[i]], g[free[j]] );
    }
  du = Unew-Uold;
  if (ens->metropolis(du)==true ) {
    rc=OK;
    utot+=du;
    dpsqr+=con->sqdist( g[0].cm, g[0].cm_trial );
    naccept++;
    nacc[iend]++;
    for (short i=0; i<iend; i++)
      g[cluster[i]].accept(*con);
    return du;
  } else rc=ENERGY;
  du=0;
  dpsqr+=0;
  for (short i=0; i<iend; i++)
    g[cluster[i]].undo(*con);
  return du;
}

bool clustertranslate::decreasing(int i, int j) {
    return (i>j);
} 

void clustertranslate::flowcheck(vector<macromolecule> &g) {
  short int iend=cluster.size();
  short int jend=free.size();
  for (short i=0; i<iend; i++)
    for (short j=0; j<jend; j++) {
      if(coll.overlap( con->p, g[cluster[i]], g[free[j]], sep))
         FLOW++;
  }
}

void clustertranslate::sort(vector<macromolecule> &g) {
  int nprot=g.size();
  free.clear();
  cluster.clear();
  for (short i=0; i<nprot; i++) {
    free.push_back(i);
  }

  int seed=slp.random_one()*nprot;              //pick seed
  cluster.push_back(seed);                    //place it in cluster
  free.erase(free.begin() + seed);            //and remove it from the free
  for (short s=0; s<cluster.size(); s++)
    for (int i=0; i<free.size(); i++) {
      if (coll.overlap(con->p, g[cluster[s]] ,g[free[i]] ,sep)) {
        cluster.push_back(free[i]);
        free.erase(free.begin() + i);
        i--;  
      }
    }
  if (nprot != cluster.size()+free.size())
    cout << "CLUSTER SORT HAS FAILED"<<endl;
}

string clustertranslate::info() {
  ostringstream o;
  o << markovmove::info();
  o << "#   CLUSTER SIZES "<<endl
    << "# --------------------------------------------------------------" <<endl;
  for (short int i=0; i<30; i++) {
    if (ntrial[i]!=0)
      o << "#  "<<i<< "-particle cluster  : "<< 100*nacc[i]/ntrial[i] << " % acc., "<<ntrial[i]<< " attempts." <<endl
        << "# -------------------------------------------------------------- " <<endl;
  }
  return o.str();
}

//---------- ROTATE CLUSTER OF MACROMOLECULE ----------------
clusterrotate::clusterrotate( ensemble &e,
    container &c, interaction<T_pairpot> &i ) : markovmove(e,c,i) {
  name = "MACROMOLECULAR CLUSTER ROTATION";
  runfraction=0.8;
  deltadp=1.;
  dp=0.5;
  nacc.assign(30,0);
  ntrial.assign(30,0);
};

double clusterrotate::move(vector<macromolecule> &g) {
  du=0;
  uold=0;
  unew=0;
  FLOW=0;
  cnt++;
  point p,cr;
  p.ranunit(slp);
  sort(g);
  angle=slp.random_half()*dp;
  cr=g[cluster[0]].cm;

  short int iend=cluster.size();
  short int jend=free.size();

  ntrial[iend]++;

  for (short i=0; i<iend; i++)
    g[cluster[i]].rotate(*con, cr, p, angle );
  flowcheck(g);
  for (short i=0; i<iend; i++) {
    if (con->collision(g[i].cm_trial)==true || FLOW!=0 || cluster.size()>3) {
      rc=HC;
      dpsqr+=0;
      for (short i=0; i<iend; i++)
        g[i].undo(*con);
      return du;
    }
  }
  double Unew, Uold;
//  #pragma omp parallel for reduction (+:Uold) reduction (+:Unew) num_threads(4)
  for (short i=0; i<iend; i++)
    for (short j=0; j<jend; j++) {
       Uold+= pot->energy( con->p, g[cluster[i]], g[free[j]] );
       Unew+= pot->energy( con->trial, g[cluster[i]], g[free[j]] );
    }
  du = Unew-Uold;
  if (ens->metropolis(du)==true) {
    rc=OK;
    utot+=du;
    dpsqr+=con->sqdist( g[0].cm, g[0].cm_trial );
    naccept++;
    nacc[iend]++;
    for (short i=0; i<iend; i++)
      g[cluster[i]].accept(*con);
    return du;
  } else rc=ENERGY;
  du=0;
  dpsqr+=0;
  for (short i=0; i<iend; i++)
    g[cluster[i]].undo(*con);
  return du;
}

bool clusterrotate::decreasing(int i, int j) {
    return (i>j);
} 

void clusterrotate::flowcheck(vector<macromolecule> &g) {
  short int iend=cluster.size();
  short int jend=free.size();
  for (short i=0; i<iend; i++)
    for (short j=0; j<jend; j++) {
      if(coll.overlap( con->p, g[cluster[i]], g[free[j]], sep))
         FLOW++;
  }
}

void clusterrotate::sort(vector<macromolecule> &g) {
  int nprot=g.size();
  free.clear();
  cluster.clear();
  for (short i=0; i<nprot; i++) {
    free.push_back(i);
  }

  int seed=slp.random_one()*nprot;              //pick seed
  cluster.push_back(seed);                    //place it in cluster
  free.erase(free.begin() + seed);            //and remove it from the free
  for (short s=0; s<cluster.size(); s++)
    for (int i=0; i<free.size(); i++) {
      if (coll.overlap(con->p, g[cluster[s]] ,g[free[i]] ,sep)) {
        cluster.push_back(free[i]);
        free.erase(free.begin() + i);
        i--;  
      }
    }
  if (nprot != cluster.size()+free.size())
    cout << "CLUSTER SORT HAS FAILED"<<endl;
}

string clusterrotate::info() {
  ostringstream o;
  o << markovmove::info();
  o << "#   CLUSTER SIZES "<<endl
    << "# --------------------------------------------------------------" <<endl;
  for (short int i=0; i<30; i++) {
    if (ntrial[i]!=0)
      o << "#  "<<i<< "-particle cluster  : "<< 100*nacc[i]/ntrial[i] << " % acc., "<<ntrial[i]<< " attempts." <<endl
        << "# -------------------------------------------------------------- " <<endl;
  }
  return o.str();
}
