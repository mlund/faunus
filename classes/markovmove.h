#ifndef _markovmove_h
#define _markovmove_h

#include "container.h"
#include "potentials.h"
#include "ensemble.h"
#include "titrate.h"
#include "slump.h"
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
 *  \code
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
/*! \brief Move one group along z-axis
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
    double r; //!< Distance between group mass centers
    dualmove( ensemble&, container&, interaction<T_pairpot>&);
    void direction(double, double, double);
    double move(macromolecule &, macromolecule &);
};

class translate : public markovmove {
  public: 
    translate( ensemble&, container&, interaction<T_pairpot>&);
    double move(macromolecule &); 
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
 */
class chargereg : public markovmove, private titrate {
  public:
    chargereg( ensemble&, container&, interaction<T_pairpot>&, group&, float);
    double titrateall();
    string info();
};
/*! \brief Grand Canonical titation of all sites
 *  \author Bjoern Persson
 */
/*
class GCchargereg : public markovmove, private GCtitrate {
  public: 
    GCchargereg( ensemble&, container&, interaction<T_pairpot>&, group&, float, float, float); //!< pH, CatPot, volume
    bool titrateall();
    string info();
*/
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
  dp=12;
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
    #pragma omp parallel
    {
      #pragma omp sections
      {
        #pragma omp section
        { uold = pot->energy(con->p, n);   }
        #pragma omp section
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
  runfraction=0.8;
  deltadp=1.;
  dp=1.;
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

  #pragma omp parallel
  {
    #pragma omp sections
    {
      #pragma omp section
      { uold = pot->energy(con->p, g);   }
      #pragma omp section
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
  #pragma omp parallel
  {
    #pragma omp sections
    {
      #pragma omp section
      { uold = pot->energy(con->p, g);   }
      #pragma omp section
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
  cite="Lund+Jonsson, Biochem. 2005, 44, 5722-5727.";
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
    #pragma omp parallel
    {
      #pragma omp sections
      {
        #pragma omp section
        { uold = pot->potential( con->p, t.site ) * con->p[t.site].charge
          + pot->potential( con->p, t.proton ) * con->p[t.proton].charge
            - con->p[t.site].potential(con->p[t.proton] )*con->p[t.proton].charge;
        }
        #pragma omp section
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
  }
  return sum;
}

/*
//-----------GRAND CANONICAL CHARGE REGULAION------------
GCchargereg::GCchargereg(
    ensemble &e, 
    container &c,
    interaction<T_pairpot> &i,
    group &g,
    float ph, float cp, float vol ) : markovmove(e,c,i), GCtitrate(c,c.p,g,ph,cp,vol)
{
  name="PROTON TITRATION";
  runfraction=0.2;
  con->trial = con->p;
}
bool GCchargereg::titrateall() {
  return true;
};
string GCchargereg::info() {
  ostringstream o;
  o <<  markovmove::info()
    << "#   pH (concentration)  = " << tit.ph << endl
    << "#   Cation potential(kT)= " << CatPot << endl
    << "#   Titrateable sites   = " << tit.sites.size() << endl
    << "#   Number of protons   = " << tit.protons.size() << endl;
  return o.str();
}  
*/

//-----------MOVE----------------------------------------
/*! \breif Class to prefom a random walk of a macromolecule
 *   in space
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
    #pragma omp parallel
    {
      #pragma omp sections
      {
        #pragma omp section
        { uold = pot->energy(con->p, g);   }
        #pragma omp section
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
    #pragma omp parallel
    {
      #pragma omp sections
      {
        #pragma omp section
        { uold = pot->energy(con->p, g);   }
        #pragma omp section
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
    container &c, interaction<T_pairpot> &i ) : markovmove(e,c,i) {
  name = "SYMMETRIC GROUP MOVE";
  cite = "Lund+Jonsson, Biophys J. 2003, 85, 2940";
  runfraction=1.0;
  deltadp=1.;
  dp=1.;
  direction(0,0,1);
}

/*! Specify unit vector that determines which coordinates
 * that will be moved. The default is (0,0,1) meaning that
 * the macromolecules will be moved in the z direction only.
 */
void dualmove::direction(double x, double y, double z) {
  v.x=x;
  v.y=y;
  v.z=z;
}
double dualmove::move(macromolecule &g1, macromolecule &g2) {
  r=con->dist(g1.cm, g2.cm_trial);
  du=0;
  cnt++;
  point p;
  group g12=g1+g2;
  p.x=v.x*dp*slp.random_half();
  p.y=v.y*dp*slp.random_half();
  p.z=v.z*dp*slp.random_half();
  g1.move(*con, p);
  g2.move(*con,-p);
  if (con->collision(g1.cm_trial)==true ||
      con->collision(g2.cm_trial)==true) {
    rc=HC;
    dpsqr+=0;
    g1.undo(*con);
    g2.undo(*con);
    return du;
  }
  #pragma omp parallel
  {
    #pragma omp sections
    {
      #pragma omp section
      { uold = pot->energy(con->p,g12) + pot->energy(con->p,g1,g2);   }
      #pragma omp section
      { unew = pot->energy(con->trial,g12) + pot->energy(con->trial,g1,g2);   }
    }
  }
  du = unew-uold;
  if (ens->metropolis(du)==true) {
    rc=OK;
    utot+=du;
    r=con->dist(g1.cm, g2.cm_trial);
    dpsqr+=r*r;
    naccept++;
    g1.accept(*con);
    g2.accept(*con);
    return du;
  } else rc=ENERGY;
  du=0;
  dpsqr+=0;
  g1.undo(*con);
  g2.undo(*con);
  return du;
}

