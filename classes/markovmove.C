#include "markovmove.h"
//--------------- MARKOV MOVE ---------------------
bool markovmove::run(float p) { return (p>slp.random_one())?true:false; }
float markovmove::accepted() { return naccept/float(cnt); }
// \param min Minimum percentage of accepted moves
// \param max Maximum percentage of accepted moves
void markovmove::adjust_dp(float min, float max) {
  float a=accepted()*100.;
  if (a>max) dp+=deltadp;
  if (a<min) dp-=deltadp;
  if (dp<=0) dp=deltadp;
}

//-------------- SALT MOVE ---------------------------------
saltmove::saltmove(
    ensemble &e, container &c, interaction<T_pairpot> &i ) : markovmove(e,c,i) {
  dp=60;
  deltadp=2;
  name="SALT DISPLACEMENTS";
}

/*! \param group Group containing mobile ions
 */
bool saltmove::move(group &g) {
  du=0;
  if (slp.runtest(runfraction)==false)
    return false;
  double sum=0;
  for (unsigned short i=0; i<g.size(); i++) {
    move(g, 1);
    sum+=du;
  }
  du=sum;
  return true;
}
bool saltmove::move(group &g, int n) {
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
      utot+=du;
      naccept++;
      con->p[n] = con->trial[n];
      return true;
    } else rc=ENERGY;
  }
  du=0;
  con->trial[n] = con->p[n];
  return false;
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
  runfraction=0.2;
  con->trial = con->p;
}

/*! \brief Exchange protons between bulk and titrateable sites.
 *
 *  This move will randomly go through the titrateable sites and
 *  try to exchange protons with the bulk. The trial energy is:
 */
bool chargereg::titrateall() {
  du=0;
  if (slp.runtest(runfraction)==false)
    return false;
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
  du=sum;
  return true;
}
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
    b->bpc(con->trial[i]);
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
/*! \breif Class to make moves of a single particle along z-axis
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

