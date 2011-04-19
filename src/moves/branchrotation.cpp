#include "faunus/moves/branchrotation.h"
#include "faunus/energy/base.h"
#include "faunus/ensemble.h"

namespace Faunus {
  branchRotation::branchRotation(ensemble &e, container &c, energybase &i, inputfile &in) : markovmove(e,c,i), rot(c) {
    name.assign("BRANCHROTATION MOVE");
    prefix.assign("branchrot_");
    runfraction=1.0;
    deltadp=0.;
    dp=in.getflt("branchrot_dp", 1.0);
    markovmove::getInput(in);
  }

  void branchRotation::setNumberOfMonomers(polymer &g) {
    len = slp.rand() % g.size();
    while ( len > g.size()/2 )
      len = slp.rand() % g.size();
  }

  bool branchRotation::findEnds(polymer &g) {
    vector<int> nb;
    v.clear();
    switch (rand() % 2) {     // random number between 0 and 1
      case 0: 
        v.push_back( g.beg );           // start with beg terminus
        v.push_back( g.beg+1 );         // add second atom
        for (int n=1; n<len; n++) {
          nb=g.neighbors(v[n]);
          if ( nb.size()==2 ) {
            for (int i=0; i<2; i++) {       // only atoms w. two bonds
              if (std::find(v.begin(), v.end(), nb[i])==v.end()) {
                v.push_back(nb[i]);
                break;
              }
            }
          } else return false; // abort of more than two bonds!
        }
        T=v.front();
        A=v.back();
        B=g.beg + (rand() % g.size());
        while ( B <= A ) 
          B=g.beg + (rand() % g.size());
        v.pop_back();
        return true;
      case 1: 
        v.push_back( g.end );
        v.push_back( g.end-1 );
        for (int n=1; n<len; n++) {
          nb=g.neighbors(v[n]);
          if ( nb.size()==2 ) {
            for (int i=0; i<2; i++) {       // only atoms w. two bonds
              if (std::find(v.begin(), v.end(), nb[i])==v.end()) {
                v.push_back(nb[i]);
                break;
              }
            }
          } else return false; // abort of more than two bonds!
        }
        T=v.front();
        A=v.back();
        B=g.beg + (rand() % g.size());
        while ( B >= A ) 
          B=g.beg + (rand() % g.size());
        v.pop_back();
        return true;
    }
    return false;
  }

  double branchRotation::move(polymer &g) {
    setNumberOfMonomers(g);
    if (slp.runtest(runfraction)==false || findEnds(g)==false)
      return 0;
    markovmove::move();

    rot.setAxis( con->p[A], con->p[B], dp*slp.random_half() );

    for (int i=0; i<v.size(); i++)
      con->trial[v[i]] = rot.rotate( con->p[v[i]] );

    bool hc=false;
    g.cm_trial = g.masscenter(*con, con->trial);
    if (con->slicecollision(g.cm_trial)==true) 
      hc=true;
    for (int i=0; i<v.size(); i++) {
      if (con->collision(con->trial[v[i]])==true) {
        hc=true;
        break;
      } else du += pot->u_monomer(con->trial,g,v[i]) - pot->u_monomer(con->p,g,v[i]);
    }

    if (hc==true) {
      rc=HC;
      du=0;
      dpsqr+=0;
      g.cm_trial=g.cm;
      for (int i=0; i<v.size(); i++)
        con->trial[v[i]] = con->p[v[i]];
      return du;
    }

    if (ens->metropolis(du)==true) {
      rc=OK;
      utot+=du;
      naccept++;
      dpsqr += g.sqmassgradius(*con, con->trial)-g.sqmassgradius(*con, con->p);
      for (int i=0; i<v.size(); i++)
        con->p[v[i]] = con->trial[v[i]];
      g.masscenter(*con);
      return du;
    } else rc=ENERGY;
    du=0;
    dpsqr+=0;
    g.cm_trial=g.cm;
    for (int i=0; i<v.size(); i++)
      con->trial[v[i]] = con->p[v[i]];
    return du;
  }

  /*!
   * This branchRotation move uses an expensive energy calculation, necessary for cm.z dependent potentials
   */
  double branchRotation::penaltymove(polymer &g) {
    setNumberOfMonomers(g);
    if (slp.runtest(runfraction)==false || findEnds(g)==false) 
      return 0;
    markovmove::move();

    rot.setAxis( con->p[A], con->p[B], dp*slp.random_half() );

    for (int i=0; i<v.size(); i++)
      con->trial[v[i]] = rot.rotate( con->p[v[i]] );

    bool hc=false;
    g.cm_trial = g.masscenter(*con, con->trial);
    if (con->slicecollision(g.cm_trial)==true) 
      hc=true;
    for (int i=0; i<v.size(); i++) {
      if (con->collision(con->trial[v[i]])==true) {
        hc=true;
        break;
      }
    }

    if (hc==true) {
      rc=HC;
      du=0;
      dpsqr+=0;
      g.cm_trial=g.cm;
      for (int i=0; i<v.size(); i++)
        con->trial[v[i]] = con->p[v[i]];
      return du;
    }

    uold = pot->energy(con->p, g) + pot->uself_polymer(con->p, g);
    unew = pot->energy(con->trial, g) + pot->uself_polymer(con->trial, g);
    du = unew-uold;

    if (ens->metropolis(du)==true) {
      rc=OK;
      utot+=du;
      naccept++;
      dpsqr += g.sqmassgradius(*con, con->trial)-g.sqmassgradius(*con, con->p);
      for (int i=0; i<v.size(); i++)
        con->p[v[i]] = con->trial[v[i]];
      g.masscenter(*con);
      return du;
    } else rc=ENERGY;
    du=0;
    dpsqr+=0;
    g.cm_trial=g.cm;
    for (int i=0; i<v.size(); i++)
      con->trial[v[i]] = con->p[v[i]];
    return du;
  }

  double branchRotation::move(polymer &g, int repeat) {
    double du=0;
    float rfbak=1.0;
    if (slp.runtest(runfraction)==false)
      return du;
    std::swap(rfbak,runfraction);
    while (repeat-- > 0)
      du+=move(g);
    std::swap(rfbak,runfraction);
    return du;
  }

  string branchRotation::info() {
    std::ostringstream o;
    if (runfraction>0) {
      o << markovmove::info();
    }
    return o.str();
  }
} // namespace
