#include "faunus/moves/crankshaft.h"
#include "faunus/energy/base.h"
#include "faunus/ensemble.h"

namespace Faunus {
  vectorRotation::vectorRotation(container &con) {
    conPtr = &con;
  }

  /*!
   * \param beg Starting point of the vector to rotate around
   * \param end Ending point of the vector to rotate around
   * \param angle How many degrees to rotate
   */
  void vectorRotation::setAxis(point beg, point end, double angle) {
    origin=beg;
    u=end-beg;
    conPtr->boundary(u);
    u=u*(1./conPtr->dist(beg,end));
    cosang=cos(angle);
    sinang=sin(angle);
    e1mcox=(1.-cosang)*u.x;
    e1mcoy=(1.-cosang)*u.y;
    e1mcoz=(1.-cosang)*u.z;
  }

  /*!
   * Rotate point around axis specified above.
   * \param p Particle to rotate
   */
  point vectorRotation::rotate(point p) {
    point b;
    b.x=p.x-origin.x;              // translate to origo...
    b.y=p.y-origin.y;
    b.z=p.z-origin.z;
    conPtr->boundary(b);                    // Apply boundary conditions
    eb=u.x*b.x + u.y*b.y + u.z*b.z;
    p.x=e1mcox*eb+cosang*b.x+sinang*(u.y*b.z - u.z*b.y) + origin.x;
    p.y=e1mcoy*eb+cosang*b.y+sinang*(u.z*b.x - u.x*b.z) + origin.y;
    p.z=e1mcoz*eb+cosang*b.z+sinang*(u.x*b.y - u.y*b.x) + origin.z;
    conPtr->boundary(p);
    return p;
  }

  void crankShaft::setNumberOfMonomers() {
    len = minMonomers + (slp.rand() % (maxMonomers-minMonomers+1));
  }
  
  void crankShaft::setNumberOfMonomers(polymer &g) {
    len = slp.rand() % g.size();
  }

  bool crankShaft::findEnds(polymer &g) {
    vector<int> nb;
    v.clear();
    v.push_back( g.random() );
    for (int n=0; n<len+2; n++) {
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
    beg=v.front();
    end=v.back();
    v.erase(v.begin());
    v.pop_back();
    return true;
  }

  crankShaft::crankShaft(ensemble &e, container &c, energybase &i, inputfile &in) : markovmove(e,c,i), rot(c) {
    name.assign("CRANKSHAFT MOVE");
    prefix.assign("crankshaft_");
    runfraction=1.0;
    deltadp=0.;
    dp=in.getflt("crankshaft_dp", 1.0);
    minMonomers=in.getint("crankshaft_min",0);
    maxMonomers=in.getint("crankshaft_max",0);
    if ( minMonomers==0 || maxMonomers==0 ) {
      useMinMaxMonomers=false;
      setNumberOfMonomers();//do we need this?
    }
    else {
      useMinMaxMonomers=true;
      setNumberOfMonomers();//do we need this?
    }
    markovmove::getInput(in);
  }

  double crankShaft::move(polymer &g) {
    if ( useMinMaxMonomers == true ) 
      setNumberOfMonomers();
    else
      setNumberOfMonomers(g);
    if (slp.runtest(runfraction)==false || findEnds(g)==false)
      return 0;
    markovmove::move();

    rot.setAxis( con->p[beg], con->p[end], dp*slp.random_half() );

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
   * This crankShaft move uses an expensive energy calculation, necessary for cm.z dependent potentials
   */
  double crankShaft::penaltymove(polymer &g) {
    if ( useMinMaxMonomers == true ) 
      setNumberOfMonomers();
    else
      setNumberOfMonomers(g);
    if (slp.runtest(runfraction)==false || findEnds(g)==false)
      return 0;
    markovmove::move();

    rot.setAxis( con->p[beg], con->p[end], dp*slp.random_half() );

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

  double crankShaft::move(polymer &g, int repeat) {
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

  string crankShaft::info() {
    std::ostringstream o;
    if (runfraction>0) {
      o << markovmove::info();
      if ( useMinMaxMonomers == true ) 
      o << "#   Min number of monomers    = " << minMonomers << endl
        << "#   Max number of monomers    = " << maxMonomers << endl;
      else
      o << "#   Min number of monomers    = " << 1 << endl
        << "#   Max number of monomers    = " << "polymer length" << endl;
    }
    return o.str();
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
        B=v.back() + 1 + ( rand() % ( g.size()-v.size() ) );
        B=rand() % g.size();        
        while ( B <= A ) 
          B=rand() % g.size();
        v.pop_back();
//cout << len << ": [" << v[0] << " " << v[1] << " " << v[2] << " (...) " << " " << v.back() << "] " << A << " " << B << endl;
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
        B=rand() % g.size();
        while ( B >= A ) 
          B=rand() % g.size();
        v.pop_back();
// cout << len << ": [" << v[0] << " " << v[1] << " " << v[2] << " (...) " << " " << v.back() << "] " << A << " " << B << endl;
        return true;
    }
  }

  branchRotation::branchRotation(ensemble &e, container &c, energybase &i, inputfile &in) : markovmove(e,c,i), rot(c) {
    name.assign("BRANCHROTATION MOVE");
    prefix.assign("branchrot_");
    runfraction=1.0;
    deltadp=0.;
    dp=in.getflt("branchrot_dp", 1.0);
    markovmove::getInput(in);
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
