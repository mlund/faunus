#include "faunus/moves/crankshaft.h"

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
    len = minMonomers + (rand() % (maxMonomers-minMonomers+1));
  }

  bool crankShaft::findEnds(polymer &g) {
    vector<unsigned short> nb;
    v.clear();
    v.push_back( g.random() );
    for (int n=0; n<len+2; n++) {
      nb=g.neighbors(v[n]);
      if ( nb.size()==2 ) {
        for (int i=0; i<2; i++) {       // only atoms w. two bonds
          if (find(v.begin(), v.end(), nb[i])==v.end()) {
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
    runfraction=1.0;
    deltadp=0.;
    dp=in.getflt("crankshaft_dp", 1.0);
    minMonomers=in.getint("crankshaft_min",1);
    maxMonomers=in.getint("crankshaft_max",10);
    setNumberOfMonomers();
  }

  double crankShaft::move(polymer &g) {
    du=0;
    setNumberOfMonomers();
    if (slp.runtest(runfraction)==false || findEnds(g)==false)
      return du;
    cnt++;

    rot.setAxis( con->p[beg], con->p[end], dp*slp.random_half() );

    for (int i=0; i<v.size(); i++)
      con->trial[v[i]] = rot.rotate( con->p[v[i]] );

    for (int i=0; i<v.size(); i++) {
      if (con->collision(con->trial[v[i]])==true) {
        du=1e3;
        break;
      } else du += pot->u_monomer(con->trial,g,v[i]) - pot->u_monomer(con->p,g,v[i]);
    }

    if (ens->metropolis(du)==true) {
      rc=OK;
      utot+=du;
      naccept++;
      for (int i=0; i<v.size(); i++)
        con->p[v[i]] = con->trial[v[i]];
      g.masscenter(con->p);
      return du;
    } else rc=ENERGY;
    du=0;
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
      o << markovmove::info()
        << "#   Min number of monomers    = " << minMonomers << endl
        << "#   Max number of monomers    = " << maxMonomers << endl;
    }
    return o.str();
  }
} // namespace
