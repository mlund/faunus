#include <faunus/move.h>
#include <faunus/space.h>
#include <faunus/slump.h>
#include <faunus/group.h>
#include <faunus/energy.h>
#include <faunus/species.h>
#include <faunus/inputfile.h>
#include <faunus/geometry.h>

namespace Faunus {

  namespace Move {

    movebase::movebase(Energy::energybase &e, space &s, string pfx) : infty(1e15) {
      pot=&e;
      spc=&s;
      prefix=pfx;
      cnt=cnt_accepted=0;
      dusum=0;
      iw=22;
      runfraction=1;
    }

    //void move::unittest(unittest&) {
    //}

    void movebase::trialmove() { cnt++; }

    void movebase::acceptmove() { cnt_accepted++; }

    double movebase::move() {
      trialmove();
      double du=energychange();
      if ( !metropolis(du) ) {
        rejectmove();
        du=0;
      }
      else {
        acceptmove();
        dusum+=du;
      }
      return du;
    }

    bool movebase::metropolis(const double &du) const {
      if (du>0)
        if ( slp_global.random_one()>std::exp(-du) )
          return false;
      return true;
    }

    bool movebase::run() const {
      if (slp_global.random_one() < runfraction)
        return true;
      return false;
    }

    void movebase::pad(std::ostringstream& o) { o << "#   " << setw(iw) << std::left; }

    string movebase::info() {
      std::ostringstream o;
      o << endl << "# MARKOV MOVE: " << title << endl;
      if (!cite.empty()) {
        pad(o); o << "More information:" << cite << endl;
      }
      pad(o); o << "Runfraction" << runfraction << endl;
      if (cnt>0) {
        pad(o); o << "Number of trials" << cnt << endl;
        pad(o); o << "Acceptance" << double(cnt_accepted)/cnt*100 << "\ufe6a" << endl;
        pad(o); o << "Total energy change" << dusum << " kT" << endl;
      }
      return o.str();
    }

    // TRANSLATE

    /*
       translate::translate(string pfx, energybase &e, space &s) : energybase(pfx,e,s) {
       title="Molecular Translation";
       }

       void translate::trialmove() {
       }

       void translate::acceptmove() {
       cnt_accepted++;
       }

       void translate::rejectmove() {
       }

       double translate::energychange() {return 0;}

       double translate::move() {
       return 0;
       }

       string translate::info() {
       std::ostringstream o;
       o << movebase::info();
       pad(o); o << "Displacement vector" << dp << endl; 
       return o.str();
       }
       */

    translate_particle::translate_particle(inputfile &in,Energy::energybase &e, space &s, string pfx) : movebase(e,s,pfx) {
      title="Single Particle Displacement";
      igroup=iparticle=-1;
      dir.x=dir.y=dir.z=1;
      iw=25;
      in.getflt(prefix+"runfraction",1);
    }

    void translate_particle::trialmove() {
      movebase::trialmove();
      if (igroup>-1)
        iparticle=spc->g[igroup]->random();
      if (iparticle>-1) {
        double dp = atom[ spc->p[iparticle].id ].dp;
        spc->trial[iparticle].x += dir.x * dp * slp_global.random_half();
        spc->trial[iparticle].y += dir.y * dp * slp_global.random_half();
        spc->trial[iparticle].z += dir.z * dp * slp_global.random_half();
        spc->geo->boundary( spc->trial[iparticle] );
      }
    }

    void translate_particle::acceptmove() {
      movebase::acceptmove();
      double r2=spc->geo->_sqdist( spc->p[iparticle], spc->trial[iparticle] );
      sqrmap[ spc->p[iparticle].id ] += r2;
      accmap[ spc->p[iparticle].id ] += 1;
      spc->p[iparticle] = spc->trial[iparticle];
    }

    void translate_particle::rejectmove() {
      spc->trial[iparticle] = spc->p[iparticle];
      sqrmap[ spc->p[iparticle].id ] += 0;
      accmap[ spc->p[iparticle].id ] += 0;
    }

    double translate_particle::energychange() {
      if (iparticle>-1) {
        if ( spc->geo->collision( spc->trial[iparticle], Geometry::geometrybase::BOUNDARY ) ) {
          //cout << "!";
          return infty;
        }
        else
          return pot->i2all(spc->trial, iparticle) - pot->i2all(spc->p, iparticle);
      }
      return 0;
    }

    double translate_particle::move() {
      if (!run()) return 0;
      if (igroup>-1) {
        double du=0;
        for (int i=0; i<spc->g[igroup]->size(); i++) {
          iparticle = spc->g[igroup]->random();
          du+=movebase::move();
        }
        iparticle=-1;
        return du;
      } else return movebase::move();
    }

    string translate_particle::info() {
      std::ostringstream o;
      o << movebase::info();
      pad(o); o << "Displacement directions" << dir << endl;
      pad(o); o << "Mean-square displacement:" << endl;
      for (mapiter it=sqrmap.begin(); it!=sqrmap.end(); ++it) {
        pad(o); o << "" << atom[it->first].name << " " << it->second << endl;
      }
      pad(o); o << "Particle acceptance:" << endl;
      for (mapiter it=accmap.begin(); it!=accmap.end(); ++it) {
        pad(o); o << "" << atom[it->first].name << " " << (it->second).avg()*100 << "\ufe6a" << endl;
      }
      return o.str();
    }

  }//namespace
}//namespace
