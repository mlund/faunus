#include <faunus/move.h>
#include <faunus/space.h>
#include <faunus/slump.h>
#include <faunus/group.h>
#include <faunus/energy.h>
#include <faunus/species.h>
#include <faunus/inputfile.h>
#include <faunus/geometry.h>
#include <faunus/faunus.h>
#include <faunus/textio.h>

namespace Faunus {

  namespace Move {

    using namespace textio;

    Movebase::Movebase(Energy::Energybase &e, Space &s, string pfx) : infty(1e15) {
      pot=&e;
      spc=&s;
      prefix=pfx;
      cnt=cnt_accepted=0;
      dusum=0;
      w=22;
      runfraction=1;
    }

    //void move::unittest(unittest&) {
    //}

    void Movebase::trialMove() {
      assert(spc->geo!=NULL); // space geometry MUST be set before moving!
      cnt++;
    }

    void Movebase::acceptMove() { cnt_accepted++; }

    double Movebase::move() {
      trialMove();
      double du=energyChange();
      if ( !metropolis(du) ) {
        rejectMove();
        du=0;
      }
      else {
        acceptMove();
        dusum+=du;
      }
      return du;
    }

    bool Movebase::metropolis(const double &du) const {
      if (du>0)
        if ( slp_global.random_one()>std::exp(-du) )
          return false;
      return true;
    }

    bool Movebase::run() const {
      if (slp_global.random_one() < runfraction)
        return true;
      return false;
    }

    double Movebase::totalEnergy(){
      return 0;
    }

    string Movebase::info() {
      std::ostringstream o;
      o << header("Markov Move: " + title);
      if (!cite.empty()) {
        o << pad(SUB,w,"More information:") << cite << endl;
      }
      o << pad(SUB,w,"Runfraction") << runfraction*100 << percent << endl;
      if (cnt>0) {
        o << pad(SUB,w,"Number of trials") << cnt << endl;
        o << pad(SUB,w,"Acceptance") << double(cnt_accepted)/cnt*100 << percent << endl;
        o << pad(SUB,w,"Total energy change") << dusum << kT << endl;
      }
      return o.str();
    }

    // TRANSLATE

    /*
       translate::translate(string pfx, Energybase &e, space &s) : Energybase(pfx,e,s) {
       title="Molecular Translation";
       }

       void translate::trialmove() {
       }

       void translate::acceptMove() {
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
       o << Movebase::info();
       pad(o); o << "Displacement vector" << dp << endl; 
       return o.str();
       }
       */

    ParticleTranslation::ParticleTranslation(InputMap &in,Energy::Energybase &e, Space &s, string pfx) : Movebase(e,s,pfx) {
      title="Single Particle Translation";
      iparticle=-1;
      igroup=NULL;
      dir.x=dir.y=dir.z=1;
      w=25;
      in.get(prefix+"runfraction",0.0);
    }

    void ParticleTranslation::setGroup(Group &g) {
      igroup=&g;
      iparticle=-1;
    }

    void ParticleTranslation::setParticle(int i) {
      iparticle=i;
      igroup=NULL;
    }

    void ParticleTranslation::trialMove() {
      Movebase::trialMove();
      if (igroup!=NULL)
        iparticle=igroup->random();
      if (iparticle>-1) {
        double dp = atom[ spc->p[iparticle].id ].dp;
        spc->trial[iparticle].x += dir.x * dp * slp_global.random_half();
        spc->trial[iparticle].y += dir.y * dp * slp_global.random_half();
        spc->trial[iparticle].z += dir.z * dp * slp_global.random_half();
        spc->geo->boundary( spc->trial[iparticle] );
      }
    }

    void ParticleTranslation::acceptMove() {
      Movebase::acceptMove();
      double r2=spc->geo->sqdist( spc->p[iparticle], spc->trial[iparticle] );
      sqrmap[ spc->p[iparticle].id ] += r2;
      accmap[ spc->p[iparticle].id ] += 1;
      spc->p[iparticle] = spc->trial[iparticle];
    }

    void ParticleTranslation::rejectMove() {
      spc->trial[iparticle] = spc->p[iparticle];
      sqrmap[ spc->p[iparticle].id ] += 0;
      accmap[ spc->p[iparticle].id ] += 0;
    }

    double ParticleTranslation::energyChange() {
      if (iparticle>-1) {
        if ( spc->geo->collision( spc->trial[iparticle], Geometry::Geometrybase::BOUNDARY ) )
          return infty;
        return
          pot->i_total(spc->trial, iparticle) - pot->i_total(spc->p, iparticle);
          //pot->i2all(spc->trial, iparticle) - pot->i2all(spc->p, iparticle);
      }
      return 0;
    }

    double ParticleTranslation::move() {
      if (!run())
        return 0;
      if (igroup!=NULL) {
        double du=0;
        for (int i=0; i<igroup->size(); i++) {
          iparticle = igroup->random(); // set random particle for trialmove()
          if ( atom[spc->p[iparticle].id].dp > 1e-5 )
            du+=Movebase::move();
        }
        iparticle=-1;
        return du;
      } else return Movebase::move();
    }

    double ParticleTranslation::totalEnergy() {
      if (iparticle>-1)
        return pot->i_total(spc->p, iparticle);
      double u=0;
      if (igroup!=NULL) {
        u = pot->g2all(spc->p, *igroup) + pot->g_internal(spc->p, *igroup);
        for (int i=igroup->beg; i<=igroup->end; ++i)
          u += pot->i_external(spc->p, i);
      }
      return u;
    }

    string ParticleTranslation::info() {
      char l=12;
      std::ostringstream o;
      o << Movebase::info()
        << pad(SUB,w,"Displacement vector") << dir << endl;
      if (cnt>0) {
        o << endl
          << indent(SUB) << "Individual particle movement:" << endl << endl
          << indent(SUBSUB) << std::left << string(7,' ')
          << setw(l-6) << "dp"
          << setw(l+1) << "Acc. "+percent
          << setw(l+5) << bracket("r"+squared) // "\u27e8\u0394r\u00b2\u27e9" 
          << rootof+bracket("r"+squared) << endl; //"\u221a\u27e8\u0394r\u00b2\u27e9" << endl;
        for (auto m : sqrmap) {
          short id=m.first;
          o << indent(SUBSUB) << std::left << setw(7) << atom[id].name
            << setw(l-6) << atom[id].dp;
          o.precision(3);
          o << setw(l) << accmap[id].avg()*100
            << setw(l) << sqrmap[id].avg()
            << setw(l) << sqrt(sqrmap[id].avg()) << endl;
        }
      }
      return o.str();
    }

  }//namespace
}//namespace
