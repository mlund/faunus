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
#include <faunus/physconst.h>

namespace Faunus {

  namespace Move {

    using namespace textio;

    Movebase::Movebase(Energy::Energybase &e, Space &s, string pfx) {
      pot=&e;
      spc=&s;
      prefix=pfx;
      cnt=cnt_accepted=0;
      dusum=0;
      w=22;
      runfraction=1;
    }

    Movebase::~Movebase() {
    }

    //void move::unittest(unittest&) {
    //}

    void Movebase::trialMove() {
      assert(spc->geo!=NULL); // space geometry MUST be set before moving!
      cnt++;
      _trialMove();
    }

    void Movebase::acceptMove() {
      cnt_accepted++;
      _acceptMove();
    }
    
    void Movebase::rejectMove() {
      _rejectMove();
    }
    
    double Movebase::energyChange() {
      return _energyChange();
    }
    
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
        if ( slp_global.randOne()>std::exp(-du) )
          return false;
      return true;
    }

    bool Movebase::run() const {
      if (slp_global.randOne() < runfraction)
        return true;
      return false;
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
      o << _info();
      return o.str();
    }

    // TRANSLATE

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

    void ParticleTranslation::_trialMove() {
      if (igroup!=NULL)
        iparticle=igroup->random();
      if (iparticle>-1) {
        double dp = atom[ spc->p[iparticle].id ].dp;
        spc->trial[iparticle].x += dir.x * dp * slp_global.randHalf();
        spc->trial[iparticle].y += dir.y * dp * slp_global.randHalf();
        spc->trial[iparticle].z += dir.z * dp * slp_global.randHalf();
        spc->geo->boundary( spc->trial[iparticle] );
      }
    }

    void ParticleTranslation::_acceptMove() {
      double r2=spc->geo->sqdist( spc->p[iparticle], spc->trial[iparticle] );
      sqrmap[ spc->p[iparticle].id ] += r2;
      accmap[ spc->p[iparticle].id ] += 1;
      spc->p[iparticle] = spc->trial[iparticle];
    }

    void ParticleTranslation::_rejectMove() {
      spc->trial[iparticle] = spc->p[iparticle];
      sqrmap[ spc->p[iparticle].id ] += 0;
      accmap[ spc->p[iparticle].id ] += 0;
    }

    double ParticleTranslation::_energyChange() {
      if (iparticle>-1) {
        if ( spc->geo->collision( spc->trial[iparticle], Geometry::Geometrybase::BOUNDARY ) )
          return pc::infty;
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

    string ParticleTranslation::_info() {
      char l=12;
      std::ostringstream o;
      o << pad(SUB,w,"Displacement vector") << dir << endl;
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

    RotateGroup::RotateGroup(InputMap &in,Energy::Energybase &e, Space &s, string pfx) : Movebase(e,s,pfx) {
      title="Group Rotation/Translation";
      igroup=NULL;
      w=30;
      dir.x=dir.y=dir.z=1;
      groupWiseEnergy=false;
      runfraction = in.get<double>(prefix+"_runfraction",1.0);
      dp_trans = in.get<double>(prefix+"_transdp", 2);
      dp_rot   = in.get<double>(prefix+"_rotdp", 100);
      if (dp_rot>4*pc::pi) // no need to rotate more than
        dp_rot=4*pc::pi;   // +/- 2 pi.
    }
    
    void RotateGroup::setGroup(Group &g) {
      assert(&g!=NULL);
      igroup=&g;
    }

    void RotateGroup::_trialMove() {
      assert(igroup!=NULL);
      angle=dp_rot*slp_global.randHalf();
      Point p;
      spc->geo->randompos(p);
      igroup->rotate(*spc, p, angle);
      p.x=dir.x * dp_trans * slp_global.randHalf();
      p.y=dir.y * dp_trans * slp_global.randHalf();
      p.z=dir.z * dp_trans * slp_global.randHalf();
      igroup->translate(*spc, p);
    }

    void RotateGroup::_acceptMove() {
      double r2 = spc->geo->sqdist( igroup->cm, igroup->cm_trial );
      sqrmap_t[ igroup->name ] += r2;
      sqrmap_r[ igroup->name ] += angle*angle;
      accmap[ igroup->name ] += 1;
      igroup->accept(*spc);
    }

    void RotateGroup::_rejectMove() {
      sqrmap_t[ igroup->name ] += 0;
      sqrmap_r[ igroup->name ] += 0;
      accmap[ igroup->name ] += 0;
      igroup->undo(*spc);
     }

    double RotateGroup::_energyChange() {
      for (int i=(*igroup).beg; i<=(*igroup).end; i++)
        if ( spc->geo->collision( spc->trial[i], Geometry::Geometrybase::BOUNDARY ) )
          return pc::infty;
      double uold = pot->g2all(spc->p, *igroup) + pot->g_external(spc->p, *igroup);
      double unew = pot->g2all(spc->trial, *igroup) + pot->g_external(spc->trial, *igroup);
      return unew-uold;
    }

    string RotateGroup::_info() {
      char l=12;
      std::ostringstream o;
      o << pad(SUB,w,"Displacement vector") << dir << endl
        << pad(SUB,w,"Max. translation") << pm << dp_trans/2 << textio::_angstrom << endl
        << pad(SUB,w,"Max. rotation") << pm << dp_rot/2*180/pc::pi << textio::degrees << endl;
      if (cnt>0) {
        o << endl
          << indent(SUB) << "Move Statistics:" << endl << endl
          << indent(SUBSUB) << std::left << setw(20) << "Group name" //<< string(20,' ')
          << setw(l+1) << "Acc. "+percent
          << setw(l+7) << rootof+bracket("dR"+squared)
          << setw(l+5) << rootof+bracket("d"+theta+squared) << endl;
        for (auto m : accmap) {
          string id=m.first;
          o << indent(SUBSUB) << std::left << setw(20) << id;
          o.precision(3);
          o << setw(l) << accmap[id].avg()*100
            << setw(l) << sqrt(sqrmap_t[id].avg())
            << setw(l) << sqrt(sqrmap_r[id].avg()) << endl;
        }
      }
      return o.str();
    }

  }//namespace
}//namespace
