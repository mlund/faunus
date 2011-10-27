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
      assert(spc->geo!=nullptr); // space geometry MUST be set before moving!
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
      igroup=nullptr;
      dir.x=dir.y=dir.z=1;
      w=25;
      in.get(prefix+"_runfraction",0.);
    }

    void ParticleTranslation::setGroup(Group &g) {
      igroup=&g;
      iparticle=-1;
    }

    void ParticleTranslation::setParticle(int i) {
      iparticle=i;
      igroup=nullptr;
    }

    void ParticleTranslation::_trialMove() {
      if (igroup!=nullptr)
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
      if (igroup!=nullptr) {
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
          << setw(l+7) << bracket("r"+squared)+"/"+angstrom+squared
          << rootof+bracket("r"+squared)+"/"+angstrom << endl;
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
      igroup=nullptr;
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
      assert(&g!=nullptr);
      igroup=&g;
    }

    //!< \todo Check random angle generation!
    void RotateGroup::_trialMove() {
      assert(igroup!=nullptr);
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
      sqrmap_r[ igroup->name ] += pow(angle*180/pc::pi, 2);
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
          << setw(l+9) << rootof+bracket("dR"+squared)+"/"+angstrom+squared
          << setw(l+5) << rootof+bracket("d"+theta+squared)+"/"+degrees << endl;
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

    Isobaric::Isobaric(InputMap &in, Energy::Hamiltonian &e, Space &s, string pfx) : Movebase(e,s,pfx) {
      title="Isobaric Volume Fluctuations";
      w=30;
      dV = in.get<double>(prefix+"_dV", 0.);
      P = in.get<double>(prefix+"_P", 0.)/1e30*pc::Nav; //pressure mM -> 1/A^3
      runfraction = in.get(prefix+"_runfraction",1.0);
      hamiltonian = &e;
      e.create( Energy::ExternalPressure( e.getGeometry(), P ) );
    }

    string Isobaric::_info() {
      using namespace textio;
      std::ostringstream o;
      const double tomM=1e30/pc::Nav;
      int N,Natom=0, Nmol=0;
      for (auto g : spc->g)
        if (g->id==Group::ATOMIC)
          Natom += g->size();
        else
          Nmol++;
      N = Natom + Nmol;
      double Pascal = P*pc::kB*pc::T*1e30;
      o << pad(SUB,w, "Displacement parameter") << dV << endl
        << pad(SUB,w, "Number of molecules") <<N<< " (" <<Nmol<< " molecular + " <<Natom<< " atomic)" << endl 
        << pad(SUB,w, "Pressure") << P*tomM << " mM" << " = " << Pascal << " Pa = " << Pascal/0.980665e5 << " atm" << endl
        << pad(SUB,w, "Temperature") << pc::T << " K" << endl;
      if (cnt>0) {
        char l=14;
        o << pad(SUB,w, "Mean displacement") << "\u221b"+rootof+bracket("dV"+squared) << " = " << pow(sqrV.avg(), 1/6.) << _angstrom << endl
          << pad(SUB,w, "Osmotic coefficient") << P / (N/V.avg()) << endl
          << endl
          << indent(SUBSUB) << std::right << setw(10) << ""
          << setw(l+5) << bracket("V")
          << setw(l+8) << "\u221b"+bracket("V")
          << setw(l+8) << bracket("N/V") << endl
          << indent(SUB) << setw(10) << "Averages"
          << setw(l) << V.avg() << _angstrom << cubed
          << setw(l) << pow(V.avg(),1/3.) << _angstrom
          << setw(l) << N/V.avg()*tomM << " mM" << endl;
      }
      return o.str();
    }

    void Isobaric::_setVolume(double V) {
      for (auto ebase : hamiltonian->baselist )
        if (&ebase->getGeometry()!=nullptr)
          ebase->getGeometry().setVolume(V);
    }

    void Isobaric::_trialMove() {
      assert(spc->g.size()>0 && "Space has empty group vector - NPT move not possible.");
      oldV = spc->geo->getVolume();
      newV = exp( log(oldV) + slp_global.randHalf()*dV );
      for (auto g : spc->g)
        g->scale(*spc, newV);
    }

    void Isobaric::_acceptMove() {
      V += newV;
      sqrV += pow( oldV-newV, 2 );
      _setVolume(newV);
      for (auto g : spc->g )
        g->accept(*spc);
    }

    void Isobaric::_rejectMove() {
      sqrV += 0;
      V += oldV;
      _setVolume(oldV);
      for (auto g : spc->g )
        g->undo(*spc);
    }

    double Isobaric::_energy(const p_vec &p) {
      double u=0;
      for (size_t i=0; i<spc->g.size()-1; ++i)      // group-group
        for (size_t j=i+1; j<spc->g.size(); ++j)
          u += pot->g2g(p, *spc->g[i], *spc->g[j]);
      for (auto g : spc->g) {
        u += pot->g_external(p, *g);
        if (g->id==Group::ATOMIC)
          u+=pot->g_internal(p, *g);
      }
      return u + pot->external();
    }

    double Isobaric::_energyChange() {
      double uold,unew;
      _setVolume( oldV );
      uold = _energy(spc->p);
      _setVolume( newV );
      unew = _energy(spc->trial);
      return unew-uold;
    }

  }//namespace
}//namespace
