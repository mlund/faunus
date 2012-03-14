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

    void Movebase::trialMove() {
      assert(spc->geo!=nullptr && "Space geometry MUST be set before moving!");
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
   
    /*! \return Energy change in units of kT */
    double Movebase::energyChange() {
      return _energyChange();
    }
  
    /*!
     * This function will perform a trial move and accept/reject using the standard
     * Metropolis criteria. That is, it will
     * perform the following actions \c n times:
     * \li Perform a trial move with \c _trialMove()
     * \li Calulate the energy change, \f$\beta\Delta U\f$ with \c _energyChange()
     * \li Accept with probability \f$ \min(1,e^{-\beta\Delta U}) \f$
     * \li Call either \c _acceptMove() or \c _rejectMove()
     *
     * \param n Perform move \c n times
     */
    double Movebase::move(int n) {
      if (!run())
        return 0;
      double utot=0;
      while (n-->0) {
        trialMove();
        double du=energyChange();
        if ( !metropolis(du) ) {
          rejectMove();
          du=0;
        }
        else {
          acceptMove();
          dusum+=du;
          utot+=du;
        }
      }
      return utot;
    }

    bool Movebase::metropolis(const double &du) const {
      if (du>0)
        if ( slp_global.randOne()>std::exp(-du) ) // core of MC!
          return false;
      return true;
    }

    bool Movebase::run() const {
      if (slp_global.randOne() < runfraction)
        return true;
      return false;
    }

    void Movebase::test(UnitTest &t) {
      if (runfraction<1e-6 || cnt==0)
        return;
      t(prefix+"_acceptance", double(cnt_accepted)/cnt*100 );
      _test(t);
    }

    void Movebase::_test(UnitTest&) {
    }

    double Movebase::getAcceptance() {
      if (cnt>0)
        return double(cnt_accepted) / cnt;
      return 0;
    }

    /*!
     * This will return a formatted multi-line information string about the move and
     * will as a minimum contain:
     * \li Name of move
     * \li Runfraction
     * \li Number of times the move has been called
     * \li Acceptance
     * \li Total energy change
     *
     * Typically additional information will be provided as well.
     */
    string Movebase::info() {
      assert(!title.empty() && "Markov Moves must have a title");
      std::ostringstream o;
      if (runfraction<1e-10)
        return o.str();
      o << header("Markov Move: " + title);
      if (!cite.empty())
        o << pad(SUB,w,"More information:") << cite << endl;
      if (cnt>0)
        o << pad(SUB,w,"Number of trials") << cnt << endl
          << pad(SUB,w,"Acceptance") << getAcceptance()*100 << percent << endl
          << pad(SUB,w,"Runfraction") << runfraction*100 << percent << endl
          << pad(SUB,w,"Total energy change") << dusum << kT << endl;
      o << _info();
      return o.str();
    }

    AtomicTranslation::AtomicTranslation(InputMap &in,Energy::Energybase &e, Space &s, string pfx) : Movebase(e,s,pfx) {
      title="Single Particle Translation";
      iparticle=-1;
      igroup=nullptr;
      dir.x=dir.y=dir.z=1;
      w=30; //width of output
      in.get(prefix+"_runfraction",0.);
    }

    void AtomicTranslation::setGroup(Group &g) {
      igroup=&g;
      iparticle=-1;
    }

    void AtomicTranslation::setParticle(int i) {
      iparticle=i;
      igroup=nullptr;
    }

    void AtomicTranslation::_trialMove() {
      if (igroup!=nullptr) {
        iparticle=igroup->random();
        gsize += igroup->size();
      }
      if (iparticle>-1) {
        double dp = atom[ spc->p[iparticle].id ].dp;
        assert(iparticle<(int)spc->p.size() && "Trial particle out of range");
        assert(dp>1e-9 && "Atomic displacement parameter seems to be zero");
        spc->trial[iparticle].x += dir.x * dp * slp_global.randHalf();
        spc->trial[iparticle].y += dir.y * dp * slp_global.randHalf();
        spc->trial[iparticle].z += dir.z * dp * slp_global.randHalf();
        spc->geo->boundary( spc->trial[iparticle] );
      }
    }

    void AtomicTranslation::_acceptMove() {
      double r2=spc->geo->sqdist( spc->p[iparticle], spc->trial[iparticle] );
      sqrmap[ spc->p[iparticle].id ] += r2;
      accmap[ spc->p[iparticle].id ] += 1;
      spc->p[iparticle] = spc->trial[iparticle];
    }

    void AtomicTranslation::_rejectMove() {
      spc->trial[iparticle] = spc->p[iparticle];
      sqrmap[ spc->p[iparticle].id ] += 0;
      accmap[ spc->p[iparticle].id ] += 0;
    }

    double AtomicTranslation::_energyChange() {
      if (iparticle>-1) {
        assert( spc->geo->collision(spc->p[iparticle])==false && "An untouched particle collides with simulation container.");
        if ( spc->geo->collision( spc->trial[iparticle], Geometry::Geometrybase::BOUNDARY ) )
          return pc::infty;
        return
          pot->i_total(spc->trial, iparticle) - pot->i_total(spc->p, iparticle);
      }
      return 0;
    }

    string AtomicTranslation::_info() {
      char l=12;
      std::ostringstream o;
      if (gsize.cnt>0)
        o << pad(SUB,w,"Average moves/particle") << cnt / gsize.avg() << endl;
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

    TranslateRotate::TranslateRotate(InputMap &in,Energy::Energybase &e, Space &s, string pfx) : Movebase(e,s,pfx) {
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
      if (dp_rot<1e-6 && dp_trans<1e-6)
        runfraction=0;
    }

    void TranslateRotate::setGroup(Group &g) {
      assert(&g!=nullptr);
      igroup=&g;
    }

    void TranslateRotate::_trialMove() {
      assert(igroup!=nullptr);
      Point p;
      if (dp_rot>1e-6) {
        angle=dp_rot*slp_global.randHalf();
        double r=2;
        while (r>1) {
          p.x=2*slp_global.randHalf(); // random vector
          p.y=2*slp_global.randHalf(); // inside a sphere
          p.z=2*slp_global.randHalf();
          r=sqrt(p.x*p.x+p.y*p.y+p.z*p.z);
        }
        p=igroup->cm+p; // set endpoint for rotation 
        igroup->rotate(*spc, p, angle);
      }
      if (dp_trans>1e-6) {
        p.x=dir.x * dp_trans * slp_global.randHalf();
        p.y=dir.y * dp_trans * slp_global.randHalf();
        p.z=dir.z * dp_trans * slp_global.randHalf();
        igroup->translate(*spc, p);
      }
    }

    void TranslateRotate::_acceptMove() {
      double r2 = spc->geo->sqdist( igroup->cm, igroup->cm_trial );
      sqrmap_t[ igroup->name ] += r2;
      sqrmap_r[ igroup->name ] += pow(angle*180/pc::pi, 2);
      accmap[ igroup->name ] += 1;
      igroup->accept(*spc);
    }

    void TranslateRotate::_rejectMove() {
      sqrmap_t[ igroup->name ] += 0;
      sqrmap_r[ igroup->name ] += 0;
      accmap[ igroup->name ] += 0;
      igroup->undo(*spc);
    }

    double TranslateRotate::_energyChange() {
      if (dp_rot<1e-6 && dp_trans<1e-6)
        return 0;

      for (auto i : *igroup)
        if ( spc->geo->collision( spc->trial[i], Geometry::Geometrybase::BOUNDARY ) )
          return pc::infty;
      double unew = pot->g_external(spc->trial, *igroup);
      if (unew==pc::infty)
        return pc::infty;   // early rejection!
      unew += pot->g2all(spc->trial, *igroup);
      double uold = pot->g2all(spc->p, *igroup) + pot->g_external(spc->p, *igroup);
      return unew-uold;
    }

    string TranslateRotate::_info() {
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
          << setw(l+9) << rootof+bracket("dR"+squared)+"/"+angstrom
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

    void TranslateRotate::_test(UnitTest &t) {
      for (auto m : accmap) {                                                                                   
        string id=m.first,
               idtrim="_"+textio::trim(id)+"_";
        t(prefix+idtrim+"acceptance", accmap[id].avg()*100);
        t(prefix+idtrim+"dRot", sqrt(sqrmap_r[id].avg()));
        t(prefix+idtrim+"dTrans", sqrt(sqrmap_t[id].avg()));
      }
    }

    /*!
     * The cluster threshold is set with the InputMap keyword \c transrot_clustersize and
     * gives specifies the SURFACE distance between macromolecular particles and atomic particles
     * to be clustered.
     */
    TranslateRotateCluster::TranslateRotateCluster(InputMap &in,Energy::Energybase &e, Space &s, string pfx) : TranslateRotate(in,e,s,pfx) {
      title="Cluster "+title;
      threshold = in.get<double>(prefix+"_clustersize",0);
      gmobile=nullptr;
    }

    void TranslateRotateCluster::setMobile(Group &g) {
      assert(&g!=nullptr);
      gmobile=&g;
    }

    string TranslateRotateCluster::_info() {
      using namespace textio;
      std::ostringstream o;
      o << TranslateRotate::_info() << endl;
      o << pad(SUB,w,"Cluster threshold") << threshold << _angstrom << endl;
      if (cnt>0) {
        o << pad(SUB,w,"Average cluster size") << avgsize.avg() << endl;
        if (threshold>1e-9)
          o << pad(SUB,w,"Average bias") << avgbias.avg() << " (0=reject, 1=accept)" << endl; 
      }
      return o.str();
    }

    void TranslateRotateCluster::_trialMove() {
      assert(gmobile!=nullptr && "Cluster group not defined");
      assert(igroup!=nullptr && "Group to move not defined");
      Point p;

      // find clustered particles
      cindex.clear();
      for (auto i : *gmobile)
        if (ClusterProbability(spc->p, i) > slp_global.randOne() )
          cindex.push_back(i); // generate cluster list
      avgsize += cindex.size();

      // rotation
      if (dp_rot>1e-6) {
        angle=dp_rot*slp_global.randHalf();
        double r=2;
        while (r>1) {
          p.x=2*slp_global.randHalf(); // random vector
          p.y=2*slp_global.randHalf(); // inside a sphere
          p.z=2*slp_global.randHalf();
          r=sqrt(p.x*p.x+p.y*p.y+p.z*p.z);
        }
        p=igroup->cm+p; // set endpoint for rotation 
        igroup->rotate(*spc, p, angle);
        vrot.setAxis(*spc->geo, igroup->cm, p, angle); // rot. around line between CM and point
        for (auto i : cindex)
          spc->trial[i] = vrot.rotate( *spc->geo, spc->p[i] ); // (boundaries are accounted for)
      }

      // translation
      if (dp_trans>1e-6) {
        p.x=dir.x * dp_trans * slp_global.randHalf();
        p.y=dir.y * dp_trans * slp_global.randHalf();
        p.z=dir.z * dp_trans * slp_global.randHalf();
        igroup->translate(*spc, p);
        for (auto i : cindex) {
          spc->trial[i] += p;
          spc->geo->boundary(spc->trial[i]);
        }
      }
    }

    void TranslateRotateCluster::_acceptMove() {
      TranslateRotate::_acceptMove();
      for (auto i : cindex)
        spc->p[i] = spc->trial[i];
    }

    void TranslateRotateCluster::_rejectMove() {
      TranslateRotate::_rejectMove();
      for (auto i : cindex)
        spc->trial[i] = spc->p[i];
    }

    double TranslateRotateCluster::_energyChange() {
      double bias=1;             // cluster bias -- see Frenkel 2nd ed, p.405
      vector<int> imoved=cindex; // index of moved particles
      for (auto l : *gmobile)
        if (std::find(cindex.begin(), cindex.end(), l)==cindex.end())
          bias *= ( 1-ClusterProbability(spc->trial, l) ) / ( 1-ClusterProbability(spc->p, l) );
      avgbias += bias;
      if (bias<1e-7)
        return pc::infty;        // don't bother to continue with energy calculation

      if (dp_rot<1e-6 && dp_trans<1e-6)
        return 0;

      for (auto i : *igroup)     // Add macromolecule to list of moved particle index
        imoved.push_back(i);

      // container boundary collision?
      for (auto i : imoved)
        if ( spc->geo->collision( spc->trial[i], Geometry::Geometrybase::BOUNDARY ) )
          return pc::infty;

      // external potential on macromolecule
      double unew = pot->g_external(spc->trial, *igroup);
      if (unew==pc::infty)
        return pc::infty; //early rejection!
      double uold = pot->g_external(spc->p, *igroup);

      // external potentia on clustered atomic species
      for (auto i : cindex) {
        uold += pot->i_external(spc->p, i);
        unew += pot->i_external(spc->trial, i);
      }

      // pair energy between static and moved particles
      double du=0;
#pragma omp parallel for reduction (+:du)
      for (int j=0; j<(int)spc->p.size(); j++)
        if ( std::find(imoved.begin(), imoved.end(), j )==imoved.end() )
          for (auto i : imoved)
            du += pot->i2i(spc->trial, i, j) - pot->i2i(spc->p, i, j);
      return unew - uold + du - log(bias);
    }

    double TranslateRotateCluster::ClusterProbability(p_vec &p, int i) {
      for (auto j : *igroup)
        if (i!=j) {
          double r=threshold+p[i].radius+p[j].radius;
          if (spc->geo->sqdist(p[i],p[j])<r*r )
            return 1;
        }
      return 0;
    }

    Isobaric::Isobaric(InputMap &in, Energy::Hamiltonian &e, Space &s, string pfx) : Movebase(e,s,pfx) {
      title="Isobaric Volume Fluctuations";
      w=30;
      dV = in.get<double>(prefix+"_dV", 0.);
      P = in.get<double>(prefix+"_P", 0.)/1e30*pc::Nav; //pressure mM -> 1/A^3
      runfraction = in.get(prefix+"_runfraction",1.0);
      if (dV<1e-6)
        runfraction=0;
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

    void Isobaric::_test(UnitTest &t) {
      t(prefix+"_averageSideLength", pow( V.avg(), 1/3.) );
      t(prefix+"_MSQDisplacement", pow(sqrV.avg(), 1/6.) );
    }

    void Isobaric::_trialMove() {
      assert(spc->g.size()>0 && "Space has empty group vector - NPT move not possible.");
      oldV = spc->geo->getVolume();
      newV = std::exp( std::log(oldV) + slp_global.randHalf()*dV );
      for (auto g : spc->g)
        g->scale(*spc, newV);
    }

    void Isobaric::_acceptMove() {
      V += newV;
      sqrV += pow( oldV-newV, 2 );
      hamiltonian->setVolume(newV);
      for (auto g : spc->g )
        g->accept(*spc);
    }

    void Isobaric::_rejectMove() {
      sqrV += 0;
      V += oldV;
      hamiltonian->setVolume(oldV);
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
      uold = _energy(spc->p);
      hamiltonian->setVolume( newV );
      for (auto g : spc->g) // In spherical geometries molecules may collide with cell boundary upon scaling mass center.
        for (auto i : *g)
          if ( spc->geo->collision( spc->trial[i], Geometry::Geometrybase::BOUNDARY ) )
            return pc::infty;
      unew = _energy(spc->trial);
      return unew-uold;
    }

    bool AtomTracker::empty() {
      return map.empty();
    }

    AtomTracker::Tid AtomTracker::randomAtomType() const {
      assert(!map.empty() && "No atom types have been added yet");
      vector<Tid> vid;
      vid.reserve( map.size() );
      for (auto &m : map) 
        vid.push_back(m.first);
      std::random_shuffle(vid.begin(), vid.end());
      return vid.front();
    }

    void AtomTracker::clear() {
      map.clear();
    }

    AtomTracker::AtomTracker(Space &s) { spc=&s; }

    AtomTracker::Tindex AtomTracker::data::random() {
      std::random_shuffle( index.begin(), index.end() );
      return index.back();
    }

    AtomTracker::data& AtomTracker::operator[](AtomTracker::Tid id) { return map[id]; }

    bool AtomTracker::insert(const particle &a, Tindex index) {
      assert( a.id == spc->p[ map[a.id].index.back() ].id && "Id mismatch");
      spc->insert(a, index);
      for (auto &m : map)
        for (auto &i : m.second.index)
          if (i>=index) i++;
      map[a.id].index.push_back(index);
      return true;
    }

    bool AtomTracker::erase(AtomTracker::Tindex index) {
      spc->erase(index);
      bool deleted=false;
      for (auto &m : map) {
        auto f=std::find(m.second.index.begin(), m.second.index.end(), index);
        if (f!=m.second.index.end()) {
          m.second.index.erase(f);
          deleted=true;
          break;
        }
      }
      if (deleted)
        for (auto &m : map)
          for (auto &i : m.second.index)
            if (i>index) i--;

#ifndef NDEBUG
      assert(deleted && "Could not delete specified index");
      for (auto &m : map) {
        for (auto &i : m.second.index)
          assert( m.first == spc->p[i].id && "Particle id mismatch");
      }
#endif
      return deleted;
    }

    GrandCanonicalSalt::GrandCanonicalSalt(InputMap &in, Energy::Hamiltonian &e, Space &s, Group &g, string pfx) :
      Movebase(e,s,pfx), tracker(s) {
        title="Grand Canonical Salt";
        w=30;
        runfraction = in.get(prefix+"_runfraction",1.0);
        e.add(Urest);
        saltPtr=&g;
        add(*saltPtr);
      }

    void GrandCanonicalSalt::add(Group &g) {
      if (g.property.find(Group::ATOMIC)==g.property.end() ) {
        assert( !"Salt group must be atomic" );
        return;
      }
      g.property.insert(Group::GRANDCANONICAL);  // mark given group as grand canonical
      spc->enroll(g);
      tracker.clear();
      for (auto i : g) {
        short id=spc->p[i].id;
        if ( atom[id].activity>1e-10 && abs(atom[id].charge)>1e-10 ) {
          map[id].p=atom[id];
          map[id].chempot=log( atom[id].activity*pc::Nav*1e-27); // beta mu
          tracker[id].index.push_back(i);
        }
      }
      assert(!tracker.empty() && "No GC ions found!");
    }

    void GrandCanonicalSalt::randomIonPair(short &id_cation, short &id_anion) {
      do id_anion  = tracker.randomAtomType(); while ( map[id_anion].p.charge>=0);
      do id_cation = tracker.randomAtomType(); while ( map[id_cation].p.charge<=0  );
      assert( !tracker[id_anion].index.empty() && "Ion list is empty");
      assert( !tracker[id_cation].index.empty() && "Ion list is empty");
    }

    void GrandCanonicalSalt::_trialMove() {
      trial_insert.clear();
      trial_delete.clear();
      randomIonPair(ida, idb);
      assert(ida>0 && idb>0 && "Ion pair id is zero (UNK). Is this really what you want?");
      int Na = abs(map[idb].p.charge);
      int Nb = abs(map[ida].p.charge);
      switch ( rand() % 2) {
        case 0:
          trial_insert.reserve(Na+Nb);
          do trial_insert.push_back( map[ida].p ); while (--Na>0);
          do trial_insert.push_back( map[idb].p ); while (--Nb>0);
          for (auto &p : trial_insert)
            spc->geo->randompos(p);
          break;
        case 1:
          trial_delete.reserve(Na+Nb);
          while ( (int)trial_delete.size()!=Na ) {
            int i=tracker[ida].random();
            assert( ida==spc->p[i].id && "id mismatch");
            if (std::find(trial_delete.begin(), trial_delete.end(), i)==trial_delete.end())
              trial_delete.push_back(i);
          }
          while ( (int)trial_delete.size()!=Na+Nb ) {
            int i=tracker[idb].random();
            assert( idb==spc->p[i].id && "id mismatch");
            if (std::find(trial_delete.begin(), trial_delete.end(), i)==trial_delete.end())
              trial_delete.push_back(i);
          }
          assert( (int)trial_delete.size()==Na+Nb );
          break;
      }
    }

    double GrandCanonicalSalt::_energyChange() {
      du_rest=0;
      int Na=0, Nb=0;                     // number of added or deleted ions
      double idfactor=1;
      double uold=0, unew=0, V=spc->geo->getVolume();
      if ( !trial_insert.empty() ) {
        for (auto &t : trial_insert)     // count added ions
          if (t.id==map[ida].p.id) Na++; else Nb++;
        for (int n=0; n<Na; n++)
          idfactor *= (tracker[ida].index.size()+1+n)/V;
        for (int n=0; n<Nb; n++)
          idfactor *= (tracker[idb].index.size()+1+n)/V;

        unew = log(idfactor) - Na*map[ida].chempot - Nb*map[idb].chempot;
        du_rest=unew;

        unew += pot->v2v(spc->p, trial_insert);
        for (auto i=trial_insert.begin(); i!=trial_insert.end()-1; i++)
          for (auto j=i+1; j!=trial_insert.end(); j++)
            unew+=pot->p2p(*i,*j);
        for (auto i=trial_insert.begin(); i!=trial_insert.end(); i++)
          unew+=pot->p_external(*i);
      }
      else if ( !trial_delete.empty() ) {
        for (auto i : trial_delete) {
          if (spc->p[i].id==map[ida].p.id) Na++;
          else if (spc->p[i].id==map[idb].p.id) Nb++;
        }
        for (int n=0; n<Na; n++)
          idfactor *= (tracker[ida].index.size()-Na+1+n)/V;
        for (int n=0; n<Nb; n++)
          idfactor *= (tracker[idb].index.size()-Nb+1+n)/V;

        unew = -log(idfactor) + Na*map[ida].chempot + Nb*map[idb].chempot;
        du_rest=unew;

        for (auto &i : trial_delete)
          uold+=pot->i_total(spc->p, i);
        for (auto i=trial_delete.begin(); i!=trial_delete.end()-1; i++)
          for (auto j=i+1; j!=trial_delete.end(); j++)
            uold-=pot->i2i(spc->p, *i, *j);
      } else {
        assert(!"No salt to insert or delete!");
        cerr << "!! No salt to insert or delete !!";
      }
      return unew-uold;
    }

    void GrandCanonicalSalt::_acceptMove() {
      if ( !trial_insert.empty() ) {
        for (auto &p : trial_insert)
          tracker.insert(p, saltPtr->back());
      }
      else if ( !trial_delete.empty() ) {
        std::sort(trial_delete.rbegin(), trial_delete.rend()); //reverse sort
        for (auto i : trial_delete)
          tracker.erase(i);
      }
      Urest.add(du_rest);
      double V = spc->geo->getVolume();
      map[ida].rho += tracker[ida].index.size() / V;
      map[idb].rho += tracker[idb].index.size() / V;
    }

    void GrandCanonicalSalt::_rejectMove() {
      double V = spc->geo->getVolume();
      map[ida].rho += tracker[ida].index.size() / V;
      map[idb].rho += tracker[idb].index.size() / V;
    }

    string GrandCanonicalSalt::_info() {
      char s=10;
      using namespace textio;
      std::ostringstream o;
      o << pad(SUB,w,"Number of GC species") 
        << endl << endl;
      o << setw(4) << "" << std::left << setw(s) << "Ion" << setw(s)
        << "activity" << setw(s+4) << bracket("c/M") << setw(s)
        << bracket( gamma+pm ) << endl;
      for (auto &m : map) {
        short id=m.first;
        o.precision(5);
        o << setw(4) << "" << setw(s) << atom[id].name
          << setw(s) << atom[id].activity << setw(s) << m.second.rho.avg()/pc::Nav/1e-27
          << setw(s) << atom[id].activity / (m.second.rho.avg()/pc::Nav/1e-27)
          << endl;
      }
      return o.str();
    }

  }//namespace
}//namespace
