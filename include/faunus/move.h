#ifndef FAUNUS_MOVE_H
#define FAUNUS_MOVE_H

#ifndef SWIG
#include <faunus/common.h>
#include <faunus/point.h>
#include <faunus/average.h>
#include <faunus/textio.h>
#include <faunus/geometry.h>
#include <faunus/energy.h>
#include <faunus/textio.h>

#ifdef ENABLE_MPI
#include <faunus/mpi.h>
#endif

#endif

namespace Faunus {

  /** @brief Monte Carlo move related classes */
  namespace Move {

    template<typename Tkey=std::string>
      class AcceptanceMap {
        private:
          typedef std::map<Tkey, Average<double> > map_type;
          map_type accmap;   //!< move acceptance ratio
          map_type sqrmap;   //!< mean square displacement map
        public:
          void accept(Tkey k, double msq) {
            accmap[k]+=1;
            sqrmap[k]+=msq;
          }
          void reject(Tkey k) {
            accmap[k]+=0;
          }
          string info(char l=10) {
            using namespace textio;
            std::ostringstream o;
            o << indent(SUB) << "Move Statistics:" << endl
              << indent(SUBSUB) << std::left << setw(20) << "Id"
              << setw(l+1) << "Acc. "+percent
              << setw(l+9) << rootof+bracket("msq"+squared)+"/"+angstrom << endl;
            for (auto m : accmap) {
              string id=m.first;
              o << indent(SUBSUB) << std::left << setw(20) << id;
              o.precision(3);
              o << setw(l) << accmap[id].avg()*100
                << setw(l) << sqrt(sqrmap[id].avg()) << endl;
            }
            return o.str();
          }
          void _test(UnitTest &t, const string &prefix) {
            for (auto &m : accmap)  {
              std::ostringstream o;
              o << m.first;
              t(prefix+"_Acceptance"+o.str(), m.second.avg());
            }
          }
      };

    /**
     * @brief Add polarization step to an arbitrary move
     *
     * This class will modify any MC move to account for polarization
     * using an iterative procedure.
     * An electric field calculation is inserted
     * after the original trial move whereafter it will iteratively
     * calculate induced dipole moments on all particles.
     * The energy change function will evaluate the *total*
     * system energy as all dipoles in the system may have changed.
     * This is thus an expensive computation and is best used with
     * MC moves that propagate many or all particles.
     *
     * @todo Unfinished - fix polarization catastrophy!
     */
    template<class Tmove>
      class PolarizeMove : public Tmove {
        private:
          using Tmove::spc;
          using Tmove::pot;
          int max_iter;                   // max numbr of iterations
          double threshold;       	  // threshold for iteration
          Eigen::MatrixXd field;  	  // field on each particle
          Average<int> numIter;           // average number of iterations per move
          bool broke_loop;
          bool groupBasedField;

          /**
           *  @brief Replaces dipole moment with permanent dipole moment plus induced dipole moment
           *  @param pot Hamiltonian
           *  @param p Particles to update
           */
          template<typename Tenergy,typename Tparticles>
            void induceDipoles(Tenergy &pot, Tparticles &p) { 
              Eigen::VectorXd mu_err_norm((int)p.size());
              int cnt=0;
              do {
                cnt++;
                mu_err_norm.setZero();
                field.setZero();
                pot.field(p,field);
                for (size_t i=0; i<p.size(); i++) {
                  Point E = field.col(i); // field on i, in e/Å
                  Point mu_trial = p[i].alpha*E + p[i].mup; // New tot dipole
                  Point mu_err = mu_trial - p[i].mu*p[i].muscalar;     // Difference between former and current state
                  mu_err_norm[i] = mu_err.norm();// Get norm of previous row
                  p[i].muscalar = mu_trial.norm();// Update dip scalar in particle
                  if (p[i].muscalar > 1e-6)
                    p[i].mu = mu_trial/p[i].muscalar;// Update article dip.
                }
                if(cnt > max_iter) {
                  cout << "Reached " << max_iter << " iterations. Breaking loop!" << endl;
                  broke_loop = true;
                  break;
                }
              } while (mu_err_norm.maxCoeff() > threshold);                 // Check if threshold is ok
              numIter+=cnt; // average number of iterations
            }

          void _trialMove() FOVERRIDE {
            Tmove::_trialMove();                     // base class MC move
            field.resize(3,Tmove::spc->trial.size());// match sizes
            induceDipoles(*Tmove::pot,Tmove::spc->trial);
          }

          double _energyChange() FOVERRIDE {
            return Energy::systemEnergy(*spc,*pot,spc->trial) - Energy::systemEnergy(*spc,*pot,spc->p);
          }

          void _rejectMove() FOVERRIDE {
            Tmove::_rejectMove();
            Tmove::spc->trial = Tmove::spc->p;
          }

          void _acceptMove() FOVERRIDE {
            Tmove::_acceptMove();
            Tmove::spc->p = Tmove::spc->trial;
          }

          string _info() FOVERRIDE {
            std::ostringstream o;
            using namespace textio;
            o << pad(SUB,Tmove::w,"Polarization iterations") << numIter.avg() << endl;
            if(broke_loop)
              o << "Maximum number of iterations reached. Loop was broken!" << endl;
            o  << Tmove::_info();
            return o.str();
          }

        public:
          template<class Tspace>
            PolarizeMove(InputMap &in, Energy::Energybase<Tspace> &e, Tspace &s) :
              Tmove(in,e,s) {
                broke_loop = false;
                threshold = in.get<double>("pol_threshold", 0.001, "Iterative polarization precision");
                max_iter = in.get<int>("max_iterations", 40, "Maximum number of iteratins");
                groupBasedField = in.get<bool>("pol_g2g", false, "Group based field calculation");
              }
      };

    /**
     * @brief Base class for Monte Carlo moves
     *
     * The is a base class that handles Monte Carlo moves and derived classes
     * are required to implement the following pure virtual (and private)
     * functions:
     *
     * - `_trialMove()`
     * - `_energyChange()`
     * - `_acceptMove()`
     * - `_rejectMove()`
     * - `_info()`
     *
     * These functions should be pretty self-explanatory and are - via wrapper
     * functions - called by move(). It is important that the _energyChange() function
     * returns the full energy associated with the move. For example, for NPT
     * moves the pV term should be included and so on. Try not to override
     * the move() function as this should be generic to all MC moves.
     *
     * @date Lund, 2007-2011
     */
    template<class Tspace>
      class Movebase {
        private:
          unsigned long int cnt_accepted;  //!< number of accepted moves
          double dusum;                    //!< Sum of all energy changes

          virtual void _test(UnitTest&);   //!< Unit testing
          virtual void _trialMove()=0;     //!< Do a trial move
          virtual void _acceptMove()=0;    //!< Accept move and config
          virtual void _rejectMove()=0;    //!< Reject move and config
          virtual double _energyChange()=0;//!< Energy change of move (kT)

          void acceptMove();               //!< Accept move (wrapper)
          void rejectMove();               //!< Reject move (wrapper)
          double energyChange();           //!< Energy (wrapper)
          bool metropolis(const double&) const;//!< Metropolis criteria

        protected:
          virtual string _info()=0;        //!< info for derived moves
          void trialMove();                //!< Do a trial move (wrapper)
          Energy::Energybase<Tspace>* pot; //!< Pointer to energy functions
          Tspace* spc;                     //!< Pointer to Space
          string title;                    //!< Title of move (mandatory!)
          string cite;                     //!< Reference, url, DOI etc.
          string prefix;                   //!< inputmap prefix
          char w;                          //!< info text width. Adjust this in constructor if needed.
          unsigned long int cnt;           //!< total number of trial moves
          virtual bool run() const;        //!< Runfraction test

          bool useAlternateReturnEnergy;   //!< Return a different energy than returned by _energyChange(). [false]
          double alternateReturnEnergy;    //!< Alternative return energy

        public:
          Movebase(Energy::Energybase<Tspace>&, Tspace&, string);//!< Constructor
          virtual ~Movebase();
          double runfraction;                //!< Fraction of times calling move() should result in an actual move. 0=never, 1=always.
          double move(int=1);                //!< Attempt \c n moves and return energy change (kT)
          string info();                     //!< Returns information string
          void test(UnitTest&);              //!< Perform unit test
          double getAcceptance();            //!< Get acceptance [0:1]
      };

    template<class Tspace>
      Movebase<Tspace>::Movebase(Energy::Energybase<Tspace> &e, Tspace &s, string pfx) {
        e.setSpace(s);
        pot=&e;
        spc=&s;
        prefix=pfx;
        cnt=cnt_accepted=0;
        dusum=0;
        w=22;
        runfraction=1;
        useAlternateReturnEnergy=false; //this has no influence on metropolis sampling!
      }

    template<class Tspace>
      Movebase<Tspace>::~Movebase() {}

    template<class Tspace>
      void Movebase<Tspace>::trialMove() {
        if (cnt==0)
          for (auto i : spc->groupList())
            i->setMassCenter(*spc);
        cnt++;
        _trialMove();
      }

    template<class Tspace>
      void Movebase<Tspace>::acceptMove() {
        cnt_accepted++;
        _acceptMove();
      }

    template<class Tspace>
      void Movebase<Tspace>::rejectMove() {
        _rejectMove();
      }

    /** @return Energy change in units of kT */
    template<class Tspace>
      double Movebase<Tspace>::energyChange() {
        double du = _energyChange();
        if (std::isnan(du)) {
          std::cerr << "Error: energy change returns not-a-number (NaN)\n";
          std::exit(1);
        }
        return du;
      }

    /**
     * This function will perform a trial move and accept/reject using the
     * Metropolis criteria (doi:10/ds736f).
     * That is, it will perform the following `n` times:
     *
     * - Perform a trial move with `_trialMove()`
     * - Calulate the energy change, \f$\beta\Delta U\f$ with `_energyChange()`
     * - Accept with probability \f$ \min(1,e^{-\beta\Delta U}) \f$
     * - Call either `_acceptMove()` or `_rejectMove()`
     *
     * @note Do not override this function in derived classes.
     * @param n Perform move `n` times (default=1)
     */
    template<class Tspace>
      double Movebase<Tspace>::move(int n) {
        double utot=0;
        if (run()) {
          while (n-->0) {
            trialMove();
            double du=energyChange();
            if ( !metropolis(du) )
              rejectMove();
            else {
              acceptMove();
              if (useAlternateReturnEnergy)
                du=alternateReturnEnergy;
              dusum+=du;
              utot+=du;
            }
          }
        }
        assert(spc->p == spc->trial && "Trial particle vector out of sync!");
        return utot;
      }

    /**
     * @param du Energy change for MC move (kT)
     * @return True if move should be accepted; false if not.
     * @note
     * One could put in `if (du>0)` before the first line, but
     * certain MPI communications require the random number
     * generator to be in sync, i.e. each rank must call
     * `slp_global()` equal number of times, independent of
     * dU.
     */
    template<class Tspace>
      bool Movebase<Tspace>::metropolis(const double &du) const {
        if ( slp_global()>std::exp(-du) ) // core of MC!
          return false;
        return true;
      }

    template<class Tspace>
      bool Movebase<Tspace>::run() const {
        if (slp_global() < runfraction)
          return true;
        return false;
      }

    template<class Tspace>
      void Movebase<Tspace>::test(UnitTest &t) {
        if (runfraction<1e-6 || cnt==0)
          return;
        t(prefix+"_acceptance", double(cnt_accepted)/cnt*100 );
        _test(t);
      }

    template<class Tspace>
      void Movebase<Tspace>::_test(UnitTest&) {
      }

    template<class Tspace>
      double Movebase<Tspace>::getAcceptance() {
        if (cnt>0)
          return double(cnt_accepted) / cnt;
        return 0;
      }

    /**
     * This will return a formatted multi-line information string about the move and
     * will as a minimum contain:
     *
     * - Name of move
     * - Runfraction
     * - Number of times the move has been called
     * - Acceptance
     * - Total energy change
     *
     * Typically, additional information will be provided as well.
     *
     * @note Do not override in derived classes - use _info().
     */
    template<class Tspace>
      string Movebase<Tspace>::info() {
        using namespace textio;
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

    /**
     * @brief Translation of atomic particles
     *
     * This Markov move can work in two modes:
     * - Move a single particle in space set by setParticle()
     * - Move single particles randomly selected in a Group set by setGroup().
     *
     * The move directions can be controlled with the dir vector - for instance if you wish
     * to translate only in the `z` direction, set `dir.x()=dir.y()=0`.
     *
     * @date Lund, 2011
     */
    template<class Tspace>
      class AtomicTranslation : public Movebase<Tspace> {
        private:
          typedef Movebase<Tspace> base;
          typedef std::map<short, Average<double> > map_type;
          bool run() const;                //!< Runfraction test
        protected:
          string _info();
          void _acceptMove() FOVERRIDE;
          void _rejectMove() FOVERRIDE;
          double _energyChange() FOVERRIDE;
          void _trialMove() FOVERRIDE;
          using base::spc;
          map_type accmap; //!< Single particle acceptance map
          map_type sqrmap; //!< Single particle mean square displacement map

          int iparticle;   //!< Select single particle to move (-1 if none, default)
          Group* igroup;   //!< Group pointer in which particles are moved randomly (NULL if none, default)
          double genericdp;//!< Generic atom displacement parameter - ignores individual dps
          Average<unsigned long long int> gsize; //!< Average size of igroup;

        public:
          AtomicTranslation(InputMap&, Energy::Energybase<Tspace>&, Tspace&, string="mv_particle");
          void setGroup(Group&); //!< Select group in which to randomly pick particles from
          void setParticle(int); //!< Select single particle in Space::p to move
          void setGenericDisplacement(double); //!< Set single displacement for all atoms
          Point dir;             //!< Translation directions (default: x=y=z=1)
      };

    /**
     * The InputMap is searched for the following keywords:
     *
     * Key                  | Description
     * :------------------- | :-------------------------------------------------------------
     * `prefix_runfraction` | Chance of running (default=1)
     * `prefix_genericdp`   | Fallback displacement paraemter if `dp` is defined in AtomData.
     *
     * The standard prefix is `mv_particle`.
     */
    template<class Tspace>
      AtomicTranslation<Tspace>::AtomicTranslation(InputMap &in,Energy::Energybase<Tspace> &e,
          Tspace &s, string pfx) : Movebase<Tspace>(e,s,pfx) {
        base::title="Single Particle Translation";
        iparticle=-1;
        igroup=nullptr;
        dir.x()=dir.y()=dir.z()=1;
        this->w=30; //width of output
        this->runfraction = in.get<double>(pfx+"_runfraction",1.);
        setGenericDisplacement( in.get<double>(pfx+"_genericdp",0) );
      }

    /**
     * The generic displacement parameter will be used only if the specific
     * atomic dp is zero.
     */
    template<class Tspace>
      void AtomicTranslation<Tspace>::setGenericDisplacement(double dp) {
        genericdp=dp;
      }

    template<class Tspace>
      void AtomicTranslation<Tspace>::setGroup(Group &g) {
        igroup=&g;
        iparticle=-1;
      }

    template<class Tspace>
      void AtomicTranslation<Tspace>::setParticle(int i) {
        iparticle=i;
        igroup=nullptr;
      }

    template<class Tspace>
      bool AtomicTranslation<Tspace>::run() const {
        if ( igroup->empty() )
          return false;
        return base::run();
      }

    template<class Tspace>
      void AtomicTranslation<Tspace>::_trialMove() {
        if (igroup!=nullptr) {
          iparticle=igroup->random();
          gsize += igroup->size();
        }
        if (iparticle>-1) {
          double dp = atom[ spc->p.at(iparticle).id ].dp;
          if (dp<1e-6) dp = genericdp;
          assert(iparticle<(int)spc->p.size()
              && "Trial particle out of range");
          Point t = dir*dp;
          t.x() *= slp_global()-0.5;
          t.y() *= slp_global()-0.5;
          t.z() *= slp_global()-0.5;
          spc->trial[iparticle].translate(spc->geo, t);

          // make sure trial mass center is updated for molecular groups
          // (certain energy functions may rely on up-to-date mass centra)
          auto gi = spc->findGroup(iparticle);
          assert(gi!=nullptr);
          assert((gi->cm - gi->cm_trial).squaredNorm()<1e-6);
          if (gi->isMolecular())
            gi->cm_trial = Geometry::massCenter(spc->geo,spc->trial,*gi);

#ifndef NDEBUG
          // are untouched particles in group synched?
          for (auto j : *gi)
            if (j!=iparticle)
              assert((base::spc->p[j] - base::spc->trial[j]).squaredNorm()<1e-6);
#endif
        }
      }

    template<class Tspace>
      void AtomicTranslation<Tspace>::_acceptMove() {
        double r2=spc->geo.sqdist( spc->p[iparticle], spc->trial[iparticle] );
        sqrmap[ spc->p[iparticle].id ] += r2;
        accmap[ spc->p[iparticle].id ] += 1;
        spc->p[iparticle] = spc->trial[iparticle];
        auto gi = spc->findGroup(iparticle);
        assert(gi!=nullptr);
        if (gi->isMolecular())
          gi->cm = gi->cm_trial;
      }

    template<class Tspace>
      void AtomicTranslation<Tspace>::_rejectMove() {
        spc->trial[iparticle] = spc->p[iparticle];
        sqrmap[ spc->p[iparticle].id ] += 0;
        accmap[ spc->p[iparticle].id ] += 0;
        auto gi = spc->findGroup(iparticle);
        assert(gi!=nullptr);
        if (gi->isMolecular())
          gi->cm_trial = gi->cm;
      }

    template<class Tspace>
      double AtomicTranslation<Tspace>::_energyChange() {
        if (iparticle>-1) {
          assert( spc->geo.collision(spc->p[iparticle])==false
              && "An untouched particle collides with simulation container.");
          if ( spc->geo.collision(
                spc->trial[iparticle], Geometry::Geometrybase::BOUNDARY ) )
            return pc::infty;
          return
            (base::pot->i_total(spc->trial, iparticle)
             + base::pot->external(spc->trial))
            - (base::pot->i_total(spc->p, iparticle)
                + base::pot->external(spc->p));
        }
        return 0;
      }

    template<class Tspace>
      string AtomicTranslation<Tspace>::_info() {
        using namespace textio;
        std::ostringstream o;
        if (gsize.cnt>0)
          o << pad(SUB,base::w,"Average moves/particle")
            << base::cnt / gsize.avg() << endl;
        o << pad(SUB,base::w,"Displacement vector")
          << dir.transpose() << endl;
        if (genericdp>1e-6)
          o << pad(SUB,base::w,"Generic displacement")
            << genericdp << _angstrom << endl;
        if (base::cnt>0) {
          char l=12;
          o << endl
            << indent(SUB) << "Individual particle movement:" << endl << endl
            << indent(SUBSUB) << std::left << string(7,' ')
            << setw(l-6) << "dp"
            << setw(l+1) << "Acc. "+percent
            << setw(l+7) << bracket("r"+squared)+"/"+angstrom+squared
            << rootof+bracket("r"+squared)+"/"+angstrom << endl;
          for (auto m : sqrmap) {
            auto id=m.first;
            o << indent(SUBSUB) << std::left << setw(7) << atom[id].name
              << setw(l-6) << ( (atom[id].dp<1e-6) ? genericdp : atom[id].dp);
            o.precision(3);
            o << setw(l) << accmap[id].avg()*100
              << setw(l) << sqrmap[id].avg()
              << setw(l) << sqrt(sqrmap[id].avg()) << endl;
          }
        }
        return o.str();
      }

    /**
     * @brief Rotate single particles
     *
     * This move works in the same way as AtomicTranslation but does
     * rotations of non-isotropic particles instead of translation. This move
     * has no effect on isotropic particles such as Faunus::PointParticle.
     */
    template<class Tspace>
      class AtomicRotation : public AtomicTranslation<Tspace> {
        protected:
          typedef AtomicTranslation<Tspace> base;
          using base::spc;
          using base::iparticle;
          using base::igroup;
          using base::w;
          using base::gsize;
          using base::genericdp;
          Geometry::QuaternionRotate rot;
          string _info();
          void _trialMove();
        public:
          AtomicRotation(InputMap&, Energy::Energybase<Tspace>&, Tspace&, string="rot_particle");
      };

    template<class Tspace>
      AtomicRotation<Tspace>::AtomicRotation(InputMap &in,Energy::Energybase<Tspace> &e, Tspace &s, string pfx) : base(in,e,s,pfx) {
        base::title="Single Particle Rotation";
      }

    template<class Tspace>
      void AtomicRotation<Tspace>::_trialMove() {
        if (igroup!=nullptr) {
          iparticle=igroup->random();
          gsize += igroup->size();
        }
        if (iparticle>-1) {
          assert(iparticle<(int)spc->p.size() && "Trial particle out of range");
          double dprot = atom[spc->p[iparticle].id ].dprot;
          if (dprot<1e-6)
            dprot = base::genericdp;

          Point u;
          u.ranunit(slp_global);
          rot.setAxis(spc->geo, Point(0,0,0), u, dprot*slp_global.randHalf() );
          spc->trial[iparticle].rotate(rot);
        }
      }

    template<class Tspace>
      string AtomicRotation<Tspace>::_info() {
        using namespace textio;
        std::ostringstream o;
        if (gsize.cnt>0)
          o << pad(SUB,w,"Average moves/particle") << base::cnt / gsize.avg() << endl;
        if (base::genericdp>1e-6)
          o << pad(SUB,w,"Generic displacement") << genericdp
            << _angstrom << endl;
        if (base::cnt>0) {
          char l=12;
          o << endl
            << indent(SUB) << "Individual particle rotation:" << endl << endl
            << indent(SUBSUB) << std::left << string(7,' ')
            << setw(l-6) << "dp"
            << setw(l+1) << "Acc. "+percent
            << setw(l+7) << bracket("r"+squared)+"/"+angstrom+squared
            << rootof+bracket("r"+squared)+"/"+angstrom << endl;
          for (auto m : base::sqrmap) {
            particle::Tid id=m.first;
            o << indent(SUBSUB) << std::left << setw(7) << atom[id].name
              << setw(l-6) << ( (atom[id].dprot<1e-6) ? genericdp : atom[id].dp);
            o.precision(3);
            o << setw(l) << base::accmap[id].avg()*100
              << setw(l) << base::sqrmap[id].avg()
              << setw(l) << sqrt(base::sqrmap[id].avg()) << endl;
          }
        }
        return o.str();
      }

    /**
     * @brief Combined rotation and rotation of groups
     *
     * This will translate and rotate groups and collect averages based on group name.
     *
     * Example:
     *
     * ~~~
     * ...
     * Group g;
     * g.name="mygroup";
     * Move::TranslateRotate tr(in, pot, spc);
     * tr.directions[g.name].z()=0; // move only on xy plane - do this only once and before calling setGroup()
     * tr.setGroup(g);              // specify which group to move
     * tr.move();                   // do the move
     * ~~~
     */
    template<class Tspace>
      class TranslateRotate : public Movebase<Tspace> {
        private:
          typedef Movebase<Tspace> base;
        protected:
          using base::spc;
          using base::pot;
          using base::w;
          using base::cnt;
          using base::prefix;
          void _test(UnitTest&);
          void _trialMove();
          void _acceptMove();
          void _rejectMove();
          double _energyChange();
          string _info();
          typedef std::map<string, Average<double> > map_type;
          map_type accmap;   //!< Group particle acceptance map
          map_type sqrmap_t; //!< Group mean square displacement map (translation)
          map_type sqrmap_r; //!< Group mean square displacement map (rotation)
          Group* igroup;
          double dp_rot;     //!< Rotational displament parameter
          double dp_trans;   //!< Translational displacement parameter
          double angle;      //!< Temporary storage for current angle
          Point dir;         //!< Translation directions (default: x=y=z=1). This will be set by setGroup()
        public:
          TranslateRotate(InputMap&, Energy::Energybase<Tspace>&, Tspace&, string="transrot");
          void setGroup(Group&); //!< Select Group to move
          bool groupWiseEnergy;  //!< Attempt to evaluate energy over groups from vector in Space (default=false)
          std::map<string,Point> directions; //!< Specify special group translation directions (default: x=y=z=1)
#ifdef ENABLE_MPI
          Faunus::MPI::MPIController* mpi;
#endif
      };

    /**
     * The InputMap is scanned for the following keys:
     *
     * Key               | Description
     * :---------------- | :-------------------------------------
     * `pfx_transdp`     | Translational displacement [angstrom]
     * `pfx_rotdp`       | Rotational displacement [radians]
     */
    template<class Tspace>
      TranslateRotate<Tspace>::TranslateRotate(InputMap &in,Energy::Energybase<Tspace> &e, Tspace &s, string pfx) : base(e,s,pfx) {
        base::title="Group Rotation/Translation";
        base::w=30;
        igroup=nullptr;
        dir.x()=dir.y()=dir.z()=1;
        groupWiseEnergy=false;
        dp_trans = in.get<double>(base::prefix+"_transdp", 2,
            "Group translationsal displacement (AA)");
        dp_rot = in.get<double>(base::prefix+"_rotdp", 3,
            "Group rotational displacement (rad)");
        if (dp_rot>4*pc::pi) // no need to rotate more than
          dp_rot=4*pc::pi;   // +/- 2 pi.
        this->runfraction = in.get<double>(base::prefix+"_runfraction",1.0);
        if (dp_rot<1e-6 && dp_trans<1e-6)
          this->runfraction=0;
#ifdef ENABLE_MPI
        mpi=nullptr;
#endif
      }

    template<class Tspace>
      void TranslateRotate<Tspace>::setGroup(Group &g) {
        assert(&g!=nullptr);
        assert(!g.name.empty() && "Group should have a name.");
        assert(g.isMolecular());
        assert(spc->geo.sqdist(g.cm,g.cm_trial)<1e-6 && "Trial CM mismatch");
        igroup=&g;
        if ( directions.find(g.name) != directions.end() )
          dir = directions[g.name];
        else
          dir.x() = dir.y() = dir.z() = 1;
      }

    template<class Tspace>
      void TranslateRotate<Tspace>::_trialMove() {
        assert(igroup!=nullptr);
        Point p;
        if (dp_rot>1e-6) {
          p.ranunit(slp_global);             // random unit vector
          p=igroup->cm+p;                    // set endpoint for rotation
          angle=dp_rot*slp_global.randHalf();
          igroup->rotate(*spc, p, angle);
        }
        if (dp_trans>1e-6) {
          p.x()=dir.x() * dp_trans * slp_global.randHalf();
          p.y()=dir.y() * dp_trans * slp_global.randHalf();
          p.z()=dir.z() * dp_trans * slp_global.randHalf();
          igroup->translate(*spc, p);
        }
      }

    template<class Tspace>
      void TranslateRotate<Tspace>::_acceptMove() {
        double r2 = spc->geo.sqdist( igroup->cm, igroup->cm_trial );
        sqrmap_t[ igroup->name ] += r2;
        sqrmap_r[ igroup->name ] += pow(angle*180/pc::pi, 2);
        accmap[ igroup->name ] += 1;
        igroup->accept(*spc);
      }

    template<class Tspace>
      void TranslateRotate<Tspace>::_rejectMove() {
        sqrmap_t[ igroup->name ] += 0;
        sqrmap_r[ igroup->name ] += 0;
        accmap[ igroup->name ] += 0;
        igroup->undo(*spc);
      }

    template<class Tspace>
      double TranslateRotate<Tspace>::_energyChange() {
        if (dp_rot<1e-6 && dp_trans<1e-6)
          return 0;

        for (auto i : *igroup)
          if ( spc->geo.collision( spc->trial[i], Geometry::Geometrybase::BOUNDARY ) )
            return pc::infty;

        double unew = pot->external(spc->trial) +pot->g_external(spc->trial, *igroup);
        if (unew==pc::infty)
          return pc::infty;       // early rejection
        double uold = pot->external(spc->p) + pot->g_external(spc->p, *igroup);

#ifdef ENABLE_MPI
        if (mpi!=nullptr) {
          double du=0;
          auto s = Faunus::MPI::splitEven(*mpi, spc->groupList().size());
          for (auto i=s.first; i<=s.second; ++i) {
            auto gi=spc->groupList()[i];
            if (gi!=igroup)
              du += pot->g2g(spc->trial, *gi, *igroup) - pot->g2g(spc->p, *gi, *igroup);
          }
          return (unew-uold) + Faunus::MPI::reduceDouble(*mpi, du);
        }
#endif

        for (auto g : spc->groupList()) {
          if (g!=igroup) {
            unew += pot->g2g(spc->trial, *g, *igroup);
            if (unew==pc::infty)
              return pc::infty;   // early rejection
            uold += pot->g2g(spc->p, *g, *igroup);
          }
        }
        return unew-uold;
      }

    template<class Tspace>
      string TranslateRotate<Tspace>::_info() {
        using namespace textio;
        std::ostringstream o;
        o << pad(SUB,w,"Max. translation") << pm << dp_trans/2 << textio::_angstrom << endl
          << pad(SUB,w,"Max. rotation") << pm << dp_rot/2*180/pc::pi << textio::degrees << endl;
        if ( !directions.empty() ) {
          o << indent(SUB) << "Group Move directions:" << endl;
          for (auto &m : directions)
            o << pad(SUBSUB,w-2,m.first)
              << m.second.transpose() << endl;
        }
        if (cnt>0) {
          char l=12;
          o << indent(SUB) << "Move Statistics:" << endl
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

    template<class Tspace>
      void TranslateRotate<Tspace>::_test(UnitTest &t) {
        for (auto m : accmap) {
          string id=m.first,
                 idtrim="_"+textio::trim(id)+"_";
          t(prefix+idtrim+"acceptance", accmap[id].avg()*100);
          t(prefix+idtrim+"dRot", sqrt(sqrmap_r[id].avg()));
          t(prefix+idtrim+"dTrans", sqrt(sqrmap_t[id].avg()));
        }
      }

    /**
      @brief Translates/rotates many groups simultaneously
      */
    template<class Tspace>
      class TranslateRotateNbody : public TranslateRotate<Tspace> {
        private:
          typedef TranslateRotate<Tspace> base;
          typedef opair<Group*> Tpair;
          std::vector<Tpair> pairlist; // interacting groups

          typename base::map_type angle2; //!< Temporary storage for angular movement
          vector<Group*> gVec;   //!< Vector of groups to move

          void _trialMove() FOVERRIDE {
            angle2.clear();
            for (auto g : gVec) {
              if (g->isMolecular()) {
                Point p;
                if (base::dp_rot>1e-6) {
                  p.ranunit(slp_global);        // random unit vector
                  p=g->cm+p;                    // set endpoint for rotation
                  double angle=base::dp_rot*slp_global.randHalf();
                  g->rotate(*base::spc, p, angle);
                  angle2[g->name] += pow(angle*180/pc::pi, 2); // sum angular movement^2
                }
                if (base::dp_trans>1e-6) {
                  p.ranunit(slp_global);
                  p=base::dp_trans*p.cwiseProduct(base::dir);
                  g->translate(*base::spc, p);
                }
              }
            }
          }

          void _acceptMove() FOVERRIDE {
            std::map<string,double> r2;
            for (auto g : gVec) {
              r2[g->name] += base::spc->geo.sqdist(g->cm, g->cm_trial);
              g->accept(*base::spc);
              base::accmap[g->name] += 1;
            }
            for (auto &i : r2)
              base::sqrmap_t[i.first] += i.second;
            for (auto &i : angle2)
              base::sqrmap_r[i.first] += i.second;
          }

          void _rejectMove() FOVERRIDE {
            std::set<string> names; // unique group names
            for (auto g : gVec) {
              names.insert(g->name);
              g->undo(*base::spc);
              base::accmap[g->name]+=0;
            }
            for (auto n : names) {
              base::sqrmap_t[n]+=0;
              base::sqrmap_r[n]+=0;
            }
          }

          string _info() {
            std::ostringstream o;
            o << textio::pad(textio::SUB,base::w,"Number of groups") << gVec.size() << endl
              << base::_info();
            return o.str();
          }

          double _energyChange() FOVERRIDE {
            int N=(int)base::spc->groupList().size();
            double du=0;

#ifdef ENABLE_MPI
            if (!pairlist.empty()) {

              // group <-> group
              auto s = Faunus::MPI::splitEven(*mpi, pairlist.size());
              for (size_t i=s.first; i<=s.second; ++i)
                du += base::pot->g2g(base::spc->trial,*pairlist[i].first,*pairlist[i].second)
                  - base::pot->g2g(base::spc->p,*pairlist[i].first,*pairlist[i].second);

              // group <-> external potential
              s = Faunus::MPI::splitEven(*mpi, base::spc->groupList().size());
              for (size_t i=s.first; i<=s.second; ++i) {
                auto gi=base::spc->groupList()[i];
                du += base::pot->g_external(base::spc->trial, *gi) - base::pot->g_external(base::spc->p, *gi);
              }

              return Faunus::MPI::reduceDouble(*mpi, du);
            }
#endif

            if (!pairlist.empty()) {
#pragma omp parallel for reduction (+:du)
              for (int i=0; i<(int)pairlist.size(); i++)
                du+=base::pot->g2g(base::spc->trial,*pairlist[i].first,*pairlist[i].second)
                  - base::pot->g2g(base::spc->p,*pairlist[i].first,*pairlist[i].second);
              for (auto g : base::spc->groupList())
                du += base::pot->g_external(base::spc->trial, *g) - base::pot->g_external(base::spc->p, *g);
              return du;
            }

            if (!gVec.empty()) {
#pragma omp parallel for reduction (+:du) schedule (dynamic)
              for (int i=0; i<N-1; i++)
                for (int j=i+1; j<N; j++) {
                  auto gi=base::spc->groupList()[i];
                  auto gj=base::spc->groupList()[j];
                  du += base::pot->g2g(base::spc->trial,*gi,*gj) - base::pot->g2g(base::spc->p,*gi,*gj);
                }
              for (auto g : base::spc->groupList())
                du += base::pot->g_external(base::spc->trial, *g) - base::pot->g_external(base::spc->p, *g);
            }
            return du;
          }
          void setGroup(std::vector<Group*> &v) {
            gVec.clear();
            pairlist.clear();
            std::set<Tpair> l; // set of all unique g2g pairs
            for (auto i : v) {
              if (i->isMolecular())
                gVec.push_back(i);
              for (auto j : v)
                l.insert(Tpair(i,j));
            }
            for (auto i : l)
              pairlist.push_back(i); // set->vector (faster)

            //auto N=gVec.size();
            //base::dp_trans /= N;
            //base::dp_rot /= N;
          }

        public:
          TranslateRotateNbody(InputMap &in,Energy::Energybase<Tspace> &e, Tspace &s, string pfx="transrot") : base(in,e,s,pfx) {
            base::title+=" (N-body)";
            setGroup(s.groupList());
          }
#ifdef ENABLE_MPI
          Faunus::MPI::MPIController* mpi;
#endif
      };

    /**
     * @brief Combined rotation and rotation of groups and mobile species around it
     *
     * This class will do a combined translational and rotational move of a group
     * along with atomic particles surrounding it.
     * To specify where to look for clustered particles, use the `setMobile()`
     * function. Whether particles are considered part of the cluster is
     * determined by the private virtual function `ClusterProbability()`.
     * By default this is a simple step function with P=1 when an atomic particle
     * in the group set by `setMobile()` is within a certain threshold to a
     * particle in the main group; P=0 otherwise.
     *
     * The implemented cluster algorithm is general - see Frenkel&Smith,
     * 2nd ed, p405 - and derived classes can re-implement `ClusterProbability()`
     * for arbitrary probability functions.
     *
     * Upon construction, the `InputMap` is scanned for the
     * following keywords,
     *
     * Keyword                | Description
     * :--------------------- | :----------------
     * `transrot_clustersize` | Surface threshold from mobile ion to particle in group (angstrom)
     *
     * @todo Energy evaluation puts all moved particles in an index vector used
     * to sum the interaction energy with static particles. This could be optimized
     * by only adding mobile ions and calculate i.e. group-group interactions in a
     * more traditional (and faster) way.
     */
    template<class Tspace>
      class TranslateRotateCluster : public TranslateRotate<Tspace> {
        protected:
          typedef TranslateRotate<Tspace> base;
          using base::pot;
          using base::w;
          using base::cnt;
          using base::igroup;
          using base::dp_trans;
          using base::dp_rot;
          using base::dir;
          Geometry::QuaternionRotate vrot;
          vector<int> cindex; //!< index of mobile ions to move with group
          void _trialMove();
          void _acceptMove();
          void _rejectMove();
          double _energyChange();
          string _info();
          Average<double> avgsize; //!< Average number of ions in cluster
          Average<double> avgbias; //!< Average bias
          Group* gmobile;          //!< Pointer to group with potential cluster particles
          virtual double ClusterProbability(p_vec&,int); //!< Probability that particle index belongs to cluster
        public:
          using base::spc;
          TranslateRotateCluster(InputMap&, Energy::Energybase<Tspace>&, Tspace&, string="transrot");
          virtual ~TranslateRotateCluster();
          void setMobile(Group&);  //!< Pool of atomic species to move with the main group
          double threshold;        //!< Distance between particles to define a cluster
      };

    template<class Tspace>
      TranslateRotateCluster<Tspace>::TranslateRotateCluster(InputMap &in,Energy::Energybase<Tspace> &e, Tspace &s, string pfx) : base(in,e,s,pfx) {
        base::title="Cluster "+base::title;
        base::cite="doi:10/cj9gnn";
        threshold = in.get<double>(pfx+"_clustersize",0);
        gmobile=nullptr;
      }

    template<class Tspace>
      TranslateRotateCluster<Tspace>::~TranslateRotateCluster() {}

    template<class Tspace>
      void TranslateRotateCluster<Tspace>::setMobile(Group &g) {
        assert(&g!=nullptr);
        gmobile=&g;
      }

    template<class Tspace>
      string TranslateRotateCluster<Tspace>::_info() {
        using namespace textio;
        std::ostringstream o;
        o << base::_info() << endl;
        o << pad(SUB,w,"Cluster threshold") << threshold << _angstrom << endl;
        if (cnt>0) {
          o << pad(SUB,w,"Average cluster size") << avgsize.avg() << endl;
          if (threshold>1e-9)
            o << pad(SUB,w,"Average bias") << avgbias.avg() << " (0=reject, 1=accept)\n";
        }
        return o.str();
      }

    template<class Tspace>
      void TranslateRotateCluster<Tspace>::_trialMove() {
        assert(gmobile!=nullptr && "Cluster group not defined");
        assert(igroup!=nullptr && "Group to move not defined");
        Point p;

        // find clustered particles
        cindex.clear();
        for (auto i : *gmobile)
          if (ClusterProbability(spc->p, i) > slp_global() )
            cindex.push_back(i); // generate cluster list
        avgsize += cindex.size();

        // rotation
        if (dp_rot>1e-6) {
          base::angle=dp_rot*slp_global.randHalf();
          p.ranunit(slp_global);
          p=igroup->cm+p; // set endpoint for rotation
          igroup->rotate(*spc, p, base::angle);
          vrot.setAxis(spc->geo, igroup->cm, p, base::angle); // rot. around line CM->p
          for (auto i : cindex)
            spc->trial[i] = vrot(spc->trial[i]); // rotate
        }

        // translation
        if (dp_trans>1e-6) {
          p.x()=dir.x() * dp_trans * slp_global.randHalf();
          p.y()=dir.y() * dp_trans * slp_global.randHalf();
          p.z()=dir.z() * dp_trans * slp_global.randHalf();
          igroup->translate(*spc, p);
          for (auto i : cindex)
            spc->trial[i].translate(spc->geo,p);
        }
      }

    template<class Tspace>
      void TranslateRotateCluster<Tspace>::_acceptMove() {
        base::_acceptMove();
        for (auto i : cindex)
          spc->p[i] = spc->trial[i];
      }

    template<class Tspace>
      void TranslateRotateCluster<Tspace>::_rejectMove() {
        base::_rejectMove();
        for (auto i : cindex)
          spc->trial[i] = spc->p[i];
      }

    template<class Tspace>
      double TranslateRotateCluster<Tspace>::_energyChange() {
        double bias=1;             // cluster bias -- see Frenkel 2nd ed, p.405
        vector<int> imoved=cindex; // index of moved particles
        for (auto l : *gmobile)    // mobile index, "l", NOT in cluster (Frenkel's "k" is the main group)
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
          if ( spc->geo.collision( spc->trial[i], Geometry::Geometrybase::BOUNDARY ) )
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
        // note: this could be optimized!
        double du=0;
#pragma omp parallel for reduction (+:du)
        for (int j=0; j<(int)spc->p.size(); j++)
          if ( std::find(imoved.begin(), imoved.end(), j )==imoved.end() )
            for (auto i : imoved)
              du += pot->i2i(spc->trial, i, j) - pot->i2i(spc->p, i, j);
        return unew - uold + du - log(bias); // exp[ -( dU-log(bias) ) ] = exp(-dU)*bias
      }

    /**
     * This is the default function for determining the probability, P,
     * that a mobile particle is considered part of the cluster. This
     * is here a simple distance critera but derived classes can reimplement
     * this (virtual) function to arbitrary probability functions.
     */
    template<class Tspace>
      double TranslateRotateCluster<Tspace>::ClusterProbability(p_vec &p, int i) {
        for (auto j : *igroup) // loop over main group
          if (i!=j) {
            double r=threshold+p[i].radius+p[j].radius;
            if (spc->geo.sqdist(p[i],p[j])<r*r )
              return 1;
          }
        return 0;
      }

    /**
     * @brief Rotate/translate group along with an extra group
     *
     * This will rotate/translate a group A around its mass center and, if
     * defined, also an extra group, B. This can be useful for sampling groups
     * joined together with springs, for example a polymer (B) joined to a
     * protein (A). The group B can consist of many molecules/chains as
     * long as these are continuous in the particle vector.
     *
     * @date Malmo 2014
     */
    template<class Tspace>
      class TranslateRotateGroupCluster : public TranslateRotateCluster<Tspace> {
        private:
          typedef TranslateRotateCluster<Tspace> base;
          void _acceptMove() {
            base::_acceptMove();
            for (auto i : base::spc->groupList())
              i->setMassCenter(*base::spc);
          }
          string _info() FOVERRIDE { return base::base::_info(); }
          double ClusterProbability(p_vec &p, int i) FOVERRIDE { return 1; }
        public:
          TranslateRotateGroupCluster(InputMap &in, Energy::Energybase<Tspace> &e,
              Tspace &s, string pfx="transrot") : base(in,e,s,pfx) {
            base::title = "Translate-Rotate w. extra group";
          }
      };

    /**
     * @brief Non-rejective cluster translation.
     *
     * This type of move will attempt to translate collective sets of macromolecules that
     * with a symmetric transition matrix (no flow through the clusters).
     * See detailed description [here](http://dx.doi.org/10/fthw8k).
     *
     * Setting the boolean `skipEnergyUpdate` to true (default is false) updates of the
     * total energy are skipped to speed up the move.
     * While this has no influence on the Markov chain it will cause an apparent energy
     * drift. It is recommended that this is enabled only for long production runs after
     * having properly checked that no drifts occur with `skipEnergyUpdate=false`.
     *
     * Upon construction the following keywords are read from `InputMap`,
     *
     * Keyword                | Description
     * :--------------------  | :-------------------------------
     * `ctransnr_dp`          | Displacement parameter
     * `ctransnr_skipenergy`  | Skip energy update, see above (default: false)
     * `ctransnr_runfraction` | Runfraction (default: 1.0)
     *
     * @note Requirements for usage:
     * - Compatible only with purely molecular systems
     * - Works only with periodic containers
     * - External potentials are ignored
     *
     * @author Bjoern Persson
     * @date Lund 2009-2010
     * @todo Energy calc. before and after can be optimized by only looping over `moved` with
     * `remaining`
     */
    template<class Tspace>
      class ClusterTranslateNR : public Movebase<Tspace> {
        private:
          typedef Movebase<Tspace> base;
          using base::w;
          using base::spc;
          using base::pot;
          using base::prefix;
          vector<int> moved, remaining;
          void _trialMove();
          void _acceptMove();
          void _rejectMove();
          double _energyChange();
          string _info();
          Average<double> movefrac; //!< Fraction of particles moved
          double dp;                //!< Displacement parameter [aa]
          vector<Group*> g;         //!< Group of molecules to move. Currently this needs to be ALL groups in the system!!
        public:
          ClusterTranslateNR(InputMap&, Energy::Energybase<Tspace>&, Tspace&, string="ctransnr");
          bool skipEnergyUpdate;    //!< True if energy updates should be skipped (faster evaluation!)
      };

    /** @brief Constructor */
    template<class Tspace>
      ClusterTranslateNR<Tspace>::ClusterTranslateNR(InputMap &in, Energy::Energybase<Tspace> &e, Tspace &s, string pfx) : base(e,s,pfx) {
        base::title="Rejection Free Cluster Translation";
        base::cite="doi:10/fthw8k";
        base::useAlternateReturnEnergy=true;
        base::runfraction=in.get<double>(prefix+"_runfraction", 1.0);
        skipEnergyUpdate=in.get<bool>(prefix+"_skipenergy", false);
        dp=in.get<double>(prefix+"_dp", 0);
        if (dp<1e-6)
          base::runfraction=0;
        g=spc->groupList(); // currently ALL groups in the system will be moved!
      }

    template<class Tspace>
      string ClusterTranslateNR<Tspace>::_info() {
        using namespace textio;
        std::ostringstream o;
        o << pad(SUB,w,"Displacement") << dp << _angstrom << endl
          << pad(SUB,w,"Skip energy update") << std::boolalpha
          << skipEnergyUpdate << endl;
        if (movefrac.cnt>0) {
          o << pad(SUB,w,"Move fraction") << movefrac.avg()*100 << percent << endl
            << pad(SUB,w,"Avg. moved groups") << movefrac.avg()*spc->groupList().size() << endl;
        }
        return o.str();
      }

    template<class Tspace>
      void ClusterTranslateNR<Tspace>::_trialMove() {
        double du=0;
        g=spc->groupList();
        moved.clear();
        remaining.resize( g.size() );

        for (size_t i=0; i<g.size(); i++) {
          remaining[i] = i;
          if (base::cnt<=1)
            g[i]->setMassCenter(*spc);
        }

        if (skipEnergyUpdate==false)
#pragma omp parallel for reduction (+:du) schedule (dynamic)
          for (size_t i=0; i<g.size()-1; i++)
            for (size_t j=i+1; j<g.size(); j++)
              du-=pot->g2g(spc->p, *g[i], *g[j]);

        Point ip(dp,dp,dp);
        ip.x()*=slp_global.randHalf();
        ip.y()*=slp_global.randHalf();
        ip.z()*=slp_global.randHalf();

        int f=slp_global()*remaining.size();
        moved.push_back(remaining[f]);
        remaining.erase(remaining.begin()+f);    // Pick first index in m to move

        for (size_t i=0; i<moved.size(); i++) {
          g[moved[i]]->translate(*spc, ip);
          for (size_t j=0; j<remaining.size(); j++) {
            double uo=pot->g2g(spc->p,     *g[moved[i]], *g[remaining[j]]);
            double un=pot->g2g(spc->trial, *g[moved[i]], *g[remaining[j]]);
            double udiff=un-uo;
            if (slp_global() < (1.-std::exp(-udiff)) ) {
              moved.push_back(remaining[j]);
              remaining.erase(remaining.begin()+j);
              j--;
            }
          }
          g[moved[i]]->accept(*spc);
        }

        if (skipEnergyUpdate==false)
#pragma omp parallel for reduction (+:du) schedule (dynamic)
          for (size_t i=0; i<g.size()-1; i++)
            for (size_t j=i+1; j<g.size(); j++)
              du+=pot->g2g(spc->p, *g[i], *g[j]);

        base::alternateReturnEnergy=du;
        movefrac+=double(moved.size()) / (moved.size()+remaining.size());

        assert( moved.size() >= 1);
        assert( spc->groupList().size() == moved.size()+remaining.size() );
      }

    template<class Tspace>
      double ClusterTranslateNR<Tspace>::_energyChange() { return 0; }

    template<class Tspace>
      void ClusterTranslateNR<Tspace>::_acceptMove() {}

    template<class Tspace>
      void ClusterTranslateNR<Tspace>::_rejectMove() {}

    /**
     * @brief Crank shaft move of linear polymers
     *
     * This will perform a crank shaft move of a linear polymer molecule.
     * Two monomers are picked at random and a rotation axis is drawn
     * between them. The particles in between are rotated around that
     * axis. By setting `minlen` and `maxlen` one can control the maximum
     * number particles to rotate. For example, for a crankshaft
     * move spanning only one bond, set `minlen=maxlen=1`.
     * The behavoir for branched molecules is currently undefined.
     *
     * ![Figure: Various polymer moves.](polymerdisplacements.jpg)
     *
     * @date Lund 2012
     */
    template<class Tspace>
      class CrankShaft : public Movebase<Tspace> {
        private:
          typedef Movebase<Tspace> base;
          void _test(UnitTest&);
          void _trialMove();
          void _acceptMove();
          void _rejectMove();
          double _energyChange();
          string _info();
          virtual bool findParticles(); //!< This will set the end points and find particles to rotate
        protected:
          using base::spc;
          using base::pot;
          using base::w;
          using base::prefix;
          Group* gPtr;       //!< Pointer to group where move is to be performed. Set by setGroup().
          double dp;         //!< Rotational displacement parameter
          double angle;      //!< Current rotation angle
          vector<int> index; //!< Index of particles to rotate
          //Geometry::VectorRotate vrot;
          Geometry::QuaternionRotate vrot;
          AcceptanceMap<string> accmap;
        public:
          CrankShaft(InputMap&, Energy::Energybase<Tspace>&, Tspace&, string="crank");
          virtual ~CrankShaft();
          void setGroup(Group&); //!< Select Group to of the polymer to move
          int minlen;            //!< Minimum number of particles to rotate (default = 1)
          int maxlen;            //!< Maximin number of particles to rotate (default = 10)
      };

    /**
     * The InputMap is searched for:
     *
     * Key            | Description
     * -------------- | --------------------------------------
     * `crank_minlen` | Minimum number of particles to rotate
     * `crank_maxlen` | Maximum number of particles to rotate
     * `crank_dp`     | Rotational displacement parameter (radians)
     */
    template<class Tspace>
      CrankShaft<Tspace>::CrankShaft(InputMap &in, Energy::Energybase<Tspace> &e, Tspace &s, string pfx) : base(e,s,pfx) {
        base::title="CrankShaft";
        w=30;
        gPtr=nullptr;
        minlen = in.get<int>(prefix+"_minlen", 1, "Min. particles to rotate");
        maxlen = in.get<int>(prefix+"_maxlen", 4, "Max. particles to rotate");
        assert(minlen<=maxlen);
        dp=in.get<double>(prefix+"_dp", 3.);
        base::runfraction = in.get<double>(prefix+"_runfraction",1.);
        if (dp<1e-6)
          base::runfraction=0;
      }

    template<class Tspace>
      CrankShaft<Tspace>::~CrankShaft() {}

    template<class Tspace>
      void CrankShaft<Tspace>::_trialMove() {
        assert(gPtr!=nullptr && "No group to perform crankshaft on.");
        if (gPtr->size()<3)
          return;
        index.clear();   // clear previous particle list to rotate
        findParticles();
        assert(!index.empty() && "No particles to rotate.");
        for (auto i : index)
          spc->trial[i] = vrot(spc->p[i]); // (boundaries are accounted for)
        gPtr->cm_trial = Geometry::massCenter(spc->geo, spc->trial, *gPtr);
      }

    template<class Tspace>
      void CrankShaft<Tspace>::_acceptMove() {
        double msq=0;
        for (auto i : index) {
          msq+=spc->geo.sqdist( spc->p[i], spc->trial[i] );
          spc->p[i] = spc->trial[i];
        }
        accmap.accept(gPtr->name, msq ) ;
        gPtr->cm = gPtr->cm_trial;
      }

    template<class Tspace>
      void CrankShaft<Tspace>::_rejectMove() {
        accmap.reject(gPtr->name);
        for (auto i : index)
          spc->trial[i] = spc->p[i];
        gPtr->cm_trial = gPtr->cm;
      }

    /**
     * @todo g_internal is not really needed - index<->g would be faster
     */
    template<class Tspace>
      double CrankShaft<Tspace>::_energyChange() {
        double du=0;
        for (auto i : index)
          if ( spc->geo.collision( spc->trial[i], Geometry::Geometrybase::BOUNDARY ) )
            return pc::infty;
        du=pot->g_internal(spc->trial, *gPtr) - pot->g_internal(spc->p, *gPtr);
        for (auto i : index)
          du+=pot->i_external(spc->trial,i) - pot->i_external(spc->p,i);
        for (auto g : spc->groupList())
          if (g!=gPtr)
            du+=pot->g2g(spc->trial, *g, *gPtr) - pot->g2g(spc->p, *g, *gPtr);
        du+=pot->external(spc->trial) - pot->external(spc->p);
        //for (auto i : index)
        //  du += pot->i2all(spc->trial, i) - pot->i2all(spc->p, i);
        return du;
      }

    /**
     * This will define the particles to be rotated (stored in index vector) and
     * also set the axis to rotate around, defined by two points.
     */
    template<class Tspace>
      bool CrankShaft<Tspace>::findParticles() {
        assert( minlen <= gPtr->size()-2 && "Minlen too big for molecule!");

        int beg,end,len;
        do {
          beg=gPtr->random();             // generate random vector to
          end=gPtr->random();             // rotate around
          len = std::abs(beg-end) - 1;    // number of particles between end points
        } while ( len<minlen || len>maxlen  );

        angle = dp*slp_global.randHalf();  // random angle
        vrot.setAxis(spc->geo, spc->p[beg], spc->p[end], angle );

        index.clear();
        if (beg>end)
          std::swap(beg,end);
        for (int i=beg+1; i<end; i++)
          index.push_back(i);             // store particle index to rotate
        assert(index.size()==size_t(len));

        return true;
      }

    template<class Tspace>
      void CrankShaft<Tspace>::setGroup(Group &g) { gPtr=&g; }

    template<class Tspace>
      string CrankShaft<Tspace>::_info() {
        using namespace textio;
        std::ostringstream o;
        o << pad(SUB,w, "Displacement parameter") << dp << endl
          << pad(SUB,w, "Min/max length to move") << minlen << " " << maxlen << endl;
        if (base::cnt>0)
          o << accmap.info();
        return o.str();
      }

    template<class Tspace>
      void CrankShaft<Tspace>::_test(UnitTest &t) {
        accmap._test(t, prefix);
      }

    /**
     * @brief Pivot move for linear polymers
     *
     * This will perform a pivot rotation of a linear polymer by the following steps:
     *
     * - Select rotation axis by two random monomers, spanning `minlen` to `maxlen` bonds
     * - Rotate monomers before or after end points of the above axis
     *
     * ![Figure: Various polymer moves.](polymerdisplacements.jpg)
     *
     * @date Asljunga 2012
     */
    template<class Tspace>
      class Pivot : public CrankShaft<Tspace> {
        protected:
          typedef CrankShaft<Tspace> base;
          using base::index;
          using base::gPtr;
          using base::spc;
          bool findParticles();
        public:
          Pivot(InputMap&, Energy::Energybase<Tspace>&, Tspace&, string="pivot");
      };

    template<class Tspace>
      Pivot<Tspace>::Pivot(InputMap &in, Energy::Energybase<Tspace> &e, Tspace &s, string pfx) : base(in,e,s,pfx) {
        base::title="Polymer Pivot Move";
        base::minlen=1; // minimum bond length to rotate around
      }

    template<class Tspace>
      bool Pivot<Tspace>::findParticles() {
        int beg(0),end(0),len;
        index.clear();
        while (index.empty()) {
          do {
            beg = gPtr->random(); // define the
            end = gPtr->random(); // axis to rotate around
            len = std::abs(beg-end);
          } while ( len<base::minlen || len>base::maxlen );

          if (slp_global.randHalf() > 0)
            for (int i=end+1; i<=gPtr->back(); i++)
              index.push_back(i);
          else
            for (int i=gPtr->front(); i<end; i++)
              index.push_back(i);
        }
        base::angle = base::dp*slp_global.randHalf();
        base::vrot.setAxis(spc->geo, spc->p[beg], spc->p[end], base::angle );
        return true;
      }

    /**
     * @brief Reptation move for linear polymers
     *
     * This will perform a reptation move of a linear, non-uniform polymer chain.
     * During construction, the InputMap is searched for the following keywords:
     *
     * Key                     | Description
     * :---------------------- | :---------------------------
     * `reptation_runfraction` | Probability to perform a move (defaults=1)
     * `reptation_bondlength`  | The bond length while moving head groups. Use -1 to use existing bondlength.
     *
     * @date Lund 2012
     */
    template<class Tspace>
      class Reptation : public Movebase<Tspace> {
        private:
          typedef Movebase<Tspace> base;
          AcceptanceMap<string> accmap;
          void _test(UnitTest&);
          void _trialMove();
          void _acceptMove();
          void _rejectMove();
          double _energyChange();
          string _info();
          Group* gPtr;
          double bondlength; //!< Reptation length used when generating new head group position
        protected:
          using base::pot;
          using base::spc;
          using base::prefix;
        public:
          Reptation(InputMap&, Energy::Energybase<Tspace>&, Tspace&, string="reptation");
          void setGroup(Group&); //!< Select Group to move
      };

    template<class Tspace>
      Reptation<Tspace>::Reptation(InputMap &in, Energy::Energybase<Tspace> &e, Tspace &s, string pfx) : base(e,s,pfx) {
        base::title="Linear Polymer Reptation";
        base::runfraction = in.get<double>(prefix+"_runfraction",1.0);
        bondlength = in.get<double>(prefix+"_bondlength", -1);
        gPtr=nullptr;
      }

    template<class Tspace>
      void Reptation<Tspace>::setGroup(Group &g) { gPtr=&g; }

    template<class Tspace>
      void Reptation<Tspace>::_test(UnitTest &t) {
        accmap._test(t, prefix);
      }

    template<class Tspace>
      void Reptation<Tspace>::_trialMove() {
        assert(gPtr!=nullptr && "Did you forget to call setGroup?");
        if (gPtr->size()<2)
          return;

        int first, second; // "first" is end point, "second" is the neighbor
        if (slp_global.randHalf()>0) {
          first=gPtr->front();
          second=first+1;
        } else {
          first=gPtr->back();
          second=first-1;
        }

        double bond;
        if (bondlength>0)
          bond=bondlength;
        else
          bond=spc->geo.dist(spc->p[first], spc->p[second]); // bond length of first or last particle

        // shift particles up or down
        for (int i=gPtr->front(); i<gPtr->back(); i++)
          if (first<second)
            spc->trial[i+1]=Point( spc->p[i] );
          else
            spc->trial[i]=Point( spc->p[i+1] );

        // generate new position for end point ("first")
        Point u;
        u.ranunit(slp_global);                          // generate random unit vector
        spc->trial[first].translate(spc->geo, u*bond); // trans. 1st w. scaled unit vector
        assert( std::abs( spc->geo.dist(spc->p[first],spc->trial[first])-bond ) < 1e-7  );

        for (auto i : *gPtr)
          spc->geo.boundary( spc->trial[i] );  // respect boundary conditions

        gPtr->cm_trial = Geometry::massCenter(spc->geo, spc->trial, *gPtr);
      }

    template<class Tspace>
      void Reptation<Tspace>::_acceptMove() {
        accmap.accept(gPtr->name, spc->geo.sqdist(gPtr->cm, gPtr->cm_trial) );
        gPtr->accept(*spc);
      }

    template<class Tspace>
      void Reptation<Tspace>::_rejectMove() {
        accmap.reject(gPtr->name);
        gPtr->undo(*spc);
      }

    template<class Tspace>
      double Reptation<Tspace>::_energyChange() {
        for (auto i : *gPtr)
          if ( spc->geo.collision( spc->trial[i], Geometry::Geometrybase::BOUNDARY ) )
            return pc::infty;

        double unew = pot->g_external(spc->trial, *gPtr) + pot->g_internal(spc->trial,*gPtr);
        if (unew==pc::infty)
          return pc::infty;       // early rejection
        double uold = pot->g_external(spc->p, *gPtr) + pot->g_internal(spc->p, *gPtr);

        for (auto g : spc->groupList()) {
          if (g!=gPtr) {
            unew += pot->g2g(spc->trial, *g, *gPtr);
            if (unew==pc::infty)
              return pc::infty;   // early rejection
            uold += pot->g2g(spc->p, *g, *gPtr);
          }
        }
        return unew-uold;
      }

    template<class Tspace>
      string Reptation<Tspace>::_info() {
        using namespace textio;
        std::ostringstream o;
        o << pad(SUB,base::w, "Bondlength") << bondlength << _angstrom + " (-1 = automatic)\n";
        if (base::cnt>0)
          o << accmap.info();
        return o.str();
      }

    /** 
     * @brief Grand canonical move for molecular species
     * @warning Under construction!
     * @date Lund, 2014
     */
    template<class Tspace>
    class GCMolecular : public Movebase<Tspace> {
    private:
      typedef Movebase<Tspace> base;
      typedef std::map<short, Average<double> > map_type;
      bool use2D;
      Point dir;
      int inmol; // molecule to be inserted
      int outmol; // molecule to be deleted
      int imol;
      int idel;
      
      // structure for molecular data
      struct MolData {
        string name;
        double activity; // in 1/angstrom^3 (3D) or 1/angstrom^2 (2D)
        typename Tspace::ParticleVector p;
      };
      
      std::vector<MolData> molvec;
      
      template<typename pvec>
      void randomRotate(pvec &p) {
        Geometry::QuaternionRotate rot;
      }
      
      template<typename pvec>
      void randomState(pvec &p) {
        Point a;
        base::spc->randomPos(a);
        a=a*dir+offset;
        auto dp = a-Geometry::massCenter(base::spc->geo, p);
        for (auto i : &p) {
          i=i+dp;
          base::spc->geo.boundary(i);
        }
      }
      
      int randomMol() {
        return slp_global.rand() % molvec.size();
      }
      
      bool insert() {
        inmol = randomMol();
        randomState( molvec[inmol] );
        return true;
      }
      
    protected:
      bool deleteBool;
      string _info();
      void _acceptMove() FOVERRIDE;
      void _rejectMove() FOVERRIDE;
      
      double _energyChange() FOVERRIDE {
        double uold=0;//, unew=0;
        
        // count molecules
        int N=0;
        for (auto g : base::spc->groupList())
          if (g->isMolecular())
            if ( g->name==molvec[imol] )
              N++;
        
        // calc. volume or area
        double V = base::spc->geo.getVolume();
        if (use2D)
          V = V / base::spc->geo.len.z();
        
        // in case of deletion
        if (deleteBool) {
          auto gi = base::spc->findGroup(idel);
          for (auto gj : base::spc->groupList())
            if (gi!=gj)
              uold += base::pot->g2g( base::spc->p, *gi, *gj )
                + base::pot->g_external( base::spc->p, *gi );
        } else {
          // in case of insertion
          
        }

        return 0;
          
      }
      
      void _trialMove() {
        
        // insert
        if (slp_global()>0.5) {
          inmol = randomMol();
          outmol=-1;
          randomState( molvec[inmol] );
        }
        // delete
        else {
          inmol=-1;
          outmol=randomMol();
          //for (auto i : base::spc->groupList())
          //  if ()
          
        }
        
      }
      using base::spc;
      map_type accmap; //!< acceptance map
      
    public:
      Point offset;
      GCMolecular(InputMap &in, Energy::Energybase<Tspace> &e, Tspace &s, string pfx="gcmol_") : base(e,s,pfx) {
        dir=Point(1,1,1);
        offset=Point(0,0,0);
        use2D=in(pfx+"2d", true);
        if (use2D)
          dir.z()=1;
      }
    };
    
    /**
     * @brief Isobaric volume move
     *
     * @details This class will perform a volume displacement and scale atomic
     * as well as molecular groups as long as these are known to Space -
     * see Space.enroll().
     * The InputMap class is scanned for the following keys:
     *
     * Key              | Description
     * :--------------- | :-----------------------------
     * `npt_dV`         | Volume displacement parameter
     * `npt_P`          | Pressure [mM]
     * `npt_runfraction`| Runfraction [default=1]
     *
     * Note that new volumes are generated according to
     * \f$ V^{\prime} = \exp\left ( \log V \pm \delta dV \right ) \f$
     * where \f$\delta\f$ is a random number between zero and one half.
     */
    template<class Tspace>
      class Isobaric : public Movebase<Tspace> {
        private:
          typedef Movebase<Tspace> base;
          using base::spc;
          using base::pot;
          using base::w;
          string _info();
          void _test(UnitTest&);
          void _trialMove();
          void _acceptMove();
          void _rejectMove();
          template<class Tpvec> double _energy(const Tpvec&);
          double _energyChange();
          double dV; //!< Volume displacement parameter
          double oldV;
          double newV;
          double P; //!< Pressure
          Average<double> sqrV;       //!< Mean squared volume displacement
          Average<double> V;          //!< Average volume
          Average<double> rV;         //!< Average 1/volume
        public:
          template<class Tenergy>
            Isobaric(InputMap&, Tenergy&, Tspace&, string="npt");
      };

    template<class Tspace>
      template<class Tenergy>
      Isobaric<Tspace>::Isobaric(InputMap &in, Tenergy &e, Tspace &s, string pfx) : base(e,s,pfx) {
        this->title="Isobaric Volume Fluctuations";
        this->w=30;
        dV=in.get<double>(pfx+"_dV", 0.,
            "NPT volume displacement parameter");
        P=in.get<double>(pfx+"_P", 0.,
            "NPT external pressure P/kT (mM)")/1e30*pc::Nav; //mM -> 1/A^3
        this->runfraction = in.get<double>(pfx+"_runfraction",1.0);
        if (dV<1e-6)
          base::runfraction=0;
      }

    template<class Tspace>
      string Isobaric<Tspace>::_info() {
        using namespace textio;
        std::ostringstream o;
        const double tomM=1e30/pc::Nav;
        int N,Natom=0, Nmol=0;
        for (auto g : spc->groupList())
          if (g->isAtomic())
            Natom += g->size();
          else
            Nmol+=g->numMolecules();
        N = Natom + Nmol;
        double Pascal = P*pc::kB*pc::T()*1e30;
        o << pad(SUB,w, "Displacement parameter") << dV << endl
          << pad(SUB,w, "Number of molecules")
          << N << " (" <<Nmol<< " molecular + " <<Natom<< " atomic)\n"
          << pad(SUB,w, "Pressure")
          << P*tomM << " mM = " << Pascal << " Pa = "
          << Pascal/0.980665e5 << " atm\n"
          << pad(SUB,w, "Temperature") << pc::T() << " K\n";
        if (base::cnt>0) {
          char l=14;
          o << pad(SUB,w, "Mean displacement")
            << cuberoot+rootof+bracket("dV"+squared)
            << " = " << pow(sqrV.avg(), 1/6.) << _angstrom << endl
            << pad(SUB,w, "Osmotic coefficient") << P / (N*rV.avg()) << endl
            << endl
            << indent(SUBSUB) << std::right << setw(10) << ""
            << setw(l+5) << bracket("V")
            << setw(l+8) << cuberoot+bracket("V")
            << setw(l+8) << bracket("1/V")
            << setw(l+8) << bracket("N/V") << endl
            << indent(SUB) << setw(10) << "Averages"
            << setw(l) << V.avg() << _angstrom + cubed
            << setw(l) << std::cbrt(V.avg()) << _angstrom
            << setw(l) << rV.avg() << " 1/" + _angstrom + cubed
            << setw(l) << N*rV.avg()*tomM << " mM\n";
        }
        return o.str();
      }

    template<class Tspace>
      void Isobaric<Tspace>::_test(UnitTest &t) {
        t(this->prefix+"_averageSideLength", std::cbrt(V.avg()) );
        t(this->prefix+"_MSQDisplacement", pow(sqrV.avg(), 1/6.) );
      }

    template<class Tspace>
      void Isobaric<Tspace>::_trialMove() {
        assert(spc->groupList().size()>0
            && "Space has empty group vector - NPT move not possible.");
        oldV = spc->geo.getVolume();
        newV = std::exp( std::log(oldV) + slp_global.randHalf()*dV );
        for (auto g : spc->groupList()) {
          g->setMassCenter(*spc);
          g->scale(*spc, newV); // scale trial coordinates to new volume
        }
      }

    template<class Tspace>
      void Isobaric<Tspace>::_acceptMove() {
        V += newV;
        sqrV += pow( oldV-newV, 2 );
        rV += 1./newV;
        spc->geo.setVolume(newV);
        pot->setSpace(*spc);
        for (auto g : spc->groupList() )
          g->accept(*spc);
      }

    template<class Tspace>
      void Isobaric<Tspace>::_rejectMove() {
        sqrV += 0;
        V += oldV;
        rV += 1./oldV;
        spc->geo.setVolume(oldV);
        pot->setSpace(*spc);
        for (auto g : spc->groupList() )
          g->undo(*spc);
      }

    /**
     * This will calculate the total energy of the configuration
     * associated with the current Hamiltonian volume
     */
    template<class Tspace>
      template<class Tpvec>
      double Isobaric<Tspace>::_energy(const Tpvec &p) {
        double u=0;
        if (dV<1e-6)
          return u;
        size_t n=spc->groupList().size();  // number of groups
        for (size_t i=0; i<n-1; ++i)      // group-group
          for (size_t j=i+1; j<n; ++j)
            u += pot->g2g(p, *spc->groupList()[i], *spc->groupList()[j]);

        for (auto g : spc->groupList()) {
          u += pot->g_external(p, *g);
          if (g->numMolecules()>1)
            u+=pot->g_internal(p, *g);
        }
        return u + pot->external(p);
      }

    /**
     * @todo Early rejection could be implemented
     *       - not relevant for geometries with periodicity, though.
     */
    template<class Tspace>
      double Isobaric<Tspace>::_energyChange() {
        double uold = _energy(spc->p);
        spc->geo.setVolume(newV);
        pot->setSpace(*spc); // potential must know about volume, too

        // In spherical geometries, molecules may collide with
        // cell boundary upon mass center scaling:
        for (auto g : spc->groupList())
          for (auto i : *g)
            if ( spc->geo.collision( spc->trial[i], Geometry::Geometrybase::BOUNDARY ) )
              return pc::infty;
        double unew = _energy(spc->trial);
        return unew-uold;
      }


    /**
     * @brief Auxillary class for tracking atomic species
     * @date Malmo 2011
     *
     * This class keeps track of individual particle positions based on particle type (id).
     * It contains functions to insert and erase particles while automatically moving particles
     * above the deletion or insertion point in sync.
     *
     * Example:
     *
     *     AtomTracker track(myspace);
     *     track.insert( myparticle, 20 );        // insert particle into Space at position 20
     *     ...
     *     int i=track[ myparticle.id ].random(); // pick a random particle of type myparticle.id
     */
    template<class Tspace>
      class AtomTracker {
        public:
          typedef int Tindex; // particle index type
        private:
          Tspace* spc;
          class data {
            public:
              vector<Tindex> index;
              Tindex random();                  //!< Pick random particle index
          };
          std::map<particle::Tid,data> map;
        public:
          AtomTracker(Tspace&);
          particle::Tid randomAtomType() const; //!< Select a random atomtype from the list
          bool insert(const particle&, Tindex); //!< Insert particle into Space and track position
          bool erase(Tindex);                   //!< Delete particle from Space at specific particle index
          data& operator[] (particle::Tid);     //!< Access operator to atomtype data
          void clear();                         //!< Clear all atom lists (does not touch Space)
          bool empty();                         //!< Test if atom list is empty
      };

    template<class Tspace>
      bool AtomTracker<Tspace>::empty() {
        return map.empty();
      }

    template<class Tspace>
      particle::Tid AtomTracker<Tspace>::randomAtomType() const {
        assert(!map.empty() && "No atom types have been added yet");
        vector<particle::Tid> vid;
        vid.reserve( map.size() );
        for (auto &m : map)
          vid.push_back(m.first);
        std::random_shuffle(vid.begin(), vid.end());
        return vid.front();
      }

    template<class Tspace>
      void AtomTracker<Tspace>::clear() {
        map.clear();
      }

    template<class Tspace>
      AtomTracker<Tspace>::AtomTracker(Tspace &s) { spc=&s; }

    template<class Tspace>
      typename AtomTracker<Tspace>::Tindex AtomTracker<Tspace>::data::random() {
        std::random_shuffle( index.begin(), index.end() );
        return index.back();
      }

    template<class Tspace>
      typename AtomTracker<Tspace>::data& AtomTracker<Tspace>::operator[](particle::Tid id) {
        return map[id];
      }

    /**
     * This will insert a particle into Space and at the same time make sure
     * that all other particles are correctly tracked.
     */
    template<class Tspace>
      bool AtomTracker<Tspace>::insert(const particle &a, Tindex index) {
        assert( a.id == spc->p[ map[a.id].index.back() ].id && "Id mismatch");
        spc->insert(a, index); // insert into Space
        for (auto &m : map)    // loop over all maps
          for (auto &i : m.second.index) // and their particle index
            if (i>=index) i++; // push forward particles beyond inserted particle
        map[a.id].index.push_back(index); // finally, add particle to appripriate map
        return true;
      }

    template<class Tspace>
      bool AtomTracker<Tspace>::erase(AtomTracker<Tspace>::Tindex index) {
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

    /**
     * @brief Grand Canonical insertion of arbitrary M:X salt pairs
     * @author Bjorn Persson and Mikael Lund
     * @date Lund 2010-2011
     * @warning Untested for asymmetric salt in this branch
     */
    template<class Tspace>
      class GrandCanonicalSalt : public Movebase<Tspace> {
        private:
          typedef Movebase<Tspace> base;
          using base::pot;
          using base::spc;
          using base::w;
          string _info();
          void _trialMove();
          void _acceptMove();
          void _rejectMove();
          double _energyChange();
          void add(Group&);       // scan group for ions with non-zero activities

          AtomTracker<Tspace> tracker;
          struct ionprop {
            particle p;
            double chempot;       // chemical potential log(1/A3)
            Average<double> rho;  // average density
          };
          std::map<particle::Tid,ionprop> map;
          void randomIonPair(particle::Tid&, particle::Tid&);  // Generate random ion pair
          p_vec trial_insert;
          vector<int> trial_delete;
          particle::Tid ida, idb;     // particle id's of current salt pair (a=cation, b=anion)

          Group* saltPtr;  // GC ions *must* be in this group

          // unit testing
          void _test(UnitTest &t) {
            for (auto &m : map) {
              auto s=base::prefix+"_"+atom[m.first].name;
              t(s+"_activity", atom[m.first].activity);
              t(s+"_conc", m.second.rho.avg()/pc::Nav/1e-27);
            }
          }

        public:
          GrandCanonicalSalt(InputMap&, Energy::Energybase<Tspace>&, Tspace&, Group&, string="saltbath");
      };

    template<class Tspace>
      GrandCanonicalSalt<Tspace>::GrandCanonicalSalt(
          InputMap &in, Energy::Energybase<Tspace> &e, Tspace &s, Group &g, string pfx) :
        base(e,s,pfx), tracker(s) {
          base::title="Grand Canonical Salt";
          base::useAlternateReturnEnergy=true;
          w=30;
          base::runfraction = in.get<double>(pfx+"_runfraction",1.0);
          saltPtr=&g;
          add(*saltPtr);
        }

    template<class Tspace>
      void GrandCanonicalSalt<Tspace>::add(Group &g) {
        assert( g.isAtomic() && "Salt group must be atomic" );
        spc->enroll(g);
        tracker.clear();
        for (auto i : g) {
          auto id=spc->p[i].id;
          if ( atom[id].activity>1e-10 && abs(atom[id].charge)>1e-10 ) {
            map[id].p=atom[id];
            map[id].chempot=log( atom[id].activity*pc::Nav*1e-27); // beta mu
            tracker[id].index.push_back(i);
          }
        }
        assert(!tracker.empty() && "No GC ions found!");
      }

    template<class Tspace>
      void GrandCanonicalSalt<Tspace>::randomIonPair(particle::Tid &id_cation, particle::Tid &id_anion) {
        do id_anion  = tracker.randomAtomType(); while ( map[id_anion].p.charge>=0);
        do id_cation = tracker.randomAtomType(); while ( map[id_cation].p.charge<=0);
        assert( !tracker[id_anion].index.empty() && "Ion list is empty");
        assert( !tracker[id_cation].index.empty() && "Ion list is empty");
      }

    template<class Tspace>
      void GrandCanonicalSalt<Tspace>::_trialMove() {
        trial_insert.clear();
        trial_delete.clear();
        randomIonPair(ida, idb);
        assert(ida>0 && idb>0 &&
            "Ion pair id is zero (UNK). Is this really what you want?");
        int Na = (int)abs(map[idb].p.charge);
        int Nb = (int)abs(map[ida].p.charge);
        switch ( rand() % 2) {
          case 0:
            trial_insert.reserve(Na+Nb);
            do trial_insert.push_back( map[ida].p ); while (--Na>0);
            do trial_insert.push_back( map[idb].p ); while (--Nb>0);
            for (auto &p : trial_insert)
              spc->geo.randompos(p);
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

    template<class Tspace>
      double GrandCanonicalSalt<Tspace>::_energyChange() {
        int Na=0, Nb=0;            // number of added or deleted ions
        double idfactor=1;
        double uold=0, unew=0, V=spc->geo.getVolume();
        double potold=0, potnew=0; // energy change due to interactions
        if ( !trial_insert.empty() ) {
          for (auto &t : trial_insert)     // count added ions
            if (t.id==map[ida].p.id) Na++; else Nb++;
          for (int n=0; n<Na; n++)
            idfactor *= (tracker[ida].index.size()+1+n)/V;
          for (int n=0; n<Nb; n++)
            idfactor *= (tracker[idb].index.size()+1+n)/V;

          unew = log(idfactor) - Na*map[ida].chempot - Nb*map[idb].chempot;

          potnew += pot->v2v(spc->p, trial_insert);
          for (auto i=trial_insert.begin(); i!=trial_insert.end()-1; i++)
            for (auto j=i+1; j!=trial_insert.end(); j++)
              potnew+=pot->p2p(*i,*j);
          for (auto i=trial_insert.begin(); i!=trial_insert.end(); i++)
            potnew+=pot->p_external(*i);
          unew+=potnew;
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

          for (auto &i : trial_delete)
            potold+=pot->i_total(spc->p, i);
          for (auto i=trial_delete.begin(); i!=trial_delete.end()-1; i++)
            for (auto j=i+1; j!=trial_delete.end(); j++)
              potold-=pot->i2i(spc->p, *i, *j);
          uold+=potold;
        } else {
          assert(!"No salt to insert or delete!");
          std::cerr << "!! No salt to insert or delete !!";
        }
        base::alternateReturnEnergy=potnew-potold; // track only pot. energy
        return unew-uold;
      }

    template<class Tspace>
      void GrandCanonicalSalt<Tspace>::_acceptMove() {
        if ( !trial_insert.empty() ) {
          for (auto &p : trial_insert)
            tracker.insert(p, saltPtr->back());
        }
        else if ( !trial_delete.empty() ) {
          std::sort(trial_delete.rbegin(), trial_delete.rend()); //reverse sort
          for (auto i : trial_delete)
            tracker.erase(i);
        }
        double V = spc->geo.getVolume();
        map[ida].rho += tracker[ida].index.size() / V;
        map[idb].rho += tracker[idb].index.size() / V;
      }

    template<class Tspace>
      void GrandCanonicalSalt<Tspace>::_rejectMove() {
        double V = spc->geo.getVolume();
        map[ida].rho += tracker[ida].index.size() / V;
        map[idb].rho += tracker[idb].index.size() / V;
      }

    template<class Tspace>
      string GrandCanonicalSalt<Tspace>::_info() {
        char s=10;
        using namespace textio;
        std::ostringstream o;
        o << pad(SUB,w,"Number of GC species") << endl << endl;
        o << setw(4) << "" << std::left
          << setw(s) << "Ion" << setw(s) << "activity"
          << setw(s+4) << bracket("c/M") << setw(s+6) << bracket( gamma+pm )
          << setw(s+4) << bracket("N") << "\n";
        for (auto &m : map) {
          particle::Tid id=m.first;
          o.precision(5);
          o << setw(4) << "" << setw(s) << atom[id].name
            << setw(s) << atom[id].activity << setw(s) << m.second.rho.avg()/pc::Nav/1e-27
            << setw(s) << atom[id].activity / (m.second.rho.avg()/pc::Nav/1e-27)
            << setw(s) << m.second.rho.avg()*spc->geo.getVolume()
            << "\n";
        }
        return o.str();
      }

#ifdef ENABLE_MPI
    /**
     * @brief Class for parallel tempering (aka replica exchange) using MPI
     *
     * This will perform replica exchange moves by the following steps:
     *
     * -# Randomly find an exchange partner with rank above/under current rank
     * -# Exchange full particle configuration with partner
     * -# Calculate energy change using Energy::systemEnergy. Note that this
     *    energy function can be replaced by setting the `ParallelTempering::usys`
     *    variable to another function with the same signature (functor wrapper).
     * -# Send/receive energy change to/from partner
     * -# Accept or reject based on *total* energy change
     *
     * Although not completely correct, the recommended way of performing a temper move
     * is to do `N` Monte Carlo passes with regular moves and then do a tempering move.
     * This is because the MPI nodes must be in sync and if you have a system where
     * the random number generator calls are influenced by the Hamiltonian we could 
     * end up in a deadlock.
     *
     * @date Lund 2012
     */
    template<class Tspace>
      class ParallelTempering : public Movebase<Tspace> {
        private:
          typedef Movebase<Tspace> base;
          using base::pot;
          using base::spc;
          using base::w;
          using base::runfraction;
          enum extradata {VOLUME=0};    //!< Structure of extra data to send
          typedef std::map<string, Average<double> > map_type;
          map_type accmap;              //!< Acceptance map
          int partner;                  //!< Exchange replica (partner)
          virtual void findPartner();   //!< Find replica to exchange with
          bool goodPartner();           //!< Is partned valid?
          double exchangeEnergy(double);//!< Exchange energy with partner
          string id();                  //!< Unique string to identify set of partners

          double currentEnergy;         //!< Energy of configuration before move (uold)
          bool haveCurrentEnergy;       //!< True if currentEnergy has been set

          string _info();
          void _trialMove();
          void _acceptMove();
          void _rejectMove();
          double _energyChange();
          std::ofstream temperPath;

          Faunus::MPI::MPIController *mpiPtr; //!< Controller class for MPI calls
          Faunus::MPI::FloatTransmitter ft;   //!< Class for transmitting floats over MPI
          Faunus::MPI::ParticleTransmitter pt;//!< Class for transmitting particles over MPI
          std::function<double (Tspace&, Energy::Energybase<Tspace>&, const p_vec&)> usys; //!< Defaults to Energy::systemEnergy but can be replaced!

        public:
          ParallelTempering(InputMap&, Energy::Energybase<Tspace>&, Tspace&, Faunus::MPI::MPIController &mpi, string="temper");
          virtual ~ParallelTempering();
          void setCurrentEnergy(double); //!< Set energy of configuration before move (for increased speed)
          void setEnergyFunction( std::function<double (Tspace&,Energy::Energybase<Tspace>&,const p_vec&)> );
      };

    template<class Tspace>
      ParallelTempering<Tspace>::ParallelTempering(
          InputMap &in,
          Energy::Energybase<Tspace> &e,
          Tspace &s,
          Faunus::MPI::MPIController &mpi,
          string pfx) : base(e,s,pfx), mpiPtr(&mpi) {
        this->title="Parallel Tempering";
        partner=-1;
        this->useAlternateReturnEnergy=true; //we don't want to return dU from partner replica (=drift)
        this->runfraction = in.get<double>(pfx+"_runfraction",1);
        pt.recvExtra.resize(1);
        pt.sendExtra.resize(1);
        pt.setFormat( in.get<string>(pfx+"_format", "XYZQI") );
        setEnergyFunction(
            Energy::systemEnergy<Tspace,Energy::Energybase<Tspace>,p_vec> );
        this->haveCurrentEnergy=false;
        //temperPath.open(textio::prefix+"temperpath.dat");
      }

    template<class Tspace>
      ParallelTempering<Tspace>::~ParallelTempering() {}

    template<class Tspace>
      void ParallelTempering<Tspace>::setEnergyFunction( std::function<double (Tspace&,Energy::Energybase<Tspace>&,const p_vec&)> f ) {
        usys = f;
      }

    template<class Tspace>
      void ParallelTempering<Tspace>::findPartner() {
        int dr=0;
        partner = mpiPtr->rank();
        if (mpiPtr->random()>0.5)
          dr++;
        else
          dr--;
        if (mpiPtr->rank() % 2 == 0)
          partner+=dr;
        else
          partner-=dr;
      }

    template<class Tspace>
      bool ParallelTempering<Tspace>::goodPartner() {
        assert(partner!=mpiPtr->rank() && "Selfpartner! This is not supposed to happen.");
        if (partner>=0)
          if ( partner<mpiPtr->nproc() )
            if ( partner!=mpiPtr->rank() )
              return true;
        return false;
      }

    template<class Tspace>
      string ParallelTempering<Tspace>::_info() {
        using namespace textio;
        std::ostringstream o;
        o << pad(SUB,w,"Process rank") << mpiPtr->rank() << endl
          << pad(SUB,w,"Number of replicas") << mpiPtr->nproc() << endl
          << pad(SUB,w,"Data size format") << short(pt.getFormat()) << endl
          << indent(SUB) << "Acceptance:"
          << endl;
        if (this->cnt>0) {
          o.precision(3);
          for (auto &m : accmap)
            o << indent(SUBSUB) << std::left << setw(12)
              << m.first << setw(8) << m.second.cnt << m.second.avg()*100
              << percent << endl;
        }
        return o.str();
      }

    template<class Tspace>
      void ParallelTempering<Tspace>::_trialMove() {
        findPartner();
        if (goodPartner()) {

          pt.sendExtra[VOLUME]=spc->geo.getVolume();  // copy current volume for sending

          pt.recv(*mpiPtr, partner, spc->trial); // receive particles
          pt.send(*mpiPtr, spc->p, partner);     // send everything
          pt.waitrecv();
          pt.waitsend();

          // update group trial mass-centers. Needed if energy calc. uses
          // cm_trial for cut-offs, for example
          for (auto g : spc->groupList())
            g->cm_trial = Geometry::massCenter(spc->geo, spc->trial, *g);

          // debug assertions
          assert(pt.recvExtra[VOLUME]>1e-6 && "Invalid partner volume received.");
          assert(spc->p.size() == spc->trial.size() && "Particle vectors messed up by MPI");

          // release assertions
          if (pt.recvExtra[VOLUME]<1e-6 || spc->p.size() != spc->trial.size())
            MPI_Abort(mpiPtr->comm, 1);
        }
      }

    /**
     * If the system energy is already known it may be specified with this
     * function to speed up the calculation. If not set, it will be calculated.
     */
    template<class Tspace>
      void ParallelTempering<Tspace>::setCurrentEnergy(double uold) {
        currentEnergy=uold;
        haveCurrentEnergy=true;
      }

    template<class Tspace>
      double ParallelTempering<Tspace>::_energyChange() {
        this->alternateReturnEnergy=0;
        if ( !goodPartner() )
          return pc::infty;
        double uold, du_partner;

        if (haveCurrentEnergy)   // do we already know the energy?
          uold = currentEnergy;
        else
          uold = usys(*spc,*pot,spc->p);

        spc->geo.setVolume( pt.recvExtra[VOLUME] ); // set new volume
        pot->setSpace(*spc);

        double unew = usys(*spc,*pot,spc->trial);

        du_partner = exchangeEnergy(unew-uold); // Exchange dU with partner (MPI)

        haveCurrentEnergy=false;                // Make sure user call setCurrentEnergy() before next move
        this->alternateReturnEnergy=unew-uold;        // Avoid energy drift (no effect on sampling!)
        return (unew-uold)+du_partner;          // final Metropolis trial energy
      }

    /**
     * This will exchange energies with replica partner
     * @todo Use FloatTransmitter::swapf() instead.
     *       Use C++11 initializer list for vectors, i.e. vector<floatp> v={mydu};
     */
    template<class Tspace>
      double ParallelTempering<Tspace>::exchangeEnergy(double mydu) {
        vector<MPI::FloatTransmitter::floatp> duSelf(1), duPartner;
        duSelf[0]=mydu;
        duPartner = ft.swapf(*mpiPtr, duSelf, partner);
        return duPartner.at(0);               // return partner energy change
      }

    template<class Tspace>
      string ParallelTempering<Tspace>::id() {
        std::ostringstream o;
        if (mpiPtr->rank() < partner)
          o << mpiPtr->rank() << " <-> " << partner;
        else
          o << partner << " <-> " << mpiPtr->rank();
        return o.str();
      }

    template<class Tspace>
      void ParallelTempering<Tspace>::_acceptMove(){
        if ( goodPartner() ) {
          //temperPath << cnt << " " << partner << endl;
          accmap[ id() ] += 1;
          for (size_t i=0; i<spc->p.size(); i++)
            spc->p[i] = spc->trial[i];  // copy new configuration
          for (auto g : spc->groupList())
            g->cm = g->cm_trial;
        }
      }

    template<class Tspace>
      void ParallelTempering<Tspace>::_rejectMove() {
        if ( goodPartner() ) {
          spc->geo.setVolume( pt.sendExtra[VOLUME] ); // restore old volume
          pot->setSpace(*spc);
          accmap[ id() ] += 0;
          for (size_t i=0; i<spc->p.size(); i++)
            spc->trial[i] = spc->p[i];   // restore old configuration
          for (auto g : spc->groupList())
            g->cm_trial = g->cm;
        }
      }
#endif

    /**
     * @brief Swap atom charges
     *
     * This move selects two particle index from a user-defined list and swaps
     * their charges.
     *
     * @date Lund, 2013
     */
    template<class Tspace>
      class SwapCharge : public Movebase<Tspace> {
        private:
          typedef Movebase<Tspace> base;
          typedef std::map<short, Average<double> > map_type;
          string _info();
        protected:
          void _acceptMove();
          void _rejectMove();
          double _energyChange();
          void _trialMove();
          using base::spc;
          using base::pot;
          map_type accmap; //!< Single particle acceptance map        
          int ip, jp;

          //int iparticle;   //!< Select single particle to move (-1 if none, default)

        public:
          SwapCharge(InputMap&, Energy::Energybase<Tspace>&, Tspace&, string="swapcharge");
          std::set<int> swappableParticles;  //!< Particle index that can be swapped
      };  

    template<class Tspace>
      SwapCharge<Tspace>::SwapCharge(InputMap &in,Energy::Energybase<Tspace> &e, Tspace &s, string pfx) : base(e,s,pfx) {
        base::title="Swap head groups of different charges";
      }

    template<class Tspace>
      void SwapCharge<Tspace>::_trialMove() {
        assert(!swappableParticles.empty());
        auto vi=swappableParticles.begin();
        auto vj=swappableParticles.begin();
        std::advance(vi, slp_global.rand() % swappableParticles.size());
        std::advance(vj, slp_global.rand() % swappableParticles.size());
        ip=*(vi);
        jp=*(vj);
        if ( spc->trial[ip].charge != spc->trial[jp].charge ) {
          std::swap( spc->trial[ip].charge, spc->trial[jp].charge );
        }
      }
    template<class Tspace>
      double SwapCharge<Tspace>::_energyChange() {
        return base::pot->i_total(spc->trial, jp) + base::pot->i_total(spc->trial, ip) 
          - base::pot->i_total(spc->p, jp) - base::pot->i_total(spc->p, ip);
      }
    template<class Tspace>
      void SwapCharge<Tspace>::_acceptMove() {
        accmap[ spc->p[ip].id ] += 1;
        spc->p[ip].charge = spc->trial[ip].charge;
        spc->p[jp].charge = spc->trial[jp].charge;
      }

    template<class Tspace>
      void SwapCharge<Tspace>::_rejectMove() {
        accmap[ spc->p[ip].id ] += 0;
        spc->trial[ip].charge = spc->p[ip].charge;
        spc->trial[jp].charge = spc->p[jp].charge;
      }
    template<class Tspace>
      string SwapCharge<Tspace>::_info() {
        using namespace textio;
        std::ostringstream o;
        o << pad(SUB,base::w,"Average moves/particle")
          << base::cnt / swappableParticles.size() << endl;
        if (base::cnt>0) {
          char l=12;
          o << endl
            << indent(SUB) << "Individual particle movement:" << endl << endl
            << indent(SUBSUB) << std::left << string(7,' ')
            << setw(l+1) << "Acc. "+percent;
          for (auto m : accmap) {
            auto id=m.first;
            o << indent(SUBSUB) << std::left << setw(7) << atom[id].name;
            o.precision(3);
            o << setw(l) << accmap[id].avg()*100;
          }
        }
        return o.str();
      }

    /**
     * @brief Make a flip-flip move on lipids
     *
     * @date Lund, 2013
     */
    template<class Tspace>
      class FlipFlop : public Movebase<Tspace> {
        private:
          typedef Movebase<Tspace> base;
        protected:
          using base::spc;
          using base::pot;
          using base::w;
          using base::cnt;
          using base::prefix;
          void _trialMove();
          void _acceptMove();
          void _rejectMove();
          double _energyChange();
          string _info();
          typedef std::map<string, Average<double> > map_type;
          map_type accmap;   //!< Group particle acceptance map
          Group* igroup;
          Point* cntr;
        public:
          FlipFlop(InputMap&, Energy::Energybase<Tspace>&, Tspace&, string="flipflop");
          void setGroup(Group&); //!< Select Group to move
          void setCenter(Point&); //!< Select Center of Mass of the vescicle
      };

    template<class Tspace>
      FlipFlop<Tspace>::FlipFlop(InputMap &in,Energy::Energybase<Tspace> &e, Tspace &s, string pfx) : base(e,s,pfx) {
        base::title="Group Flip-Flop Move";
        base::w=30;
        igroup=nullptr;
        cntr=nullptr;
        this->runfraction = in.get<double>(base::prefix+"_runfraction",1.0);
      }

    template<class Tspace>
      void FlipFlop<Tspace>::setGroup(Group &g) {
        assert(&g!=nullptr);
        assert(g.isMolecular());
        igroup=&g;
      }

    template<class Tspace>
      void FlipFlop<Tspace>::setCenter(Point &center) {
        assert(&center!=nullptr);
        cntr=&center;
      }

    template<class Tspace>
      void FlipFlop<Tspace>::_trialMove() {
        assert(igroup!=nullptr);
        assert(cntr!=nullptr);
        Point startpoint=spc->p[igroup->back()];
        Point head=spc->p[igroup->front()];
        cntr->z()=head.z()=startpoint.z();
        Point dir = spc->geo.vdist(*cntr, startpoint)
          / sqrt(spc->geo.sqdist(*cntr, startpoint)) * 1.1*spc->p[igroup->back()].radius;
        if (spc->geo.sqdist(*cntr, startpoint) > spc->geo.sqdist(*cntr, head))
          startpoint.translate(spc->geo,-dir);      // set startpoint for rotation
        else startpoint.translate(spc->geo, dir);
        double x1=cntr->x();
        double y1=cntr->y();
        double x2=startpoint.x();
        double y2=startpoint.y();
        Point endpoint; // rot endpoint for on axis ⊥ to line connecting cm of cylinder with 2nd TL
        endpoint.x()=x2+1;
        endpoint.y()=-(x2-x1)/(y2-y1)+y2;
        endpoint.z()=startpoint.z();
        double angle=pc::pi;
        Geometry::QuaternionRotate vrot;
        vrot.setAxis(spc->geo, startpoint, endpoint, angle); //rot around startpoint->endpoint vec
        for (auto i : *igroup)
          spc->trial[i] = vrot(spc->trial[i]);
        igroup->cm_trial = vrot(igroup->cm_trial);
      }

    template<class Tspace>
      void FlipFlop<Tspace>::_acceptMove() {
        accmap[ igroup->name ] += 1;
        igroup->accept(*spc);
      }

    template<class Tspace>
      void FlipFlop<Tspace>::_rejectMove() {
        accmap[ igroup->name ] += 0;
        igroup->undo(*spc);
      }

    template<class Tspace>
      double FlipFlop<Tspace>::_energyChange() {

        for (auto i : *igroup)
          if ( spc->geo.collision( spc->trial[i], Geometry::Geometrybase::BOUNDARY ) )
            return pc::infty;

        double unew = pot->g_external(spc->trial, *igroup);
        if (unew==pc::infty)
          return pc::infty;       // early rejection
        double uold = pot->g_external(spc->p, *igroup);

        for (auto g : spc->groupList()) {
          if (g!=igroup) {
            unew += pot->g2g(spc->trial, *g, *igroup);
            if (unew==pc::infty)
              return pc::infty;   // early rejection
            uold += pot->g2g(spc->p, *g, *igroup);
          }
        }
        return unew-uold;
      }

    template<class Tspace>
      string FlipFlop<Tspace>::_info() {
        using namespace textio;
        std::ostringstream o;
        if (cnt>0) {
          char l=12;
          o << indent(SUB) << "Move Statistics:" << endl
            << indent(SUBSUB) << std::left << setw(20) << "Group name" //<< string(20,' ')
            << setw(l+1) << "Acc. "+percent << endl;
          for (auto m : accmap) {
            string id=m.first;
            o << indent(SUBSUB) << std::left << setw(20) << id;
            o.precision(3);
            o << setw(l) << accmap[id].avg()*100 << endl;
          }
        }
        return o.str();
      }

    /** @brief Atomic translation with dipolar polarizability */
    //typedef PolarizeMove<AtomicTranslation> AtomicTranslationPol;

    /** @brief Atomic rotation with dipolar polarizability */
    //typedef PolarizeMove<AtomicRotation> AtomicRotationPol;

  }//namespace
}//namespace
#endif
