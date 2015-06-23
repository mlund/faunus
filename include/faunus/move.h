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
#include <faunus/json.h>
#include <faunus/titrate.h>

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
              << setw(l) << "Nmoves"
              << setw(l+9) << rootof+bracket("msq"+squared)+"/"+angstrom << endl;
            for (auto m : accmap) {
              string id=m.first;
              o << indent(SUBSUB) << std::left << setw(20) << id;
              o.precision(3);
              o << setw(l) << accmap[id].avg()*100
                << setw(l) << accmap[id].cnt
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
     * @brief Add polarisation step to an arbitrary move
     *
     * This class will modify any MC move to account for polarization
     * using an iterative procedure.
     * An electric field calculation is inserted
     * after the original trial move whereafter it will iteratively
     * calculate induced dipole moments on all particles.
     * The energy change function will evaluate the *total*
     * system energy as all dipoles in the system may have changed.
     * This is thus an expensive computation and is best used with
     * MC moves that propagate  all particles.
     *
     * ### Update frequency ###
     *
     * Updating induced moments is an iterative N*N operation and
     * very inefficient for MC moves that update only a subset of the system.
     * In liquid systems that propagate only slowly as a function of MC steps
     * one may attempt to update induced dipoles less frequently at the
     * expense of accuracy.
     * For repeating moves -- i.e. molecular translate/rotate or atomic
     * translation -- polarisation is updated only after all moves have
     * been carried out.
     *
     * @note Will currently not work for Grand Caninical moves
     */
    template<class Tmove>
      class PolarizeMove : public Tmove {
        private:
          using Tmove::spc;
          using Tmove::pot;
          int Ntrials;                    // Number of repeats within move
          int max_iter;                   // max numbr of iterations
          double threshold;       	  // threshold for iteration
          bool updateDip;                 // true if ind. dipoles should be updated
          Eigen::MatrixXd field;  	  // field on each particle
          Average<int> numIter;           // average number of iterations per move

          /**
           *  @brief Updates dipole moment w. permanent plus induced dipole moment
           *  @param pot Hamiltonian
           *  @param p Particles to update
           */
          template<typename Tenergy,typename Tparticles>
            void induceDipoles(Tenergy &pot, Tparticles &p) { 

              int cnt=0;
              Eigen::VectorXd mu_err_norm( (int)p.size() );

              do {
                cnt++;
                mu_err_norm.setZero();
                field.setZero();
                pot.field(p,field);
                for ( size_t i=0; i<p.size(); i++ ) {
                  Point E = field.col(i);                  // field on i
                  Point mu_trial = p[i].alpha*E + p[i].mup;// new tot. dipole
                  Point mu_err = mu_trial - p[i].mu*p[i].muscalar;// mu difference
                  mu_err_norm[i] = mu_err.norm();          // norm of previous row
                  p[i].muscalar = mu_trial.norm();         // update dip scalar in particle
                  if (p[i].muscalar > 1e-6)
                    p[i].mu = mu_trial/p[i].muscalar;      // update article dip.
                }
                if ( cnt > max_iter )
                  throw std::runtime_error("Field induction reached maximum number of iterations.");

              } while ( mu_err_norm.maxCoeff() > threshold ); // is threshold OK?

              numIter += cnt; // average number of iterations
            }

          void _trialMove() FOVERRIDE {
            Tmove::_trialMove();

            Ntrials++;
            int updateAt = 1;  // default: dipoles are always updated
            if ( ! Tmove::mollist.empty() )
              updateAt = Tmove::mollist[ Tmove::currentMolId ].repeat;
            else
              Ntrials = 1;      // in case move(n) is called w. n>1

            updateDip = ( Ntrials == updateAt ); 

            if ( updateDip ) {
              field.resize( 3, spc->trial.size() );
              induceDipoles( *pot, spc->trial );
            }
          }

          double _energyChange() FOVERRIDE {
            if ( updateDip )
              return Energy::systemEnergy( *spc, *pot, spc->trial )
                - Energy::systemEnergy( *spc, *pot, spc->p );
            else
              return Tmove::_energyChange();
          }

          void _rejectMove() FOVERRIDE {
            Tmove::_rejectMove();
            if (updateDip)
              Tmove::spc->trial = Tmove::spc->p;
          }

          void _acceptMove() FOVERRIDE {
            Tmove::_acceptMove();
            if (updateDip)
              Tmove::spc->p = Tmove::spc->trial;
          }

          string _info() FOVERRIDE {
            std::ostringstream o;
            using namespace textio;
            o << pad(SUB,Tmove::w,"Polarisation updates") << numIter.cnt << "\n"
              << pad(SUB,Tmove::w,"Polarisation threshold") << threshold << "\n"
              << pad(SUB,Tmove::w,"Polarisation iterations") << numIter.avg()
              << " (max. " << max_iter << ")" << "\n"
              << Tmove::_info();
            return o.str();
          }

        public:
          template<class Tspace>
            PolarizeMove(InputMap &in, Energy::Energybase<Tspace> &e, Tspace &s) :
              Tmove(in,e,s) {
                threshold = in.get<double>("pol_threshold", 0.001, "Iteration precision");
                max_iter = in.get<int>("max_iterations", 40, "Max. iterations");
              }

          template<class Tspace>
            PolarizeMove(Energy::Energybase<Tspace> &e, Tspace &s, Tmjson &j) :
              Tmove(e,s,j) {
                threshold = j["pol_threshold"] | 0.001;
                max_iter  = j["max_iterations"] | 40;
              }

          PolarizeMove(const Tmove &m) : max_iter(40), threshold(0.001), Tmove(m) {};

          double move(int n) FOVERRIDE {
            Ntrials = 0;
            return Tmove::move(n);
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
    template<class Tspace=Space<class Tgeometry,class Tparticle> >
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

          TimeRelativeOfTotal<std::chrono::microseconds> timer;

        protected:
          virtual string _info()=0;        //!< info for derived moves
          void trialMove();                //!< Do a trial move (wrapper)
          Energy::Energybase<Tspace>* pot; //!< Pointer to energy functions
          Tspace* spc;                     //!< Pointer to Space
          string title;                    //!< Title of move (mandatory!)
          string cite;                     //!< Reference, url, DOI etc.
          string jsondir;                  //!< inputmap section
          char w;                          //!< info text width. Adjust this in constructor if needed.
          unsigned long int cnt;           //!< total number of trial moves
          virtual bool run() const;        //!< Runfraction test

          bool useAlternateReturnEnergy;   //!< Return a different energy than returned by _energyChange(). [false]
          double alternateReturnEnergy;    //!< Alternative return energy

          struct MolListData {
            double prob;  // probability of performing a move
            bool perAtom; // repeat move for each molecule?
            bool perMol;  // repeat move for atom in molecules?
            int repeat;   // total number of repeats
            unsigned long Nattempts;  // # of attempted moves
            unsigned long Naccepted;  // # of accepted moves
            Point dir;    // translational move directions
            double dp1;   // displacement parameter 1
            double dp2;   // displacement parameter 2

            MolListData() : prob(1.0), perAtom(false), perMol(false),
            repeat(1), Nattempts(0), Naccepted(0), dir(1,1,1), dp1(0), dp2(0) {}

            MolListData( Tmjson &j ) {
              *this = MolListData();
              prob = j["prob"] | 1.0;
              perMol = j["permol"] | false;
              perAtom = j["peratom"] | false;
              dir << ( j["dir"] | std::string("1 1 1") );
            }
          };

          std::map<int,MolListData> mollist;    //!< Move acts on these molecule id's

          /**
           * @brief Iterate over json object where each key is a molecule
           *        name and the value is read as `MolListData`.
           */
          void fillMolList( Tmjson &j ) {
            for (auto it=j.begin(); it!=j.end(); ++it) {  // iterate over molecules
              auto mol = spc->molList().find( it.key() ); // is molecule defined?
              if ( mol == spc->molList().end() ) {
                std::cerr << "Error: molecule '" << it.key() << "' not defined.\n";
                exit(1); 
              } else
                addMol( mol->id, MolListData( it.value() ) );
            }
          }

        public:
          Movebase(Energy::Energybase<Tspace>&, Tspace&, string="moves");//!< Constructor
          virtual ~Movebase();
          double runfraction;                //!< Fraction of times calling move() should result in an actual move. 0=never, 1=always.
          virtual double move(int=1);                //!< Attempt `n` moves and return energy change (kT)
          string info();                     //!< Returns information string
          void test(UnitTest&);              //!< Perform unit test
          double getAcceptance();            //!< Get acceptance [0:1]

          void addMol(int, const MolListData&d=MolListData()); //!< Specify molecule id to act upon
          Group* randomMol();
          int randomMolId();                 //!< Random mol id from mollist
          int currentMolId;                  //!< Current molid to act upon

#ifdef ENABLE_MPI
          Faunus::MPI::MPIController* mpiPtr;
#endif
      };

    /**
     * @brief Constructor
     * @param e Energy class
     * @param s Space
     * @param rootsec Name of section in input file to search for parameters
     */
    template<class Tspace>
      Movebase<Tspace>::Movebase(Energy::Energybase<Tspace> &e, Tspace &s, string rootsec) {
        jsondir = rootsec;
        e.setSpace(s);
        pot=&e;
        spc=&s;
        cnt=cnt_accepted=0;
        dusum=0;
        w=30;
        runfraction=1;
        useAlternateReturnEnergy=false; //this has no influence on metropolis sampling!
#ifdef ENABLE_MPI
        mpiPtr=nullptr;
#endif
      }

    template<class Tspace>
      Movebase<Tspace>::~Movebase() {}

    template<class Tspace>
      void Movebase<Tspace>::addMol(int molid, const MolListData &d) {
        mollist[ molid ] = d;
      }

    template<class Tspace>
      int Movebase<Tspace>::randomMolId() {
        if ( !mollist.empty() ) {
          auto it = propagation_slump.element( mollist.begin(), mollist.end() );
          if (it != mollist.end() ) {
            it->second.repeat = 1;
            if ( it->second.perMol )
              it->second.repeat *= spc->numMolecules( it->first );
            if ( it->second.perAtom )
              it->second.repeat *= spc->findMolecules( it->first ).front()->size(); 
            return it->first;
          }
        }
        return -1;
      }

    /**
     * Returns pointer to a random group matching a molecule id
     * in `mollist`
     */
    template<class Tspace>
      Group* Movebase<Tspace>::randomMol() {
        Group* gPtr=nullptr;
        if ( !mollist.empty() ) {
          auto it = propagation_slump.element( mollist.begin(), mollist.end() );
          auto g = spc->findMolecules( it->first ); // vector of group pointers
          if ( !g.empty() )
            gPtr = *propagation_slump.element( g.begin(), g.end() );
        }
        return gPtr;
      }

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
        if (std::isnan(du))
          std::cerr << "Warning: energy change from move returns not-a-number (NaN)" << endl;
        return du;
      }

    /**
     * This function performs trial move and accept/reject using 
     * the Metropolis criteria.
     * It carries out the following `n` times:
     *
     * - Perform a trial move with `_trialMove()`
     * - Calulate the energy change, \f$\beta\Delta U\f$ with `_energyChange()`
     * - Accept with probability \f$ \min(1,e^{-\beta\Delta U}) \f$
     * - Call either `_acceptMove()` or `_rejectMove()`
     *
     * @note Do not override this function in derived classes.
     * @param n Perform move `n` times (default=1)
     *
     * [More info](http://dx.doi.org/10.1063/1.1699114)
     */
    template<class Tspace>
      double Movebase<Tspace>::move(int n) {
        timer.start();
        double utot=0;

        if ( ! mollist.empty() ) {
          currentMolId = randomMolId();
          n = mollist[ currentMolId ].repeat;
          runfraction = mollist[ currentMolId ].prob;
        }

        if ( run() ) {
          while ( n-->0 ) {
            trialMove();
            double du = energyChange();
            bool outcome = metropolis(du);
            if ( !outcome ) {
              rejectMove();
              utot += pot->penalty_update(outcome);
            }
            else {
              acceptMove();
              if ( useAlternateReturnEnergy )
                du = alternateReturnEnergy;
              utot += pot->penalty_update(outcome);
              dusum += du;
              utot += du;
            }
          }
        }
        assert(spc->p == spc->trial && "Trial particle vector out of sync!");
        timer.stop();
        return utot;
      }

    /**
     * @param du Energy change for MC move (kT)
     * @return True if move should be accepted; false if not.
     * @note
     * One could put in `if (du>0)` before the first line, but
     * certain MPI communications require the random number
     * generator to be in sync, i.e. each rank must call
     * `slump()` equal number of times, independent of
     * dU.
     */
    template<class Tspace>
      bool Movebase<Tspace>::metropolis(const double &du) const {
        if ( slump()>std::exp(-du) ) // core of MC!
          return false;
        return true;
      }

    template<class Tspace>
      bool Movebase<Tspace>::run() const {
        if (propagation_slump() < runfraction)
          return true;
        return false;
      }

    template<class Tspace>
      void Movebase<Tspace>::test(UnitTest &t) {
        if (runfraction<1e-6 || cnt==0)
          return;
        t(jsondir+"_acceptance", double(cnt_accepted)/cnt*100 );
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
            << pad(SUB,w,"Relative time consumption") << timer.result() << endl
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

          AtomicTranslation(Energy::Energybase<Tspace>&, Tspace&, Tmjson&, string="atomtranslate");

          void setGenericDisplacement(double); //!< Set single displacement for all atoms

          Point dir;             //!< Translation directions (default: x=y=z=1)
      };

    /**
     * @brief Constructor
     *
     * By default input is read from json section `atomtranslate`
     * with each element being the molecule name with the following
     * properties:
     *
     * Value                | Description
     * :------------------- | :-------------------------------------------------------------
     * `dir`                | Move directions (default: "1 1 1" = xyz)
     * `peratom`            | Repeat move for each atom in molecule (default: false)
     * `permol`             | Repeat move for each molecule in system (default: false)
     * `prob`               | Probability of performing the move (default: 1)
     *
     * Example:
     *
     *     atomictranslation {
     *       "salt" : { "dir":"1 1 0", "peratom":true }
     *     }
     *
     * Atomic displacement parameters are read from `Faunus::AtomData`.
     */
    template<class Tspace>
      AtomicTranslation<Tspace>::AtomicTranslation(
          Energy::Energybase<Tspace> &e,
          Tspace &s,
          Tmjson &j,
          string sec) : Movebase<Tspace>( e, s ) {

        base::title="Single Particle Translation";
        base::jsondir="moves/"+sec; // temp. fix to be removed
        iparticle=-1;
        igroup=nullptr;
        dir={1,1,1};
        genericdp = 0;

        base::fillMolList( j["moves"][sec] );
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
      bool AtomicTranslation<Tspace>::run() const {
        if ( igroup != nullptr )
          if ( igroup->empty() )
            return false;
        return base::run();
      }

    template<class Tspace>
      void AtomicTranslation<Tspace>::_trialMove() {

        if ( ! this->mollist.empty() ) {
          auto gvec = spc->findMolecules( this->currentMolId );
          assert( !gvec.empty() );
          igroup = *slump.element( gvec.begin(), gvec.end() );
          assert( ! igroup->empty() );
          dir = this->mollist[ this->currentMolId ].dir;
        }

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
          t.x() *= slump()-0.5;
          t.y() *= slump()-0.5;
          t.z() *= slump()-0.5;
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
          assert( spc->geo.collision(spc->p[iparticle], spc->p[iparticle].radius)==false
              && "An untouched particle collides with simulation container.");
          if ( spc->geo.collision(
                spc->trial[iparticle], spc->trial[iparticle].radius, Geometry::Geometrybase::BOUNDARY ) )
            return pc::infty;
          return
            (base::pot->i_total(spc->trial, iparticle)
             + base::pot->external(spc->trial) + base::pot->penalty(spc->trial))
            - (base::pot->i_total(spc->p, iparticle)
                + base::pot->external(spc->p) + base::pot->penalty(spc->p));
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
            << setw(l) << "Nmoves"
            << setw(l+7) << bracket("r"+squared)+"/"+angstrom+squared
            << rootof+bracket("r"+squared)+"/"+angstrom << endl;
          for (auto m : sqrmap) {
            auto id=m.first;
            o << indent(SUBSUB) << std::left << setw(7) << atom[id].name
              << setw(l-6) << ( (atom[id].dp<1e-6) ? genericdp : atom[id].dp);
            o.precision(3);
            o << setw(l) << accmap[id].avg()*100
              << setw(l) << accmap[id].cnt
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
          using base::accmap;
          using base::sqrmap;
          Geometry::QuaternionRotate rot;
          string _info();
          void _trialMove();
          void _acceptMove();
          void _rejectMove();
          double dprot;      //!< Temporary storage for current angle

        public:
          AtomicRotation(Energy::Energybase<Tspace>&, Tspace&,
              Tmjson&, string="atomrotate");
      };

    template<class Tspace>
      AtomicRotation<Tspace>::AtomicRotation(
          Energy::Energybase<Tspace> &e,
          Tspace &s,
          Tmjson &j,
          string sec) : base(e, s, j, sec) {
        base::title="Single Particle Rotation";
      }
    template<class Tspace>
      void AtomicRotation<Tspace>::_trialMove() {
        if ( ! this->mollist.empty() ) {
          igroup = spc->randomMol( this->currentMolId );
          if ( igroup != nullptr ) {
            iparticle = igroup->random();
            gsize += igroup->size();
          } else return;
        } else return;

        if (iparticle>-1) {
          assert( iparticle<(int)spc->p.size() && "Trial particle out of range");
          dprot = atom[spc->p[iparticle].id ].dprot;
          if (dprot<1e-6)
            dprot = base::genericdp;

          Point u;
          u.ranunit(slump);
          rot.setAxis(spc->geo, Point(0,0,0), u, dprot* slump.half() );
          spc->trial[iparticle].rotate(rot);
        }
      }

    template<class Tspace>
      void AtomicRotation<Tspace>::_acceptMove() {
        sqrmap[ spc->p[iparticle].id ] += pow(dprot*180/pc::pi, 2);
        accmap[ spc->p[iparticle].id ] += 1;
        spc->p[iparticle] = spc->trial[iparticle];
      }

    template<class Tspace>
      void AtomicRotation<Tspace>::_rejectMove() {
        spc->trial[iparticle] = spc->p[iparticle];
        sqrmap[ spc->p[iparticle].id ] += 0;
        accmap[ spc->p[iparticle].id ] += 0;
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
            << setw(l+7) << bracket("d"+theta+squared)+"/"+degrees
            << rootof+bracket("d"+theta+squared)+"/"+degrees << endl;
          for (auto m : sqrmap) {
            auto id=m.first;
            o << indent(SUBSUB) << std::left << setw(7) << atom[id].name
              << setw(l-6) << ( (atom[id].dprot<1e-6) ? genericdp : atom[id].dprot*180/pc::pi);
            o.precision(3);
            o << setw(l) << accmap[id].avg()*100
              << setw(l) << sqrmap[id].avg()
              << setw(l) << sqrt(sqrmap[id].avg()) << endl;
          }
        }
        return o.str();
      }

    /**
     * @brief Combined rotation and rotation of groups
     *
     * This will translate and rotate groups and collect averages based on group name.
     * See constructor for usage.
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
          using base::jsondir;
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

          TranslateRotate(Energy::Energybase<Tspace>&, Tspace&, Tmjson&,
              string="moltransrot");

          void setGroup(Group&); //!< Select Group to move
          bool groupWiseEnergy;  //!< Attempt to evaluate energy over groups from vector in Space (default=false)
          std::map<string,Point> directions; //!< Specify special group translation directions (default: x=y=z=1)
      };

    /**
     * @brief Constructor
     *
     * The default JSON entry is read from section `moltransrot`
     * with each element being the molecule name with the following
     * values:
     *
     * Value      | Description
     * :--------- | :-------------------------------------------------------------
     * `dir`      | Move directions (default: "1 1 1" = xyz)
     * `permol`   | Repeat move for each molecule in system (default: true) 
     * `prob`     | Probability of performing the move (default: 1)
     * `dp`       | Translational displacement parameter (angstrom, default: 0)
     * `dprot`    | Angular displacement parameter (radians, default: 0)
     *
     * Example:
     *
     *     moltransrot {
     *       "water"   : { "dp":0.5, "dprot":0.5 },
     *       "polymer" : { ... }
     *     } 
     *
     * Atomic displacement parameters are read from `AtomData`.
     *
     * @param e Energy function
     * @param s Space
     * @param j JSON object - typically `moves`.
     * @param sec JSON section -- default `moltransrot`.
     */
    template<class Tspace>
      TranslateRotate<Tspace>::TranslateRotate(
          Energy::Energybase<Tspace> &e,
          Tspace &s, Tmjson &j, string sec ) : base( e, s ) {

        base::title="Group Rotation/Translation";
        base::w=30;
        igroup=nullptr;
        groupWiseEnergy=false;
        base::jsondir = "moves/"+sec; // temp. fix to be removed (for _test())

        auto m = j["moves"][sec];
        base::fillMolList( m );// find molecules to be moved

        for (auto &i : this->mollist) { // loop over molecules to be moved
          string molname = spc->molList()[ i.first ].name;
          i.second.dp1 = m[molname]["dp"] | 0.0;
          i.second.dp2 = m[molname]["dprot"] | 0.0;
          if (i.second.dp2>4*pc::pi)    // no need to rotate more than
            i.second.dp2=4*pc::pi;      // +/- 2 pi.
         }
      }

    template<class Tspace>
      void TranslateRotate<Tspace>::setGroup(Group &g) {
        assert( this->mollist.empty() && "Use either JSON data or setGroup");
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

        // if `mollist` has data, favor this over `setGroup()`
        // Note that `currentMolId` is set by Movebase::move()
        if ( ! this->mollist.empty() ) {
          auto gvec = spc->findMolecules( this->currentMolId );
          assert( !gvec.empty() );
          igroup = *slump.element( gvec.begin(), gvec.end() );
          assert( ! igroup->empty() );
          auto it = this->mollist.find( this->currentMolId );
          if ( it != this->mollist.end() ) {
            dp_trans = it->second.dp1;
            dp_rot = it->second.dp2;
            dir = it->second.dir;
          }
        }

        assert(igroup!=nullptr);
        Point p;
        if (dp_rot>1e-6) {
          p.ranunit(slump);             // random unit vector
          p=igroup->cm+p;                    // set endpoint for rotation
          angle=dp_rot* slump.half();
          igroup->rotate(*spc, p, angle);
        }
        if (dp_trans>1e-6) {
          p.x()=dir.x() * dp_trans * slump.half();
          p.y()=dir.y() * dp_trans * slump.half();
          p.z()=dir.z() * dp_trans * slump.half();
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
          if ( spc->geo.collision( spc->trial[i], spc->trial[i].radius, Geometry::Geometrybase::BOUNDARY ) )
            return pc::infty;

        double unew = pot->external(spc->trial) + pot->penalty(spc->trial) + pot->g_external(spc->trial, *igroup);
        if (unew==pc::infty)
          return pc::infty;       // early rejection
        double uold = pot->external(spc->p) + pot->penalty(spc->p) + pot->g_external(spc->p, *igroup);

/*#ifdef ENABLE_MPI
        if (base::mpiPtr!=nullptr) {
          double du=0;
          auto s = Faunus::MPI::splitEven(*base::mpiPtr, spc->groupList().size());
          for (auto i=s.first; i<=s.second; ++i) {
            auto gi=spc->groupList()[i];
            if (gi!=igroup)
              du += pot->g2g(spc->trial, *gi, *igroup) - pot->g2g(spc->p, *gi, *igroup);
          }
          return (unew-uold) + Faunus::MPI::reduceDouble(*base::mpiPtr, du);
        }
#endif*/

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
          t(jsondir+idtrim+"acceptance", accmap[id].avg()*100);
          t(jsondir+idtrim+"dRot", sqrt(sqrmap_r[id].avg()));
          t(jsondir+idtrim+"dTrans", sqrt(sqrmap_t[id].avg()));
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
                  p.ranunit(slump);        // random unit vector
                  p=g->cm+p;                    // set endpoint for rotation
                  double angle=base::dp_rot* slump.half();
                  g->rotate(*base::spc, p, angle);
                  angle2[g->name] += pow(angle*180/pc::pi, 2); // sum angular movement^2
                }
                if (base::dp_trans>1e-6) {
                  p.ranunit(slump);
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
            if (!pairlist.empty() && base::mpiPtr!=nullptr) {

              // group <-> group
              auto s = Faunus::MPI::splitEven(*base::mpiPtr, pairlist.size());
              for (size_t i=s.first; i<=s.second; ++i)
                du += base::pot->g2g(base::spc->trial,*pairlist[i].first,*pairlist[i].second)
                  - base::pot->g2g(base::spc->p,*pairlist[i].first,*pairlist[i].second);

              // group <-> external potential
              s = Faunus::MPI::splitEven(*base::mpiPtr, base::spc->groupList().size());
              for (size_t i=s.first; i<=s.second; ++i) {
                auto gi=base::spc->groupList()[i];
                du += base::pot->g_external(base::spc->trial, *gi) - base::pot->g_external(base::spc->p, *gi);
              }

              return Faunus::MPI::reduceDouble(*base::mpiPtr, du);
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
     * following keywords in the section `moves/moltransrot`,
     *
     * Keyword         | Description
     * :---------------| :----------------
     * `clusterradius` | Surface threshold from mobile ion to particle in group (angstrom)
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
          typedef typename Tspace::ParticleVector Tpvec;
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
          virtual double ClusterProbability(Tpvec&,int); //!< Probability that particle index belongs to cluster
        public:
          using base::spc;
          TranslateRotateCluster(InputMap&, Energy::Energybase<Tspace>&, Tspace&, string="");
          virtual ~TranslateRotateCluster();
          void setMobile(Group&);  //!< Pool of atomic species to move with the main group
          double threshold;        //!< Distance between particles to define a cluster
      };

    template<class Tspace>
      TranslateRotateCluster<Tspace>::TranslateRotateCluster(InputMap &in,Energy::Energybase<Tspace> &e, Tspace &s, string pfx) : base(in,e,s,pfx) {
        base::title="Cluster "+base::title;
        base::cite="doi:10/cj9gnn";
        threshold = in.get( "clusterradius", 0.0 );
        gmobile=nullptr;
      }

    template<class Tspace>
      TranslateRotateCluster<Tspace>::~TranslateRotateCluster() {}

    template<class Tspace>
      void TranslateRotateCluster<Tspace>::setMobile(Group &g) {
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
          if (ClusterProbability(spc->p, i) > slump() )
            cindex.push_back(i); // generate cluster list
        avgsize += cindex.size();

        // rotation
        if (dp_rot>1e-6) {
          base::angle=dp_rot* slump.half();
          p.ranunit(slump);
          p=igroup->cm+p; // set endpoint for rotation
          igroup->rotate(*spc, p, base::angle);
          vrot.setAxis(spc->geo, igroup->cm, p, base::angle); // rot. around line CM->p
          for (auto i : cindex)
            spc->trial[i] = vrot(spc->trial[i]); // rotate
        }

        // translation
        if (dp_trans>1e-6) {
          p.x()=dir.x() * dp_trans * slump.half();
          p.y()=dir.y() * dp_trans * slump.half();
          p.z()=dir.z() * dp_trans * slump.half();
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
          if ( spc->geo.collision( spc->trial[i], spc->trial[i].radius, Geometry::Geometrybase::BOUNDARY ) )
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
      double TranslateRotateCluster<Tspace>::ClusterProbability(Tpvec &p, int i) {
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
          double ClusterProbability(typename base::Tpvec &p, int i) FOVERRIDE { return 1; }
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
          using base::jsondir;
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
        base::runfraction=in.get<double>("_runfraction", 1.0);
        skipEnergyUpdate=in.get<bool>("_skipenergy", false);
        dp=in.get<double>("_dp", 0);
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
        ip.x()*= slump.half();
        ip.y()*= slump.half();
        ip.z()*= slump.half();

        int f= slump()*remaining.size();
        moved.push_back(remaining[f]);
        remaining.erase(remaining.begin()+f);    // Pick first index in m to move

        for (size_t i=0; i<moved.size(); i++) {
          g[moved[i]]->translate(*spc, ip);
          for (size_t j=0; j<remaining.size(); j++) {
            double uo=pot->g2g(spc->p,     *g[moved[i]], *g[remaining[j]]);
            double un=pot->g2g(spc->trial, *g[moved[i]], *g[remaining[j]]);
            double udiff=un-uo;
            if (slump() < (1.-std::exp(-udiff)) ) {
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
          std::map<int,int> _minlen, _maxlen; 
          using base::spc;
          using base::pot;
          using base::w;
          using base::jsondir;
          Group* gPtr;       //!< Pointer to group where move is to be performed. Set by setGroup().
          double dp;         //!< Rotational displacement parameter
          double angle;      //!< Current rotation angle
          vector<int> index; //!< Index of particles to rotate
          //Geometry::VectorRotate vrot;
          Geometry::QuaternionRotate vrot;
          AcceptanceMap<string> accmap;
        public:
          CrankShaft(Energy::Energybase<Tspace>&,
              Tspace&, Tmjson&, string="crankshaft");
          virtual ~CrankShaft();
          void setGroup(Group&); //!< Select Group to of the polymer to move
          int minlen;            //!< Minimum number of particles to rotate (default = 1)
          int maxlen;            //!< Maximin number of particles to rotate (default = 10)
      };

    /**
     * The section `moves/crankshaft` is searched for:
     *
     * Key      | Description
     * -------- | -------------------------------------------
     * `minlen` | Minimum number of particles to rotate (default: 1)
     * `maxlen` | Maximum number of particles to rotate (default: 4)
     * `dp`     | Rotational displacement parameter (radians)
     */
    template<class Tspace>
      CrankShaft<Tspace>::CrankShaft(Energy::Energybase<Tspace> &e,
          Tspace &s, Tmjson &j, string sec) : base(e,s) {
        base::title = "CrankShaft";
        base::jsondir = "moves/"+sec;
        w=30;
        gPtr=nullptr;

        auto m = j["moves"][sec];
        base::fillMolList( m );
        for (auto &i : this->mollist) {
          string name = spc->molList()[ i.first ].name;
          i.second.dp1       = m[name]["dp"] | 0.0;
          _minlen[ i.first ] = m[name]["minlen"] | 1.0;
          _maxlen[ i.first ] = m[name]["maxlen"] | 4.0;
        }
      }

    template<class Tspace>
      CrankShaft<Tspace>::~CrankShaft() {}

    template<class Tspace>
      void CrankShaft<Tspace>::_trialMove() {

        if ( ! this->mollist.empty() ) {
          auto gvec = spc->findMolecules( this->currentMolId );
          assert( !gvec.empty() );
          gPtr = *slump.element( gvec.begin(), gvec.end() );
          assert( ! gPtr->empty() );
          dp = this->mollist[ this->currentMolId ].dp1;
          minlen = _minlen[ this->currentMolId ];
          maxlen = _maxlen[ this->currentMolId ];
        }

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
          if ( spc->geo.collision( spc->trial[i], spc->trial[i].radius, Geometry::Geometrybase::BOUNDARY ) )
            return pc::infty;
        du=pot->g_internal(spc->trial, *gPtr) - pot->g_internal(spc->p, *gPtr);
        for (auto i : index)
          du+=pot->i_external(spc->trial,i) - pot->i_external(spc->p,i);
        for (auto g : spc->groupList())
          if (g!=gPtr)
            du+=pot->g2g(spc->trial, *g, *gPtr) - pot->g2g(spc->p, *g, *gPtr);
        du+=pot->external(spc->trial) + pot->penalty(spc->trial) 
          - pot->external(spc->p) - pot->penalty(spc->p);
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

        angle = dp* slump.half();  // random angle
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
        accmap._test(t, jsondir);
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
          Pivot(Energy::Energybase<Tspace>&, Tspace&, Tmjson&, string="pivot");
      };

    template<class Tspace>
      Pivot<Tspace>::Pivot( Energy::Energybase<Tspace> &e, Tspace &s, Tmjson &j, string sec)
      : base( e, s, j, sec ) {
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

          if (slump.half() > 0)
            for (int i=end+1; i<=gPtr->back(); i++)
              index.push_back(i);
          else
            for (int i=gPtr->front(); i<end; i++)
              index.push_back(i);
        }
        base::angle = base::dp* slump.half();
        base::vrot.setAxis(spc->geo, spc->p[beg], spc->p[end], base::angle );
        return true;
      }

    /**
     * @brief Reptation move for linear polymers
     *
     * This will perform a reptation move of a linear, non-uniform polymer chain.
     * During construction, the json object is searched, molecule-wise, for the following keywords
     * in the section `[moves][reptation]`:
     *
     * Key           | Description
     * :------------ | :---------------------------------------------------------------------------
     * `prob`        | Probability to perform a move (defaults=1)
     * `bondlength`  | The bond length while moving head groups. Use -1 to use existing bondlength.
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
          using base::jsondir;
        public:
          Reptation(Energy::Energybase<Tspace>&, Tspace&, Tmjson&, string="reptation");
          void setGroup(Group&); //!< Select Group to move
      };

    template<class Tspace>
      Reptation<Tspace>::Reptation(Energy::Energybase<Tspace> &e,
          Tspace &s, Tmjson &j, string sec) : base(e,s) {

        base::jsondir = "moves/reptation";
        base::title="Linear Polymer Reptation";
        gPtr=nullptr;

        auto m = j["moves"][sec];
        base::fillMolList( m );         // find molecules to be moved
        for (auto &i : this->mollist) { // loop over molecules to be moved
          string molname = spc->molList()[ i.first ].name;
          i.second.dp1 = m[molname]["bondlength"] | -1.0;
         }
      }

    template<class Tspace>
      void Reptation<Tspace>::_test(UnitTest &t) {
        accmap._test(t, jsondir);
      }

    template<class Tspace>
      void Reptation<Tspace>::_trialMove() {

        gPtr=nullptr;
        if ( ! this->mollist.empty() ) {
          auto gvec = spc->findMolecules( this->currentMolId );
          if ( !gvec.empty() ) {
            gPtr = *slump.element( gvec.begin(), gvec.end() );
            bondlength = this->mollist[ this->currentMolId ].dp1;
          }
        }

        if (gPtr==nullptr)
          throw std::runtime_error("Molecule "+gPtr->name+" not found in space");
        if ( gPtr->size() < 2 )
          throw std::runtime_error("Molecule "+gPtr->name+" too short for reptation.");

        int first, second; // "first" is end point, "second" is the neighbor
        if (slump.half()>0) {
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
        u.ranunit(slump);                          // generate random unit vector
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
          if ( spc->geo.collision( spc->trial[i], spc->trial[i].radius, Geometry::Geometrybase::BOUNDARY ) )
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
     * @brief Isobaric volume move
     *
     * @details This class will perform a volume displacement and scale atomic
     * as well as molecular groups as long as these are known to Space -
     * see Space.enroll().
     * The json object is scanned for the following keys in `moves/isobaric`:
     *
     * Key     | Description
     * :-------| :-----------------------------
     * `dV`    | Volume displacement parameter
     * `P`     | Pressure [mM]
     * `prob`  | Runfraction [default=1]
     *
     * Note that new volumes are generated according to
     * \f$ V^{\prime} = \exp\left ( \log V \pm \delta dp \right ) \f$
     * where \f$\delta\f$ is a random number between zero and one half.
     */
    template<class Tspace>
      class Isobaric : public Movebase<Tspace> {
        private:
          typedef Movebase<Tspace> base;
        protected:
          string _info();
          void _test(UnitTest&);
          void _trialMove();
          void _acceptMove();
          void _rejectMove();
          template<class Tpvec> double _energy(const Tpvec&);
          double _energyChange();
          using base::spc;
          using base::pot;
          using base::w;
          double P; //!< Pressure
          double dp; //!< Volume displacement parameter
          double oldval;
          double newval;
          Point oldlen;
          Point newlen;
          Average<double> msd;       //!< Mean squared volume displacement
          Average<double> val;          //!< Average volume
          Average<double> rval;         //!< Average 1/volume
        public:
          template<typename Tenergy>
            Isobaric(Tenergy&, Tspace&, Tmjson&, string="isobaric");
      };

    template<class Tspace>
      template<class Tenergy> Isobaric<Tspace>::Isobaric(
          Tenergy &e, Tspace &s, Tmjson &j, string sec) : base(e,s) {

        this->title="Isobaric Volume Fluctuations";
        this->jsondir = "moves/isobaric";
        this->w=30;
        auto m = j["moves"][sec];
        dp = m["dp"] | 0.0;
        P  = ( m["pressure"] | 0.0 ) * 1.0_mM;
        base::runfraction = m["prob"] | 1.0;
        if (dp<1e-6)
          base::runfraction=0;

        auto t = e.tuple(); // tuple w. pointers to all energy terms
        auto ptr = TupleFindType::get< Energy::ExternalPressure<Tspace>* >( t );
        if ( ptr != nullptr )
          (*ptr)->setPressure( P ); 
        else {
          std::cerr << "Error: Volume move requires pressure term in Hamiltonian." << endl;
          exit(1);
        }
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
        o << pad(SUB,w, "Displacement parameter") << dp << endl
          << pad(SUB,w, "Number of molecules")
          << N << " (" <<Nmol<< " molecular + " <<Natom<< " atomic)\n"
          << pad(SUB,w, "Pressure")
          << P*tomM << " mM = " << Pascal << " Pa = "
          << Pascal/0.980665e5 << " atm\n"
          << pad(SUB,w, "Temperature") << pc::T() << " K\n";
        if (base::cnt>0) {
          char l=14;
          o << pad(SUB,w, "Mean displacement")
            << cuberoot+rootof+bracket("dp"+squared)
            << " = " << pow(msd.avg(), 1/6.) << _angstrom << endl
            << pad(SUB,w, "Osmotic coefficient") << P / (N*rval.avg()) << endl
            << endl
            << indent(SUBSUB) << std::right << setw(10) << ""
            << setw(l+5) << bracket("V")
            << setw(l+8) << cuberoot+bracket("V")
            << setw(l+8) << bracket("1/V")
            << setw(l+8) << bracket("N/V") << endl
            << indent(SUB) << setw(10) << "Averages"
            << setw(l) << val.avg() << _angstrom + cubed
            << setw(l) << std::cbrt(val.avg()) << _angstrom
            << setw(l) << rval.avg() << " 1/" + _angstrom + cubed
            << setw(l) << N*rval.avg()*tomM << " mM\n";
        }
        return o.str();
      }

    template<class Tspace>
      void Isobaric<Tspace>::_test(UnitTest &t) {
        t(this->jsondir+"_averageSideLength", std::cbrt(val.avg()) );
        t(this->jsondir+"_MSQDisplacement", pow(msd.avg(), 1/6.) );
      }

    template<class Tspace>
      void Isobaric<Tspace>::_trialMove() {
        assert(spc->groupList().size()>0
            && "Space has empty group vector - NPT move not possible.");
        oldval = spc->geo.getVolume();
        oldlen = newlen = spc->geo.len;
        newval = std::exp( std::log(oldval) + slump.half()*dp );
        Point s = Point(1,1,1);
        double xyz = cbrt(newval/oldval);
        double xy = sqrt(newval/oldval);
        newlen.scale(spc->geo,s,xyz,xy);
        for (auto g : spc->groupList()) {
          g->setMassCenter(*spc);
          g->scale(*spc, s, xyz, xy); // scale trial coordinates to new volume
        }
      }

    template<class Tspace>
      void Isobaric<Tspace>::_acceptMove() {
        val += newval;
        msd += pow( oldval-newval, 2 );
        rval += 1./newval;
        spc->geo.setlen(newlen);
        pot->setSpace(*spc);
        for (auto g : spc->groupList() )
          g->accept(*spc);
      }

    template<class Tspace>
      void Isobaric<Tspace>::_rejectMove() {
        msd += 0;
        val += oldval;
        rval += 1./oldval;
        spc->geo.setlen(oldlen);
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
        if (dp<1e-6)
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
        return u + pot->external(p) + pot->penalty(p);
      }

    /**
     * @todo Early rejection could be implemented
     *       - not relevant for geometries with periodicity, though.
     */
    template<class Tspace>
      double Isobaric<Tspace>::_energyChange() {
        double uold = _energy(spc->p);
        spc->geo.setlen(newlen);
        pot->setSpace(*spc); // potential must know about volume, too

        // In spherical geometries, molecules may collide with
        // cell boundary upon mass center scaling:
        for (auto g : spc->groupList())
          for (auto i : *g)
            if ( spc->geo.collision( spc->trial[i], spc->trial[i].radius, Geometry::Geometrybase::BOUNDARY ) )
              return pc::infty;
        double unew = _energy(spc->trial);
        return unew-uold;
      }

    /**
     * @brief Isochoric scaling move in Cuboid geometry
     *
     * @details This class will expand/shrink along the z-axis
     * and shrink/expand in the xy-plane atomic as well as molecular groups 
     * as long as these are known to Space - see Space.enroll().
     * The json object class is scanned for the following keys:
     *
     * Key                | Description
     * :----------------- | :-----------------------------
     * `nvt_dz`           | Length displacement parameter
     * `nvt_runfraction`  | Runfraction [default=1]
     *
     */
    template<class Tspace>
       class Isochoric : public Isobaric<Tspace> {
         protected:
           typedef Isobaric<Tspace> base;
           using base::spc;
           using base::pot;
           using base::w;
           using base::dp;
           using base::oldval;
           using base::newval;
           using base::oldlen;
           using base::newlen;
           using base::msd;
           using base::val;
           void _trialMove();
           string _info();
         public:
           template<typename Tenergy>
             Isochoric(Tenergy&, Tspace&, Tmjson&, string="isochoric");
       };

    template<class Tspace>
      template<class Tenergy> Isochoric<Tspace>::Isochoric(Tenergy &e, Tspace &s, Tmjson &j, string sec) : base(e,s,j,"isochoric") {
        this->title="Isochoric Side Lengths Fluctuations";
        this->w=30;
        auto m = j["moves"][sec];
        dp = m["dp"] | 0.0;
        base::runfraction = m["prob"] | 1.0;
        if (dp<1e-6)
          base::runfraction=0;
      }

    template<class Tspace>
      string Isochoric<Tspace>::_info() {
        using namespace textio;
        std::ostringstream o;
        o << pad(SUB,w, "Displacement parameter") << dp << endl
          << pad(SUB,w, "Temperature") << pc::T() << " K\n";
        if (base::cnt>0) {
          char l=14;
          o << pad(SUB,w, "Mean displacement")
            << rootof+bracket("dz"+squared)
            << " = " << pow(msd.avg(), 1/2.) << _angstrom << endl
            << endl
            << indent(SUBSUB) << std::right << setw(10) << ""
            << setw(l+5) << bracket("Lz") << endl
            << indent(SUB) << setw(10) << "Averages"
            << setw(l) << val.avg() << _angstrom + cubed;
        }
        return o.str();
      }

    template<class Tspace>
      void Isochoric<Tspace>::_trialMove() {
        assert(spc->groupList().size()>0
            && "Space has empty group vector - Isochoric scaling move not possible.");
        oldlen = spc->geo.len;
        newlen = oldlen;
        oldval = spc->geo.len.z();
        newval = std::exp( std::log(oldval) + slump.half()*dp );
        //newval = oldval+ slump.half()*dp;
        Point s;
        s.z() = newval / oldval;
        s.x() = s.y() = 1 / std::sqrt(s.z());
        newlen.scale(spc->geo,s);
        for (auto g : spc->groupList()) {
          g->scale(*spc,s); // scale trial coordinates to new coordinates
        }
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
          typedef typename Tspace::ParticleType::Tid Tid;
          typedef typename Tspace::ParticleType Tparticle;
          Tspace* spc;
          class data {
            public:
              vector<Tindex> index;
              Tindex random();                  //!< Pick random particle index
          };
          std::map<Tid,data> map;
        public:
          AtomTracker(Tspace&);
          Tid randomAtomType() const;           //!< Select a random atomtype from the list
          bool insert(const Tparticle&, Tindex);//!< Insert particle into Space and track position
          bool erase(Tindex);                   //!< Delete particle from Space at specific particle index
          data& operator[] (Tid);               //!< Access operator to atomtype data
          void clear();                         //!< Clear all atom lists (does not touch Space)
          bool empty();                         //!< Test if atom list is empty
      };

    template<class Tspace>
      bool AtomTracker<Tspace>::empty() {
        return map.empty();
      }

    template<class Tspace>
      typename AtomTracker<Tspace>::Tid AtomTracker<Tspace>::randomAtomType() const {
        assert(!map.empty() && "No atom types have been added yet");
        vector<Tid> vid;
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
      typename AtomTracker<Tspace>::data& AtomTracker<Tspace>::operator[](Tid id) {
        return map[id];
      }

    /**
     * This will insert a particle into Space and at the same time make sure
     * that all other particles are correctly tracked.
     */
    template<class Tspace>
      bool AtomTracker<Tspace>::insert(const Tparticle &a, Tindex index) {
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
     *
     * This will do GC moves of salt pairs, automatically combined
     * from their valencies. JSON input:
     *
     * ~~~~
     * "moves" : {
     *   "atomgc" : { "molecule":"mysalt", "prob":1.0 }
     * }
     * ~~~~
     *
     * where `mysalt` must be an atomic molecule. Only atom types with
     * non-zero activities will be considered.
     *
     * @date Lund 2010-2011
     * @warning Untested for asymmetric salt in this branch
     */
    template<class Tspace>
      class GrandCanonicalSalt : public Movebase<Tspace> {
        protected:
          typedef Movebase<Tspace> base;
          typedef typename Tspace::ParticleType Tparticle;
          typedef typename Tspace::ParticleVector Tpvec;
          typedef typename Tparticle::Tid Tid;
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
            Tparticle p;
            double chempot;       // chemical potential log(1/A3)
            Average<double> rho;  // average density
          };
          std::map<typename Tparticle::Tid,ionprop> map;
          void randomIonPair(typename Tparticle::Tid&, typename Tparticle::Tid&);  // Generate random ion pair
          Tpvec trial_insert;
          vector<int> trial_delete;
          typename Tparticle::Tid ida, idb;     // particle id's of current salt pair (a=cation, b=anion)

          Group* saltPtr;  // GC ions *must* be in this group

          // unit testing
          void _test(UnitTest &t) {
            for (auto &m : map) {
              auto s=base::jsondir+"_"+atom[m.first].name;
              t(s+"_activity", atom[m.first].activity);
              t(s+"_conc", m.second.rho.avg()/pc::Nav/1e-27);
            }
          }

        public:
          GrandCanonicalSalt(Energy::Energybase<Tspace>&, Tspace&, Tmjson&, string="atomgc");
      };

    template<class Tspace>
      GrandCanonicalSalt<Tspace>::GrandCanonicalSalt(
          Energy::Energybase<Tspace> &e, Tspace &s, Tmjson &j,
          string sec) : base(e,s), tracker(s) {

          base::title="Grand Canonical Salt";
          base::useAlternateReturnEnergy=true;
          base::jsondir = "moves/" + sec;

          string saltname = j["moves"][sec]["molecule"] | string();
          auto v = spc->findMolecules( saltname );
          if ( v.empty() ) { // insert if no atomic species found
            auto it = spc->molList().find( saltname ); 
            if ( it != spc->molList().end() )
              saltPtr = spc->insert( it->id, it->getRandomConformation(spc->geo, spc->p) );
          } else {
            if ( v.size() != 1 )
              throw std::runtime_error( "Number of atomic GC groups must be exactly ONE." );
            if ( v.front()->isMolecular() )
              throw std::runtime_error( "Atomic GC group must be atomic.");
            saltPtr=v.front();
          }
          add(*saltPtr);
        }

    template<class Tspace>
      void GrandCanonicalSalt<Tspace>::add( Group &g ) {
        assert( g.isAtomic() && "Salt group must be atomic" );
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
      void GrandCanonicalSalt<Tspace>::randomIonPair(typename Tparticle::Tid &id_cation, typename Tparticle::Tid &id_anion) {
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
          std::cerr << "!! No salt to insert or delete !!" << endl;
          exit(1);
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
          Tid id=m.first;
          o.precision(5);
          o << setw(4) << "" << setw(s) << atom[id].name
            << setw(s) << atom[id].activity << setw(s) << m.second.rho.avg()/pc::Nav/1e-27
            << setw(s) << atom[id].activity / (m.second.rho.avg()/pc::Nav/1e-27)
            << setw(s) << m.second.rho.avg()*spc->geo.getVolume()
            << "\n";
        }
        return o.str();
      }

    /**
     * @brief Grand Canonical Titration derived from Grand Canonical Salt
     *
     * @date Lund 2015
     * @warning Untested
     *
     * TODO: contains lots of redundant code from SwapMove, could inherit from there 
     * as well
     */
    template<class Tspace>
      class GrandCanonicalTitration : public GrandCanonicalSalt<Tspace> {
        protected:
          typedef GrandCanonicalSalt<Tspace> base;
          typedef typename Tspace::ParticleType Tparticle;
          typedef typename Tspace::ParticleVector Tpvec;
          typedef typename Tparticle::Tid Tid;
          Energy::EquilibriumEnergy<Tspace>* eqpot;
          using base::spc;
          using base::pot;
          using base::w;
          void _trialMove();
          void _acceptMove();
          void _rejectMove();
          string _info();
          double _energyChange();
          void add(Group&) {};       // scan group for ions with non-zero activities

          unsigned long int cnt_tit, cnt_salt, cnt_tit_acc, cnt_salt_acc;
          Tid pid;     // particle id's of current salt pair (a=cation, b=anion)
	  int  N, isite=-1;                    // switch with a built in message
          int k;
          bool protonation;          // if yes, the process shuld lead to protonation
          bool gcyes;

          std::map<int, Average<double> > accmap; //!< Site acceptance map
          std::map<int, std::map<int, Average<double> >> molCharge;

          void updateMolCharge( int pindex ) {
            auto g = spc->findGroup( pindex );
            molCharge[ g->molId ][ pindex - g->front() ] += spc->p[pindex].charge; 
          }

        public:
          template<class Tenergy>
            GrandCanonicalTitration(Tenergy&, Tspace&, Tmjson&, string="gctit");

          template<class Tpvec>
            int findSites( Tpvec &p ){
              accmap.clear();
              return eqpot->findSites(p);
            }

          template<class Tpvec>
            void applycharges( const Tpvec &p ) {
              eqpot->eq.applycharges(p);
            }  
      };

      template<class Tspace>
        template<class Tenergy>

        GrandCanonicalTitration<Tspace>::GrandCanonicalTitration(
            Tenergy &e,
            Tspace &s,
            Tmjson &j,
            string sec) : base(e,s,j,sec) {

            base::title += " Titration";
            base::useAlternateReturnEnergy=true;
            auto t = e.tuple();
            auto ptr = TupleFindType::get< Energy::EquilibriumEnergy<Tspace>* >( t );
            if ( ptr != nullptr )
              eqpot = *ptr;
            else {
              std::cerr << "Error: Equilibrium energy required in Hamiltonian\
                for Grand Canonical Titration moves." << endl;
              exit(1);
            }
            eqpot->eq = Energy::EquilibriumController( j );
            findSites( spc->p );
            cnt_tit = cnt_salt = cnt_tit_acc = cnt_salt_acc = 0;

            /* Sync particle charges with `AtomMap` */
            for (auto i : eqpot->eq.sites)
              spc->trial[i].charge = spc->p[i].charge = atom[ spc->p[i].id ].charge;

            // neutralise system, if needed, using GC ions
            if ( j["moves"][sec]["neutralize"] | true ) {
              double Z = netCharge( s.p, Group(0,s.p.size()-1) );
              if ( fabs(Z) > 1e-9 ) {
                Tid id;
                double z = 0;
                int maxtry = 1000;
                cout << "# Neutralizing system with GC ions. Initial charge = "
                  << Z << "e." << endl;
                do {
                  id = base::tracker.randomAtomType();
                  z = atom[id].charge;
                  if ( --maxtry==0 ) {
                    std::cerr << "Error: Failed to find GC ions capable of "
                      << "neutralizing system." << endl;
                    exit(1);
                  }
                } while (
                    ( (z<0 && Z>0) || (Z<0 && z>0) )
                    && ( fabs( fmod(Z,z) ) < 1e-9 ) );

                int n = round(-Z/z);
                assert( n>0 && fabs(n*z+Z) < 1e-9 );

                typename Tspace::ParticleType a;
                a = atom[id];
                for ( int i=0; i<n; i++ ) {
                  s.geo.randompos(a);
                  base::tracker.insert( a, base::saltPtr->back() );
                }
                Z = netCharge( s.p, Group(0,s.p.size()-1) );
                cout << "# Final charge = " << Z << "e." << endl;
                assert( fabs(Z)<1e-9 ) ;
              }
            }
        }

      template<class Tspace> 
        void GrandCanonicalTitration<Tspace>::_trialMove()  { 
          gcyes=false;
          int switcher = slump.range(0,1);
          if (eqpot->eq.number_of_sites()==0){ // If no sites associated with processess
            gcyes=true, switcher=0;            // fall back to plain gc
          } 
          switch ( switcher ) {
            case 0:   // Go for the inheritance
              cnt_salt++;
              gcyes=true;
              base::_trialMove(); 
              break; 
            case 1:   // Brand new deal
              cnt_tit++;
              base::trial_insert.clear(); // First some cleaning in the attic
              base::trial_delete.clear();
              do {              // Pick a monovalent ion
                pid = base::tracker.randomAtomType();
              } while (atom[pid].charge*atom[pid].charge != 1);
              if (!eqpot->eq.sites.empty()) {
                int i = slump.range( 0, eqpot->eq.sites.size()-1); // pick random site (local in *eq)
                isite = eqpot->eq.sites.at(i); // and corresponding particle (index in spc->p)
                do {
                  k = slump.range( 0, eqpot->eq.process.size()-1 );// pick random process..
                } while (!eqpot->eq.process[k].one_of_us( this->spc->p[isite].id )); //that match particle isite

                eqpot->eq.process[k].swap( this->spc->trial[isite] ); // change state and get intrinsic energy change
              }
              if ( !eqpot->eq.process[k].bound(this->spc->trial[isite].id) ) {// have action lead to deprotonation?
                protonation = false;
              } else {
                protonation = true;
              }
              N = -1;
              // The following section is hardcoded for monovalent salt
              if (base::map[pid].p.charge>0) {  // Determine weather cat-/anion
                N = 0;            // N==0 cation, N==1 anion
              } else if ( base::map[pid].p.charge<0 ) {      
                N = 1;
              } else {  
                std::cerr << " Error, something fails !"<<std::endl, exit(0);
              }
              int iIon=-1;
              if ( protonation == true ) {
                if ( N == 0) { // Protonation and deletion of cation
                  iIon=base::tracker[pid].random();
                  assert( pid==spc->p[iIon].id && "id mismatch");
                  base::trial_delete.push_back(iIon);
                } else if ( N == 1 ) { // protonation and addition of anion
                  base::trial_insert.push_back( base::map[pid].p );
                  base::spc->geo.randompos(base::trial_insert[0]);
                  assert( pid==base::trial_insert[0].id );
                } else {
                  std::cerr << " Process error !"<< std::endl; exit(1);
                }
              } else {
                if ( N == 0) { // Deprotonation and addition of cation
                  base::trial_insert.push_back( base::map[pid].p );
                  base::spc->geo.randompos(base::trial_insert[0]);
                  assert( pid==base::trial_insert[0].id );
                } else if ( N == 1 ) { // Deprotonation and deletion of anion
                  iIon=base::tracker[pid].random();
                  assert( pid==spc->p[iIon].id && "id mismatch");
                  base::trial_delete.push_back(iIon);
                } else {
                  std::cerr << " Process error !"<< std::endl; exit(1);
                }
              }
              break;
          }
        };

      template<class Tspace>
        double GrandCanonicalTitration<Tspace>::_energyChange() {
          if (gcyes == true) // Go about the old habbit
            return base::_energyChange();
          double idfactor=1;
          double uold=0, unew=0, V=spc->geo.getVolume();
          double potold=0, potnew=0; // energy change due to interactions
          double site_old=0, salt_old=0;
          double site_new=0, salt_new=0;
          potnew = pot->i_internal(spc->trial, isite);  // Intrinsic energies for new
          potold = pot->i_internal(spc->p, isite);      // and old state
          if (protonation == true && N==1) { // Protonate and ins. anion
            idfactor *= (base::tracker[pid].index.size()+1)/V;
            unew = log(idfactor) - base::map[pid].chempot;
            salt_new += pot->all2p(spc->trial, base::trial_insert[0]);

            site_new+=pot->i2all(spc->trial, isite);

            site_old+=pot->i2all(spc->p, isite);
          }
          if (protonation == true && N==0) { // Protonate and del. cat
            idfactor *= V/base::tracker[pid].index.size();
            unew = log(idfactor) + base::map[pid].chempot;

            salt_old+=pot->i2all(spc->p, base::trial_delete[0]);

            site_new+=pot->i2all(spc->trial, isite);
            site_new-=pot->i2i(spc->trial, base::trial_delete[0], isite);

            site_old+=pot->i2all(spc->p, isite);
            site_old-=pot->i2i(spc->p, base::trial_delete[0], isite); //Subtracted from previuous double count
          }
          if (protonation == false && N==0) { // Deprotonate and ins.
            idfactor *= (base::tracker[pid].index.size()+1)/V;   // cation
            unew = log(idfactor) - base::map[pid].chempot;
            salt_new += pot->all2p(spc->trial, base::trial_insert[0]);
            site_new+=pot->i2all(spc->trial, isite);

            site_old+=pot->i2all(spc->p, isite);
          }
          if (protonation == false && N==1) { // Deprotonate and del.
            idfactor *= V/base::tracker[pid].index.size();       // anion
            unew = log(idfactor) + base::map[pid].chempot;

            salt_old+=pot->i2all(spc->p, base::trial_delete[0]);

            site_new+=pot->i2all(spc->trial, isite);
            site_new-=pot->i2i(spc->trial, base::trial_delete[0], isite);

            site_old+=pot->i2all(spc->p, isite);
            site_old-=pot->i2i(spc->p, base::trial_delete[0], isite);
          }
          unew+=potnew+salt_new+site_new;
          uold+=potold+salt_old+site_old;
          potnew+=salt_new+site_new;
          potold+=salt_old+site_old;
          base::alternateReturnEnergy=potnew-potold; // track only pot. energy
          return unew-uold;
        };

      template<class Tspace>
        void GrandCanonicalTitration<Tspace>::_acceptMove() { 
          if (gcyes==true) {
            base::_acceptMove();
            cnt_salt_acc++;
          } else{
            assert( spc->p[isite].id != spc->trial[isite].id );
            spc->p[isite] = spc->trial[isite];
            if ( !base::trial_insert.empty() ) {
              for (auto &p : base::trial_insert)
                base::tracker.insert(p, base::saltPtr->back());
            } else if ( !base::trial_delete.empty() ) {
              std::sort(base::trial_delete.rbegin(), base::trial_delete.rend()); //reverse sort
              for (auto i : base::trial_delete)
                base::tracker.erase(i);
            }
            double V = spc->geo.getVolume();
            base::map[pid].rho += base::tracker[pid].index.size() / V;
            accmap[isite] += 1;
            updateMolCharge( isite );
            cnt_tit_acc++;
          } 
        };

      template<class Tspace>
        void GrandCanonicalTitration<Tspace>::_rejectMove() { 
          if(gcyes==true) {
            base::_rejectMove();
          } else {
            assert( spc->p[isite].id != spc->trial[isite].id );
            spc->trial[isite] = spc->p[isite];
            accmap[isite] += 0;
            updateMolCharge( isite );
            double V = spc->geo.getVolume();
            base::map[pid].rho += base::tracker[pid].index.size() / V;
          }
        };

      template<class Tspace>
        string GrandCanonicalTitration<Tspace>::_info() {
          char s=10;
          using namespace textio;
          std::ostringstream o;
          o << pad(SUB,w,"Number of GC species") << endl << endl;
          o << setw(4) << "" << std::left
            << setw(s) << "Ion" << setw(s) << "activity"
            << setw(s+4) << bracket("c/M") << setw(s+6) << bracket( gamma+pm )
            << setw(s+4) << bracket("N") << "\n";
          for (auto &m : base::map) {
            Tid id=m.first;
            o.precision(5);
            o << setw(4) << "" << setw(s) << atom[id].name
              << setw(s) << atom[id].activity << setw(s) << m.second.rho.avg()/pc::Nav/1e-27
              << setw(s) << atom[id].activity / (m.second.rho.avg()/pc::Nav/1e-27)
              << setw(s) << m.second.rho.avg()*spc->geo.getVolume()
              << "\n";
          }
          for (auto &m : molCharge) {
            int molid = m.first;
            o << "\n" << indent(SUB) << "Molecule: " << spc->molList()[ molid ].name << "\n\n"
              << std::left << "    " << setw(8) << "index" << setw(12) << "name"
              << setw(12) << "Z" << "\n";
            for (auto &i : m.second)
              o << "    " << setw(8) << i.first
                << setw(12) << atom[ spc->molList()[molid].atoms[i.first] ].name
                << setw(12) << i.second << "\n"; 
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
            typedef typename Tspace::ParticleVector Tpvec;
            typedef Movebase<Tspace> base;
            using base::pot;
            using base::spc;
            using base::w;
            using base::runfraction;
            using base::mpiPtr;
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

            Faunus::MPI::FloatTransmitter ft;   //!< Class for transmitting floats over MPI
            Faunus::MPI::ParticleTransmitter<Tpvec> pt;//!< Class for transmitting particles over MPI

            typedef std::function<double(Tspace&, Energy::Energybase<Tspace>&, const Tpvec&)> Tenergyfunc;
            Tenergyfunc usys; //!< Defaults to Energy::systemEnergy but can be replaced!

          public:
            ParallelTempering(
                Energy::Energybase<Tspace>&, Tspace&, Tmjson&, MPI::MPIController&, string="temper");

            virtual ~ParallelTempering();

            void setCurrentEnergy(double); //!< Set energy before move (for increased speed)

            void setEnergyFunction( Tenergyfunc );
        };

      template<class Tspace>
        ParallelTempering<Tspace>::ParallelTempering(
            Energy::Energybase<Tspace> &e,
            Tspace &s,
            Tmjson &j,
            MPI::MPIController &mpi,
            string sec) : base( e, s ) {

          this->title   = "Parallel Tempering";
          this->jsondir = "moves/"+sec; // to be changed - compatibility w. tests
          this->mpiPtr  = &mpi;
          partner=-1;
          this->useAlternateReturnEnergy=true; //dont return dU from partner replica (=drift)
          this->runfraction = j["moves"][sec]["prob"] | 1.0;
          pt.recvExtra.resize(1);
          pt.sendExtra.resize(1);
          pt.setFormat( j["moves"][sec]["format"] | string("XYZQI") );

          setEnergyFunction(
              Energy::systemEnergy<Tspace,Energy::Energybase<Tspace>,Tpvec> );

          this->haveCurrentEnergy=false;
          assert( this->mpiPtr != nullptr );
        }

      template<class Tspace>
        ParallelTempering<Tspace>::~ParallelTempering() {}

      template<class Tspace>
        void ParallelTempering<Tspace>::setEnergyFunction( Tenergyfunc f ) {
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
          //std::advance(vi, slump.rand() % swappableParticles.size());
          //std::advance(vj, slump.rand() % swappableParticles.size());
          //ip=*(vi);
          //jp=*(vj);
          ip=*( slump.element(swappableParticles.begin(), swappableParticles.end() ) );
          jp=*( slump.element(swappableParticles.begin(), swappableParticles.end() ) );
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
       * @brief Flip-flop move of lipids in planar and cylindrical geometry
       *
       * Key                    | Description
       * :--------------------- | :-----------------------------
       * `flipflop_geometry`    | Geometry of the bilayer [planar(default) or cylindrical]
       * `flipflop_runfraction` | Runfraction [default=1]
       *
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
            using base::jsondir;
            void _trialMove();
            void _acceptMove();
            void _rejectMove();
            double _energyChange();
            string _info();
            typedef std::map<string, Average<double> > map_type;
            map_type accmap;   //!< Group particle acceptance map
            Group* igroup;
            Point* cntr;
            string geometry;
          public:
            FlipFlop(InputMap&, Energy::Energybase<Tspace>&, Tspace&, string="flipflop"); // if cylindrical geometry, string=cylinder
            void setGroup(Group&); //!< Select Group to move
            void setCenter(Point&); //!< Select Center of Mass of the bilayer
        };

      template<class Tspace>
        FlipFlop<Tspace>::FlipFlop(InputMap &in,Energy::Energybase<Tspace> &e, Tspace &s, string pfx) : base(e,s,pfx) {
          base::title="Group Flip-Flop Move";
          base::w=30;
          igroup=nullptr;
          cntr=nullptr;
          geometry = in( "geometry", string("planar") );
          this->runfraction=in( "prob", 1.0 );
        }

      template<class Tspace>
        void FlipFlop<Tspace>::setGroup(Group &g) {
          assert(g.isMolecular());
          igroup=&g;
        }

      template<class Tspace>
        void FlipFlop<Tspace>::setCenter(Point &center) {
          cntr=&center;
        }

      template<class Tspace>
        void FlipFlop<Tspace>::_trialMove() {
          assert(igroup!=nullptr);
          assert(cntr!=nullptr);
          Point startpoint=spc->p[igroup->back()];
          Point endpoint=*cntr;
          startpoint.z()=cntr->z();
          if (geometry.compare("cylindrical") == 0) { // MC move in case of cylindrical geometry
            startpoint=spc->p[igroup->back()];
            Point head=spc->p[igroup->front()];
            cntr->z()=head.z()=startpoint.z();
            Point dir = spc->geo.vdist(*cntr, startpoint)
              / sqrt(spc->geo.sqdist(*cntr, startpoint)) * 1.1*spc->p[igroup->back()].radius;
            if (spc->geo.sqdist(*cntr, startpoint) > spc->geo.sqdist(*cntr, head))
              startpoint.translate(spc->geo,-dir); // set startpoint for rotation
            else startpoint.translate(spc->geo, dir);
            double x1=cntr->x();
            double y1=cntr->y();
            double x2=startpoint.x();
            double y2=startpoint.y();
            endpoint.x()=x2+1; // rot endpoint of axis ⊥ to line connecting cm of cylinder with 2nd TL
            endpoint.y()=-(x2-x1)/(y2-y1)+y2;
            endpoint.z()=startpoint.z();
          }
          double angle=pc::pi; // MC move in case of planar geometry
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
            if ( spc->geo.collision( spc->trial[i], spc->trial[i].radius, Geometry::Geometrybase::BOUNDARY ) )
              return pc::infty;

          double unew = pot->external(spc->trial) + pot->g_external(spc->trial, *igroup);
          if (unew==pc::infty)
            return pc::infty;       // early rejection
          double uold = pot->external(spc->p) + pot->g_external(spc->p, *igroup);

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

      /**
       * @brief Grand Canonical Monte Carlo Move
       *
       * This is a general class for GCMC that can handle both
       * atomic and molecular species at constant chemical potential.
       *
       * @todo Currently tested only with rigid, molecular species.
       * @date Brno/Lund 2014
       * @author Lukas Sukenik and Mikael Lund
       */
      template<class Tspace, class base=Movebase<Tspace> >
        class GreenGC : public base {
          private:

            typedef typename Tspace::ParticleVector Tpvec;
            using base::spc;
            using base::pot;
            using base::w;

            Faunus::Tracker<Group*> molTrack;    // tracker for molecules
            Faunus::Tracker<int> atomTrack;      // tracker for atoms
            vector<Group*> molDel;               // groups to delete
            vector<int> atomDel;                 // atom index to delete
            MoleculeCombinationMap<Tpvec> comb;  // map of combinations to insert
            std::map<int,int> molcnt, atomcnt;   // id's and number of inserted/deleted mols and atoms
            std::multimap<int, Tpvec> pmap;      // coordinates of mols and atoms to be inserted
            unsigned int Ndeleted, Ninserted;    // Number of accepted deletions and insertions
            bool insertBool;                     // current status - either insert or delete
            typename MoleculeCombinationMap<Tpvec>::iterator it; // current combination

            /** @brief Perform an insertion or deletion trial move */
            void _trialMove() {

              // pick random combination and count mols and atoms
              base::alternateReturnEnergy=0;
              molcnt.clear();
              atomcnt.clear();
              it = comb.random();                 // random combination
              for ( auto id : it->molComb ) {     // loop over molecules in combination
                if ( spc->molecule[id].isAtomic() )
                  for ( auto i : spc->molecule[id].atoms )
                    atomcnt[i]++;                 // count number of atoms per type
                else
                  molcnt[id]++;                   // count number of molecules per type
              }

              insertBool = slump.range(0,1) == 1;

              // try delete move
              if (!insertBool) {
                molDel.clear();
                atomDel.clear();

                // find atom index and group pointers
                bool empty = false;           // true if too few atoms/mols are present
                for ( auto &a : atomcnt )     // (first=type, second=count)
                  if ( !atomTrack.find( a.first, a.second, atomDel ) )
                    empty = true;
                for ( auto &m : molcnt )
                  if ( !spc->molecule[m.first].isAtomic() )
                    if ( !molTrack.find( m.first, m.second, molDel ) )
                      empty = true;
                if (empty) {        // nothing left to delete
                  molDel.clear();
                  atomDel.clear();
                  pmap.clear();
                } else assert( !molDel.empty() || !atomDel.empty() );
              }

              // try insert move (nothing is actually inserter - just a proposed configuration)
              if (insertBool) {
                pmap.clear();
                for ( auto molid : it->molComb ) // loop over molecules in combination
                  pmap.insert(
                      { molid, spc->molecule[molid].getRandomConformation( base::spc->geo, base::spc->p ) } );
                assert( !pmap.empty() );
              }
            }

            // contribution from external chemical potential
            double externalEnergy() {
              double u    = 0;
              double V    = spc->geo.getVolume();
              int    bit  = ( insertBool ) ?  1 : 0;
              double sign = ( insertBool ) ?  1 : -1;

              for ( auto i : molcnt )                  // loop over molecule types
                if (!spc->molecule[i.first].isAtomic())
                  for (int n=0; n<i.second; n++)       // loop over n number of molecules
                    u += log(( molTrack.size(i.first) + bit ) / V) - spc->molecule[i.first].chemPot;

              for ( auto i : atomcnt )                 // loop over atom types
                for (int n=0; n<i.second; n++)         // loop over n number of atoms
                  u += log(( atomTrack.size(i.first) + bit ) / V) - atom[i.first].chemPot;

              return sign * u;
            }

            double _energyChange() {

              double u=0;         // change in potential energy (kT)
              double uinternal=0; // change in internal, molecular energy (kT)

              // energy if insertion move
              if ( insertBool ) {
                for ( auto &p : pmap ) {                         // loop over molecules
                  Group g( 0, p.second.size()-1 );               // (first=id, second=pvec)
                  g.molId = p.first;
                  g.setMolSize( p.second.size() );

                  u += pot->g_external( p.second, g );           // ...atoms/mols with external pot

                  if ( spc->molecule[g.molId].isAtomic() ) {
                    u += pot->g_internal( p.second, g );         // ...between inserted atoms
                    for ( auto &pi : p.second )                  // ...atoms with all particles
                      u += pot->all2p( spc->p, pi );
                  }
                  else {
                    for ( auto g2 : spc->groupList() )           // ...molecules with all groups
                      u += pot->g1g2(p.second, g, spc->p, *g2);
                    uinternal += pot->g_internal( p.second, g ); // ...internal mol energy (dummy)
                  }
                }

                for ( auto i=pmap.begin(); i!=pmap.end(); ++i )  //...between inserted molecules
                  for ( auto j=i; ++j!=pmap.end(); ) {
                    Group gi( 0, i->second.size()-1 );
                    Group gj( 0, i->second.size()-1 );
                    gi.molId = i->first;
                    gj.molId = j->first;
                    u += pot->g1g2( i->second, gi, j->second, gj);
                  }

                assert( !pmap.empty() );
                base::alternateReturnEnergy = u + uinternal;
                return u + externalEnergy();
              }

              // energy if deletion move
              else {
                if ( ! molDel.empty() || ! atomDel.empty() ) {
                  for ( auto i : molDel ) {                     // loop over molecules/atoms
                    u += pot->g_external( spc->p, *i );         // molecule w. external pot.

                    if ( ! spc->molecule[i->molId].isAtomic() ) {
                      for ( auto j : spc->groupList() )         // molecule w. all groups
                        if ( find(molDel.begin(), molDel.end(), j) == molDel.end() ) // slow!
                          u += pot->g2g( spc->p, *i, *j );
                      uinternal += pot->g_internal( spc->p, *i);// internal mol energy (dummy)
                    }
                  }

                  // energy between deleted molecules
                  for ( auto i=molDel.begin(); i!=molDel.end(); ++i )
                    for ( auto j=i; ++j!=molDel.end(); )
                      u += pot->g2g( spc->p, **i, **j );

                  for ( auto i : atomDel )                        // atoms w. all particles
                    u += pot->i_total( spc->p, i );

                  for ( int i=0; i<(int)atomDel.size()-1; i++ )   // subtract double counted 
                    for ( int j=i+1; j<(int)atomDel.size(); j++ ) // internal energy (atoms)
                      u -= pot->i2i(spc->p, i, j);

                  base::alternateReturnEnergy = -u - uinternal;
                  return -u + externalEnergy(); // ...add activity terms
                }
              }

              // if we reach here, we're out of particles -> reject

              assert(!insertBool);
              assert(fabs(u)<1e-10);
              base::alternateReturnEnergy = 0;
              return pc::infty;
            }

            void _acceptMove() {

              // accept a deletion move
              if ( !insertBool ) {
                Ndeleted++;
                for ( auto m : molDel ) { // loop over Group pointers
                  molTrack.erase(m->molId, m);
                  base::spc->eraseGroup( spc->findIndex(m) );
                }
                for ( auto i : atomDel ) {// loop over particle index
                  assert(1==2 && "Under construction");
                  base::spc->erase(i);
                }

                atomTrack.update( spc->p );
                return;
              }

              // accept an insertion move
              if ( insertBool ) {
                Ninserted++;
                for ( auto &p : pmap ) { // loop over sets of new coordinates
                  auto molid = p.first;
                  if ( spc->molecule[molid].isAtomic() ) {
                    assert(1==2 && "Under construction");
                    spc->insert( molid, p.second );
                    atomTrack.update( spc->p );
                  }
                  else {
                    assert( !p.second.empty() );
                    auto g = spc->insert( molid, p.second ); // auto gen. group
                    molTrack.insert( molid, g );
                    assert( molTrack.size(molid) > 0 );
                  }
                }
              }
              molTrack.updateAvg();   // update average number of molecules
              atomTrack.updateAvg();  // ...and atoms
              pot->setSpace( *spc );
            }

            void _rejectMove() {
              molTrack.updateAvg();   // update average number of molecules
              atomTrack.updateAvg();  // ...and atoms
            }

            string _info() {
              using namespace textio;
              std::ostringstream o;

              o << pad( SUB,base::w,"Accepted insertions" ) << Ninserted << "\n"
                << pad( SUB,base::w,"Accepted deletions" ) << Ndeleted << "\n"
                << pad( SUB,base::w,"Flux (Nins/Ndel)" ) << Ninserted / double(Ndeleted) << "\n"
                << "\n";

              double V=spc->geo.getVolume();
              o << std::left
                << setw(w+5) << "  Molecule/Atom"
                << setw(w) << "a (mol/l)"
                << setw(w) << "c (mol/l)"
                << setw(w) << textio::gamma+"=a/c" << "\n"
                << "  " << string(4*w,'-') << "\n";

              for (auto &m : spc->molecule) {
                if ( m.activity > 1e-10 )
                  if ( molTrack.getAvg(m.id).cnt>0 )
                    o << setw(w+5) << ("  "+m.name) << setw(w) << m.activity
                      << setw(w) << molTrack.getAvg(m.id)/V/1.0_molar
                      << setw(w) << m.activity / (molTrack.getAvg(m.id)/V/1.0_molar) << "\n";
              }

              o << "\n";

              for (auto &m : atom) {
                if ( m.activity > 1e-6 )
                  if ( atomTrack.getAvg(m.id).cnt>0 )
                    o << setw(w+5) << ("  "+m.name) << setw(w) << m.activity
                      << setw(w) << atomTrack.getAvg(m.id)/V/1.0_molar
                      << setw(w) << m.activity / (atomTrack.getAvg(m.id)/V/1.0_molar) << "\n";
              }

              return o.str() + spc->molecule.info() + comb.info();
            }

            void _test( UnitTest &t ) {
              double V=spc->geo.getVolume();
              t( this->jsondir + "_flux", Ninserted / double(Ndeleted) );
              for (auto &m : spc->molecule)
                if ( m.activity > 1e-6 )
                  if ( molTrack.getAvg(m.id).cnt>0 )
                    if ( !m.name.empty() )
                      t( this->jsondir + "_mol_" + m.name + "_gamma", m.activity / (molTrack.getAvg(m.id)/V/1.0_molar));
              for (auto &m : atom)
                if ( m.activity > 1e-6 )
                  if ( !m.name.empty() )
                    if ( atomTrack.getAvg(m.id).cnt>0 )
                      t( this->jsondir + "_atom_" + m.name + "_gamma", m.activity / (atomTrack.getAvg(m.id)/V/1.0_molar));
            }

            void init() { // call this upon construction
              Ninserted = 0;
              Ndeleted  = 0;
              base::title = "Grand Canonical";
              base::useAlternateReturnEnergy = true;

              // update tracker with GC molecules and atoms
              for ( auto g : spc->groupList() )
                if ( spc->molecule[g->molId].isAtomic() ) {
                  for ( auto i : *g )
                    if ( atom[ spc->p[i].id ].activity > 1e-9 ) 
                      atomTrack.insert( spc->p[i].id, i );
                }
                else
                  if ( spc->molecule[g->molId].activity > 1e-9 )
                    molTrack.insert( g->molId, g );
            }

          public:

            /** @brief Constructor -- load combinations, initialize trackers */
            GreenGC(
                Energy::Energybase<Tspace> &e,
                Tspace &s,
                Tmjson &j,
                const string &sec="gc") : base( e, s ), comb( s.molecule ) {

              init();
              base::jsondir = "moves/" + sec;
              base::runfraction = j["moves"][sec]["prob"] | 1.0;
              comb.include( j ); // load combinations
            }
        };

      /**
       * @brief Move for swapping species types - i.e. implicit titration
       *
       * Upon construction this class will add an instance of
       * Energy::EquilibriumEnergy to the Hamiltonian. For details
       * about the titration procedure see `Energy::EquilibriumController`.
       *
       * Upon construction the following are read from input section
       * `moves/titrate`:
       *
       *  Keyword       |  Description
       * :------------- |  :---------------------------------
       * `processfile`  |  json file name with processes
       * `prob`         |  probability of running (default: 1)
       */
      template<class Tspace>
        class SwapMove : public Movebase<Tspace> {
          private:
            typedef Movebase<Tspace> base;
            std::map<int, Average<double> > accmap; //!< Site acceptance map
            string _info();
            void _trialMove();
            void _acceptMove();
            void _rejectMove();

            std::map<int, std::map<int, Average<double> >> molCharge;

            void updateMolCharge( int pindex ) {
              auto g = spc->findGroup( pindex );
              molCharge[ g->molId ][ pindex - g->front() ] += spc->p[pindex].charge; 
            }

          protected:
            using base::spc;
            using base::pot;

            double _energyChange();
            int ipart;                              //!< Particle to be swapped
            Energy::EquilibriumEnergy<Tspace>* eqpot;

          public:
            template<class Tenergy>
              SwapMove(Tenergy&, Tspace&, Tmjson&, string="titrate"); //!< Constructor

            template<class Tpvec>
              int findSites(const Tpvec&); //!< Search for titratable sites (old ones are discarded)

            double move(int n=1) override {
              double du=0;
              if (this->run()) {
                eqpot->findSites( this->spc->p );
                size_t i = eqpot->eq.sites.size();
                while (i-->0)
                  du += base::move();
                eqpot->eq.sampleCharge(spc->p);
              }
              return du;
            }

            template<class Tpvec>
              void applycharges(Tpvec &);
        };

      template<class Tspace>
        template<class Tenergy> SwapMove<Tspace>::SwapMove(
            Tenergy &e, Tspace &spc, Tmjson &j, string sec) : base(e,spc) {

          base::title="Site Titration - Swap Move";
          base::runfraction = j["moves"][sec]["prob"] | 1.0;
          base::w = 30;
          ipart=-1;

          auto t = e.tuple();
          auto ptr = TupleFindType::get< Energy::EquilibriumEnergy<Tspace>* >( t );
          if ( ptr != nullptr )
            eqpot = *ptr; 
          else {
            std::cerr << "Error: Equilibrium energy required in\
              Hamiltonian for titration swap moves." << endl;
            exit(1);
          }

          eqpot->eq = Energy::EquilibriumController( j );
          findSites( spc.p );

          /* Sync particle charges with `AtomMap` */
          for (auto i : eqpot->eq.sites)
            spc.trial[i].charge = spc.p[i].charge = atom[ spc.p[i].id ].charge;
        }

      /**
       * @brief Search for titratable sites and store internally
       *
       * Use this to re-scan for titratable sites. Called by default
       * in the constructor
       */
      template<class Tspace>
        template<class Tpvec>
        int SwapMove<Tspace>::findSites(const Tpvec &p) {
          accmap.clear();
          return eqpot->findSites(p);
        }

      template<class Tspace>
        void SwapMove<Tspace>::_trialMove() {
          if (!eqpot->eq.sites.empty()) {
            int i = slump.range( 0, eqpot->eq.sites.size()-1); // pick random site
            ipart = eqpot->eq.sites.at(i);                      // and corresponding particle
            int k;
            do {
              k = slump.range( 0, eqpot->eq.process.size()-1 );// pick random process..
            } while (!eqpot->eq.process[k].one_of_us( this->spc->p[ipart].id )); //that match particle j

            eqpot->eq.process[k].swap( this->spc->trial[ipart] ); // change state and get intrinsic energy change
          }
        }

      template<class Tspace>
        double SwapMove<Tspace>::_energyChange() {
          assert( spc->geo.collision( spc->p[ipart], spc->p[ipart].radius )==false
              && "Accepted particle collides with container");

          if (spc->geo.collision(spc->trial[ipart], spc->trial[ipart].radius))  // trial<->container collision?
            return pc::infty;
          double uold = pot->external(spc->p) + pot->i_total(spc->p,ipart);
          double unew = pot->external(spc->trial) + pot->i_total(spc->trial,ipart);
#ifdef ENABLE_MPI
          if ( base::mpiPtr != nullptr ) {
            double sum=0;
            auto r = Faunus::MPI::splitEven(*base::mpiPtr, (int)spc->p.size());
            for (int i=r.first; i<=r.second; ++i)
              if (i!=ipart)
                sum+=pot->i2i(spc->trial,i,ipart) - pot->i2i(spc->p,i,ipart);

            sum = Faunus::MPI::reduceDouble(*base::mpiPtr, sum);

            return sum + pot->i_external(spc->trial, ipart) - pot->i_external(spc->p, ipart)
              + pot->i_internal(spc->trial, ipart) - pot->i_internal(spc->p, ipart);
          }
#endif

          return unew - uold;
        }

      template<class Tspace>
        void SwapMove<Tspace>::_acceptMove() {
          accmap[ipart] += 1;
          spc->p[ipart] = spc->trial[ipart];
          updateMolCharge( ipart );
        }

      template<class Tspace>
        void SwapMove<Tspace>::_rejectMove() {
          accmap[ipart] += 0;
          spc->trial[ipart] = spc->p[ipart];
          updateMolCharge( ipart );
        }

      template<class Tspace>
        template<class Tpvec>
        void SwapMove<Tspace>::applycharges(Tpvec &p){
          eqpot->eq.applycharges(p);
        }

      template<class Tspace>
        string SwapMove<Tspace>::_info() {
          using namespace textio;
          std::ostringstream o;
          for (auto &m : molCharge) {
            int molid = m.first;
            o << "\n" << indent(SUB) << "Molecule: " << spc->molList()[ molid ].name << "\n\n"
              << std::left << "    " << setw(8) << "index" << setw(12) << "name"
              << setw(12) << "Z" << "\n";
            for (auto &i : m.second)
              o << "    " << setw(8) << i.first
                << setw(12) << atom[ spc->molList()[molid].atoms[i.first] ].name
                << setw(12) << i.second << "\n"; 
          }
          return o.str();
        }

      /**
       * @brief As SwapMove but Minimizes Short Ranged interactions
       *        within a molecule upon swapping
       *
       * Before calculating dU of an attempted swap move, radii on
       * particles within the SAME group are set to minus radius of
       * the swapped particle and hydrophobicity is set to false.
       * This to minimize large interactions in molecules with overlapping
       * particles - i.e LJ will be zero. It can also be used to avoid
       * internal hydrophobic interactions in rigid groups upon swapping
       * between hydrophobic and non-hydrophobic species.
       */
      template<class Tspace>
        class SwapMoveMSR : public SwapMove<Tspace> {
          private:
            using SwapMove<Tspace>::spc;
            using SwapMove<Tspace>::pot;
            std::map<int, double> radiusbak;    // backup for radii
            std::map<int, bool> hydrophobicbak; // backup for hydrophobic state

            void modify() {
              radiusbak.clear();
              hydrophobicbak.clear();
              for (auto g : spc->groupList() )   // loop over all groups
                if (g->find(this->ipart)) {  //   is ipart part of a group?
                  for (auto i : *g)    //     if so, loop over that group
                    if (i!=this->ipart) {    //       and ignore ipart
                      assert( abs(spc->p[i].radius-spc->trial[i].radius)<1e-9);
                      assert( spc->p[i].hydrophobic==spc->trial[i].hydrophobic);

                      //radiusbak[i]         = spc->p[i].radius;
                      //spc->p[i].radius     = -spc->p[ipart].radius;
                      //spc->trial[i].radius = -spc->p[ipart].radius;

                      hydrophobicbak[i]         = spc->p[i].hydrophobic;
                      spc->p[i].hydrophobic     = false;
                      spc->trial[i].hydrophobic = false;
                    }
                  return; // a particle can be part of a single group, only
                }
            }

            void restore() {
              for (auto &m : radiusbak) {
                spc->p[m.first].radius = m.second;
                spc->trial[m.first].radius = m.second;
              }
              for (auto &m : hydrophobicbak) {
                spc->p[m.first].hydrophobic = m.second;
                spc->trial[m.first].hydrophobic = m.second;
              }
            }

            double _energyChange() {
              double du_orig = SwapMove<Tspace>::_energyChange();
              modify();
              double du = SwapMove<Tspace>::_energyChange();
              restore();
              this->alternateReturnEnergy=du_orig;
              return du;
            }

          public:
            SwapMoveMSR(
                InputMap &in, Energy::Energybase<Tspace> &ham, Tspace &spc,
                string pfx="swapmv_") : SwapMove<Tspace>(in,ham,spc,pfx)
            {
              this->title+=" (min. shortrange)";
              this->useAlternateReturnEnergy=true;
            }
        };


      /**
       * @brief Multiple moves controlled via JSON input
       *
       * This is a move class that randomly picks between a number of
       * moves as defined in a JSON file in the section `moves`.
       * The available moves are shown
       * in the table below; each can occur only once and are picked
       * with uniform weight.
       *
       * Keyword         | Move class                | Description
       * :-------------- | :------------------------ | :----------------
       * `atomtranslate` | `Move::AtomicTranslation` | Translate atoms
       * `atomrotate`    | `Move::AtomicRotation`    | Rotate atoms
       * `atomgc`        | `Move::GrandCanonicalSalt`| GC salt move (muVT ensemble)
       * `crankshaft`    | `Move::CrankShaft`        | Crank shaft polymer move
       * `gc`            | `Move::GreenGC`           | Grand canonical move (muVT ensemble)
       * `isobaric`      | `Move::Isobaric`          | Volume move (NPT ensemple)
       * `moltransrot`   | `Move::TranslateRotate`   | Translate/rotate molecules
       * `pivot`         | `Move::Pivot`             | Pivot polymer move
       * `titrate`       | `Move::SwapMove`          | Particle swap move
       */
      template<typename Tspace, bool polarise=false, typename base=Movebase<Tspace>>
        class Propagator : public base {
          private:
            typedef std::shared_ptr<base> basePtr;
            std::vector<basePtr> mPtr; 

            string _info() override {
              string o;
              for ( auto i : mPtr )
                o += i->info();
              return o;
            }

            void _acceptMove() override { assert(1==2); }
            void _rejectMove() override { assert(1==2); }
            void _trialMove()  override { assert(1==2); }
            double _energyChange() override { assert(1==2); return 0;}

            template<typename Tmove>
              basePtr toPtr(Tmove m) {
                typedef typename std::conditional<polarise, PolarizeMove<Tmove>, Tmove>::type T;
                return basePtr( new T(m) );
              }

          public:
            template<typename Tenergy>
              Propagator( InputMap &in, Tenergy &e, Tspace &s ) : base( e, s )
            {
              this->title = "P R O P A G A T O R S";

              auto m = in["moves"];
              for ( auto i=m.begin(); i!=m.end(); ++i) {
                if (i.key()=="atomtranslate")
                  mPtr.push_back( toPtr( AtomicTranslation<Tspace>(e,s,in) ) );
                if (i.key()=="atomrotate")
                  mPtr.push_back( toPtr( AtomicRotation<Tspace>(e,s,in) ) );
                if (i.key()=="atomgc")
                  mPtr.push_back( toPtr( GrandCanonicalSalt<Tspace>(e,s,in) ) );
                if (i.key()=="gctit")
                  mPtr.push_back( toPtr( GrandCanonicalTitration<Tspace>(e,s,in) ) );
                if (i.key()=="moltransrot")
                  mPtr.push_back( toPtr( TranslateRotate<Tspace>(e,s,in) ) );
                if (i.key()=="isobaric")
                  mPtr.push_back( toPtr( Isobaric<Tspace>(e,s,in) ) );
                if (i.key()=="isochoric")
                  mPtr.push_back( toPtr( Isochoric<Tspace>(e,s,in) ) );
                if (i.key()=="gc")
                  mPtr.push_back( toPtr( GreenGC<Tspace>(e,s,in) ) );
                if (i.key()=="titrate")
                  mPtr.push_back( toPtr( SwapMove<Tspace>(e,s,in) ) );
                if (i.key()=="crankshaft")
                  mPtr.push_back( toPtr( CrankShaft<Tspace>(e,s,in) ) );
                if (i.key()=="pivot")
                  mPtr.push_back( toPtr( Pivot<Tspace>(e,s,in) ) );
              }
              if ( mPtr.empty() )
                throw std::runtime_error("No moves defined - check JSON file.");
            }

            double move(int n=1) override {
              return ( mPtr.empty() ) ?
                0 : (*propagation_slump.element( mPtr.begin(), mPtr.end() ))->move();
            }

            void test(UnitTest &t) { for (auto i : mPtr) i->test(t); }

#ifdef ENABLE_MPI
            void setMPI( Faunus::MPI::MPIController* mpi ) {
              base::mpiPtr = mpi;
              for ( auto i : mPtr )
                i->mpiPtr = mpi;
            }
#endif

        };

      /** @brief Atomic translation with dipolar polarizability */
      //typedef PolarizeMove<AtomicTranslation> AtomicTranslationPol;

      /** @brief Atomic rotation with dipolar polarizability */
      //typedef PolarizeMove<AtomicRotation> AtomicRotationPol;

      }//namespace
  }//namespace
#endif
