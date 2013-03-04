#ifndef FAUNUS_MCMOVE_H
#define FAUNUS_MCMOVE_H

#ifndef SWIG
#include <faunus/common.h>
#include <faunus/point.h>
#include <faunus/average.h>
#include <faunus/textio.h>
#include <faunus/geometry.h>
#include <faunus/energy.h>

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
     * @brief Add polarization step to a move
     *
     * This will insert an electric field calculation
     * after the original trial move and iteratively
     * calculate induced dipole moments on all particles.
     */
    template<class Tmove>
      class PolarizeMove : public Tmove {
        protected:
          void _trialMove() {
            Tmove::_trialMove();
            // ... update induced moments
          }
        public:
          PolarizeMove(InputMap &in, Energy::Energybase &e, Space &s) :
            Tmove(in,e,s) {}
      };

    template<class Tmove>
      class EwaldMove : public Tmove {
        protected:
          void _trialMove() {
            Tmove::_trialMove();
            // ... update induced moments
          }
        public:
          EwaldMove(InputMap &in, Energy::Energybase &e, Space &s) :
            Tmove(in,e,s) {}
      };

    /*!
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
     * moves the pV term should be included and so on. Please do stride not try override
     * the move() function as this should be generic to all MC moves.
     *
     * @date Lund, 2007-2011
     */

    class Movebase {
      private:
        unsigned long int cnt_accepted;        //!< number of accepted moves
        double dusum;                          //!< Sum of all energy changes made by this move

        virtual void _test(UnitTest&);         //!< Unit testing
        virtual string _info()=0;              //!< Specific info for derived moves
        virtual void _trialMove()=0;           //!< Do a trial move
        virtual void _acceptMove()=0;          //!< Accept move, store new coordinates.
        virtual void _rejectMove()=0;          //!< Reject move, revert to old coordinates.
        virtual double _energyChange()=0;      //!< Returns energy change of trialMove (kT)

        void acceptMove();                     //!< Accept move, store new coordinates etc. (wrapper)
        void rejectMove();                     //!< Reject move, revert to old coordinates etc. (wrapper)
        double energyChange();                 //!< Returns energy change of trialMove i kT (wrapper)
        bool metropolis(const double&) const;  //!< Metropolis criteria

      protected:
        void trialMove();                      //!< Do a trial move (wrapper)
        Energy::Energybase* pot;         //!< Pointer to energy functions
        Space* spc;                      //!< Pointer to Space (particles and groups are stored there)
        string title;                    //!< title of move (mandatory!)
        string cite;                     //!< litterature reference, url, DOI etc.
        string prefix;                   //!< inputmap prefix
        char w;                          //!< info string text width. Adjust this in constructor if needed.
        unsigned long int cnt;           //!< total number of trial moves
        virtual bool run() const;        //!< Runfraction test

        bool useAlternateReturnEnergy;   //!< Return a different energy than returned by _energyChange(). [false]
        double alternateReturnEnergy;    //!< Alternative return energy (kT).

      public:
        Movebase(Energy::Energybase&, Space&, string);//!< Constructor
        virtual ~Movebase();
        double runfraction;                //!< Fraction of times calling move() should result in an actual move. 0=never, 1=always.
        double move(int=1);                //!< Attempt \c n moves and return energy change (kT)
        string info();                     //!< Returns information string
        void test(UnitTest&);              //!< Perform unit test
        double getAcceptance();            //!< Get acceptance [0:1]
    };

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
    class AtomicTranslation : public Movebase {
      private:
        typedef std::map<short, Average<double> > map_type;
        string _info();
        void _trialMove();
        void _acceptMove();
        void _rejectMove();
        double _energyChange();
        bool run() const;                //!< Runfraction test
      protected:
        map_type accmap; //!< Single particle acceptance map
        map_type sqrmap; //!< Single particle mean square displacement map

        int iparticle;   //!< Select single particle to move (-1 if none, default)
        Group* igroup;   //!< Group pointer in which particles are moved randomly (NULL if none, default)
        double genericdp;//!< Generic atom displacement parameter - ignores individual dps
        Average<unsigned long long int> gsize; //!< Average size of igroup;

      public:
        AtomicTranslation(InputMap&, Energy::Energybase&, Space&, string="mv_particle");
        void setGroup(Group&); //!< Select group in which to randomly pick particles from
        void setParticle(int); //!< Select single particle in Space::p to move
        void setGenericDisplacement(double); //!< Set single displacement for all atoms
        Point dir;             //!< Translation directions (default: x=y=z=1)
    };

    /**
     * @brief Rotate single particles
     *
     * This move works in the same way as AtomicTranslation but does
     * rotations of non-isotropic particles instead of translation. This move
     * has no effect on isotropic particles such as Faunus::PointParticle.
     */
    class AtomicRotation : public AtomicTranslation {
      private:
        void _trialMove();
        string _info();
        Geometry::QuaternionRotate rot;
      public:
        AtomicRotation(InputMap&, Energy::Energybase&, Space&, string="rot_particle");
    };

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
    class TranslateRotate : public Movebase {
      protected:
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
        TranslateRotate(InputMap&, Energy::Energybase&, Space&, string="transrot");
        void setGroup(Group&); //!< Select Group to move
        bool groupWiseEnergy;  //!< Attempt to evaluate energy over groups from vector in Space (default=false)
        std::map<string,Point> directions; //!< Specify special group translation directions (default: x=y=z=1)
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
     * in the group set by setMobile is closer than a certain threshold to a
     * particle in the main group; P=0 otherwise.
     *
     * The implemented cluster algorithm is general - see Frenkel&Smith,
     * 2nd ed, p405 - and derived classes can re-implement `ClusterProbability()`
     * for arbitrary probability functions.
     *
     */
    class TranslateRotateCluster : public TranslateRotate {
      private:
        //Geometry::VectorRotate vrot;
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
        TranslateRotateCluster(InputMap&, Energy::Energybase&, Space&, string="transrot");
        virtual ~TranslateRotateCluster();
        void setMobile(Group&); //!< Select atomic species to move with the main group
        double threshold;  //!< Distance between particles to define a cluster
    };

    /**
     * @brief Non-rejective cluster translation.
     *
     * This type of move will attempt to move collective sets of macromolecules that
     * obeys some criteria (here a hardcore overlap(?)) with a symmetric transition
     * matrix (no flow through the clusters).
     *
     * Setting the boolen `skipEnergyUpdate` to true (default is false) updates of the
     * total energy are skipped to speed up the move.
     * While this has no influence on the Markov chain it will cause an apparent energy
     * drift. It is recommended that this is enabled only for long production runs after
     * having properly checked that no drifts occur with `skipEnergyUpdate=false`.
     *
     * @author Bjoern Persson
     * @date Lund 2009-2010
     * @note Requirements for usage:
     * - Compatible only with purely molecular systems
     * - Works only with periodic containers
     * - External potentials are ignored
     */
    class ClusterTranslateNR : public Movebase {
      private:
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
        ClusterTranslateNR(InputMap&, Energy::Energybase&, Space&, string="ctransnr");
        bool skipEnergyUpdate;    //!< True if energy updates should be skipped (faster evaluation!)
    };

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
    class CrankShaft : public Movebase {
      private:
        void _test(UnitTest&);
        void _trialMove();
        void _acceptMove();
        void _rejectMove();
        double _energyChange();
        string _info();
        virtual bool findParticles(); //!< This will set the end points and find particles to rotate
      protected:
        Group* gPtr;       //!< Pointer to group where move is to be performed. Set by setGroup().
        double dp;         //!< Rotational displacement parameter
        double angle;      //!< Current rotation angle
        vector<int> index; //!< Index of particles to rotate
        //Geometry::VectorRotate vrot;
        Geometry::QuaternionRotate vrot;
        AcceptanceMap<string> accmap;
      public:
        CrankShaft(InputMap&, Energy::Energybase&, Space&, string="crank");
        virtual ~CrankShaft();
        void setGroup(Group&); //!< Select Group to of the polymer to move
        int minlen;            //!< Minimum number of particles to rotate (default = 1)
        int maxlen;            //!< Maximin number of particles to rotate (default = 10)
    };

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
    class Pivot : public CrankShaft {
      private:
        bool findParticles();
      public:
        Pivot(InputMap&, Energy::Energybase&, Space&, string="pivot");
    };

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
    class Reptation : public Movebase {
      private:
        AcceptanceMap<string> accmap;
        void _test(UnitTest&);
        void _trialMove();
        void _acceptMove();
        void _rejectMove();
        double _energyChange();
        string _info();
        Group* gPtr;
        double bondlength; //!< Reptation length used when generating new head group position
      public:
        Reptation(InputMap&, Energy::Energybase&, Space&, string="reptation");
        void setGroup(Group&); //!< Select Group to move
    };

    /**
     * @brief Isobaric volume move
     *
     * @details This class will perform a volume displacement and scale atomic
     * as well as molecular groups as long as these are known to Space -
     * see Space.enroll().
     * The constructor will automatically add an instance of Energy::ExternalPressure
     * to the Hamiltonian. The InputMap class is scanned for the following keys:
     *
     * Key           | Description
     * :------------ | :-----------------------------
     * `npt_dV`      | Volume displacement parameter
     * `npt_P`       | Pressure [mM]
     *
     * Note that new volumes are generated according to
     * \f$ V^{\prime} = \exp\left ( \log V \pm \delta dV \right ) \f$
     * where \f$\delta\f$ is a random number between zero and one half.
     *
     * Example:
     *
     *     Energy::Hamiltonian pot;        // we need a hamiltonian
     *     ...                             // insert interactions...
     *     Space spc( pot.getGeometry() ); // ...and a space
     *     spc.enroll( myprotein );        // register some groups
     *     spc.enroll( mysolvent );
     *     Move::Isobaric(in, pot, spc);   // setup and add pressure to hamiltonian
     *
     */
    class Isobaric : public Movebase {
      private:
        Energy::Hamiltonian* hamiltonian;
        string _info();
        void _test(UnitTest&);
        void _trialMove();
        void _acceptMove();
        void _rejectMove();
        double _energy(const p_vec&);
        double _energyChange();
        double dV; //!< Volume displacement parameter
        double oldV;
        double newV;
        double P; //!< Pressure
        Average<double> sqrV;       //!< Mean squared volume displacement
        Average<double> V;          //!< Average volume
        Average<double> rV;         //!< Average 1/volume
      public:
        Isobaric(InputMap&, Energy::Hamiltonian&, Space&, string="npt");
    };

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
    class AtomTracker {
      public:
        typedef int Tindex; // particle index type
      private:
        Space* spc;
        class data {
          public:
            vector<Tindex> index;
            Tindex random();                  //!< Pick random particle index
        };
        std::map<particle::Tid,data> map;
      public:
        AtomTracker(Space&);
        particle::Tid randomAtomType() const; //!< Select a random atomtype from the list
        bool insert(const particle&, Tindex); //!< Insert particle into Space and track position
        bool erase(Tindex);                   //!< Delete particle from Space at specific particle index
        data& operator[] (particle::Tid);     //!< Access operator to atomtype data
        void clear();                         //!< Clear all atom lists (does not touch Space)
        bool empty();                         //!< Test if atom list is empty
    };

    /**
     * @brief Grand Canonical insertion of arbitrary M:X salt pairs
     * @author Bjorn Persson and Mikael Lund
     * @date Lund 2010-2011
     * @warning Untested in this branch
     */
    class GrandCanonicalSalt : public Movebase {
      private:
        string _info();
        void _trialMove();
        void _acceptMove();
        void _rejectMove();
        double _energyChange();
        void add(Group&);       // add salt group and scan for ions with non-zero activities

        AtomTracker tracker;
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

        Energy::EnergyRest Urest;   // store non-Hamiltonian energy here to correctly calculate energy drift
        double du_rest;
        Group* saltPtr;  // GC ions *must* be in this group

      public:
        GrandCanonicalSalt(InputMap&, Energy::Hamiltonian&, Space&, Group&, string="saltbath");
    };

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
    class ParallelTempering : public Movebase {
      private:
        enum extradata {VOLUME=0};    //!< Indicates structure of extra data to send
        typedef std::map<string, Average<double> > map_type;
        map_type accmap;              //!< Acceptance map
        int partner;                  //!< Other replica (partner) to exchange with
        virtual void findPartner();   //!< Find which replica to exchange with
        bool goodPartner();           //!< Detemine if found partner is valid
        double exchangeEnergy(double);//!< Exchange energy with partner
        string id();                  //!< Get unique string to identify set of partners

        double currentEnergy;         //!< Energy of configuration before move (uold)
        bool haveCurrentEnergy;       //!< True if currentEnergy has been set

        string _info();
        void _trialMove();
        void _acceptMove();
        void _rejectMove();
        double _energyChange();
        std::ofstream temperPath;

        Energy::Energybase* hamiltonian;   //!< Hamiltonian class needed for volume displacement
        Faunus::MPI::MPIController *mpiPtr; //!< Controller class for MPI calls
        Faunus::MPI::FloatTransmitter ft;   //!< Class for transmitting floats over MPI
        Faunus::MPI::ParticleTransmitter pt;//!< Class for transmitting particles over MPI
        std::function<double (Space&, Energy::Energybase&, const p_vec&)> usys; //!< Defaults to Energy::systemEnergy but can be replaced!

      public:
        ParallelTempering(InputMap&, Energy::Energybase&, Space&, Faunus::MPI::MPIController &mpi, string="temper");
        virtual ~ParallelTempering();
        void setCurrentEnergy(double); //!< Set energy of configuration before move (for increased speed)
        void setEnergyFunction( std::function<double (Space&,Energy::Energybase&,const p_vec&)> );
    };

#endif

  }//namespace
}//namespace
#endif
