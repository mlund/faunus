#ifndef FAUNUS_MCMOVE_H
#define FAUNUS_MCMOVE_H

#include <faunus/common.h>
#include <faunus/point.h>
#include <faunus/average.h>
#include <faunus/energy.h>

namespace Faunus {
  //class average;
  class unittest;
  class Energybase;

  /*!
   * \brief Monte Carlo move related classes
   */
  namespace Move {

    /*
     * \brief Optimize Monte Carlo displacement parameters for optimal mean square displament
     * \date Lund, 2011
     * \author Mikael Lund
     */
    class DisplacementOptimizer {
      private:
        double* dpPtr;   //!< Pointer to displacement parameter
        //average<double>* msqPtr; //!< Pointer to mean-square displacement
        double dp_min;
        double dp_max;

      public:
        //void set_target(double&, average&);
    };

    /*!
     * \brief Base class for Monte Carlo moves
     * \author Mikael Lund
     * \date Lund, 2007-2011
     *
     * The is a base class that handles Monte Carlo moves and derived classes
     * are required to implement the following pure virtual (and private)
     * functions:
     *
     * \li \c _trialMove()
     * \li \c _energyChange()
     * \li \c _acceptMove()
     * \li \c _rejectMove()
     * \li \c _info()
     *
     * These functions should be pretty self-explanatory and are - via wrapper
     * functions - called by the move(). It is important that the _energyChange() function
     * returns the full energy associated with the move. For example, for NPT
     * moves the pV term should be included and so on. Please do stride not try override
     * the move() function as this should be generic to all MC moves.
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
        virtual double _energyChange()=0;      //!< Returns energy change of trialMove

        void trialMove();                      //!< Do a trial move (wrapper)
        void acceptMove();                     //!< Accept move, store new coordinates etc. (wrapper)
        void rejectMove();                     //!< Reject move, revert to old coordinates etc. (wrapper)
        double energyChange();                 //!< Returns energy change of trialMove (wrapper)
        bool metropolis(const double&) const;  //!< Metropolis criteria

      protected:
        Energy::Energybase* pot;         //!< Pointer to energy functions
        Space* spc;
        string title;                    //!< title of move
        string cite;                     //!< litterature reference, url, DOI etc.
        string prefix;                   //!< inputmap prefix
        char w;                          //!< info string text width
        unsigned long int cnt;           //!< total number of trial moves
        bool run() const;                //!< Runfraction test

      public:
        Movebase(Energy::Energybase&, Space&, string);//!< Constructor
        virtual ~Movebase();
        double runfraction;                //!< Fraction of times calling move() should result in an actual move. 0=never, 1=always.
        virtual double move(int=1) final;  //!< Attempt \c n moves and return energy change
        virtual string info() final;       //!< Returns information string (wrapper)
        virtual void test(UnitTest&) final;//!< Perform unit test
        double getAcceptance();            //!< Get acceptance [0:1]
    };

    /*!
     * \brief Translation of atomic particles
     * \author Mikael Lund
     * \date Lund, 2011
     *
     * This Markov move can work in two modes:
     * \li Move a single particle in space set by setParticle()
     * \li Move single particles randomly selected in a Group set by setGroup().
     *
     * The move directions can be controlled with the dir vector - for instance if you wish
     * to translate only in the \c z direction, set \c dir.x=dir.y=0.
     */
    class AtomicTranslation : public Movebase {
      private:
        typedef std::map<short, Average<double> > map_type;
        map_type accmap; //!< Single particle acceptance map
        map_type sqrmap; //!< Single particle mean square displacement map
        Group* igroup;   //!< Group pointer in which particles are moved randomly (NULL if none, default)
        int iparticle;   //!< Select single particle to move (-1 if none, default)
        Average<unsigned long long int> gsize; //!< Average size of igroup;

        string _info();
        void _trialMove();
        void _acceptMove();
        void _rejectMove();
        double _energyChange();
      public:
        AtomicTranslation(InputMap&, Energy::Energybase&, Space&, string="mv_particle");
        void setGroup(Group&); //!< Select group in which to randomly pick particles from
        void setParticle(int); //!< Select single particle in Space::p to move
        Point dir;             //!< Translation directions (default: x=y=z=1)
    };

    /*!
     * \brief Combined rotation and rotation of groups
     * \author Mikael Lund
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
      public:
        TranslateRotate(InputMap&, Energy::Energybase&, Space&, string="transrot");
        void setGroup(Group&); //!< Select Group to move
        Point dir;             //!< Translation directions (default: x=y=z=1)
        bool groupWiseEnergy;  //!< Attempt to evaluate energy over groups from vector in Space (default=false)
    };

    /*!
     * \brief Combined rotation and rotation of groups and mobile species around it
     * \author Mikael Lund
     *
     * This class will do a combined translational and rotational move of a group along with
     * atomic particles surrounding it. To specify where to look for clustered particles, use
     * the setMobile() function. Whether particles are considered part of the cluster is
     * determined by the private virtual function ClusterProbability(). By default this is a simple
     * step function with P=1 when an atomic particle in the group set by setMobile is closer
     * than a certain threshold to a particle in the main group; P=0 otherwise.
     * The implemented cluster algorithm is general - see Frenkel and Smith, 2nd ed, p405 - and derived classes
     * can re-implement ClusterProbability() for arbitrary probability functions.
     */
     class TranslateRotateCluster : public TranslateRotate {
      private:
        Geometry::VectorRotate vrot;
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
        void setMobile(Group&); //!< Select atomic species to move with the main group
        double threshold;  //!< Distance between particles to define a cluster
    };

    /*!
     * \brief Isobaric volume move
     *
     * This class will perform a volume displacement and scale atomic as well as molecular
     * groups as long as these are known to Space -- see Space.enroll().
     * The constructor will automatically add an instance of Energy::ExternalPressure
     * to the Hamiltonian. The InputMap class is scanned for the following keys:
     * \li \c npt_dV \n Volume displacement parameter
     * \li \c npt_P \n Pressure
     * \li \c npt_Punit \n Pressure unit: mM [default] or 1/A3
     *
     * Note that the volume displacement is done by:
     *
     * \f$ V^{\prime} = \exp\left ( \log V \pm \delta dV \right ) \f$ where \f$\delta\f$ is a random number
     * between zero and one half.
     *
     * Example:
     * \code
     * Energy::Hamiltonian pot;        // we need a hamiltonian
     * ...                             // insert interactions...
     * Space spc( pot.getGeometry() ); // ...and a space
     * spc.enroll( myprotein );        // register some groups
     * spc.enroll( mysolvent );
     * Move::Isobaric(in, pot, spc);   // setup and add pressure to hamiltonian
     * \endcode
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
      public:
        Isobaric(InputMap&, Energy::Hamiltonian&, Space&, string="npt");
    };

    /*!
     * \brief Auxillary class for tracking atomic species
     * \author Mikael Lund
     * \date Malmo 2011
     *
     * This class keeps track of individual particle positions based on particle type (id).
     * It contains functions to insert and erase particles while automatically moving particles
     * above the deletion or insertion point in sync.
     *
     * Example:
     * \code
     * AtomTracker track(myspace);
     * track.insert( myparticle, 20 );        // insert particle into Space at position 20
     * ...
     * int i=track[ myparticle.id ].random(); // pick a random particle of type myparticle.id
     * \endcode
     */
    class AtomTracker {
      private:
        typedef short Tid;  // particle id type
        typedef int Tindex; // particle index type
        Space* spc;
        class data {
          public:
            vector<Tindex> index;
            Tindex random();                  //!< Pick random particle index
        };
        std::map<Tid,data> map; 
      public:
        AtomTracker(Space&);
        Tid randomAtomType() const;           //!< Select a random atomtype from the list
        bool insert(const particle&, Tindex); //!< Insert particle into Space and track position
        bool erase(Tindex);                   //!< Delete particle from Space at specific particle index
        data& operator[] (Tid);               //!< Access operator to atomtype data
        void clear();                         //!< Clear all atom lists (does not touch Space)
        bool empty();                         //!< Test if atom list is empty
    };

    /*!
     * \brief Grand Canonical insertion of arbitrary M:X salt pairs
     * \author Bjorn Persson and Mikael Lund
     * \date Lund 2010-2011
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
        std::map<short,ionprop> map;
        void randomIonPair(short&,short&);  // Generate random ion pair
        p_vec trial_insert;
        vector<int> trial_delete;
        short ida, idb;                    // particle id's of current salt pair (a=cation, b=anion)

        Energy::EnergyRest Urest;   // store non-Hamiltonian energy here to correctly calculate energy drift
        double du_rest;
        Group* saltPtr;  // GC ions *must* be in this group

      public:
        GrandCanonicalSalt(InputMap&, Energy::Hamiltonian&, Space&, Group&, string="saltbath");
    };

  }//namespace
}//namespace
#endif
