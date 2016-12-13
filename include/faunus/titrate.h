#ifndef FAUNUS_TITRATE_H
#define FAUNUS_TITRATE_H

#ifndef SWIG
#include <faunus/common.h>
#include <faunus/energy.h>
#include <faunus/average.h>
#include <faunus/json.h>

#endif

namespace Faunus
{

  namespace Energy
  {

    /**
     * @brief  Class for implicit titration of species with fixed chemical potential.
     *
     * Consider the dissociation process AX=A+X. This class will locate
     * all species of type AX and A and make a MC swap move between them.
     * X is implicit, meaning that it enters only with its chemical potential
     * (activity). The titrating species, their dissociation constants
     * and the chemical potential of the titrant are read from a
     * `processes` JSON object.
     * For example, for proton titration of phosphate one would
     * use the following JSON input (pH 7.0):
     *
     *     {
     *       "processes" :
     *       {
     *         "K1" : { "bound":"H3PO4", "free":"H2PO4", "pKd":2.12,  "pX":7.0 },
     *         "K2" : { "bound":"H2PO4", "free":"HPO4",  "pKd":7.21,  "pX":7.0 },
     *         "K3" : { "bound":"HPO4",  "free":"PO4",   "pKd":12.67, "pX":7.0 }
     *       }
     *     }
     *
     * All species and their properties must be defined in `AtomMap` before
     * initializing this class.
     *
     * @date Malmo, October 2010
     * @author Mikael Lund and Chris Evers
     */
    class EquilibriumController
    {
    private:
        typedef int Tid;                 //!< Particle type id
    public:
        class processdata
        {
        public:
            typedef int Tid;             //!< Particle type id
            double mu_AX;                //!< chemical potential of AX
            double mu_A;                 //!< chemical potential of A
            double mu_X;                 //!< chemical potential of X (this is the titrant)
            double ddG;                  //!< ddG = mu_A + mu_X - mu_AX
            int cnt;                     //!< number of sites for this process

        public:
            Tid id_AX, id_A;   //!< Particle id's for AX and A

            bool one_of_us( const Tid & );//!< Does the particle belong to this process?

            double energy( const Tid & ); //!< Returns intrinsic energy of particle id

            template<class Tparticle>
            double swap( Tparticle & );   //!< Swap AX<->A and return intrinsic energy change

            void set( double, double );     //!< Set activity of X and the pKd value

            void set_mu_AX( double );      //!< Set chemical potential of species AX - mu_A then follows.

            void set_mu_A( double );       //!< Set chemical potential of species A  - mu_AX then follows.
            bool bound( const Tid & );      //!< Returns true if state is bound for given process
        };

        std::map<int, Average < double> >
        q;       //!< Map of average charges per site
        std::vector<processdata> process;        //!< Vector of processes.

        EquilibriumController( Tmjson & );

        template<class Tpvec>
        void findSites( const Tpvec & );            //!< Locate all titratable sites

        double intrinsicEnergy( const Tid & );        //!< Intrinsic energy of particle id (kT)
        string info( char= 25 );                      //!< Get information string

        template<class Tpvec>
        processdata &random( const Tpvec &, int & ); //!< Random titratable particle and assiciated random process

        std::vector<int> sites;                    //!< List of titratable sites

        template<class Tpvec>
        void sampleCharge( const Tpvec & );         //!< Updates the average charge vector titrate::q

        template<class Tpvec>
        void copyAvgCharge( Tpvec & );             //!< Copy average charges to particles in the particle vector

        template<class Tpvec>
        double avgcharge( const Tpvec &, int & );    //!< Print average charges of process i

        int number_of_sites() const { return sites.size(); }
    };

    /**
     * @param pX  Negative logarithm of the X activity (titrant)
     * @param pKd Negative logarithm of dissociation constant.
     */
    void EquilibriumController::processdata::set( double pX, double pKd )
    {
        ddG = -log(pow(10., -pKd));
        mu_X = -log(pow(10., -pX));
        set_mu_AX(0);
    }

    void EquilibriumController::processdata::set_mu_AX( double mu )
    {
        mu_AX = mu;
        mu_A = ddG + mu_AX - mu_X;
    }

    void EquilibriumController::processdata::set_mu_A( double mu )
    {
        mu_A = mu;
        mu_AX = mu_A + mu_X - ddG;
    }

    /**
     * Returns `true` if the particle either matches AX or A.
     */
    bool EquilibriumController::processdata::one_of_us( const Tid &id )
    {
        if ( id == id_AX || id == id_A )
            return true;
        return false;
    }

    bool EquilibriumController::processdata::bound( const Tid &id )
    {
        return (id == id_AX) ? true : false;
    }

    /*!
     * Returns the intrinsic energy of the given particle. Intrinsic
     * means the energy stemming from the equilibrium expression when
     * no external interactions are accounted for (activity factors unity).
     */
    double EquilibriumController::processdata::energy( const Tid &id )
    {
        if ( id == id_AX )
            return mu_AX;
        if ( id == id_A )
            return mu_A;
        return 0;
    }

    /**
     * This will swap the state of given particle from AX<->A and
     * return the energy change associated with the process.
     *
     * @note This will *NOT* swap particle radii nor masses!
     */
    template<class Tparticle>
    double EquilibriumController::processdata::swap( Tparticle &p )
    {
        double uold = energy(p.id);
        double oldradius = p.radius;
        double oldmw = p.mw;
        Point pos = p;           // backup coordinate
        if ( p.id == id_AX )
            p = atom[id_A];
        else if ( p.id == id_A )
            p = atom[id_AX];
        p = pos;                    // apply old coordinates
        p.radius = oldradius;       // apply old radius
        p.mw = oldmw;
        return energy(p.id) - uold; // return intrinsic energy change
    }

    /**
     * @brief Construct from `InputMap`.
     *
     * Call this *after* particles have been loaded into `Space`, i.e.
     * typically just before starting the Markov chain. Also make
     * sure that `AtomMap` has been loaded with all atomic properties
     * as these will be used to reset the charge, radii, weight etc.
     * on all particles in the system.
     */
    EquilibriumController::EquilibriumController( Tmjson &j )
    {
        auto m = j["processes"];
        for ( auto p = m.begin(); p != m.end(); ++p )
        {
            cout << "# Reading process " << p.key() << " ... ";
            string bound = p.value()["bound"] | string();
            string free = p.value()["free"] | string();
            double pKd = p.value()["pKd"] | 0.0;
            double pX = p.value()["pX"] | 0.0;
            processdata d;
            d.id_AX = atom[bound].id;
            d.id_A = atom[free].id;
            d.set(pX, pKd);
            if ( d.id_AX != 0 && d.id_A != 0 )
            {
                process.push_back(d);
                cout << "OK!\n";
            }
            else
                cout << "ignored.\n";
        }

        // update reference states
        if ( !process.empty())
            for ( size_t i = 0; i < process.size() - 1; i++ )
            {
                for ( size_t j = i + 1; j < process.size(); j++ )
                {
                    Tid i_AX = process[i].id_AX;
                    Tid i_A = process[i].id_A;
                    Tid j_AX = process[j].id_AX;
                    Tid j_A = process[j].id_A;
                    if ( j_A == i_A )
                        process[j].set_mu_A(process[i].mu_A);
                    if ( j_A == i_AX )
                        process[j].set_mu_A(process[i].mu_AX);
                    if ( j_AX == i_A )
                        process[j].set_mu_AX(process[i].mu_A);
                    if ( j_AX == i_AX )
                        process[j].set_mu_AX(process[i].mu_AX);
                }
            }
    }

    /*!
     * This will go through the specified particle vector and
     * locate titratable sites. Their indexes will be stored
     * in the sites vector.
     */
    template<class Tpvec>
    void EquilibriumController::findSites( const Tpvec &p )
    {
        q.clear(); // clear average charge vector
        sites.clear(); // empty sites vector
        sites.reserve(p.size());

        for ( auto &prs : process )
            prs.cnt = 0;

        for ( size_t i = 0; i < p.size(); i++ )
            for ( size_t j = 0; j < process.size(); j++ )
                if ( process[j].one_of_us(p[i].id))
                {
                    sites.push_back(i);
                    process[j].cnt++;
                    break; // no double counting of sites
                }
    }

    /**
     * Returns the intrinsic energy, i.e. the ideal free energy connected with -log(Kd)
     * and the current state of the site. Explicit interactions with the surroundings
     * are not included.
     */
    double EquilibriumController::intrinsicEnergy( const Tid &id )
    {
        for ( auto &p : process )
            if ( p.one_of_us(id))
                return p.energy(id);
        return 0;
    }

    template<class Tpvec>
    void EquilibriumController::sampleCharge( const Tpvec &p )
    {
        for ( auto i : sites )
            q[i] += p[i].charge;
    }

    /**
     * This function will take the average charges from `q`
     * and apply these to the specified particle vector.
     *
     * @param p Particle vector to modify
     * @return Average net charge
     *
     * @warning After this function you can no longer perform any
     *          titration steps. It is meant to be called before saving coordinates
     *          to disk, so as to include the partial charges of the system.
     */
    template<class Tpvec>
    void EquilibriumController::copyAvgCharge( Tpvec &p )
    {
        for ( auto i : sites )
        {
            auto it = q.find(i);
            if ( it != q.end())
                p.at(i).charge = it->second.avg();
        }
    }

    /**
     * This function gives the average charge for all particles which
     * are titrated by process i, or -nan if no particles are part of process i
     */
    template<class Tpvec>
    double EquilibriumController::avgcharge( const Tpvec &p, int &k )
    {
        Average<double> qavg;
        for ( auto i : sites )
            if ( process[k].one_of_us(p[i].id))
                qavg += q[i].avg();
        return qavg.avg();
    }

    template<class Tpvec>
    EquilibriumController::processdata &EquilibriumController::random( const Tpvec &p, int &j )
    {
        int i = slump.range(0, sites.size() - 1);// pick random titratable site
        //int i=slump.rand() % sites.size();     // pick random titratable site
        j = sites[i];                                 // corresponding particle
        int k;
        do
            k = slump.range(0, process.size() - 1);// pick random process..
        while ( !process[k].one_of_us(p[j].id));   // ..that involves particle j
        return process[k];
    }

    string EquilibriumController::info( char w )
    {
        using namespace Faunus::textio;
        std::ostringstream o;
        o << pad(SUB, w, "Number of sites") << sites.size() << endl
          << indent(SUB) << "Processes:" << endl << endl;
        w = 8;
        if ( !process.empty())
        {
            o << indent(SUB)
              << setw(5) << "AX" << setw(5) << "<-> "
              << std::left << setw(5) << "A" << std::right << setw(w) << "pKd"
              << setw(w) << "pX"
              << setw(w) << "pAX"
              << setw(w) << "pA"
              << setw(w) << "N"
              << setw(12) << bracket("Z") << endl
              << indent(SUB) << string(70, '-') << endl;
            for ( size_t i = 0; i < process.size(); i++ )
            {
                o << indent(SUB)
                  << setw(5) << atom[process[i].id_AX].name
                  << setw(5) << "<-> " << setw(5) << std::left
                  << atom[process[i].id_A].name << std::right
                  << setw(w) << -log10(exp(-process[i].ddG))
                  << setw(w) << -log10(exp(-process[i].mu_X))
                  << setw(w) << -log10(exp(-process[i].mu_AX))
                  << setw(w) << -log10(exp(-process[i].mu_A))
                  << setw(w) << process[i].cnt;
                o << endl;
            }
        }
        return o.str();
    }

    /**
     * @brief Energy class for implicit titration of species
     *        used with Move::SwapMove.
     *
     *  This is a Hamiltonian for swapping atomic species according
     *  to their chemical potential and equilibrium constant as
     *  explained in `EquilibriumController`.
     */
    template<class Tspace>
    class EquilibriumEnergy : public Energybase<Tspace>
    {

    private:
        string _info() override { return eq.info(); }

    protected:
        std::map<int, double> energymap;//!< Intrinsic site energy

    public:
        EquilibriumController eq;

        EquilibriumEnergy( Tmjson &j ) : eq(j)
        {
            this->name = "Equilibrium State Energy";
        }

        auto tuple() -> decltype(std::make_tuple(this)) { return std::make_tuple(this); }

        template<class Tpvec>
        int findSites( const Tpvec &p )
        {
            eq.findSites(p);
            for ( auto &s : atom )
                for ( auto &process : eq.process )
                    if ( process.one_of_us(s.id))
                        energymap[s.id] = process.energy(s.id);
                    else
                        energymap[s.id] = 0;
            return eq.sites.size();
        }

        double i_internal( const typename Tspace::p_vec &p, int i ) override
        {
            return eq.intrinsicEnergy(p[i].id);
        }

        double g_internal( const typename Tspace::p_vec &p, Group &g ) override
        {
            double u = 0;
            for ( auto i : g )
                u += i_internal(p, i);
            return u;
        }

        void setSpace( Tspace &s ) override
        {
            Energybase<Tspace>::setSpace(s);
            findSites(s.p);
        }

    };

  }//Energy namespace 
}//Faunus namespace
#endif
