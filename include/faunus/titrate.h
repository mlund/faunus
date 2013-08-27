#ifndef FAUNUS_TITRATE_H
#define FAUNUS_TITRATE_H

#ifndef SWIG
#include <faunus/common.h>
#include <faunus/move.h>
#include <faunus/average.h>
#include <faunus/json.h>
#endif

namespace Faunus {

  namespace Energy {

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
    class EquilibriumController {
      private:
        bool includeJSON(const string&); //!< Read equilibrium processes
      public:
        //friend class EquilibriumEnergy;
        class processdata {
          //friend class EquilibriumController;
          public:
            double mu_AX;                //!< chemical potential of AX
            double mu_A;                 //!< chemical potential of A
            double mu_X;                 //!< chemical potential of X (this is the titrant)
            double ddG;                  //!< ddG = mu_A + mu_X - mu_AX
            int cnt;                     //!< number of sites for this process
          public:
            particle::Tid id_AX, id_A;   //!< Particle id's for AX and A
            bool one_of_us(const particle::Tid&);//!< Does the particle belong to this process?
            double energy(const particle::Tid&); //!< Returns intrinsic energy of particle id
            template<class Tparticle>
              double swap(Tparticle&);      //!< Swap AX<->A and return intrinsic energy change
            void set(double,double);     //!< Set activity of X and the pKd value
            void set_mu_AX(double);      //!< Set chemical potential of species AX - mu_A then follows.
            void set_mu_A(double);       //!< Set chemical potential of species A  - mu_AX then follows.
        };

        std::map<int, Average<double> > q;       //!< Map of average charges per site
        std::vector<processdata> process;        //!< Vector of processes.

        EquilibriumController(InputMap&, string="eq_");
        bool include(string);                    //!< Read equilibrium processes
        template<class Tpvec>
          void findSites(const Tpvec&);            //!< Locate all titratable sites
        double intrinsicEnergy(const particle::Tid&);    //!< Intrinsic energy of particle id (kT)
        string info(char=25);                    //!< Get information string
        template<class Tpvec>
          processdata& random(const Tpvec&, int&); //!< Random titratable particle and assiciated random process

        std::vector<int> sites;                  //!< List of titratable sites

        template<class Tpvec>
          void sampleCharge(const Tpvec&);         //!< Updates the average charge vector titrate::q
        template<class Tpvec>
          double applycharges(Tpvec &);            //!< Copy average charges to particles in the particle vector
        template<class Tpvec>
          double avgcharge(const Tpvec&, int&);    //!< Print average charges of process i
    };

    /*!
     * \param pX  Negative logarithm of the X activity (titrant)
     * \param pKd Negative logarithm of dissociation constant.
     */
    void EquilibriumController::processdata::set(double pX, double pKd) {
      ddG = -log(pow(10., -pKd));
      mu_X= -log(pow(10., -pX));
      set_mu_AX(0);
    }

    void EquilibriumController::processdata::set_mu_AX(double mu) {
      mu_AX = mu;
      mu_A = ddG  + mu_AX - mu_X;
    }

    void EquilibriumController::processdata::set_mu_A(double mu) {
      mu_A = mu;
      mu_AX = mu_A + mu_X - ddG;
    }

    /**
     * Returns `true` if the particle either matches AX or A.
     */
    bool EquilibriumController::processdata::one_of_us(const particle::Tid &id) {
      if (id==id_AX || id==id_A)
        return true;
      return false;
    }

    /*!
     * Returns the intrinsic energy of the given particle. Intrinsic
     * means the energy stemming from the equilibrium expression when
     * no external interactions are accounted for (activity factors unity).
     */
    double EquilibriumController::processdata::energy(const particle::Tid &id) {
      if (id==id_AX) return mu_AX;
      if (id==id_A)  return mu_A;
      return 0;
    }

    /**
     * This will swap the state of given particle from AX<->A and
     * return the energy change associated with the process.
     *
     * @note This will *NOT* swap particle radii nor masses!
     */
    template<class Tparticle>
      double EquilibriumController::processdata::swap(Tparticle &p) {
        double uold=energy(p.id);
        double oldradius=p.radius;
        double oldmw=p.mw;
        Point pos=p;           // backup coordinate
        if (p.id==id_AX)
          p=atom[id_A];
        else if (p.id==id_A)
          p=atom[id_AX];
        p=pos;                    // apply old coordinates
        p.radius=oldradius;       // apply old radius
        p.mw=oldmw;
        return energy(p.id)-uold; // return intrinsic energy change
      }

    /**
     * @brief Construct from `InputMap`.
     *
     * Call this *after* particles have been loaded into `Space`, i.e.
     * typically just before starting the Markov chain. Also make
     * sure that `AtomMap` has been loaded with all atomic properties
     * as these will be used to reset the charge, radii, weight etc.
     * on all particles in the system.
     *
     * The `InputMap` is searched for:
     *
     * - `processfile` Name of process file
     */
    EquilibriumController::EquilibriumController(InputMap &in, string pfx) {
      string prefix=pfx;
      include( in.get<string>(prefix+"processfile","eq.process"));
    }

    bool EquilibriumController::includeJSON(const string &file) {
      auto j=json::open(file);
      for (auto &p : json::object("processes", j)) {
        cout << "Reading process " << p.first << " ... ";
        string bound = json::value<string>(p.second, "bound", string());
        string free  = json::value<string>(p.second, "free", string());
        double pKd   = json::value<double>(p.second, "pKd", 0);
        double pX    = json::value<double>(p.second, "pX", 0);
        processdata d;
        d.id_AX=atom[bound].id;
        d.id_A=atom[free].id;
        d.set(pX, pKd);
        if (d.id_AX!=0 && d.id_A!=0) {
          process.push_back(d);
          cout << "OK!\n";
        }
        else
          cout << "ignored.\n";
      }
      return (process.empty() ? false : true);
    }

    bool EquilibriumController::include(string file) {
      if (file.substr(file.find_last_of(".") + 1) == "json") {
        includeJSON(file);
        // update reference states
        particle::Tid i_AX,i_A, j_AX, j_A;
        if (!process.empty())
          for (size_t i=0; i<process.size()-1; i++) {
            for (size_t j=i+1; j<process.size(); j++) {
              i_AX = process[i].id_AX;
              i_A  = process[i].id_A;
              j_AX = process[j].id_AX;
              j_A  = process[j].id_A;
              if ( j_A  == i_A  ) process[j].set_mu_A ( process[i].mu_A  );
              if ( j_A  == i_AX ) process[j].set_mu_A ( process[i].mu_AX );
              if ( j_AX == i_A  ) process[j].set_mu_AX( process[i].mu_A  );
              if ( j_AX == i_AX ) process[j].set_mu_AX( process[i].mu_AX );
            }
          }
        return true;
      }
      return false;
    }

    /*!
     * This will go through the specified particle vector and
     * locate titratable sites. Their indexes will be stored
     * in the sites vector.
     */
    template<class Tpvec>
      void EquilibriumController::findSites(const Tpvec &p) {
        q.clear(); // clear average charge vector
        sites.clear(); // empty sites vector
        for (auto &prs : process)
          prs.cnt=0;

        for (size_t i=0; i<p.size(); i++)
          for (size_t j=0; j<process.size(); j++)
            if ( process[j].one_of_us( p[i].id )) {
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
    double EquilibriumController::intrinsicEnergy(const particle::Tid &id) {
      for (auto &p : process)
        if ( p.one_of_us(id) )
          return p.energy(id);
      return 0;
    }

    template<class Tpvec>
      void EquilibriumController::sampleCharge(const Tpvec &p) {
        for (auto i : sites)
          q[i] += p[i].charge;
      }

    /**
     * This function will take the average charges from `q`
     * and apply these to the specified particle vector.
     *
     * @warning After this function you can no longer perform any
     *          titration steps. It is meant to be called before saving coordinates
     *          to disk, so as to include the partial charges of the system.
     */
    template<class Tpvec>
      double EquilibriumController::applycharges(Tpvec &p) {
        double qtot=0;
        for (auto i : sites) {
          if (q[i].cnt>0)
            p[i].charge = q[i].avg();
          qtot += q[i].avg();
        }
        return qtot;
      }

    /**
     * This function gives the average charge for all particles which
     * are titrated by process i, or -nan if no particles are part of process i
     */
    template<class Tpvec>
      double EquilibriumController::avgcharge(const Tpvec &p, int &k) {
        Average<double> qavg;
        for (auto i : sites)
          if (process[k].one_of_us( p[i].id ))
            qavg+=q[i].avg();
        return qavg.avg();
      }

    template<class Tpvec>
      EquilibriumController::processdata& EquilibriumController::random(const Tpvec &p, int &j) {
        int i=slp_global.rand() % sites.size();     // pick random titratable site
        j=sites[i];                                 // corresponding particle
        int k;
        do
          k=slp_global.rand() % process.size();     // pick random process..
        while (!process[k].one_of_us( p[j].id ));   // ..that involves particle j
        return process[k];
      }

    string EquilibriumController::info(char w) {
      using namespace Faunus::textio;
      std::ostringstream o;
      o << pad(SUB,w,"Number of sites") << sites.size() << endl
        << indent(SUB) << "Processes:" << endl << endl;
      w=8;
      if (!process.empty()) {
        o << indent(SUB)
          << setw(5) << "AX" << setw(5) << "<-> "
          << std::left << setw(5) << "A" << std::right << setw(w) << "pKd"
          << setw(w) << "pX"
          << setw(w) << "pAX"
          << setw(w) << "pA"
          << setw(w) << "N"
          << setw(12) << bracket("Z") << endl
          << indent(SUB) << string(70,'-') << endl;
        for (size_t i=0; i<process.size(); i++) {
          o << indent(SUB)
            << setw(5) << atom[process[i].id_AX].name
            << setw(5) << "<-> " << setw(5) << std::left
            << atom[process[i].id_A].name << std:: right
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
      class EquilibriumEnergy : public Energybase<Tspace> {
        
        private:
          string _info() { return eq.info(); }
        
        protected:
          std::map<particle::Tid, double> energymap;//!< Intrinsic site energy
        
        public:
          EquilibriumController eq;
        
          EquilibriumEnergy(InputMap &in) : eq(in) {
            this->name="Equilibrium State Energy";
          }

          template<class Tpvec>
            int findSites(const Tpvec &p) {
              eq.findSites(p);
              for (auto &s : atom.list )
                for (auto &process : eq.process )
                  if ( process.one_of_us( s.id ) )
                    energymap[s.id] = process.energy( s.id );
                  else energymap[s.id]=0;
              return eq.sites.size();
            }

          double i_internal(const typename Tspace::p_vec &p, int i) FOVERRIDE {
            return eq.intrinsicEnergy( p[i].id );
          }

          double g_internal(const typename Tspace::p_vec &p, Group &g) FOVERRIDE {
            double u=0;
            for (auto i : g)
              u+=i_internal(p, i);
            return u;
          }
      };

  }//Energy namespace 

  namespace Move {

    /**
     * @brief Move for swapping species types - i.e. implicit titration
     *
     * Upon construction this class will add an instance of
     * Energy::EquilibriumEnergy to the Hamiltonian. For details
     * about the titration procedure see Energy::EquilibriumController.
     */
    template<class Tspace>
      class SwapMove : public Movebase<Tspace> {
        private:
          std::map<int, Average<double> > accmap; //!< Site acceptance map
          string _info();
          void _trialMove();
          void _acceptMove();
          void _rejectMove();
        protected:
          using Movebase<Tspace>::spc;
          using Movebase<Tspace>::pot;
          double _energyChange();
          int ipart;                              //!< Particle to be swapped
          Energy::EquilibriumEnergy<Tspace>* eqpot;
        public:
          SwapMove(InputMap&, Energy::Energybase<Tspace>&, Tspace&, Energy::EquilibriumEnergy<Tspace>&, string="swapmv_"); //!< Constructor
          template<class Tpvec>
            int findSites(const Tpvec&); //!< Search for titratable sites (old ones are discarded)
          double move();
          template<class Tpvec>
            void applycharges(Tpvec &);
        
#ifdef ENABLE_MPI
          Faunus::MPI::MPIController* mpi;
#endif
      };

    /**
     * This will set up swap move routines and search for
     * titratable sites in `Space`.
     */
    template<class Tspace>
      SwapMove<Tspace>::SwapMove(
          InputMap &in,
          Energy::Energybase<Tspace> &ham,
          Tspace &spc,
          Energy::EquilibriumEnergy<Tspace> &eq,
          string pfx) : Movebase<Tspace>(ham,spc,pfx) {

        this->title="Site Titration - Swap Move";
        this->runfraction=in.get<double>(pfx+"runfraction",1);
        eqpot=&eq;
        ipart=-1;

        findSites(spc.p);

        /* Sync particle charges with `AtomMap` */
        for (auto i : eqpot->eq.sites)
          spc.trial[i].charge = spc.p[i].charge = atom[ spc.p[i].id ].charge;
#ifdef ENABLE_MPI
        mpi=nullptr;
#endif
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
          int i=slp_global.rand() % eqpot->eq.sites.size(); // pick random site
          ipart=eqpot->eq.sites.at(i);                      // and corresponding particle
          int k;
          do {
            k=slp_global.rand() % eqpot->eq.process.size(); // pick random process..
          } while (!eqpot->eq.process[k].one_of_us( this->spc->p[ipart].id )); //that match particle j
          eqpot->eq.process[k].swap( this->spc->trial[ipart] ); // change state and get intrinsic energy change
        }
      }

    template<class Tspace>
      double SwapMove<Tspace>::_energyChange() {
        assert( spc->geo.collision(this->spc->p[ipart])==false
            && "Accepted particle collides with container");
        if (spc->geo.collision(spc->trial[ipart]))  // trial<->container collision?
          return pc::infty;
        
#ifdef ENABLE_MPI
        if (mpi!=nullptr) {
          double sum=0;
          auto r = Faunus::MPI::splitEven(*mpi, (int)spc->p.size());
          for (int i=r.first; i<=r.second; ++i)
            if (i!=ipart)
              sum+=pot->i2i(spc->trial,i,ipart) - pot->i2i(spc->p,i,ipart);
          sum = Faunus::MPI::reduceDouble(*mpi, sum);
          
          return sum + pot->i_external(spc->trial, ipart) - pot->i_external(spc->p, ipart)
          + pot->i_internal(spc->trial, ipart) - pot->i_internal(spc->p, ipart);
        }
#endif
        return pot->i_total(spc->trial,ipart) - pot->i_total(spc->p,ipart);
      }

    template<class Tspace>
      void SwapMove<Tspace>::_acceptMove() {
        accmap[ipart] += 1;
        spc->p[ipart] = spc->trial[ipart];
      }

    template<class Tspace>
      void SwapMove<Tspace>::_rejectMove() {
        accmap[ipart] += 0;
        spc->trial[ipart] = spc->p[ipart];
      }

    template<class Tspace>
      double SwapMove<Tspace>::move() {
        double du=0;
        if (this->run()) {
          size_t i=eqpot->eq.sites.size();
          while (i-->0)
            du+=Movebase<Tspace>::move();
          eqpot->eq.sampleCharge(spc->p);
        }
        return du;
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
        if (this->cnt>0 && !eqpot->eq.sites.empty()) {
          o << indent(SUB) << "Site statistics:" << endl
            << indent(SUBSUB) << std::left
            << setw(16) << "Site"
            << setw(14) << bracket("z")
            << "Acceptance" << endl;
          for (auto i : eqpot->eq.sites) {
            if (accmap[i].cnt>0) {
              std::ostringstream a;
              o.precision(5);
              o.setf( std::ios::fixed, std::ios::floatfield );
              a << std::left << setw(5) << atom[ spc->p[i].id ].name << std::right << setw(5) << i;
              o << pad(SUBSUB,15, a.str())
                << setw(8) << std::right << eqpot->eq.q[i].avg()
                << setw(11) << std::right << accmap[i].avg()*100. << " " << percent
                << endl;
            }
          }
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

  }// Move namespace
}//Faunus namespace
#endif
