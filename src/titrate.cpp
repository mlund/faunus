#include "faunus/common.h"
#include "faunus/point.h"
#include "faunus/space.h"
#include "faunus/energy.h"
#include "faunus/inputfile.h"
#include <faunus/titrate.h>
#include <faunus/species.h>

namespace Faunus {

  namespace Energy {

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
    double EquilibriumController::processdata::swap(particle &p) {
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
     * sure that `AtomTypes` has been loaded with all atomic properties
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

    bool EquilibriumController::include(string file) {
      using namespace textio;
      string t;
      processdata d;
      std::ifstream f(file.c_str());
      cout << "Reading titration process data from '" << file << "'. ";
      if (!f) {
        cout << "FAILED!\n";
        return false;
      }
      cout << "OK!\n";
      while (!f.eof()) {
        t.clear();
        f >> t;
        if (t=="Process") {
          string a,ax;
          double pkd, pX;
          f >> ax >> a >> pkd >> pX;
          if (atom[ax].id==0 || atom[a].id==0) {
            cout << indent(SUB) << "Warning: species " << ax << " or " << a << " are unknown and will be ignored." << endl;
          } else {
            d.id_AX=atom[ax].id;
            d.id_A=atom[a].id;
            d.set(pX, pkd);
            process.push_back(d);
          }
        }
      }
      f.close();

      // update reference states
      particle::Tid i_AX,i_A, j_AX, j_A;
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

    /*!
     * This will go through the specified particle vector and
     * locate titratable sites. Their indexes will be stored
     * in the sites vector.
     */
    void EquilibriumController::findSites(const p_vec &p) {
      q.clear(); // clear average charge vector
      sites.clear(); // empty sites vector
      for (auto &prs : process)
        prs.cnt=0;

      int mismatch=0; // number of charge mismatches
      for (size_t i=0; i<p.size(); i++)
        for (size_t j=0; j<process.size(); j++)
          if ( process[j].one_of_us( p[i].id )) {
            if ( abs(p[i].charge-atom[p[i].id].charge)>1e-6 )
              mismatch++;
            sites.push_back(i);
            process[j].cnt++;
            break; // no double counting of sites
          }
      if (mismatch>0)
        std::cerr
          << "Warning: Found " << mismatch
          << " mismatched charge(s) while searching for titratable sites.\n";
    }

    /*!
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

    void EquilibriumController::sampleCharge(const p_vec &p) {
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
    double EquilibriumController::applycharges(p_vec &p) {
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
    double EquilibriumController::avgcharge(const p_vec &p, int &k) {
      Average<double> qavg;
      for (auto i : sites)
        if (process[k].one_of_us( p[i].id ))
          qavg+=q[i].avg();
      return qavg.avg();
    }

    EquilibriumController::processdata& EquilibriumController::random(const p_vec &p, int &j) {
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
      return o.str();
    }

    EquilibriumEnergy::EquilibriumEnergy(InputMap &in) : eq(in) {
      name="Equilibrium State Energy";
    }

    int EquilibriumEnergy::findSites(const p_vec &p) {
      eq.findSites(p);
      for (auto &s : atom.list )
        for (auto &process : eq.process )
          if ( process.one_of_us( s.id ) )
            energymap[s.id] = process.energy( s.id );
          else energymap[s.id]=0;
      return eq.sites.size();
    }

    double EquilibriumEnergy::i_internal(const p_vec &p, int i) {
      return eq.intrinsicEnergy( p[i].id );
    }

    double EquilibriumEnergy::g_internal(const p_vec &p, Group &g) {
      double u=0;
      for (auto i : g)
        u+=EquilibriumEnergy::i_internal(p, i);
      return u;
    }

    string EquilibriumEnergy::_info() {
      return eq.info();
    }

  }//namespace Energy

  namespace Move {

    /**
     * This will set up swap move routines and search for
     * titratable sites in `Space`.
     * @warning This will use properties from `AtomTypes` to
     *          override those already stored in the particle
     *          vector.
     */
    SwapMove::SwapMove(
        InputMap &in,
        Energy::Hamiltonian &ham,
        Space &spc,
        string pfx) : Movebase(ham,spc,pfx), eqpot(in) {

      title="Site Titration - Swap Move";
      runfraction=in.get<double>(pfx+"runfraction",1);
      ham.add(eqpot);
      ipart=-1;
      findSites(spc.p);

      /* Sync particles with `AtomTypes` */
      for (auto &i : eqpot.eq.sites ) {
        spc.p[i] = atom[ spc.p[i].id ]; 
        spc.trial[i] = spc.p[i];
      }
    }

    int SwapMove::findSites(const p_vec &p) {
      accmap.clear();
      return eqpot.findSites(p);
    }

    void SwapMove::_trialMove() {
      int i=slp_global.rand() % eqpot.eq.sites.size(); // pick random site
      ipart=eqpot.eq.sites.at(i);                      // and corresponding particle
      int k;
      do {
        k=slp_global.rand() % eqpot.eq.process.size(); // pick random process..
      } while (!eqpot.eq.process[k].one_of_us( spc->p[ipart].id )); //that match particle j
      eqpot.eq.process[k].swap( spc->trial[ipart] ); // change state and get intrinsic energy change
    }

    double SwapMove::_energyChange() {
      assert( spc->geo->collision(spc->p[ipart])==false
          && "Accepted particle collides with container");
      if (spc->geo->collision(spc->trial[ipart]))  // trial<->container collision?
        return pc::infty;
      return pot->i_total(spc->trial,ipart) - pot->i_total(spc->p,ipart);
    }

    void SwapMove::_acceptMove() {
      accmap[ipart] += 1;
      spc->p[ipart] = spc->trial[ipart];
    }

    void SwapMove::_rejectMove() {
      accmap[ipart] += 0;
      spc->trial[ipart] = spc->p[ipart];
    }

    double SwapMove::move() {
      double du=0;
      if (run()) {
        size_t i=eqpot.eq.sites.size();
        while (i-->0)
          du+=Movebase::move();
        eqpot.eq.sampleCharge(spc->p);
      }
      return du;
    }

    void SwapMove::applycharges(p_vec &p){
      eqpot.eq.applycharges(p);
    }

    string SwapMove::_info() {
      using namespace textio;
      std::ostringstream o;
      if (cnt>0) {
        o << indent(SUB) << "Site statistics:" << endl
          << indent(SUBSUB) << std::left
          << setw(15) << "Site"
          << setw(14) << bracket("z")
          << "Acceptance" << endl;
        for (auto i : eqpot.eq.sites) {
          if (accmap[i].cnt>0) {
            std::ostringstream a;
            o.precision(4);
            a << atom[ spc->p[i].id ].name << " " << i;
            o << pad(SUBSUB,15, a.str())
              << setw(10) << eqpot.eq.q[i].avg()
              << accmap[i].avg()*100. << percent
              << endl;
          }
        }
      }
      return o.str();
    }

    SwapMoveMSR::SwapMoveMSR(
        InputMap &in, Energy::Hamiltonian &ham, Space &spc, string pfx) : SwapMove(in,ham,spc,pfx) {
      title+=" (min. shortrange)";
      ham.add( potrest );
    }

    void SwapMoveMSR::modify() {
      radiusbak.clear();
      hydrophobicbak.clear();
      for (auto g : spc->groupList() )   // loop over all groups
        if (g->find(ipart)) {  //   is ipart part of a group?
          for (auto i : *g)    //     if so, loop over that group
            if (i!=ipart) {    //       and ignore ipart
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

    void SwapMoveMSR::restore() {
      for (auto &m : radiusbak) {
        spc->p[m.first].radius = m.second;
        spc->trial[m.first].radius = m.second;
      }
      for (auto &m : hydrophobicbak) {
        spc->p[m.first].hydrophobic = m.second;
        spc->trial[m.first].hydrophobic = m.second;
      }
    }

    double SwapMoveMSR::_energyChange() {
      double du_orig = SwapMove::_energyChange();
      modify();
      double du = SwapMove::_energyChange();
      restore();
      potrest.add( du-du_orig );
      return du;
    }

  }//Move namespace
}//Faunus namespace

