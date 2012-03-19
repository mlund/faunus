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

    /*!
     * Returns true if the particle either matches AX or A.
     */
    bool EquilibriumController::processdata::one_of_us(const short &id) {
      if (id!=id_AX)
        if (id!=id_A)
          return false;
      return true;
    }

    /*!
     * Returns the intrinsic energy of the given particle. Intrinsic
     * means the energy stemming from the equilibrium expression when
     * no external interactions are accounted for (activity factors unity).
     */
    double EquilibriumController::processdata::energy(const short &id) {
      if (id==id_AX)
        return mu_AX;
      if (id==id_A)
        return mu_A;
      return 0;
    }

    /*!
     * This will swap the state of given particle from AX<->A and
     * return the energy change associated with the process.
     * \note This will NOT swap particle radii nor masses!
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

    /*!
     * Constructor. Should be called when particles have been loaded into the container.
     * Note that the particle charges, radii, weight etc. (save positions) will be reset
     * to the default species values defined in the atoms class.
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
      short i_AX,i_A, j_AX, j_A;
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
      q.clear();
      sites.clear();
      for (size_t j=0; j<process.size(); j++)
        process[j].cnt=0;
      for (size_t i=0; i<p.size(); i++)
        for (size_t j=0; j<process.size(); j++)
          if ( process[j].one_of_us( p[i].id )) {
            if ( abs(p[i].charge-atom[p[i].id].charge)>1e-9 )
              cout << "Warning: Charge mismatch while searching for titratable sites - this can be serious.\n";
            sites.push_back(i);
            process[j].cnt++;
            break; // no double counting of sites
          }
    }

    /*!
     * Returns the intrinsic energy, i.e. the ideal free energy connected with -log(Kd)
     * and the current state of the site. Explicit interactions with the surroundings
     * are not included.
     */
    double EquilibriumController::intrinsicEnergy(const short &id) {
      for (auto &p : process) 
        if ( p.one_of_us(id) )
          return p.energy(id);
      return 0;
    }

    void EquilibriumController::sampleCharge(const p_vec &p) {
      for (auto i : sites)
        q[i] += p[i].charge;
    }

    /*!
     * This function will take the average charges from titrate::q
     * and apply these to the specified particle vector. 
     * \warning After this function you can no longer perform any
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

    /*!
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
      //if (energymap.find( p[i].id )!=energymap.end())
      //  return energymap[ p[i].id ];
      //return 0;
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

    /*!
     * This will set up the swap move routines and the particles in the space are automatically
     * searched for titratable sites. Before starting the titration the charge of titratable
     * particles must match that loaded into Faunus::atoms -- this is also handled by the
     * constructor.
     */
    SwapMove::SwapMove(
        InputMap &in,
        Energy::Hamiltonian &ham,
        Space &spc,
        string pfx) : Movebase(ham,spc,pfx), eqpot(in) {

      title="Site Titration - Swap Move";
      ham.add(eqpot);
      ipart=-1;
      findSites(spc.p);
      for (auto &i : eqpot.eq.sites ) // make sure charges are in sync w. atom list
        spc.p[i].charge = spc.trial[i].charge = atom[ spc.p[i].id ].charge;
    }

    int SwapMove::findSites(const p_vec &p) {
      accmap.clear();
      return eqpot.findSites(p);
    }

    void SwapMove::_trialMove() {
      int i=slp_global.rand() % eqpot.eq.sites.size();     // pick random titratable site
      ipart=eqpot.eq.sites.at(i);                          // corresponding particle
      int k;
      do {
        k=slp_global.rand() % eqpot.eq.process.size();               // pick random process..
      } while (!eqpot.eq.process[k].one_of_us( spc->p[ipart].id ));  //   which matches particle j
      eqpot.eq.process[k].swap( spc->trial[ipart] );       // change state and get intrinsic energy change
    }

    double SwapMove::_energyChange() {
      assert( spc->geo->collision(spc->p[ipart])==false && "An accepted particle collides with simulation container.");
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

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
    double SwapMove::move() {
      if (!run())
        return 0;
      double du=0;
      for (auto &n : eqpot.eq.sites )
        du+=Movebase::move();
      eqpot.eq.sampleCharge(spc->p);
      return du;
    }
#pragma GCC diagnostic pop

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
          std::ostringstream a;
          o.precision(4);
          a << atom[ spc->p[i].id ].name << " " << i;
          o << pad(SUBSUB,15, a.str())
            << setw(10) << eqpot.eq.q[i].avg()
            << accmap[i].avg()*100. << percent
            << endl;
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
      for (auto g : spc->g )   // loop over all groups
        if (g->find(ipart)) {  //   is ipart part of a group?
          for (auto i : *g)    //     if so, loop over that group
            if (i!=ipart) {    //       and ignore ipart
              assert( abs(spc->p[i].radius-spc->trial[i].radius)<1e-9 && "Untouched radii must be in sync!" );
              assert( spc->p[i].hydrophobic==spc->trial[i].hydrophobic && "Untouched particles changed.");

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

