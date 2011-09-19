#include "faunus/common.h"
#include "faunus/point.h"
#include "faunus/space.h"
#include "faunus/energy.h"
#include "faunus/inputfile.h"
#include <faunus/titrate.h>
#include <faunus/species.h>

namespace Faunus {

  /*!
   * \param pX  Negative logarithm of the X activity (titrant)
   * \param pKd Negative logarithm of dissociation constant.
   */
  void equilibria::data::set(double pX, double pKd) {
    ddG = -log(pow(10., -pKd));
    mu_X= -log(pow(10., -pX));
    set_mu_AX(0);
  }

  void equilibria::data::set_mu_AX(double mu) {
    mu_AX = mu;
    mu_A = ddG  + mu_AX - mu_X;
  }

  void equilibria::data::set_mu_A(double mu) {
    mu_A = mu;
    mu_AX = mu_A + mu_X - ddG;
  }

  /*!
   * Returns true if the particle either matches AX or A.
   */
  bool equilibria::data::one_of_us(const short &id) {
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
  double equilibria::data::energy(const short &id) {
    if (id==id_AX)
      return mu_AX;
    if (id==id_A)
      return mu_A;
    return 0;
  }

  /*!
   * This will swap the state of given particle from AX<->A and
   * return the energy change associated with the process.
   */
  double equilibria::data::swap(particle &p) {
    double uold=energy(p.id);
    point pos=p;           // backup coordinate
    if (p.id==id_AX) {
      p=atom[id_A];
    }
    else if (p.id==id_A) {
      p=atom[id_AX];
    }
    p=pos;                 // apply old coordinates
    return energy(p.id)-uold; // return intrinsic energy change
  }

  /*!
   * Constructor. Should be called when particles have been loaded into the container.
   * Note that the particle charges, radii, weight etc. (save positions) will be reset
   * to the default species values defined in the atoms class.
   */
  equilibria::equilibria(space &s, inputfile &in, string pfx) { 
    string prefix=pfx;
    load(in.getstr("eqtit_processes","eq.in"));
    //atom.reset_properties(spc->p); // sets the faunatom particle properties into particle vector
    //s.trial=s.p;
    findSites(s.p);
    samplesites(s.p);
  }

  bool equilibria::load(string file) {
    process.clear();
    string t;
    data d;
    std::ifstream f(file.c_str());
    if (!f)
      return false;
    while (!f.eof()) {
      t.clear();
      f >> t;
      if (t=="Process") {
        string a,ax;
        double pkd, pX;
        f >> ax >> a >> pkd >> pX;
        d.id_AX=atom[ax].id;
        d.id_A=atom[a].id;
        d.set(pX, pkd);
        process.push_back(d);
      }
    }
    f.close();

    // update reference states
    char i_AX,i_A, j_AX, j_A;
    for (int i=0; i<process.size()-1; i++) {
      for (int j=i+1; j<process.size(); j++) {
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
  void equilibria::findSites(p_vec &p) {
    sites.clear();
    for (int j=0; j<process.size(); j++)
      process[j].cnt=0;
    for (int i=0; i<p.size(); i++)
      for (int j=0; j<process.size(); j++)
        if (process[j].one_of_us( p[i].id )) {
          sites.push_back(i);
          process[j].cnt++;
          break; // no double counting of sites
        }
    q.resize( sites.size() );
  }

  /*!
   * Returns the intrinsic energy, i.e. the ideal free energy connected with -log(Kd)
   * and the current state of the site. Explicit interactions with the surroundings
   * are not included.
   */
  double equilibria::intrinsicenergy(const short &id) {
    for (int i=0; i<process.size(); i++)
      if ( process[i].one_of_us(id) )
        return process[i].energy(id);
    return 0;
  }

  /*!
   * Calls intrinsicenergy(particle) on all particles.
   */
  double equilibria::intrinsicenergy(p_vec &p) {
    double u=0;
    for (int i=0; i<sites.size(); i++)
      u+=intrinsicenergy( p[sites[i]].id );
    return u;
  }

  /*!
   * \note Why is this public?
   */
  void equilibria::samplesites(p_vec &p) {
    for (unsigned int i=0; i<sites.size(); i++) {
      q[i] += p[sites[i]].charge;
    }
  }

  /*!
   * This function will take the average charges from titrate::q
   * and apply these to the specified particle vector. 
   * WARNING! After this function you can no longer perform any
   * titration steps. It is meant to be called before saving coordinates
   * to disk, so as to include the partial charges of the system.
   */
  double equilibria::applycharges(p_vec &p) {
    double qtot=0;
    for (unsigned int i=0; i<sites.size(); i++) {
      p[sites[i]].charge = q[i].avg();
      qtot += q[i].avg();
    }
    return qtot;
  }

  /*!
   * This function gives the average charge for all particles which 
   * are titrated by process i, or -nan if no particles are part of process i
   */  
  double equilibria::avgcharge(p_vec &p, int &i) {
    double qavg=0;
    int n=0;
    for (unsigned int j=0; j<sites.size(); j++) {
      if (process[i].one_of_us( p[sites[j]].id )) {
        n++;
        qavg+=q[j].avg();
      }
    }
    return qavg/n;
  }

  string equilibria::info() {
    char w=10;
    std::ostringstream o;
    o << "#   Number of sites           = " << sites.size() << endl;
    o << "#   Processes:\n"
      << "#    " << setw(5) << "AX" << setw(5) << "<-> "
      << std::left << setw(5) << "A" << std::right << setw(w) << "pKd"
      << setw(w) << "pX"
      << setw(w) << "pAX"
      << setw(w) << "pA" << endl
      << "#  -------------------------------------------------------------------"
      << endl;
    for (int i=0; i<process.size(); i++) {
      o << "#    "
        << setw(5) << atom[process[i].id_AX].name
        << setw(5) << "<-> " << setw(5) << std::left
        << atom[process[i].id_A].name << std:: right 
        << setw(w) << -log10(exp(-process[i].ddG))
        << setw(w) << -log10(exp(-process[i].mu_X))
        << setw(w) << -log10(exp(-process[i].mu_AX))
        << setw(w) << -log10(exp(-process[i].mu_A))
        << endl;
    }
    return o.str();
  }

  string equilibria::info(p_vec &p) {
    char w=8;
    std::ostringstream o;
    o << "#   Number of sites           = " << sites.size() << endl;
    o << "#   Processes:\n"
      << "#    " << setw(5) << "AX" << setw(5) << "<-> "
      << std::left << setw(5) << "A" << std::right << setw(w) << "pKd"
      << setw(w) << "pX"
      << setw(w) << "pAX"
      << setw(w) << "pA" 
      << setw(w) << "N"
      << setw(12) << "<Z>" << endl
      << "#  -----------------------------------------------------------------------------"
      << endl;
    for (int i=0; i<process.size(); i++) {
      o << "#    "
        << setw(5) << atom[process[i].id_AX].name
        << setw(5) << "<-> " << setw(5) << std::left
        << atom[process[i].id_A].name << std:: right 
        << setw(w) << -log10(exp(-process[i].ddG))
        << setw(w) << -log10(exp(-process[i].mu_X))
        << setw(w) << -log10(exp(-process[i].mu_AX))
        << setw(w) << -log10(exp(-process[i].mu_A))
        << setw(w) << process[i].cnt;
      if (process[i].cnt>0)
        o << setw(12) << avgcharge(p,i) << std::setprecision (4);
      o << endl;
    }
    return o.str();
  }
}//namespace

