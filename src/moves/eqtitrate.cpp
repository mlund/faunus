#include "faunus/common.h"
#include "faunus/point.h"
#include "faunus/moves/eqtitrate.h"
#include "faunus/species.h"
#include "faunus/container.h"
#include "faunus/energy/base.h"
#include "faunus/moves/base.h"
#include "faunus/ensemble.h"

namespace Faunus {
    
  /*!
   * \param pX  Negative logarithm of the X activity (titrant)
   * \param pKd Negative logarithm of dissociation constant.
   */
  void eqtitrate::data::set(double pX, double pKd) {
    ddG = -log(pow(10., -pKd));
    mu_X= -log(pow(10., -pX));
    set_mu_AX(0);
  }

  void eqtitrate::data::set_mu_AX(double mu) {
    mu_AX = mu;
    mu_A = ddG  + mu_AX - mu_X;
  }

  void eqtitrate::data::set_mu_A(double mu) {
    mu_A = mu;
    mu_AX = mu_A + mu_X - ddG;
  }

  /*!
   * Returns true if the particle either matches AX or A.
   */
  bool eqtitrate::data::one_of_us(particle &p) {
    if (p.id!=id_AX)
      if (p.id!=id_A)
        return false;
    return true;
  }

  /*!
   * Returns the intrinsic energy of the given particle. Intrinsic
   * means the energy stemming from the equilibrium expression when
   * no external interactions are accounted for (activity factors unity).
   */
  double eqtitrate::data::energy(particle &p) {
    if (p.id==id_AX)
      return mu_AX;
    if (p.id==id_A)
      return mu_A;
    return 0;
  }
 
  /*!
   * This will swap the state of given particle from AX<->A and
   * return the energy change associated with the process.
   */
  double eqtitrate::data::swap(particle &p) {
    double uold=energy(p);
    point pos=p;           // backup coordinate
    if (p.id==id_AX) {
      p=atom(id_A);
    }
    else if (p.id==id_A) {
      p=atom(id_AX);
    }
    p=pos;                 // apply old coordinates
    return energy(p)-uold; // return intrinsic energy change
  }

  eqtitrate::eqtitrate(ensemble &e, container &c, energybase &p, inputfile &in, string pfx) : markovmove(e,c,p) {
    name.assign("IMPLICIT EQUILIBRIUM AX=A+X TITRATION");
    prefix=pfx;
    markovmove::getInput(in);
    load(in.getstr("eqtit_processes","eq.in"));
    findSites(c.p);
    samplesites(con->p);
  }
  
  bool eqtitrate::load(string file) {
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
  void eqtitrate::findSites(vector<particle> &p) {
    sites.clear();
    for (int j=0; j<process.size(); j++)
      process[j].cnt=0;
    for (int i=0; i<p.size(); i++)
      for (int j=0; j<process.size(); j++)
        if (process[j].one_of_us(p[i])) {
          sites.push_back(i);
          process[j].cnt++;
          break; // no double counting of sites
        }
    q.resize( sites.size() );
  }
  
  double eqtitrate::move() {
    if (slp.runtest( runfraction )==false)
      return 0;
    double usum=0;
    int n=sites.size();
    while (n-->0) {
      markovmove::move();
      int i=slp.rand() % sites.size();     // pick random titratable site
      int j=sites[i];                      // corresponding particle
      int k;
      do {
        k=slp.rand() % process.size();                  // pick random process..
      } while (!process.at(k).one_of_us( con->p[j] ));  //   which matches particle j
      du  = process[k].swap( con->trial.at(j) );        // change state and get intrinsic energy change
      du += energy( con->trial, j) - energy( con->p, j);// get interactions energy change
      if (ens->metropolis(du)) {
        naccept++;
        utot+=du;
        con->p[j]=con->trial[j]; // accept new state
      } else {
        du=0;
        con->trial[j]=con->p[j]; // revert to old state
      }
      usum+=du;
    }
    samplesites(con->p);
    return usum;
  }

  /*!
   * Returns the intrinsic energy, i.e. the ideal free energy connected with -log(Kd)
   * and the current state of the site. Explicit interactions with the surroundings
   * are not included.
   */
  double eqtitrate::intrinsicenergy(particle &p) {
    for (int i=0; i<process.size(); i++)
      if ( process[i].one_of_us(p) )
        return process[i].energy(p);
    return 0;
  }
  
  /*!
   * Calls intrinsicenergy(particle) on all particles.
   */
  double eqtitrate::intrinsicenergy(vector<particle> &p) {
    double u=0;
    for (int i=0; i<sites.size(); i++)
      u+=intrinsicenergy(p[ sites[i] ]);
    return u;
  }
  
  /*! 
   * This virtual energy function is called whenever a swap move is made on
   * a particle. The default is to just call the corresponding energy function
   * in energybase but this can be customized for special needs - for example
   * if a change in charge changes the self-energy of a molecule (implicit
   * proton titration with Debye-Huckel).
   */
  double eqtitrate::energy(vector<particle> &p, int i) {
    return pot->energy(p,i);
  }

  void eqtitrate::samplesites(vector<particle> &p) {
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
  double eqtitrate::applycharges(vector<particle> &p) {
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
  double eqtitrate::avgcharge(vector<particle> &p, int &i) {
    double qavg=0;
    int cnt=0;
    for (unsigned int j=0; j<sites.size(); j++) {
      if (process[i].one_of_us(p[sites[j]])) {
        cnt++;
        qavg+=q[j].avg();
      }
    }
      return qavg/cnt;
  }
  
  string eqtitrate::info() {
    char w=10;
    std::ostringstream o;
    o << markovmove::info();
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

  string eqtitrate::info(vector<particle> &p) {
    char w=8;
    std::ostringstream o;
    o << markovmove::info();
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

