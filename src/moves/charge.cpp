#include <faunus/moves/charge.h>
#include <faunus/titrate.h>

namespace Faunus {
 //---------- CHARGE REG ---------------------
  string chargereg::info() {
    std::ostringstream o;
    o <<  markovmove::info()
      << "#   pH (concentration)  = " << ph << endl
      << "#   Titrateable sites   = " << sites.size() << endl
      << "#   Number of protons   = " << protons.size() << endl;
    return o.str();
  }

  chargereg::chargereg(ensemble &e, container &c, energybase &i, group &g, float ph ) : markovmove(e,c,i), titrate(c.p,g,ph)
  {
    name.assign("PROTON TITRATION");
    cite.assign("Biochem. 2005, 44, 5722-5727.");
    runfraction=0.2;
    con->trial = con->p;
  }

  /*! \brief Exchange protons between bulk and titrateable sites.
   *
   *  This move will randomly go through the titrateable sites and
   *  try to exchange protons with the bulk. The trial energy is:
   */
  double chargereg::titrateall() {
    du=0;
    if (slp.runtest(runfraction)==false)
      return du;
    action t;
    double sum=0;
    for (unsigned short i=0; i<sites.size(); i++) {
      cnt++;
      t=exchange(con->trial);
      //#pragma omp parallel
      {
        //#pragma omp sections
        {
          //#pragma omp section
          {
            uold = pot->potential( con->p, t.site ) * con->p[t.site].charge
              + pot->potential( con->p, t.proton ) * con->p[t.proton].charge
              - con->p[t.site].charge * con->p[t.proton].charge * pot->tokT /
              sqrt( con->sqdist( con->p[t.site], con->p[t.proton]) );
          }
          //#pragma omp section
          {
            unew = pot->potential( con->trial, t.site ) * con->trial[t.site].charge
              + pot->potential( con->trial, t.proton ) * con->trial[t.proton].charge
              - con->trial[t.site].charge * con->trial[t.proton].charge * pot->tokT / 
              sqrt( con->sqdist( con->trial[t.site], con->trial[t.proton]) );
          }
        }
      }
      du = (unew-uold);

      if (ens->metropolis( energy(con->trial, du, t) )==true) {
        rc=OK;
        utot+=du;
        naccept++;
        con->p[t.site].charge   = con->trial[t.site].charge;
        con->p[t.proton].charge = con->trial[t.proton].charge;
        sum+=du;
      } else {
        rc=ENERGY;
        exchange(con->trial, t);
      }
      samplesites(con->p);  // Average charges on all sites
    }
    return sum;
  }

  //---------- INHERITED GCTITRATE ---------
  /*!
   * \param ph pH value
   * \param mu Proton excess chemical potential
   */
  HAchargereg::HAchargereg(ensemble &e,
      container &c,
      energybase &i,
      group &g, float ph, float mu ) : chargereg(e,c,i,g,ph)
  {
    this->name.assign("GC PROTON TITRATION...");
    this->cite.assign("Labbez+Jonsson....");
    CatPot=mu;
  }

  string HAchargereg::info() {
    std::ostringstream o;
    o << this->chargereg::info();
    o << "#   Excess chem. pot.   = " << CatPot << endl;
    return o.str();
  }

  double HAchargereg::energy( vector<particle> &p, double du, titrate::action &a ) {
    int i=p[a.site].id;
    if (a.action==this->PROTONATED)
      return du+( log(10.)*( this->ph - atom.list[i].pka ) ) + 
        CatPot - log( this->protons.size() / this->con->getvolume() ) ;
    else
      return du-( log(10.)*( this->ph - atom.list[i].pka ) ) -
        CatPot + log( (this->protons.size()+1) / this->con->getvolume() );
  }


  // ----------- DH Titration ----------------

 DHchargereg::DHchargereg(ensemble &e, container &c, energybase &i, float ph, float muH ) : markovmove(e,c,i), titrate_implicit(c.p,ph,muH)
  {
    name.assign("DH-PROTON TITRATION");
    //cite.assign("...");
    runfraction=0.2;
    con->trial = con->p;
  }
  /*! \brief Exchange protons between bulk and titrateable sites.
   *
   *  This move will randomly go through the titrateable sites and
   *  try to exchange protons with the bulk. The trial energy is:
   */
  double DHchargereg::titrateall() {
    du=0;
    if (slp.runtest(runfraction)==false)
      return du;
    double sum=0;
    for (unsigned short i=0; i<sites.size(); i++) {
      cnt++;
      int k=exchange(con->trial);

      uold = pot->energy( con->p, k );
      unew = pot->energy( con->trial, k );
      du = (unew-uold);

      if (ens->metropolis( energy(con->trial, du, k) )==true) {
        rc=OK;
        utot+=du;
        naccept++;
        con->p[k].charge = con->trial[k].charge;
        sum+=du;
      } else {
        rc=ENERGY;
        exchange(con->trial, k);
      }
      //samplesites(con->p);  // Average charges on all sites
    }
    return sum;
  }

  string DHchargereg::info() {
    std::ostringstream o;
    o <<  markovmove::info()
      << "#   pH (concentration)  = " << ph << endl
      << "#   Proton excess (kT)  = " << mu_proton << endl
      << "#   Titrateable sites   = " << sites.size() << endl;
    return o.str();
  }

}//namespace
