#include <faunus/moves/switch.h>
#include <faunus/titrate.h>
#include <faunus/physconst.h>
#include "faunus/energy/base.h"
#include "faunus/ensemble.h"

namespace Faunus {
  
 //---------- SWITCH ---------------------
  string switch_type::info() {
    std::ostringstream o;
    o <<  markovmove::info()
      << "#   pc (concentration)        = " << pc << endl
      << "#   pK                        = " << pK << endl
      << "#   Titrateable sites         = " << sites.size() << endl
      << "#   Types:  "<< titrate_switch::fType1 << "<-->" << titrate_switch::fType2 <<endl
      << "#   Fraction of "<<fType2<<"    ="<<frac_2.avg()<<endl;
    return o.str();
  }

  switch_type::switch_type(ensemble &e, container &c, energybase &i, inputfile &in) : markovmove(e,c,i), titrate_switch(c,in)
  {
    runfraction=1.0;
    prefix.assign("sw_");
    name.assign("Type switcher");
    cite.assign("None so far.");
    con->trial = con->p;
  }

  /*! \brief Exchange protons between bulk and titrateable sites.
   *
   *  This move will randomly go through the titrateable sites and
   *  try to exchange protons with the bulk. The trial energy is:
   */
  double switch_type::titrateall() {
    du=0;
    double dNmu=0;
    double nT1, nT2;                                 // number of each type in initial config.
    nT1=nT2=0.0;

    for (unsigned short i=0; i<sites.size(); i++)    // Correction to get a full hamiltonian.
      if (con->p[sites[i]].id==fType2)
        nT2++;
    dNmu+=nT2*log(10.)*(pc-pK);                       // The energy of fType1 states are choosen to 0. Subtract inital energy (note sign)
                                                     // See faunus::titrate_switch 
    double sum=0;
    for (unsigned short i=0; i<sites.size(); i++) {
      cnt++;
      exchange(con->trial);
      //#pragma omp parallel
      {
        //#pragma omp sections
        {
          //#pragma omp section
          {
            uold = pot->energy( con->p, recent );
          }
          //#pragma omp section
          {
            unew = pot->energy( con->trial, recent );
          }
        }
      }
      du = (unew-uold);

      if (ens->metropolis( energy(*con, du) )==true) {
        rc=OK;
        utot+=du;
        naccept++;
        exchange(con->p,recent);
        sum+=du;
      } else {
        rc=ENERGY;
        exchange(con->trial, recent);
      }
      samplesites(*con);  // Average charges on all sites
    }
    for (unsigned short i=0; i<sites.size(); i++)    
      if (con->p[sites[i]].id==fType2)
        nT2++;
    dNmu-=nT2*log(10.)*(pc-pK);                    // Add ideal term of final state (note sign)   

    return (full_ham==true ? sum + dNmu : sum);
  }

  void switch_type::samplesites(container &c) {
    n=0;
    for(int i=0; i<sites.size(); i++) {
      if (c.p[sites[i]].id==fType2)
        n++;
    }
    frac_2.add(double(n)/double(sites.size()));
  }
}//namespace
