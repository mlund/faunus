#include <faunus/moves/charge.h>
#include <faunus/titrate.h>
#include <faunus/physconst.h>
#include "faunus/energy/base.h"
#include "faunus/ensemble.h"

namespace Faunus {
  
  //---------- CHARGE REG GAUSSIAN ---------------------
  chargeregGaussian::chargeregGaussian(ensemble &e, container &c, energybase &i) : markovmove(e,c,i) {
    name.assign("Gaussian Charge Fluctuations");
    cite.assign("ML, unpublished");
    prefix="chargereggauss";
    dp=1.0;
    data t;
    for (unsigned int i=0; i<c.p.size(); i++) // search for titratable sites
      if ( atom[ c.p[i].id ].variance>1e-4) {
        t.n=i;
        sites.push_back(t);
      }
    if (sites.size()==0)
      runfraction=0;
  }
  
  string chargeregGaussian::info() {
    std::ostringstream o;
    if (dp>0)
      o <<  markovmove::info()
        << "#   Nr. of titratable sites   = " << sites.size() << endl;
    return o.str();
  }
  
  double chargeregGaussian::titrateall() {
    double sum=0;
    for (unsigned int i=0; i<sites.size(); i++) {
      markovmove::move();
      unsigned int j=markovmove::slp.rand() % sites.size(),         //!< Pick a random site
                   n=sites.at(j).n;                 //!< Get particle number
      con->trial[n].charge += dp*markovmove::slp.random_half(); //!< Displace charge
      
      uold = pot->energy(con->p, n) + gaussianEnergy(con->p[n]);
      unew = pot->energy(con->trial, n) + gaussianEnergy(con->trial[n]);
      du = unew-uold;
      
      if (ens->metropolis(du)==true) {
        rc=OK;
        utot+=du;
        naccept++;
        dpsqr+=pow( con->p[n].charge - con->trial[n].charge,2);
        con->p[n].charge = con->trial[n].charge;
        sum+=du;
      } else {
        rc=ENERGY;
        con->trial[n].charge = con->p[n].charge;
        dpsqr+=0;
      }
    }
    return sum;
  }
  
  double chargeregGaussian::gaussianEnergyTotal(vector<particle> &p) {
    double u=0;
    for (int i=0; i<sites.size(); i++)
      u+=gaussianEnergy(p.at( sites[i].n ));
    return u;
  }
  
 //---------- CHARGE REG ---------------------
  string chargereg::info() {
    std::ostringstream o;
    o <<  markovmove::info()
      << "#   pH (concentration)        = " << ph << endl
      << "#   Titrateable sites         = " << sites.size() << endl
      << "#   Number of protons         = " << protons.size() << endl;
    return o.str();
  }

  chargereg::chargereg(ensemble &e, container &c, energybase &i, group &g, float ph ) : markovmove(e,c,i), titrate(c.p,g,ph)
  {
    runfraction=0.2;
    prefix.assign("tit_");
    name.assign("PROTON TITRATION");
    cite.assign("Biochem. 2005, 44, 5722-5727.");
    con->trial = con->p;
  }

  /*! \brief Exchange protons between bulk and titrateable sites.
   *
   *  This move will randomly go through the titrateable sites and
   *  try to exchange protons with the bulk. The trial energy is:
   */
  double chargereg::titrateall() {
    du=0;
    if (markovmove::slp.runtest(runfraction)==false)
      return du;
    action t;
    double sum=0;
    for (int i=0; i<sites.size(); i++) {
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

  //---------- GCTITRATE ---------
  /*!
   * \param ph pH value
   * \param mu Proton excess chemical potential
   */
  GCchargereg::GCchargereg(grandcanonical &e,
      container &c,
      energybase &i,
      inputfile &in) : markovmove(e,c,i), titrate_gc(c,in,e)
  {
    prefix.assign("tit_");
    this->name.assign("GRAND CANONICAL PROTON TITRATION");
    this->cite.assign("doi:10.1007/978-3-540-75755-9_8");
    markovmove::getInput(in);
  }

  string GCchargereg::info() {
    std::ostringstream o;
    o << markovmove::info()
      << titrate_gc::info(*con);
    return o.str();
  }

  double GCchargereg::titrateall() {
    double sum=0;
    if (markovmove::slp.random_one()>runfraction)
      return sum;
    unsigned int I, s;
    particle J=atom(nameA);
    for (int i=0; i<sites.size(); i++) {
      cnt++;
      int s=exchange(con->trial);
      if (recent==PROTONATED ) {     //Take ion from bulk        
        I=gcPtr->gp[gcPtr->findgroup(nameA)]->beg + markovmove::slp.random_one()*gcPtr->gp[gcPtr->findgroup(nameA)]->size();
        unew=pot->energy(con->trial, s)- pot->energy(con->trial[s], con->trial[I]);
        uold=pot->energy(con->p, s)+ pot->energy(con->p, I)
             -pot->energy(con->p[s], con->p[I]);
      }
      if (recent==DEPROTONATED) {                      //Release ion to bulk
        con->randompos(J);
        unew=pot->energy(con->trial, s)+ pot->energy(con->trial, J);
        uold=pot->energy(con->p, s);
      }
      du=unew-uold;

      if (recent==PROTONATED && gcPtr->gp[gcPtr->findgroup(nameA)]->size()==0) {
        rc=ENERGY;                                     //Cumbersome way to make sure we don't erase
        exchange(con->trial,s);                        //what isn't there.
      } else {
        if (ens->metropolis( titrate_gc::energy(*con, du, s)) == true ) {
          rc=OK;
          utot+=du;
          sum+=du;
          naccept++;
          if(recent==PROTONATED) {
            if(gcPtr->erase(*con,I)==false)
              std::cerr<<" Couldn't errase in GCtit!!!"<<endl;
            con->p[s].charge=con->trial[s].charge;
          } 
          if(recent==DEPROTONATED) {
            if(gcPtr->insert(*con,J)==false)
              std::cerr<<" Couldn't insert in GCtit!!!"<<endl;
            con->p[s].charge=con->trial[s].charge;
          }
        } else {
          rc=ENERGY;
          exchange(con->trial,s);
        }
      }
    }    
    for (int i=0; i<sites.size(); i++)
      zavg.at(i) += con->p.at( sites[i] ).charge;
    return sum;
  }


  // ----------- DH Titration ----------------

 DHchargereg::DHchargereg(ensemble &e, container &c, energybase &i, float ph, float muH ) : markovmove(e,c,i), titrate_implicit(c.p,ph)
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
    if (markovmove::slp.runtest(runfraction)==false)
      return du;
    double sum=0;
    for (int i=0; i<sites.size(); i++) {
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
      << "#   Titrateable sites   = " << sites.size() << endl;
    return o.str();
  }

#ifdef BABEL
  //GLU3titration constuctor
  glu3corechargereg::glu3corechargereg(ensemble &e, container &c, energybase &i, inputfile &in, group &g) 
  : chargereg(e,c,i,g,in.getflt("pH",7)) {
    name.assign("CHARGEREG / PORPHYRIN DOUBLE-TITRATION");
    cite.assign("NONE");
    runfraction=1;
    con->trial = con->p;
    porphyrinpKa = in.getflt("pKa_core", 5.0);
  }
  //GLU3 info
  string glu3corechargereg::info(){
    std::ostringstream o;
    o << chargereg::info()
      << "#   pKa of double core    = "<<porphyrinpKa<<endl
      << "#   Average charges on core"<<endl
      << "#   P1 = "<<p1.avg()<<" ("<<p1.stdev()<<")"<<endl
      << "#   P2 = "<<p2.avg()<<" ("<<p2.stdev()<<")"<<endl<<endl;
    return o.str();
  }
  //GLU3 corrected weight
  double glu3corechargereg::energy(vector<particle> &p,
      double du, action &a, action &b) {
//    int i=p[a.site].id;
    if (p[a.site].charge>0.001)
      return du+( log(10.)*( 2*ph - porphyrinpKa ) );
    else
      return du-( log(10.)*( 2*ph - porphyrinpKa ) );
  }
  //GLU3 Markov step
  double glu3corechargereg::move(glu3 &g) {
    du=0;
    if (markovmove::slp.runtest(runfraction)==false)
      return du;
    action t1, t2;
    double sum=0;
    cntcore++;
    if (con->trial[g.beg+4].charge>0.1) {
      t1=titrate::moveToBulk(con->trial,g.beg+4);
      t2=titrate::moveToBulk(con->trial,g.beg+13);
    } else {
      t1=titrate::takeFromBulk(con->trial,g.beg+4);
      t2=titrate::takeFromBulk(con->trial,g.beg+13);
    }  
    if(t1.proton==t2.proton)
    std::cout<<"site "<<t1.site<<" and site "<<t2.site<<" shares the same proton, neutron!"<<std::endl;
    //#pragma omp parallel
    {
      //#pragma omp sections
      {
        //#pragma omp section
        {
          uold = pot->potential( con->p, t1.site ) * con->p[t1.site].charge
            + pot->potential( con->p, t1.proton ) * con->p[t1.proton].charge
            - con->p[t1.site].charge * con->p[t2.site].charge * pot->tokT /
            sqrt( con->sqdist( con->p[t1.site], con->p[t2.site]) )
            +pot->potential( con->p, t2.site ) * con->p[t2.site].charge
            + pot->potential( con->p, t2.proton ) * con->p[t2.proton].charge
            - con->p[t1.proton].charge * con->p[t2.proton].charge * pot->tokT /
            sqrt( con->sqdist( con->p[t1.proton], con->p[t2.proton]) );
        }
        //#pragma omp section
        {
          unew = pot->potential( con->trial, t1.site ) * con->trial[t1.site].charge
            + pot->potential( con->trial, t1.proton ) * con->trial[t1.proton].charge
            - con->trial[t1.site].charge * con->trial[t2.site].charge * pot->tokT / 
            sqrt( con->sqdist( con->trial[t1.site], con->trial[t2.site]) )
            + pot->potential( con->trial, t2.site ) * con->trial[t2.site].charge
            + pot->potential( con->trial, t2.proton ) * con->trial[t2.proton].charge
            - con->trial[t1.proton].charge * con->trial[t2.proton].charge * pot->tokT / 
            sqrt( con->sqdist( con->trial[t1.proton], con->trial[t2.proton]) );
        }
      }
    }
    du = (unew-uold);
    if (ens->metropolis( energy(con->trial, du, t1,t2) )==true) {
      rc=OK;
      utot+=du;
      naccept++;
      con->p[t1.site].charge   = con->trial[t1.site].charge;
      con->p[t1.proton].charge = con->trial[t1.proton].charge;
      con->p[t2.site].charge   = con->trial[t2.site].charge;
      con->p[t2.proton].charge = con->trial[t2.proton].charge;
      sum+=du;
    } else {
      rc=ENERGY;
      exchange(con->trial, t1);
      exchange(con->trial, t2);
      du=0.;
    }
    p1+=con->p[g.beg+ 4].charge;
    p2+=con->p[g.beg+13].charge;
    
    return du;

  }
  //GLU3GCtitration constuctor
  GCglu3corechargereg::GCglu3corechargereg(grandcanonical &e, container &c, energybase &i, inputfile &in) 
  : GCchargereg(e,c,i,in) {
    name.assign("GCCHARGEREG / PORPHYRIN DOUBLE-TITRATION");
    cite.assign("NONE");
    runfraction=1;
    con->trial = con->p;
    porphyrinpKa = in.getflt("pKa_core", 5.0);
    i1=atom(nameA), i2=atom(nameA);
  }
  //GLU3 info
  string GCglu3corechargereg::info(){
    std::ostringstream o;
    o << GCchargereg::info()
      << "#   pKa of double core    = "<<porphyrinpKa<<endl
      << "#   Average charges on core"<<endl
      << "#   P1 = "<<p1.avg()<<" ("<<p1.stdev()<<")"<<endl
      << "#   P2 = "<<p2.avg()<<" ("<<p2.stdev()<<")"<<endl<<endl;
    return o.str();
  }
  //GLU3 corrected weight
  double GCglu3corechargereg::energy(vector<particle> &p,
      double du, int i) {
    double nA=double(gcPtr->gp[gcPtr->findgroup(nameA)]->size());
    if (con->trial[i].charge>0.001)  //Did we protonate?
      return du+( log(10.)*( 2*ph - porphyrinpKa ) ) - log(nA/con->getvolume()*(nA-1.0)/con->getvolume()) + 2*atom[nameA].chempot;
    else
      return du-( log(10.)*( 2*ph - porphyrinpKa ) ) + log((nA+2)/con->getvolume()*(nA+1.0)/con->getvolume()) - 2*atom[nameA].chempot;
  }
  //GLU3 Markov step
  double GCglu3corechargereg::move(glu3 &g) {
    du=0;
    if (markovmove::slp.runtest(runfraction)==false)
      return du;
    double sum=0;
    cntcore++;
    // New configuration
    if (con->trial[g.beg+4].charge>0.1) {     //Is the core protonated?
      con->randompos(i1), con->randompos(i2); //DEPROTONATE
      con->trial[g.beg+4].charge=0.0, con->trial[g.beg+13].charge=0.0;
      unew = pot->energy(con->trial, g.beg+4) + pot->energy(con->trial, g.beg+13) -
             pot->energy(con->trial[g.beg+4], con->trial[g.beg+13]) +
             pot->energy(con->trial, i1) + pot->energy(con->trial, i2)+ pot->energy(i1,i2);
      uold = pot->energy(con->p, g.beg+4) + pot->energy(con->p, g.beg+13) -
             pot->energy(con->p[g.beg+4], con->p[g.beg+13]);
    } else {
      o1=gcPtr->gp[gcPtr->findgroup(nameA)]->beg+markovmove::slp.random_one()*gcPtr->gp[gcPtr->findgroup(nameA)]->size();
      o2=o1;
      while(o2==o1) //PROTONATE
        o2=gcPtr->gp[gcPtr->findgroup(nameA)]->beg+markovmove::slp.random_one()*gcPtr->gp[gcPtr->findgroup(nameA)]->size();
      con->trial[g.beg+4].charge=1.0, con->trial[g.beg+13].charge=1.0;
      unew = pot->energy(con->trial, g.beg+4) + pot->energy(con->trial, g.beg+13)-
             pot->energy(con->trial[g.beg+4], con->trial[g.beg+13]) - pot->energy(con->trial[o1], con->trial[g.beg+4])-
             pot->energy(con->trial[o1], con->trial[g.beg+13]) - pot->energy(con->trial[o2], con->trial[g.beg+4])-
             pot->energy(con->trial[o2], con->trial[g.beg+13]);

      uold = pot->energy(con->p,g.beg+4) + pot->energy(con->p,g.beg+13) + pot->energy(con->p,o1)+
             pot->energy(con->p, o2) - pot->energy(con->p[g.beg+4], con->p[g.beg+13])- pot->energy(con->p[o1], con->p[o2])-
             pot->energy(con->p[g.beg+4], con->p[o1]) - pot->energy(con->p[g.beg+13], con->p[o1]) - 
             pot->energy(con->p[g.beg+4], con->p[o2]) - pot->energy(con->p[g.beg+13], con->p[o2]);
    }  
    du = (unew-uold);
    // Check so that there are ions enough!
    if (con->trial[g.beg+4].charge>0.1 && gcPtr->gp[gcPtr->findgroup(nameA)]->size()<2) {//Did we protonate with insufficent number of ions?
      rc=ENERGY;
      con->trial[g.beg+4].charge  = con->p[g.beg+4].charge;
      con->trial[g.beg+13].charge = con->p[g.beg+13].charge;
      du=0.;
      return du;
    }
    // Weight of new conf.
    if (ens->metropolis( energy(con->trial, du, g.beg+4) )==true) {
      rc=OK;
      utot+=du;
      naccept++;
      con->p[g.beg+4].charge  = con->trial[g.beg+4].charge;
      con->p[g.beg+13].charge = con->trial[g.beg+13].charge;
      if (con->trial[g.beg+4].charge>0.1) { // Did we protonate?
        if(o1>o2) {                         // If we begin with the large index we will not
        gcPtr->erase(*con ,o1);             // have to adjust for the last erase
        gcPtr->erase(*con ,o2);
        } else {
        gcPtr->erase(*con ,o2);
        gcPtr->erase(*con ,o1);
        }
      } else { 
        gcPtr->insert(*con, i1);
        gcPtr->insert(*con, i2);
      }
    } else {
      rc=ENERGY;
      con->trial[g.beg+4].charge  = con->p[g.beg+4].charge;
      con->trial[g.beg+13].charge = con->p[g.beg+13].charge;
      du=0.;
    }
    p1+=con->p[g.beg+ 4].charge;
    p2+=con->p[g.beg+13].charge;
    
    return du;

  }
#endif

  // ----------- AT Titration ----------------
  ATchargereg::ATchargereg( ensemble& e, container& c, energybase& i, float ph, inputfile& in , pot_debyehuckel& pair)
  :markovmove(e,c,i), titrate_implicit(c,ph) {

    name.assign("IMPLICIT SALT TITRATION");
    cite.assign("doi:10.1021/ct1003093");
    double sytem_charge = con->charge();
    pairpot = &pair;

    // Getting config from input
    protein_conc = in.getflt("ProteinConc", 0.1);
    double salt1_conc   = in.getflt("Salt1Conc"  , 0.1);
    double salt2_conc   = in.getflt("Salt2Conc"  , 0.1);
    double salt1_charge = in.getflt("Salt1Charge", 1.0);
    double salt2_charge = in.getflt("Salt2Charge", 1.0);
    float lB = in.getflt("bjerrum", 7.12);

    // Evaluating physical quantities
    ionic_str0 = .5 * ( salt1_conc * pow(salt1_charge,2) + salt2_conc * pow(salt2_charge,2) );
    ionic_str1 = .5 * ( protein_conc * abs(ionic_str0) );
    const_kappa = sqrt(8*pyc.pi*lB*pyc.Nav/1e27);
    kold = (const_kappa * sqrt(ionic_str0+ionic_str1));
    set_kappa(kold);
  }

  double ATchargereg::titrateall( vector<macromolecule>& g ) {
    du=0;
    if (markovmove::slp.runtest(runfraction)==false)
      return du;
    double sum=0;
    int t;
    int who;
    for (int i=0; i<sites.size(); i++) {
      cnt++;
      t    = exchange(con->trial);         // Protonate/Deprotonate site
      who  = who_is_tit(g,t);              // Find which protein is titrating (in case we have many proteins)
      kold = calc_kappa(g,con->p);         // Update old kappa value
      knew = calc_kappa(g,con->trial);     // Update new kappa value

      // Evaluate energy variation of the move...
      uold = pot->potential(con->p,t,g[who],kold) * con->p[t].charge
           + prot_ion_u(g,con->p,kold);
      unew = pot->potential(con->trial,t,g[who],knew) * con->trial[t].charge
           + prot_ion_u(g,con->trial,knew);
      du = (unew-uold);

      // Accept or not the move...
      if (ens->metropolis( energy(con->trial, du, t) )==true) {
        rc=OK;
        utot+=du;
        naccept++;
        con->p[t].charge = con->trial[t].charge;
        sum+=du;
        set_kappa(knew);
      }
      else {
        rc=ENERGY;
        exchange(con->trial, t);
        set_kappa(kold);
      }
      //samplesites(con->p);  // Average charges on all sites
    }
    return sum;
  }

  /*!
   * \brief Returns information about the titration.
   * \author Andre Teixeira
   * \date Jul2009
   */
  string ATchargereg::info() {
    std::ostringstream o;
    o <<  markovmove::info()
      << "#   pH (activity scale)       = " << ph             << endl
      << "#   Titrateable sites         = " << sites.size()   << endl
      << "#   Protein concentration (M) = " << protein_conc   << endl
      << "#   Ionic strenght (M)        = " << ionic_str0     << endl
      << "#   Debye lenght (AA)         = " << 1./ (const_kappa * sqrt(ionic_str0+ionic_str1)) << endl;
    return o.str();
  }

  /*!
   * \brief Calculates the value of kappa for a given vector of particles.
   * \author Andre Teixeira
   * \date Jul 2009
   * \note Shouldn't the charge be squared?!
   */
  double ATchargereg::calc_kappa( vector<macromolecule>& g, vector<particle>& p) {
    ionic_str1 = 0;
    for (int i=0 ; i<g.size() ; i++)
      ionic_str1 += g[i].conc * abs(g[i].getcharge(p));
    ionic_str1 *= 0.5;
    return const_kappa * sqrt(ionic_str0+ionic_str1);
  }

  double ATchargereg::totalenergy(vector<macromolecule> &g, vector<particle> &p ) {
    double u=0,
           k=calc_kappa(g,p);
    for  (int ig=0 ; ig<g.size() ; ig++)
      for (int ip=g[ig].beg; ip<=g[ig].end; ip++)
        u+=p[ip].charge * pot->potential( p, ip, g[ig], k );
    return u/2 + prot_ion_u(g,p,k);
  }
  
  /*!
   * \brief Evaluates the protein-counter ion energy for all proteins for a calculated value of kappa
   * and a given vector of particles.
   * \author Andre Teixeira
   * \date Jul 2009
   */
  double ATchargereg::prot_ion_u( vector<macromolecule>& g, vector<particle>& p) {
    double k = calc_kappa(g,p);
    return prot_ion_u(g,p,k);
  }

  /*!
   * \brief Evaluate the protein-counter ion energy for all proteins for a given value of kappa
   * and a given vector of particles.
   * \author Andre Teixeira
   * \date Jul 2009
   */
  double ATchargereg::prot_ion_u( vector<macromolecule>& g, vector<particle>& p, double& k) {
    double u=0, zp;
    for (int i=0 ; i<g.size() ; i++) {
      zp = g[i].getcharge(p);
      u += zp*zp / ( 1. + k*g[i].cm.radius );
    }
    return -0.5 * pairpot->f * k * u;
  }

  /*!
   * \brief Find the macromolecule which the current titrating site belongs.
   * \author Andre Teixeira
   * \date Jul 2009
   */
  int ATchargereg::who_is_tit( vector<macromolecule>& g, int& t) {
    for (int i=0 ; i<g.size() ; i++)
      if ( t >= g[i].beg && t <= g[i].end)
        return i;
    std::cout << "\n\n#Site not found!!\n\n" << endl;
    return -1;
  }

  /*!
   * \brief Update kappa value used in intermolecular interactions calculations.
   * \authot Andre Teixeira
   * \date Jul 2009
   */
  void ATchargereg::set_kappa( double& k ) {
    pairpot->k = k;
    kold = k;
  }

}//namespace
