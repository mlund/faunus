#include <faunus/moves/charge.h>
#include <faunus/titrate.h>

namespace Faunus {
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
    if (slp.runtest(runfraction)==false)
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
#endif


    // ----------- AT Titration ----------------
  ATchargereg::ATchargereg( ensemble& e, container& c, energybase& i, float ph, inputfile& in , pot_debyehuckel& pair)
  :markovmove(e,c,i), titrate_implicit(c,ph) {

	name.assign("IMPLICIT IONS TITRATION");
	double sytem_charge = con->charge();
	pairpot = &pair;

	// Getting config from input
    protein_conc = in.getflt("ProteinConc", 0.1);
    double salt1_conc   = in.getflt("Salt1Conc"  , 0.1);
    double salt2_conc   = in.getflt("Salt2Conc"  , 0.1);
    double salt1_charge = in.getflt("Salt1Charge", 1.0);
    double salt2_charge = in.getflt("Salt2Charge", 1.0);
    float lB = in.getflt("bjerrum"  , 7.12);

    // Evaluating physical quantities
    ionic_str0 = 0.5 * ( salt1_conc * pow(salt1_charge,2) + salt2_conc * pow(salt2_charge,2) );
    ionic_str1 = 0.5 * ( protein_conc * abs(ionic_str0) );
    const_kappa = sqrt( 80 * 3.14159265358979 * 6.02214179e23 * lB )*1e-14;
    kold = (const_kappa * sqrt(ionic_str0+ionic_str1));
    set_kappa(kold);
  };




  double ATchargereg::titrateall( vector<macromolecule>& g ) {
    du=0;
    if (slp.runtest(runfraction)==false)
      return du;
    double sum=0;
    int t;
    short who;
    for (unsigned short i=0; i<sites.size(); i++) {
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
      };
      //samplesites(con->p);  // Average charges on all sites
    };
    return sum;
  };
    
    

  /*!
   * \brief Returns information about the titration.
   * \author Andre Teixeira
   * \date Jul2009
   */
  string ATchargereg::info() {
    std::ostringstream o;
    o <<  markovmove::info()
      << "#   pH (concentration)    = " << ph             << endl
      << "#   Titrateable sites     = " << sites.size()   << endl
      << "#   Protein concentration = " << protein_conc   << endl
      << "#   Salt Ionic Strenght   = " << ionic_str0     << endl
      << "#   Debye Lenght          = " << 1./ (const_kappa * sqrt(ionic_str0+ionic_str1)) << endl;
    return o.str();
  };

  /*!
   * \brief Calculates the value of kappa for a given vector of particles.
   * \author Andre Teixeira
   * \date Jul 2009
   */
  double ATchargereg::calc_kappa( vector<macromolecule>& g, vector<particle>& p) {
    ionic_str1 = 0;
    for (int i=0 ; i<g.size() ; i++)
      ionic_str1 += g[i].conc * abs(g[i].getcharge(p));
    ionic_str1 *= 0.5;
    return const_kappa * sqrt(ionic_str0+ionic_str1);
  };

  /*!
   * \brief Evaluates the protein-counter ion energy for all proteins for a calculated value of kappa
   * and a given vector of particles.
   * \author Andre Teixeira
   * \date Jul 2009
   */
  double ATchargereg::prot_ion_u( vector<macromolecule>& g, vector<particle>& p) {
    double u = 0;
    double k = calc_kappa(g,p);
    double zp;
    for (int i=0 ; i<g.size() ; i++) {
      zp = g[i].getcharge(p);
      u += -0.5 * pairpot->f * zp*zp * k / ( 1. + k*g[i].cm.radius );
    }
    return u;
  };
  
  
  /*!
   * \brief Evaluate the protein-counter ion energy for all proteins for a given value of kappa
   * and a given vector of particles.
   * \author Andre Teixeira
   * \date Jul 2009
   */
  double ATchargereg::prot_ion_u( vector<macromolecule>& g, vector<particle>& p , double& k) {
    double u = 0;
    double zp;
    for (int i=0 ; i<g.size() ; i++) {
      zp = g[i].getcharge(p);
      u += -0.5 * pairpot->f * zp*zp * k / ( 1. + k*g[i].cm.radius );
    }
    return u;
  };
 

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
  };



  /*!
   * \brief Update kappa value used in intermolecular interactions calculations.
   * \authot Andre Teixeira
   * \date Jul 2009
   */
  void ATchargereg::set_kappa( double& k ) {
    pairpot->k = k;
    kold = k;
  };

}//namespace
