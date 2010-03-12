#ifndef FAU_ISOBARIC_H
#define FAU_ISOBARIC_H

#include "faunus/moves/base.h"

namespace Faunus {
  /*!
   * \brief Fluctuate the volume against an external pressure 
   * \author Bjoern Persson and Mikael Lund
   * \todo Should take an arbritray vector of groups and implement a more efficent energycalculation.
   * Salt scaling not yet implemented
   * \note Tested and found valid using debyehuckel_P3 and point charges and ideal systems. 
   * \note To compare with ideal systems sample the INVERSE of the volume to avoid (N+1)/N error! Intrinis property of the ensemble.
   * \date Lund/Canberra, Jan-Feb 2008
   *
   * This move is intended to maintain constant pressure in a canonical
   * Monte Carlo simulation. It will perform volume fluctuations and scale
   * particle positions according to the scaling algorithm specified for the container.
   */
  template<typename T> class isobaric : public markovmove {
    public:
      isobaric( ensemble&, container&, T &, double, double, double, int, int, int);
      string info();
      void move(group &, unsigned short);     //!< Scale group mass center and move
      void move(unsigned short, group &,...);
      double move(vector<macromolecule> &);   // Markov move 
      average<double> vol, len, ilen, ivol;    // Standard analysis
      double pen;                             // Penalty to update penalty function
      double initpen;
      double scalepen;
      void loadpenaltyfunction(string);       // Load old penalty function 
      void printupdatedpenalty(string);       // Print penalty, old+new
      void printpenalty(string);              // Old penalty
      void updatepenalty();                   // Add penalty and newpenalty
      void check(checkValue &);               //!< Unit test routine 
    private:
      unsigned short N;
      T trialpot;                             // Copy of potential class for volume changes
      double P, dV, dh;                       // Pressure, volume difference and hamiltonian difference
      int maxlen, minlen, scale;              // Window parameters
      double newV;                            // New volume
      void newvolume();                       // Calculate new volume
      void accept();                          // Update volume and pair potential
      bool penalize;                          // Use penalty function
      vector<double> penalty;                 // Initial penalty function
      vector<double> newpenalty;              // New penalty function
      T *pot;                                 // Override markovmove pointer.
      double pendiff();                       // Retruns the penalty potential diff of old and new volume
      int penbin();                           // Retruns the bin pointer of the present volume
  };

  //---------- ISOBARIC ----------
  /*!
   * \note If your pair-potential depends on the geometry of the 
   *       container (minimum image, for example) make sure that
   *       the pair-potential has a correct implementation of the
   *       setvolume() function. To allow for parallization this 
   *       class will instantiate a new (private) interaction class
   *       and in this way allow for two different container sizes.
   *
   * \param pressure bulk pressure in A^-3
   * \param i Potential class. Make sure that the T_pairpot::setvolume()
   *          function is appropriate.
   */
  template<typename T> isobaric<T>::isobaric
    ( ensemble &e, container &c, T &i,
      double pressure, double PEN, double scpen, int maxsize, int minsize, int sc ) : markovmove(e,c,i), trialpot(i) 
  {
    name="ISOBARIC";
    cite="none so far";
    P=pressure;
    minlen=minsize;
    maxlen=maxsize;
    scale=sc;
    runfraction=1;
    dp=100; 
    N=0;
    // Add new fuction to set penalty fuctions (minlen-maxlen)/scale
    int penbins= (maxlen-minlen)/scale;
    if ((maxlen-minlen)%scale!=0) penbins++; //number of bins, add one if range%scale!=0 
    penalty.clear();
    penalty.resize(penbins, 1.0e-20);
    pen=PEN;
    initpen=PEN;
    scalepen=scpen;
    newpenalty.clear();
    newpenalty.resize(penbins, 1.0e-20);
    //
    penalize=false;
    pot=&i;
  }

  template<typename T> void isobaric<T>::printupdatedpenalty(string file) {
    std::ofstream f(file.c_str());
    if (f) {
      int j=penalty.size();
      f.precision(12);
      for (int i=0; i<j; i++)
        f << (i+minlen)*scale <<" "<< penalty[i]+newpenalty[i]<<endl;
      f.close();
    }
  }

  template<typename T> void isobaric<T>::printpenalty(string file) {
    std::ofstream f(file.c_str());
    if (f) {
      int j=penalty.size();
      f.precision(12);
      for (int i=0; i<j; i++)
        f << (i+minlen)*scale <<" "<< penalty[i]<<endl;
      f.close();
    }
  }

  template<typename T> void isobaric<T>::loadpenaltyfunction(string file) {
    //Set up penalty functions based on min,max and scale and enables penalty function updating
    penalize=true;    //Boolean to determine update of newpenalty
    std::cout <<"# Let's go BIASED!!!"<<endl;
    int penbins= (maxlen-minlen)/scale;
    if ((maxlen-minlen)%scale!=0) penbins++; //number of bins, add one if range%scale!=0 
    penalty.clear();
    penalty.resize(penbins,1.0e-20);
    newpenalty.clear();
    newpenalty.resize(penbins,1.0e-20);
    vector<string> v;
    string s;
    v.clear();
    std::ifstream f(file.c_str() );
    if (f) {
      while (getline(f,s))
        v.push_back(s);
      f.close();
      std::cout << "# Penalty function loaded!! " << endl;
      short c=v.size();
      if (c==penbins) {
        double val;
        double bin;
        for (short i=0;i<c;i++) {
          val=0;
          bin=0;
          std::stringstream o;
          o << v[i];
          o >> bin >> val;
          penalty[(bin-minlen)/scale]=val;
        }
      }else {
        std::cout    <<"# Loaded penalty function does not fit range" <<endl
          <<"# Wrong min/maxlen or scale?"<<endl;  
        std::ofstream fback("penaltybackup.dat");
        for (int i=0; i<v.size(); i++) {
          fback <<v[i]<<endl;
        }
        fback.close();
        std::cout    <<"# Backup of input penalty as penaltybackup.dat"<<endl;
      }
    }
    else std::cout << "# WARNING! Penalty function " << file << " NOT READ!\n";
  }

  template<typename T> string isobaric<T>::info() {
    std::ostringstream o;
    o << markovmove::info();
    o <<   "#         Min/Max box length  = "<<minlen<<" / "<<maxlen<<endl
      <<   "#         External pressure   = "<< P <<"(A^-3) = "<<P*1660<<" (M)"<< endl
      <<   "#         Average volume      = "<< vol.avg() << " ("<<vol.stdev()<<") A^3" << endl
      <<   "#         Average box length  = "<< len.avg() << " ("<<len.stdev()<<") A^3" << endl;
    if(penalize==true) {
      o << "#         Penalty             = "<< pen<<" kT per observation"<<endl
        << "#         Initial penalty     = "<< initpen<<" kT per observation"<<endl;
    }
    return o.str();
  }

  template<typename T> void isobaric<T>::check(checkValue &test) {
    markovmove::check(test);
    test.check("ISOBARIC_AvgBoxLen", len.avg() );
  }

  template<typename T> void isobaric<T>::newvolume() {
    newV = exp(log(con->getvolume())  // Random move in ln(V)
        + slp.random_half()*dp);
    trialpot.pair.setvolume(newV);    // Set volume for trial-pair-potential
    dV=newV - con->getvolume();
  }
  template<typename T> void isobaric<T>::accept() {
    naccept++;
    rc=OK;
    utot+=du;
    dpsqr+=pow(newV,1./3.);
    con->setvolume(newV);         // save the trial volume in the container
    pot->pair.setvolume(newV);    // ...AND in the pair-potential function
  }
  /*!\param n Number of groups
   * \param g List of groups
   * \note Unfinished
   */
  template<typename T> void isobaric<T>::move(unsigned short n, group &g,...) {
    /*
       group first=g; // check if first group is included!!
       double newlen=pow(newV,1./3);
       group *gPtr;
       va_list ap;
       va_start(ap, g);
       while (n>0) { 
    //gPtr=&va_arg(ap, group);
    //gPtr->isobaricmove(*con, newlen);
    //cout << gPtr->name << endl;
    //N+=gPtr->nummolecules();
    n--;
    }
    */
  }
  template<typename T> double isobaric<T>::pendiff() {
    double penval;
    penval=penalty[int( pow(newV, 1./3.) - minlen)/scale]
      -penalty[int( pow(con->getvolume(),1./3.) - minlen)/scale];
    return penval;
  }
  template<typename T> int isobaric<T>::penbin() {
    int binnr;
    binnr=int(pow(con->getvolume(),1./3.)-minlen)/scale;
    return binnr;
  }
  template<typename T> void isobaric<T>::updatepenalty() {
    for (int i=0; i<penalty.size(); i++) {
      penalty[i]+=newpenalty[i];
      newpenalty[i]=0;
    }
    pen*=scalepen;
  }
  template<typename T> double isobaric<T>::move(vector<macromolecule> &g) {
    du=0;
    if (slp.runtest(runfraction)==false)
      return du;
    cnt++;
    double uold=0,unew=0;
    newvolume();                                // Generate new trial volume - prep. trial potential
    if (newV>pow(double(minlen), 3.) && newV<pow(double(maxlen),3.)) {
      int i=g.size();
      //#pragma omp parallel for
      for (int n=0; n<i; n++)                      // Loop over macromolecules
        g[n].isobaricmove(*con, pow(newV,1./3.));  // ..and scale their mass-centers
      int n=g.size();
#pragma omp parallel for reduction (+:uold,unew) schedule (dynamic)
      for (int j=0; j<n-1; j++)
        for (int k=j+1; k<n; k++) {
          uold += pot->energy(con->p, g[j], g[k]);         // calc. old energy with original potential class
          unew += trialpot.energy(con->trial, g[j], g[k]); // calc. new energy using a COPY of the potential class
        }
      du = unew-uold;
      dh = du + P*dV-(i+1)*log( newV/con->getvolume() );
      if (penalize==true) 
        dh += pendiff();
      if (ens->metropolis(dh)==true ) {
        accept();
        for (unsigned short n=0; n<i; n++) 
          g[n].accept(*con);
        vol += con->getvolume();
        ivol+= 1./con->getvolume();
        if (pen!=0.)
          newpenalty[penbin()]+=pen;
        len += pow(con->getvolume(),1./3.);
        ilen+= 1./pow(con->getvolume(),1./3.);
        return du;
      } else rc=ENERGY;
      du=0;
      dpsqr+=0;
      for (unsigned short n=0; n<i; n++)
        g[n].undo(*con);
      vol += con->getvolume();
      ivol+= 1./con->getvolume();
      if (pen!=0.)
        newpenalty[penbin()]+=pen;
      len += pow(con->getvolume(),1./3.);
      ilen+= 1./pow(con->getvolume(),1./3.);
      return du;
    }
    dpsqr+=0;
    vol += con->getvolume();
    ivol+= 1./con->getvolume();
    len += pow(con->getvolume(),1./3.);
    ilen+= 1./pow(con->getvolume(),1./3.);
    if (pen!=0.)
      newpenalty[penbin()]+=pen;
    return du;
  } 
}
#endif
