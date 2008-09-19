#ifndef FAU_ISOBARIC_H
#define FAU_ISOBARIC_H

#include "faunus/moves/base.h"

namespace Faunus {
  /*!
   * \brief Fluctuate the volume against an external pressure 
   * \author Bjoern Persson and Mikael Lund
   * \todo Should take an arbritray vector of groups and implement a more efficent energycalculation.
   * Salt scaling not yet implemented
   * \note Tested and found valid using debyehuckel_P3 and point charges.
   * \date Lund/Canberra, Jan-Feb 2008
   *
   * This move is intended to maintain constant pressure in a canonical
   * Monte Carlo simulation. It will perform volume fluctuations and scale
   * particle positions according to the scaling algorithm specified for the container.
   */
  template<typename T> class isobaric : public markovmove {
    public:
      isobaric( ensemble&, container&, interaction<T>&, double, double, int);
      string info();
      void move(group &, unsigned short);  //!< Pressure scale group
      void move(unsigned short, group &,...);
      double move(vector<macromolecule> &);
      average<float> vol, vol2, len, len2;
      double pen;
      vector<double> penalty;
      void loadpenaltyfunction(string);
      void printpenalty(string, histogram &);
      void printpenalty(string);
    private:
      unsigned short N;
      interaction<T> trialpot; // Copy of potential class for volume changes
      double P, dV, dh; // Pressure, volume difference and hamiltonian difference
      double newV;      // New volume
      void newvolume();
      void accept();
      vector<double> newpenalty;
      interaction<T> *pot; // Override markovmove pointer.
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
   *        function is appropriate.
   */
  template<typename T> isobaric<T>::isobaric (
      ensemble &e,
      container &c,
      interaction<T> &i, double pressure, double PEN,int maxsize ) : markovmove(e,c,i), trialpot(i) {
    name="ISOBARIC";
    cite="none so far";
    P=pressure;
    runfraction=0.10;
    dp=100; 
    N=0;
    penalty.clear();
    penalty.resize(maxsize,0);
    pen=PEN;
    newpenalty.clear();
    newpenalty.resize(maxsize,0);
    pot=&i;
  }
  template<typename T> void isobaric<T>::printpenalty(string file, histogram &ld) {
    std::ofstream f(file.c_str());
    if (f) {
      int j=penalty.size();
      f.precision(12);
      for (int i=0; i<j; i++)
        f << i <<" "<< penalty[i]+newpenalty[i]/cnt+log(ld.get(i))<<endl;
      f.close();
    }
  }
  template<typename T> void isobaric<T>::printpenalty(string file) {
    std::ofstream f(file.c_str());
    if (f) {
      int j=penalty.size();
      f.precision(12);
      for (int i=0; i<j; i++)
        f << i <<" "<< penalty[i]<<endl;
      f.close();
    }
  }
  template<typename T> void isobaric<T>::loadpenaltyfunction(string file) {
    vector<string> v;
    string s;
    v.clear();
    std::ifstream f(file.c_str() );
    if (f) {
      while (getline(f,s))
        v.push_back(s);
      f.close();
      std::cout << "# Penalty function loaded!! Let's go biased" << endl;
    }
    else std::cout << "# WARNING! Penalty function " << file << " NOT READ!\n";
    short c=v.size();
    double val;
    int bin;
    for (short i=0;i<c;i++) {
      val=0;
      bin=0;
      std::stringstream o;
      o << v[i];
      o >> bin >> val;
      penalty[bin]=val;
    }
  }
  template<typename T> string isobaric<T>::info() {
    std::ostringstream o;
    o << markovmove::info();
    o << "#   External pressure   = "<< P <<"(A^-3) = "<<P*1660<<" (M)"<< endl
      << "#   Average volume      = "<< vol.avg() << " A^3" << endl
      << "#     sqrt(<V^2>-<V>^2) = "<< sqrt(vol2.avg()-pow(vol.avg(),2)) << endl
      << "#   Average box length  = "<< len.avg() << " A^3" << endl
      << "#     sqrt(<l^2>-<l>^2) = "<< sqrt(len2.avg()-pow(len.avg(),2)) << endl;
    return o.str();
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
    vol += con->getvolume();
    vol2+= pow(con->getvolume(),2);
    len += pow(con->getvolume(),1./3);
    len2+= pow(con->getvolume(),2./3);
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
  template<typename T> double isobaric<T>::move(vector<macromolecule> &g) {
    cnt++;
    du=0;
    newvolume();                                // Generate new trial volume - prep. trial potential
    unsigned short i=g.size();
    for (unsigned short n=0; n<i; n++)          // Loop over macromolecules
      g[n].isobaricmove(*con, pow(newV,1./3));  // ..and scale their mass-centers
    uold = pot->energy(con->p);                 // calc. old energy with original potential class
    unew = trialpot.energy(con->trial);         // calc. new energy using a COPY of the potential class
    du = unew-uold;
    dh = du + P*dV-(i+1)*log( newV/con->getvolume() );
    if (pen!=0) 
      dh += penalty[int(pow(newV,1./3))]-penalty[int(pow(con->getvolume(),1./3))];
    if (ens->metropolis(dh)==true) {
      accept();
      for (unsigned short n=0; n<i; n++) 
        g[n].accept(*con);
      if (pen!=0)
        newpenalty[int(pow(con->getvolume(),1./3))]+=pen;
      return du;
    } else rc=ENERGY;
    du=0;
    dpsqr+=0;
    for (unsigned short n=0; n<i; n++)
      g[n].undo(*con);
    vol += con->getvolume();
    vol2+= pow(con->getvolume(),2);
    if (pen!=0)
      newpenalty[int(pow(con->getvolume(),1./3))]+=pen;
    len += pow(con->getvolume(),1./3);
    len2+= pow(con->getvolume(),2./3);
    return du;
  } 
}
#endif
