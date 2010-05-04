#ifndef FAUNUS_COARSEGRAIN_H
#define FAUNUS_COARSEGRAIN_H

#include "faunus/energy/base.h"

namespace Faunus {
  
  /*!
   * \brief Treats far-away groups as monopoles for faster energy evaluation
   * \author Mikael Lund
   * \date Lund 2009
   * \warning Be careful when you have non-molecular groups such as salt.
   *
   * This interaction class redefines group-group and group-particle interactions
   * so that groups (i.e. molecules) far away will be seen as monopoles centered
   * at their charge-centers. I.e. the interaction between two charged proteins beyond
   * some threshold distance will interact as two point charges. Likewise a particle
   * far way from a group will see it only as a single particle
   */
  template<class T>
  class interaction_monopole : public interaction<T> {
  private:
    particle monopole(vector<particle> &p, const group &g) {
      double sum=0;
      particle cm;
      point t, o = p[g.beg]; // set origo to first particle
      for (unsigned short i=g.beg; i<=g.end; i++) {
        t = p[i]-o;        // translate to origo
        interaction<T>::pair.boundary(t);       // periodic boundary (if any)
        cm += t * p[i].mw;
        sum += p[i].mw; 
        cm.charge+=p[i].charge;
      }
      cm=cm*(1./sum) + o;
      interaction<T>::pair.boundary(cm);
      return cm;
    }
  public:
    double cut_g2g; //!< Cut-off distance for group-group interactions
    double cut_g2p; //!< Cut-off distance for group-particle interactions
    interaction_monopole(inputfile &in) : interaction<T>(in) {
      interaction<T>::name+=" w. monopole cut-offs";
      cut_g2g = in.getflt( "threshold_g2g", 1e6 );
      cut_g2p = in.getflt( "threshold_g2p", 1e6 );
    }
    double energy(vector<particle> &p ) { return interaction<T>::energy(p); } 
    double energy(vector<particle> &p, const group &g1, const group &g2) {
      particle mp1=monopole(p,g1), mp2=monopole(p,g2);
      return ( sqrt(interaction<T>::pair.sqdist(mp1, mp2)) > cut_g2g ) ?
      interaction<T>::energy( mp1, mp2 ) : interaction<T>::energy(p, g1, g2);
    }
    double energy(vector<particle> &p, const group &g) {
      double u=0;
      particle mp=monopole(p,g);
      for (int i=0; i<g.beg; i++)
        u+= (sqrt(interaction<T>::pair.sqdist(mp,p[i]))>cut_g2p) ? interaction<T>::energy(mp,p[i]) : interaction<T>::energy(p,g,i);
      for (int i=g.end+1; i<p.size(); i++)
        u+= (sqrt(interaction<T>::pair.sqdist(mp,p[i]))>cut_g2p) ? interaction<T>::energy(mp,p[i]) : interaction<T>::energy(p,g,i);
      return u;
    }
    string info() {
      std::ostringstream o;
      o << interaction<T>::info()
      << "#   Group-Group threshold    = " << cut_g2g << endl
      << "#   Group-Particle threshold = " << cut_g2p << endl;
      return o.str();
    }
  };
  
  /*!
   * \brief Treats far-away groups as (charged) dipoles for faster energy evaluation
   * \author Mikael Lund
   * \date Lund 2010
   * \warning Be careful when you have non-molecular groups such as salt.
   *
   * This interaction class redefines group-group interactions so that groups
   * (i.e. molecules) far away will be seen as dipoles with a cation and an anion centered
   * at their charge-centers. I.e. the interaction between two charged proteins beyond
   * some threshold distance will interact as dipoles. To avoid energy drifts it is
   * important that all energy evaluations are performed on a group-group basis.
   *
   * For proteins this procedure potentially moves charges from the surface to
   * two charges further in. When using Debye-Huckel potentials this leads to
   * an artifact that must be corrected for by the prefactor:
   *
   * \f$ \sinh(\kappa R_1)/\kappa R_1 \cdot \sinh(\kappa R_2)/\kappa R_2 \f$
   *
   * where 1/kappa is the screening length and R1 and R2 are the
   * radii of the two interacting spheres.
   * Here we have two particles per group and the sizes for each charge
   * center is calculated according to
   *
   * \f$ R_{\pm} = R_p - |\mathbf{R_{cm}R_{\pm}}| \f$
   *
   * where Rp is the protein radius, Rcm is the center of mass vector and
   * R+= is the center of of charge vector.
   *
   */
  
  template<class T>
  class interaction_dipole : public interaction<T> {
  private:
    struct dipole {
      point cm;
      particle anion, cation;
    };
    double kappa;
    unsigned long int cnt, cntDip; // no. of group-group interactions, no. of dipole interactions
    
    void centerOfMass(const vector<particle> &p, const group &g, point &cm) {
      point t;
      double sum=0;
      cm.clear();
      for (int i=g.beg; i<=g.end; i++) {
        t = p[i]-p[g.beg];
        interaction<T>::pair.boundary(t);
        cm += t*p[i].mw;
        sum+= p[i].mw; 
      }
      cm=cm*(1/sum) + p[g.beg];
      interaction<T>::pair.boundary(cm);
    }
    
    void calcDipole(const vector<particle> &p, const group &g, dipole &a) {
      point t;
      a.anion.clear();
      a.cation.clear();
      for (int i=g.beg; i<=g.end; i++) {
        t = p[i]-a.cm;
        interaction<T>::pair.boundary(t);
        if (p[i].charge>0) {
          a.cation.charge+=p[i].charge;
          a.cation += t*p[i].charge;
        }
        if (p[i].charge<0) {
          a.anion.charge+=p[i].charge;
          a.anion += t*p[i].charge;
        }
      }
      if (a.anion.charge!=0)  a.anion = a.anion*(1/a.anion.charge);
      if (a.cation.charge!=0) a.cation= a.cation*(1/a.cation.charge);
      a.anion+=a.cm;
      a.cation+=a.cm;
      interaction<T>::pair.boundary(a.anion);
      interaction<T>::pair.boundary(a.cation);
    }
    
    //!< Calculate Debye-Huckel prefactor
    double f(particle &p1, particle &p2) {
      return std::sinh(p1.radius*kappa) * std::sinh(p2.radius*kappa)
      / (kappa*kappa*p1.radius*p2.radius) ;
    }
    
  public:
    double cut_g2g; //!< Cut-off distance for group-group interactions
    double R1;      //!< Radius of group 1 (used for DH correction)
    double R2;      //!< Radius of group 2 (used for DH correction)
    
    interaction_dipole(inputfile &in) : interaction<T>(in) {
      interaction<T>::name="Group-group dipole cut-off";
      cut_g2g = in.getflt( "threshold_g2g", 1e6 );
      kappa = 1/in.getflt("debyelen", 2e4);
      R1 = in.getflt("groupradius_g2g",20);
      R2 = R1;
      cnt=cntDip=0;
    }
    
    double energy(vector<particle> &p, const group &g1, const group &g2) {
      cnt++;
      dipole mu1, mu2;
      centerOfMass(p, g1, mu1.cm);
      centerOfMass(p, g2, mu2.cm);
      if ( sqrt(interaction<T>::pair.sqdist(mu1.cm, mu2.cm)) < cut_g2g )
        return interaction<T>::energy(p, g1, g2);
      else {
        cntDip++;
        calcDipole(p,g1,mu1);
        calcDipole(p,g2,mu2);
        mu1.anion.radius =R1-sqrt(interaction<T>::pair.sqdist(mu1.anion, mu1.cm));
        mu1.cation.radius=R1-sqrt(interaction<T>::pair.sqdist(mu1.cation, mu1.cm));
        mu2.anion.radius =R2-sqrt(interaction<T>::pair.sqdist(mu2.anion, mu2.cm));
        mu2.cation.radius=R2-sqrt(interaction<T>::pair.sqdist(mu2.cation, mu2.cm));
        return interaction<T>::pair.f * (
                                         f(mu1.anion, mu2.anion) * interaction<T>::pair.pairpot(mu1.anion, mu2.anion) +
                                         f(mu1.anion, mu2.cation) * interaction<T>::pair.pairpot(mu1.anion, mu2.cation) +
                                         f(mu1.cation, mu2.anion) * interaction<T>::pair.pairpot(mu1.cation, mu2.anion) +
                                         f(mu1.cation, mu2.cation) * interaction<T>::pair.pairpot(mu1.cation, mu2.cation) );
      }
    }
    
    string info() {
      std::ostringstream o;
      o << interaction<T>::info()
      << "#   Cut-off:" << endl
      << "#     Group-Group threshold   = " << cut_g2g << endl
      << "#     Group radii             = " << R1 << " " << R2 << endl;
      //<< "#     Dipole 1 radii (+,-)    = " << mu1.cation.radius << " " << mu1.anion.radius << endl
      //<< "#     Dipole 2 radii (+,-)    = " << mu2.cation.radius << " " << mu2.anion.radius << endl;
      if (cnt>100)
        o << "#     Cut-off fraction        = " << double(cntDip)/cnt << endl;
      return o.str();
    }
  };
  
};//namespace
#endif
