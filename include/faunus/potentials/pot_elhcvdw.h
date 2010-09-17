#ifndef FAU_POT_ELHCVDW_H
#define FAU_POT_ELHCVDW_H
#include "faunus/potentials/base.h"
#include "faunus/species.h"

namespace Faunus {
  /*!
   * \brief Coulomb + HC - vdW potential. Minimum image
   * \author Mikael Lund
   * \date 2008
   */
  class pot_elhcvdw {
    private:
      vector< vector<double> >
        sigma2, // ij squared avg. diameter
        eps;    // ij quadratic epsilon
    
    public:
      string name;
      double f, len, lenh;
    
      pot_elhcvdw( inputfile &in ) {
        name+="Coulomb + HC - 4*eij*[(sij/r)^6] min. image";
        f=in.getflt("bjerrum",7.1);
        len=in.getflt("boxlen", 100);
        lenh=len/2.;
      }
    
      inline double pairpot(const particle &p1, const particle &p2) {
        double a=std::abs(p1.x-p2.x),
               b=std::abs(p1.y-p2.y),
               c=std::abs(p1.z-p2.z);
        if (a>lenh) a-=len;
        if (b>lenh) b-=len;
        if (c>lenh) c-=len;
        a=a*a+b*b+c*c; // squared dist.
        c=sigma2[p1.id][p2.id]; // (s1+s2)^2
        if (a<c)
          return 400.;
        b=c/a;
        return p1.charge*p2.charge/sqrt(a) - eps[p1.id][p2.id]*b*b*b;
      }
    
      inline double sqdist(const point &p1, const point &p2) { return p1.sqdist(p2); }
      string info() {
        std::ostringstream o;
        o << "#   Type              = " << name << std::endl
          << "#   Bjerrum length    = " << f << std::endl
          << "#   Box length        = " << len << std::endl;
        return o.str();
      }
    
      void setvolume(double v) {}
    
      void init(atoms &a) {
        short i,j,n=a.list.size();
        eps.resize(n);
        sigma2.resize(n);
        for (i=0; i<n; i++) {
          eps[i].resize(n);
          sigma2[i].resize(n);
          for (j=0; j<n; j++)
          {
            eps[i][j]=4/f*sqrt( a[i].eps * a[j].eps );
            sigma2[i][j]=std::pow( a[i].radius+a[j].radius, 2);
          }
        }
      }
  };
}//namespace
#endif
