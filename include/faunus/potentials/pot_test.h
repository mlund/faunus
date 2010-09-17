#ifndef FAU_POT_TEST_H
#define FAU_POT_TEST_H
#include "faunus/potentials/base.h"
#include "faunus/species.h"

namespace Faunus {
  /*!
   * \brief Coulomb + LJ potential with pair parameters
   * \author Mikael Lund
   * \date 2008
   */
  class pot_test {
    private:
      vector< vector<double> >
        sigma2, // ij squared avg. diameter
        eps;    // ij quadratic epsilon
    public:
      string name;
      double f;
      pot_test( inputfile &in ) {
        name+="Coulomb + 4*eij*[(sij/r)^12-(sij/r)^6]";
        f=in.getflt("bjerrum",7.1);
      }
      inline double pairpot(const particle &p1, const particle &p2) {
        register double a=p1.x-p2.x,
                        b=p1.y-p2.y,
                        c=p1.z-p2.z;
        a=1./(a*a+b*b+c*c);
        b=a*sigma2[p1.id][p2.id];
        c=b*b*b;
        return eps[p1.id][p2.id] * (c*c-c) + p1.charge*p2.charge*sqrt(a);
      }
      inline double sqdist(const point &p1, const point &p2) { return p1.sqdist(p2); }
      string info() {
        std::ostringstream o;
        o << "#   Type              = " << name << std::endl
          << "#   L-J parm. pairs   = " << eps.size() << std::endl
          << "#   Bjerrum length    = " << f << std::endl;
        return o.str();
      }
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

  class pot_testminim {
    private:
      vector< vector<double> >
        sigma2, // ij squared avg. diameter
        eps;    // ij quadratic epsilon
    public:
      string name;
      double f, len, lenh;
    
      void setvolume(double v) {
        len=pow(v,1/3.);
        lenh=len/2;
      }
    
      pot_testminim( inputfile &in ) {
        name+="Coulomb + 4*eij*[(sij/r)^12-(sij/r)^6] minimum image";
        f=in.getflt("bjerrum",7.1);
        len=in.getflt("boxlen", -1);
        lenh=len/2.;
      }
    
      inline double pairpot(const particle &p1, const particle &p2) {
        register double a=std::abs(p1.x-p2.x),
                        b=std::abs(p1.y-p2.y),
                        c=std::abs(p1.z-p2.z);
                        if (a>lenh) a-=len;
                        if (b>lenh) b-=len;
                        if (c>lenh) c-=len;
        a=1./(a*a+b*b+c*c);
        b=a*sigma2[p1.id][p2.id];
        c=b*b*b;
        return eps[p1.id][p2.id] * (c*c-c) + p1.charge*p2.charge*sqrt(a);
      }
    
      inline double sqdist(const point &p1, const point &p2) { return p1.sqdist(p2); }
      string info() {
        std::ostringstream o;
        o << "#   Type              = " << name << std::endl
          << "#   L-J parm. pairs   = " << eps.size() << std::endl
          << "#   Bjerrum length    = " << f << std::endl
          << "#   Box length        = " << len << std::endl;
        return o.str();
      }
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
