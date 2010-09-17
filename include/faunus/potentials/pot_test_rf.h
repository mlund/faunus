#ifndef FAU_POT_TEST_RF_H
#define FAU_POT_TEST_RF_H
#include "faunus/potentials/base.h"
#include "faunus/potentials/pot_rfield.h"
#include "faunus/species.h"

namespace Faunus {
  /*!
   * \brief Coulomb + LJ potential with pair parameters
   * \author Mikael Lund
   * \date 2008
   */
  class pot_test_rf {
    private:
      vector< vector<double> >
        sigma2, // ij squared avg. diameter
        eps;    // ij quadratic epsilon
    public:
      pot_rfield rf;
      string name;
      double f;
    
      pot_test_rf( inputfile &in ) : rf(in) {
        name+="Reactionfield + 4*eij*[(sij/r)^12-(sij/r)^6]";
        f=rf.f;
      }
    
      inline double pairpot(const particle &p1, const particle &p2) {
        register double a=p1.x-p2.x,
                        b=p1.y-p2.y,
                        c=p1.z-p2.z;
        a=1./(a*a+b*b+c*c);
        b=a*sigma2[p1.id][p2.id];
        c=b*b*b;
        //return eps[p1.id][p2.id] * (c*c-c) + p1.charge*p2.charge*sqrt(a);
        return eps[p1.id][p2.id] * (c*c-c) + rf.pairpot(p1,p2);
      }
    
      inline double sqdist(const point &p1, const point &p2) { return p1.sqdist(p2); }
      string info() {
        std::ostringstream o;
        o << "#   Type              = " << name << std::endl
          << "#   L-J parm. pairs   = " << eps.size() << std::endl
          << rf.info();
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
  
      void setvolume(double V) {}
    
  };
}//namespace
#endif
