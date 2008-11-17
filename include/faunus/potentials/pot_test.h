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
      vector< vector<double> > sigma2, eps;
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
        eps=a.eps;
        sigma2=a.sigma;
        short i,j,n=eps.size();
        for (i=0; i<n; i++)
          for (j=0; j<n; j++)
          {
            eps[i][j]=4/f*a.eps[i][j];
            sigma2[i][j]=pow(a.sigma[i][j],2);
          }
      }
  };
}//namespace
#endif
