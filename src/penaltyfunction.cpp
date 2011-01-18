#include "faunus/energy/penaltyfunction.h"

namespace Faunus {

  penaltyfunction::penaltyfunction(double xmin, double xmax, double xres, double upenalty)
    : v(xres,xmin,xmax), cnt(xres,xmin,xmax)
  {
    number_of_updates=0;
    du=upenalty;
    for (int i=0; i<v.y.size(); i++) {
      v.y.at(i)=0;
      cnt.y.at(i)=0;
    }
  }
  
  void penaltyfunction::update(double x) {
    number_of_updates++;
    cnt(x)++;
    v(x)+=du;
  }

  double penaltyfunction::energy(double x) {
    return v(x);
  }
  
  void penaltyfunction::write(string filename) {
    std::ostringstream o;
    for (int x=v.xmin; x<v.xmax(); x+=v.xres) {
      o << x << " " << -v(x) << " " << -std::log( cnt(x)+1e-4 ) << std::endl;
    }
    std::cout << o.str();
  }
  
  string penaltyfunction::info() {
    std::ostringstream o;
    o << "#  Penalty function:\n"
      << "#     penalty energy (kT) = " << du << std::endl
      << "#     number of updates   = " << number_of_updates << std::endl;
    return o.str();
  }
}//namespace

