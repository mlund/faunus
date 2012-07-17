#ifndef FAU_DRIFT_H
#define FAU_DRIFT_H

#ifndef SWIG
#include <faunus/common.h>
#include <faunus/average.h>
#endif

namespace Faunus {
  class EnergyDrift {
    private:
      double delta;
      double initial;
      Average<double> avg;
    public:
      double drift;
      EnergyDrift();
      void init(const double&);
      double current() const;
      EnergyDrift& operator+=(const double&);
      double checkDrift(const double&);
      string info();
      void test(UnitTest&);
  };
}//namespace
#endif
