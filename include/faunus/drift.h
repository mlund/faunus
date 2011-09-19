#ifndef FAU_DRIFT_H
#define FAU_DRIFT_H
#include <faunus/common.h>
#include <faunus/average.h>

namespace Faunus {
  class energydrift {
    private:
      double delta;
      double initial;
      average<double> avg;
    public:
      double drift;
      energydrift();
      void init(const double&);
      double current() const;
      energydrift& operator+=(const double&);
      double checkdrift(const double&);
      string info();
  };
}//namespace
#endif
