#ifndef FAUNUS_DRIFT_H
#define FAUNUS_DRIFT_H

#include <faunus/common.h>
#include <faunus/average.h>

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
      EnergyDrift& operator()(const std::pair<double,double>&);
      double checkDrift(const double&);
      string info();
      void test(UnitTest&);
  };

  EnergyDrift::EnergyDrift() {
    delta=initial=drift=0;
  }

  void EnergyDrift::init(const double &u0) {
    avg.reset();
    delta=drift=0;
    initial=u0;
    avg+=u0;
  }

  double EnergyDrift::current() const {
    return initial + delta;
  }

  EnergyDrift& EnergyDrift::operator+=(const double &du) {
    delta+=du;
    avg+=current();
    return *this;
  }

  double EnergyDrift::checkDrift(const double &snapshot) {
    drift = snapshot-current();
    return drift;
  }

  string EnergyDrift::info() {
    using namespace Faunus::textio;
    std::ostringstream o;
    o << header("System Energy and Drift");
    if (avg.cnt>0) {
      char w=25;
      o << textio::pad(SUB,w, "Average") << avg.avg() << kT << ", "
        << sigma << "=" << avg.stdev() << endl
        << textio::pad(SUB,w, "Initial energy") << initial << kT << endl
        << textio::pad(SUB,w, "Initial + changes") << current() << kT << endl;
      o.precision(4);
      o << pad(SUB,w, "Total energy drift") << drift << kT
        << " (" << drift/current()*100.
        << percent << ")" << endl;
    }
    return o.str();
  }

  void EnergyDrift::test(UnitTest &t) {
    //t("initialEnergy", initial, 1e-3);
    t("energyAverage", avg.avg() );
    t("relativeEnergyDrift", std::abs(drift/current()), 10.0 ); // allow 200% deviation    
  }
}//namespace
#endif
