#include <faunus/drift.h>
#include <faunus/faunus.h>
namespace Faunus {
  energydrift::energydrift() {
    delta=initial=drift=0;
  }

  void energydrift::init(const double &u0) {
    avg.reset();
    delta=drift=0;
    initial=u0;
    avg+=u0;
  }

  double energydrift::current() const {
    return initial + delta;
  }

  energydrift& energydrift::operator+=(const double &du) {
    delta+=du;
    avg+=current();
    return *this;
  }

  double energydrift::checkdrift(const double &snapshot) {
    drift = snapshot-current();
    return drift;
  }

  string energydrift::info() {
    char w=25;
    std::ostringstream o;
    o << header("System Energy (kT)");
    if (avg.cnt>0) {
      o << pad("Average <U> s",w,SUB) << avg.avg() << " " << avg.stdev() << "\n"
        << pad("Initial energy",w,SUB) << initial << endl
        << pad("Initial + changes",w,SUB) << current() << endl;
      o.precision(4);
      o << pad("Total energy drift",w,SUB) << drift << " (" << drift/current()*100. << "\ufe6a)" << endl;
    }
    return o.str();
  }
}//namespace
