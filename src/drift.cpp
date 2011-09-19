#include <faunus/drift.h>
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
    std::ostringstream o;
    o << endl << "# SYSTEM ENERGY (kT):" << endl;
    if (avg.cnt>0) {
      o << "#   Average <U> s      = " << avg.avg() << " " << avg.stdev() << "\n"
        << "#   Initial energy     = " << initial << endl
        << "#   Initial + changes  = " << current() << endl;
      o.precision(4);
      o << "#   Total energy drift = " << drift << " (" << drift/current()*100. << "\ufe6a)" << endl;
    }
    return o.str();
  }
}//namespace
