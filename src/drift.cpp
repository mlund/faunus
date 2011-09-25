#include <faunus/drift.h>
#include <faunus/textio.h>
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
    using namespace Faunus::textio;
    std::ostringstream o;
    o << header("System Energy and Drift");
    if (avg.cnt>0) {
      o << textio::pad("Average",w,SUB) << avg.avg() << " kT, \u03C3=" << avg.stdev() << endl
        << textio::pad("Initial energy",w,SUB) << initial << " kT" << endl
        << textio::pad("Initial + changes",w,SUB) << current() << " kT" << endl;
      o.precision(4);
      o << pad("Total energy drift",w,SUB) << drift << " kT (" << drift/current()*100. << "\ufe6a)" << endl;
    }
    return o.str();
  }
}//namespace
