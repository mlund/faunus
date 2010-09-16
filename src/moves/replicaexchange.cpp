#include "faunus/moves/replicaexchange.h"
#include "faunus/bottles/base.h"

namespace Faunus {

  bool replicaexchange::swap(simBottle &i, simBottle &j) {
    string id=i.prefix+"<->"+j.prefix;
    if (cnt.find(id)<0) {
      cnt[id]=0;
      accepted[id]=0;
    }
    cnt[id]++;

    i.cPtr->trial = j.cPtr->p;
    j.cPtr->trial = i.cPtr->p;

    double inew,iold,jnew,jold;
    inew = i.systemEnergy();
    iold = i.usys.sum;
    jnew = j.systemEnergy();
    jold = j.usys.sum;

    if (nvt.metropolis( (inew+jnew)-(iold+jold) )==true ) {
      accepted[id]++;
      i.cPtr->p = i.cPtr->trial;
      j.cPtr->p = j.cPtr->trial;
      i.usys+=inew-iold;
      j.usys+=jnew-jold;
      return true;
    }
    i.cPtr->trial = i.cPtr->p;
    j.cPtr->trial = j.cPtr->p;
    return false;
  }

  string replicaexchange::info() {
    std::ostringstream o;
    if (cnt.size()>0) {
      int w=12;
      o << "\n# PARALLEL TEMPERING / REPLICA EXCHANGE:" << endl
        << "#   More information:   doi:10.1039/b509983h" << endl
        << "#   Exchange statistics:" << endl
        << "#     " << std::left << setw(w) << "System"
        << setw(8) << "Trials"
        << setw(w) << "Acceptance" << endl;
      for (int i=0; i<cnt.size(); i++) {
        o << "#     " << setw(w) << cnt.id(i)
          << setw(8) << cnt.at(i)
          << setw(w) << accepted.at(i)*100./cnt.at(i) << endl;
      }
    }
    return o.str();
  }

}//namespace

