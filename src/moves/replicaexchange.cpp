#include "faunus/moves/replicaexchange.h"
#include "faunus/bottles/base.h"

namespace Faunus {
  
  bool replicaexchange::swap(bottle &i, bottle &j) {
    // Register systems
    string id=i.prefix+"<->"+j.prefix;
    if (cnt.find(id)<0) {
      cnt[id]=0;
      accepted[id]=0;
    }
    cnt[id]++;

    // Exchange particle positions
    // note: we change "p" not "trial" as usual this is because
    // the systemEnergy() function uses "p".
    for (int k=0; k<i.cPtr->p.size(); k++) {
      i.cPtr->p[k].x = j.cPtr->trial[k].x;
      i.cPtr->p[k].y = j.cPtr->trial[k].y;
      i.cPtr->p[k].z = j.cPtr->trial[k].z;
      j.cPtr->p[k].x = i.cPtr->trial[k].x;
      j.cPtr->p[k].y = i.cPtr->trial[k].y;
      j.cPtr->p[k].z = i.cPtr->trial[k].z;
    }

    // Exchange volumes
    double Vi=i.cPtr->getvolume();
    double Vj=j.cPtr->getvolume();
    double dPV=(i.P*Vj+j.P*Vi)-(i.P*Vi+j.P*Vj); // see doi:10.1063/1.1416491
    i.cPtr->setvolume(Vj); // both in container
    i.pPtr->setvolume(Vj); // and in hamiltonian...
    j.cPtr->setvolume(Vi);
    j.pPtr->setvolume(Vi);
    
    // Evaluate energies
    double inew,iold,jnew,jold;
    inew = i.systemEnergy();
    jnew = j.systemEnergy();
    iold = i.usys.sum;
    jold = j.usys.sum;
    
    if (nvt.metropolis( (inew+jnew)-(iold+jold) + dPV )==true ) {
      accepted[id]++;
      for (int k=0; k<i.cPtr->p.size(); k++) {
        i.cPtr->trial[k].x = i.cPtr->p[k].x;
        i.cPtr->trial[k].y = i.cPtr->p[k].y;
        i.cPtr->trial[k].z = i.cPtr->p[k].z;
        j.cPtr->trial[k].x = j.cPtr->p[k].x;
        j.cPtr->trial[k].y = j.cPtr->p[k].y;
        j.cPtr->trial[k].z = j.cPtr->p[k].z;
      }
      i.usys+=inew-iold;
      j.usys+=jnew-jold;
      return true;
    }
    for (int k=0; k<i.cPtr->p.size(); k++) {
      i.cPtr->p[k].x = i.cPtr->trial[k].x;
      i.cPtr->p[k].y = i.cPtr->trial[k].y;
      i.cPtr->p[k].z = i.cPtr->trial[k].z;
      j.cPtr->p[k].x = j.cPtr->trial[k].x;
      j.cPtr->p[k].y = j.cPtr->trial[k].y;
      j.cPtr->p[k].z = j.cPtr->trial[k].z;
    }
    i.cPtr->setvolume(Vi);   // Restore volumes
    i.pPtr->setvolume(Vi);   // in containers and potentials
    j.cPtr->setvolume(Vj);
    j.pPtr->setvolume(Vj);
    return false;
  }

  string replicaexchange::info() {
    std::ostringstream o;
    if (cnt.size()>0) {
      int w=12;
      o << "\n# PARALLEL TEMPERING / REPLICA EXCHANGE:" << endl
        << "#   More information:   doi:10.1039/b509983h  (review)" << endl
        << "#                       doi:10.1063/1.1416491 (npt ensemble)" << endl
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

