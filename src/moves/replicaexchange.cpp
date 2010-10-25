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
    // Due to double conventions, remember to copy both id and radius.
    for (int k=0; k<i.cPtr->p.size(); k++) {
      i.cPtr->p[k].x = j.cPtr->trial[k].x;
      i.cPtr->p[k].y = j.cPtr->trial[k].y;
      i.cPtr->p[k].z = j.cPtr->trial[k].z;
      i.cPtr->p[k].id = j.cPtr->trial[k].id;
      i.cPtr->p[k].radius = j.cPtr->trial[k].radius;
      j.cPtr->p[k].x = i.cPtr->trial[k].x;
      j.cPtr->p[k].y = i.cPtr->trial[k].y;
      j.cPtr->p[k].z = i.cPtr->trial[k].z;
      j.cPtr->p[k].id = i.cPtr->trial[k].id;
      j.cPtr->p[k].radius = i.cPtr->trial[k].radius;
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

    // Average energy differences
    dU[id] += jold-iold;
    expdU[id] += exp( -(jold-iold) ); 
    
    if (nvt.metropolis( (inew+jnew)-(iold+jold) + dPV )==true ) {
      accepted[id]++;
      for (int k=0; k<i.cPtr->p.size(); k++) {
        i.cPtr->trial[k].x = i.cPtr->p[k].x;
        i.cPtr->trial[k].y = i.cPtr->p[k].y;
        i.cPtr->trial[k].z = i.cPtr->p[k].z;
        i.cPtr->trial[k].id = i.cPtr->p[k].id;
        i.cPtr->trial[k].radius = i.cPtr->p[k].radius;
        j.cPtr->trial[k].x = j.cPtr->p[k].x;
        j.cPtr->trial[k].y = j.cPtr->p[k].y;
        j.cPtr->trial[k].z = j.cPtr->p[k].z;
        j.cPtr->trial[k].id = j.cPtr->p[k].id;
        j.cPtr->trial[k].radius = j.cPtr->p[k].radius;
      }
      i.usys+=inew-iold;
      j.usys+=jnew-jold;
      return true;
    }
    for (int k=0; k<i.cPtr->p.size(); k++) {
      i.cPtr->p[k].x = i.cPtr->trial[k].x;
      i.cPtr->p[k].y = i.cPtr->trial[k].y;
      i.cPtr->p[k].z = i.cPtr->trial[k].z;
      i.cPtr->p[k].id = i.cPtr->trial[k].id;
      i.cPtr->p[k].radius = i.cPtr->trial[k].radius;
      j.cPtr->p[k].x = j.cPtr->trial[k].x;
      j.cPtr->p[k].y = j.cPtr->trial[k].y;
      j.cPtr->p[k].z = j.cPtr->trial[k].z;
      j.cPtr->p[k].id = j.cPtr->trial[k].id;
      j.cPtr->p[k].radius = j.cPtr->trial[k].radius;
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
        << setw(w) << "Acceptance"
        << setw(w) << "<dU/kT>" << setw(w) << "dA/kT"<< endl;
      for (int i=0; i<cnt.size(); i++) {
        o << "#     " << setw(w) << cnt.id(i)
          << setw(8) << cnt.at(i)
          << setw(w) << accepted.at(i)*100./cnt.at(i);
        if (dU.at(i).cnt>0)
          o << setw(w) << dU.at(i).avg() << setw(w) << -log(expdU.at(i).avg());
        o << endl;
      }
    }
    return o.str();
  }

}//namespace

