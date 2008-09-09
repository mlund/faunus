#include "iobabel.h"
namespace Faunus {
  void iobabel::read(string filename, string format) {
    p.clear();
    mol.Clear();
    obconv.SetInFormat(format.c_str());
    bool notatend = obconv.ReadFile(&mol,filename.c_str());
    while (notatend) {
      for (unsigned int i=1; i<=mol.NumAtoms(); i++)
        p.push_back(get(i));
      mol.Clear();
      notatend = obconv.Read(&mol);
    }
  }
  particle iobabel::get(unsigned int i) {
    atom = mol.GetAtom(i);
    v=atom->GetVector();
    v.Get(c);
    a.x=c[0]; a.y=c[1]; a.z=c[2];
    a.mw=atom->GetAtomicMass();
    if (a.mw==0)
      a.mw=1;  //we don't like weightless atoms.
    a.charge=atom->GetPartialCharge();
    return a;
  }
}
