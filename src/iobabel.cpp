#include "faunus/iobabel.h"
namespace Faunus {
  void iobabel::p2atom(particle &p) {
    atom.SetVector(p.x, p.y, p.z);
    atom.SetPartialCharge(p.charge);
  }
  void iobabel::read(string filename) {
    p.clear();
    mol.Clear();
    obconv.SetInFormat(obconv.FormatFromExt(filename.c_str()));
    bool notatend = obconv.ReadFile(&mol,filename.c_str());
    while (notatend) {
      for (unsigned int i=1; i<=mol.NumAtoms(); i++)
        p.push_back(get(i));
      mol.Clear();
      notatend = obconv.Read(&mol);
    }
  }
  bool iobabel::write(string filename, const vector<particle> &) {
    mol.Clear();
    for (unsigned int i=0; i<p.size(); i++) {
      p2atom(p[i]);
      mol.AddAtom(atom);
    }
    obconv.SetOutFormat(obconv.FormatFromExt(filename.c_str()));
    return obconv.WriteFile(&mol,filename.c_str());
  }
  particle iobabel::get(unsigned int i) {
    atomPtr = mol.GetAtom(i);
    v=atomPtr->GetVector();
    v.Get(c);
    a.x=c[0]; a.y=c[1]; a.z=c[2];
    a.mw=atomPtr->GetAtomicMass();
    if (a.mw==0)
      a.mw=1;  //we don't like weightless atoms.
    a.charge=atomPtr->GetPartialCharge();
    return a;
  }
}
