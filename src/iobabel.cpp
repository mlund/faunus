#include "faunus/iobabel.h"
namespace Faunus {

  iobabel::iobabel(atoms &a) { faunatomsPtr=&a; }

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
      //mol.Clear();
      notatend = obconv.Read(&mol);
    }
  }

  vector<unsigned short> iobabel::neighbors(unsigned short i) {
    vector<unsigned short> nb;
    if (mol.NumAtoms()>0) {
      OpenBabel::OBAtom* atomPtr, nbrPtr;
      atomPtr = mol.GetAtom(i+1);  // Babel starts at 1; Faunus at 0.
      FOR_NBORS_OF_ATOM(nbrPtr, atomPtr) {  // Maybe FOR_BONDS_OF_ATOM?
        nb.push_back( nbrPtr->GetIdx()-1 );
      }
    }
    return nb;
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
    a.id=faunatomsPtr->get(faunatomsPtr->find( "NA" )).id;
    a.x=c[0]; a.y=c[1]; a.z=c[2];
    a.mw=atomPtr->GetAtomicMass();
    if (a.mw<1e-5)
      a.mw=1;   //we don't like weightless atoms.
    a.charge=atomPtr->GetPartialCharge();

    if (atomPtr->HasData("Radius")) {
      OpenBabel::OBPairData *gdat = dynamic_cast<OpenBabel::OBPairData *>( atomPtr->GetData("Radius") );
      a.radius = atof( gdat->GetValue().c_str() );
    } else
      a.radius=2;
    return a;
  }
}
