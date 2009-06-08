#include "faunus/iobabel.h"
namespace Faunus {

  iobabel::iobabel(atoms &a) {
    faunatomsPtr=&a;
    // Teach Babel some new "elements"!! (see openbabel "data.cpp")
    unsigned int n = OpenBabel::etab.GetNumberOfElements();
    for (unsigned int i=0; i<a.list.size(); i++) {
      std::ostringstream o;
      o << n++ << " "              // "atomic" number
        << a.list[i].name << " "   // symbol
        << 0 << " "                // AREneg
        << a.list[i].radius << " " // radius covalent
        << a.list[i].radius << " " // radius vdw
        << 20 << " "               // max bonds
        << a.list[i].mw << " "     // weight
        << "0 0 0 0 0 0 Faunus";
      OpenBabel::etab.ParseLine(o.str().c_str());     
    }
  }

  void iobabel::p2atom(particle &p) {
    atom.SetVector(p.x, p.y, p.z);
    atom.SetPartialCharge(p.charge);
  }

  void iobabel::read(string filename) {
    p.clear();
    mol.Clear();
    obconv.SetInFormat(obconv.FormatFromExt(filename.c_str()));
    bool notatend = obconv.ReadFile(&mol,filename.c_str());
    for (unsigned int i=1; i<=mol.NumAtoms(); i++)
      p.push_back(get(i));
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

  /*!
   * This function will convert between a babel atom and the
   * faunus particle approach. Coordinates and molecular weight
   * can be obtained from most file formats. Charge and radius
   * are currently supported only by the PQR format (at least
   * as far as we know).
   *
   * \note Since we cannot fetch the atomname from the original
   *       structure file, opened by babel, the recognition of
   *       particles is done via their molecular weight. That is,
   *       if you want to recognize a non-atomic particle one could
   *       use an exotic element and assign the same weight to a
   *       particle in the "faunatoms.dat" file. Ugly, we know but
   *       it'll have to do for now.
   */
  particle iobabel::get(unsigned int i) {
    atomPtr = mol.GetAtom(i);
    v=atomPtr->GetVector();
    v.Get(c);
    a.x=c[0]; a.y=c[1]; a.z=c[2];
    a.mw=atomPtr->GetAtomicMass();
    if (a.mw<1e-5)
      a.mw=1;   //we don't like weightless atoms.
    string name=string( OpenBabel::etab.GetSymbol( atomPtr->GetAtomicNum() ) );
    a.id=faunatomsPtr->get( faunatomsPtr->find(name) ).id;
    a.charge=atomPtr->GetPartialCharge();

    if (atomPtr->HasData("Radius")) {
      OpenBabel::OBPairData *gdat = dynamic_cast<OpenBabel::OBPairData *>( atomPtr->GetData("Radius") );
      a.radius = atof( gdat->GetValue().c_str() );
    } else
      a.radius=2;
    return a;
  }
}
