#ifdef BABEL
#include "faunus/iobabel.h"
namespace Faunus {

  iobabel::iobabel() {
    // Teach Babel some new "elements"!! (see openbabel "data.cpp")
    unsigned int n = OpenBabel::etab.GetNumberOfElements();
    for (unsigned int i=0; i<atom.list.size(); i++) {
      std::ostringstream o;
      o << n++ << " "              // "atomic" number
        << atom.list[i].name << " "   // symbol
        << 0 << " "                // AREneg
        << atom.list[i].radius << " " // radius covalent
        << atom.list[i].radius << " " // radius vdw
        << 20 << " "               // max bonds
        << atom.list[i].mw << " "     // weight
        << "0 0 0 0 0 0 Faunus";
      OpenBabel::etab.ParseLine(o.str().c_str());     
    }
  }

  void iobabel::p2atom(particle &p) {
    obatom.SetVector(p.x, p.y, p.z);
    obatom.SetPartialCharge(p.charge);
  }

  /*!
   * This will read a molecular file format into the particle
   * vector. If the structure file contains residue names and
   * if the number of pacticles is equal to the number of
   * residues, the residue name will be used to identify the
   * particle type.
   */
  void iobabel::read(string filename, bool resnaming) {
    p.clear();
    obmol.Clear();
    obconv.SetInFormat(obconv.FormatFromExt(filename.c_str()));
    bool notatend = obconv.ReadFile(&obmol,filename.c_str());
    for (unsigned int i=1; i<=obmol.NumAtoms(); i++)
      p.push_back(get(i));

    // Use residue names for particle recognition
    if ( resnaming==true && obmol.NumResidues() == p.size() ) {
      unsigned int i=0;
      OpenBabel::OBResidue* obres;
      FOR_RESIDUES_OF_MOL(obres, obmol) {
        p[i++].id = atom[ obres->GetName() ].id;
      }
    }
  }

  vector<unsigned short> iobabel::neighbors(unsigned short i) {
    vector<unsigned short> nb;
    if (obmol.NumAtoms()>0) {
      OpenBabel::OBAtom* obatomPtr, nbrPtr;
      obatomPtr = obmol.GetAtom(i+1);  // Babel starts at 1; Faunus at 0.
      FOR_NBORS_OF_ATOM(nbrPtr, obatomPtr) {  // Maybe FOR_BONDS_OF_ATOM?
        nb.push_back( nbrPtr->GetIdx()-1 );
      }
    }
    return nb;
  }

  bool iobabel::write(string filename, const vector<particle> &) {
    obmol.Clear();
    for (unsigned int i=0; i<p.size(); i++) {
      p2atom(p[i]);
      obmol.AddAtom(obatom);
    }
    obconv.SetOutFormat(obconv.FormatFromExt(filename.c_str()));
    return obconv.WriteFile(&obmol,filename.c_str());
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
    obatomPtr = obmol.GetAtom(i);
    v=obatomPtr->GetVector();
    v.Get(c);
    a.x=c[0]; a.y=c[1]; a.z=c[2];
    a.mw=obatomPtr->GetAtomicMass();
    if (a.mw<1e-5)
      a.mw=1;   //we don't like weightless atoms.
    string name=string( OpenBabel::etab.GetSymbol( obatomPtr->GetAtomicNum() ) );
    a.id=atom[name].id;
    a.charge=obatomPtr->GetPartialCharge();
//    std::cout << int(a.id) << " " << name << " " << obatomPtr->GetAtomicNum() << std::endl;

    if (obatomPtr->HasData("Radius")) {
      OpenBabel::OBPairData *gdat = dynamic_cast<OpenBabel::OBPairData *>( obatomPtr->GetData("Radius") );
      a.radius = atof( gdat->GetValue().c_str() );
    } else
      a.radius=2;
    return a;
  }
}//namespace
#endif

