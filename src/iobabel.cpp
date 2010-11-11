#ifdef BABEL
#include "faunus/iobabel.h"
namespace Faunus {

  iobabel::iobabel() {
  }

  void iobabel::p2atom(particle &p) {
    obatom.SetVector(p.x, p.y, p.z);
    obatom.SetPartialCharge(p.charge);
  }

  /*!
   * This will read a molecular file format into the particle
   * vector. Residue name is used to identify the particle and
   * one should make sure that the number of residues match the
   * number of particles.
   */
  void iobabel::read(string filename) {
    p.clear();
    obmol.Clear();
    obconv.SetInFormat(obconv.FormatFromExt(filename.c_str()));
    bool notatend = obconv.ReadFile(&obmol,filename.c_str());
    for (int i=1; i<=obmol.NumAtoms(); i++)
      p.push_back(get(i));

    // Did number of residues correspond to number of atoms?
    if (obmol.NumResidues() != p.size() )
      std::cerr << "# Warning while loading " << filename << " - Number of residues != Number of particles.\n";
  }

  vector<int> iobabel::neighbors(int i) {
    vector<int> nb;
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
    for (int i=0; i<p.size(); i++) {
      p2atom(p[i]);
      obmol.AddAtom(obatom);
    }
    obconv.SetOutFormat(obconv.FormatFromExt(filename.c_str()));
    return obconv.WriteFile(&obmol,filename.c_str());
  }

  /*!
   * This function will convert between a babel atom and the
   * faunus particle approach. In OB atoms are real atoms from
   * the periodic table which is usually not particularly useful
   * in Faunus. We therefore use the RESIDUE name to identify the
   * particle.
   */
  particle iobabel::get(int i) {
    particle a;
    double c[3];
    obatomPtr = obmol.GetAtom(i);        // Get pointer to i'th OB atom
    obresPtr  = obatomPtr->GetResidue(); // Get residue that i'th OB atom belongs to
    a=atom( obresPtr->GetName() );       // and use residue name to determine Faunus species type
    
    v=obatomPtr->GetVector();            // Get XYZ coordinates
    v.Get(c);
    a.x=c[0]; a.y=c[1]; a.z=c[2];
    
    a.charge=obatomPtr->GetPartialCharge(); // Get partial charge
    
    if (obatomPtr->HasData("Radius")) {
      OpenBabel::OBPairData *gdat = dynamic_cast<OpenBabel::OBPairData *>( obatomPtr->GetData("Radius") );
      a.radius = atof( gdat->GetValue().c_str() );
    }
        
    if (atom[ obresPtr->GetName() ].id==0)
      std::cerr << "# Warning: OpenBabel atom " << i << " is unknown\n";
    
    return a;
  }
}//namespace
#endif

