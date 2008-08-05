#include "iobabel.h"

void iobabel::read(string filename, string format) {
  obconv.SetInFormat(format.c_str());
  bool notatend = obconv.ReadFile(&mol, filename.c_str());
  while (notatend) {
    cout << "Molecular Weight: " << mol.GetMolWt() << endl
      << "Residues = " << mol.NumResidues() << endl
      << "Atoms    = " << mol.NumAtoms() << endl
      << "NumRotors= " << mol.NumRotors() << endl;
    mol.Clear();
    notatend = obconv.Read(&mol);
  }
}
