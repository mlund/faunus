#include <iostream>
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>

using namespace OpenBabel;
using namespace std;

int main()
{
  OBConversion obconversion;
  obconversion.SetInFormat("pdb");
  OBMol mol;

  bool notatend = obconversion.ReadFile(&mol,"1BRS.pdb");
  while (notatend)
  {
    cout << "Molecular Weight: " << mol.GetMolWt() << endl
         << "Residues = " << mol.NumResidues() << endl
         << "Atoms    = " << mol.NumAtoms() << endl
         << "NumRotors= " << mol.NumRotors() << endl;

    mol.Clear();
    notatend = obconversion.Read(&mol);
  }

  return(0);
};
