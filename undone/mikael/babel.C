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

  bool notatend = obconversion.ReadFile(&mol,"test.pdb");
  while (notatend)
  {
    std::cout << "Molecular Weight: " << mol.GetMolWt() << std::endl;
    cout << "Residues = " << mol.NumResidues() << endl
         << "Atoms    = " << mol.NumAtoms() << endl
         << "NumRotors= " << mol.NumRotors() << endl;

    mol.Clear();
    notatend = obconversion.Read(&mol);
  }

  return(0);
};
