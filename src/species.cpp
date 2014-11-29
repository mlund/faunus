#include <faunus/species.h>
#include <faunus/inputfile.h>

namespace Faunus {

  bool AtomMap::includefile(const string &file) {
    return base::includefile(file);
  }

  bool AtomMap::includefile(InputMap &in) {
    string file = in("atomlist", string("") );
    return base::includefile(file);
  }

  bool MoleculeMap::includefile(const string &file) {
    return base::includefile(file);
  }

  bool MoleculeMap::includefile(InputMap &in) {
    string file = in("moleculelist", string("") );
    return base::includefile(file);
  }

  bool MoleculeCombinationMap::includefile(InputMap &in) {
    return base::includefile( in("moleculecombinations", string("") ) );
  }

  AtomMap atom; // Instantiate global copy
  MoleculeMap molecule;

}//namespace
