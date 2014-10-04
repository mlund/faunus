#include <faunus/species.h>
#include <faunus/inputfile.h>

namespace Faunus {

  bool AtomMap::includefile(const string &file) {
    return base::includefile(file);
  }

  bool AtomMap::includefile(InputMap &in) {
    string file = in("atomlist", string("atom.json") );
    return base::includefile(file);
  }

  bool MoleculeMap::includefile(const string &file) {
    return base::includefile(file);
  }

  bool MoleculeMap::includefile(InputMap &in) {
    string file = in("topology", string("topo.json") );
    return base::includefile(file);
  }

  bool MoleculeComboMap::includefile(InputMap &in) {
    return base::includefile( in("moleculecombo", string("") ) );
  }

  AtomMap atom; // Instantiate global copy
  MoleculeMap molecule;

}//namespace
