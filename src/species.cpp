#include <faunus/species.h>
#include <faunus/inputfile.h>

namespace Faunus {

  bool AtomMap::includefile(const string &file) {
    return base::includefile(file, "atomlist" );
  }

  bool AtomMap::includefile(InputMap &in) {
    string file = in("atomlist", string("atom.json") );
    return includefile(file);
  }

  AtomMap atom; // Instantiate global copy

}//namespace
