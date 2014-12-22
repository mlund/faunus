#include <faunus/species.h>
#include <faunus/inputfile.h>
#include <faunus/io.h>

namespace Faunus {

  bool AtomMap::includefile(const string &file) {
    return base::includefile(file);
  }

  bool AtomMap::includefile(InputMap &in) {
    string file = in("atomlist", string("") );
    return base::includefile(file);
  }

  AtomMap atom; // Instantiate global copy

}//namespace
