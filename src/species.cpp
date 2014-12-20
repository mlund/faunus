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

  /** @brief Read data from a picojson object */
  void MoleculeData::readJSON(const MoleculeData::Tjson &molecule) {
    string::size_type pos = 0;
    string::size_type oldPos = 0;
    string token;

    name = molecule.first;

    activity = json::value<double>(molecule.second, "activity", 0);
    chemPot = log( activity * 1.0_molar );

    _isAtomic = json::value<bool>(molecule.second, "atomic", false);

    string line = json::value<string>(molecule.second, "atoms", "Error");

    // read coordinates from disk
    string structure = json::value<string>( molecule.second, "structure", "" );
    if ( !structure.empty() ) {
      p_vec v;
      if ( FormatAAM::load( structure, v ) ) {
        if ( !v.empty() ) {
          cout << "# added molecular structure: " << structure << endl;
          conformations.push_back( v );
        }
      }
      else
        std::cerr << "# error loading molecule: " << structure << endl;
    }

    // tokenize atoms string and save as atom TID
    while(pos != string::npos) {
      pos = line.find_first_of(',', oldPos);
      token = line.substr(oldPos, pos-oldPos);
      oldPos = pos+1;

      for(auto &i : atom) {
        if(i.name.compare(token) == 0) {
          atoms.push_back(i.id);
          break;
        }
      }
    }

    // for atomic groups, add particles
    if ( isAtomic() ) {
      decltype(conformations)::value_type _vec;
      decltype(_vec)::value_type _p;
      for (auto id : atoms) {
        _p = atom[id];
        _vec.push_back( _p );
      }
      conformations.push_back( _vec );
    }
  }

  AtomMap atom; // Instantiate global copy
  MoleculeMap molecule;

}//namespace
