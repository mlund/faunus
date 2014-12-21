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

    name = molecule.first;

    activity = json::value<double>(molecule.second, "activity", 0);
    chemPot = log( activity * 1.0_molar );

    _isAtomic = json::value<bool>(molecule.second, "atomic", false);

    // read conformation from disk
    string structure = json::value<string>( molecule.second, "structure", "" );
    if ( !structure.empty() ) {
      p_vec v;
      if ( FormatAAM::load( structure, v ) ) {
        if ( !v.empty() ) {
          conformations.push_back( v ); // add conformation
          for ( auto &p : v )           // add atoms to atomlist
            atoms.push_back(p.id);
          cout << "# added molecular structure: " << structure << endl;
        }
      }
      else
        std::cerr << "# error loading molecule: " << structure << endl;
    }

    // add atoms to atom list
    if ( atoms.empty() ) {
      string atomlist = json::value<string>(molecule.second, "atoms", "");
      for ( auto &a : textio::words2vec<string>(atomlist) )
        if ( atom[a].id > 0 )
          atoms.push_back( atom[a].id );
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
    
    //assert( id==0 || !atoms.empty() );
    
  }

  AtomMap atom; // Instantiate global copy
  MoleculeMap molecule;

}//namespace
