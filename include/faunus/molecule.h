#ifndef faunus_molecule_h
#define faunus_molecule_h

#include <faunus/io.h>
#include <faunus/inputfile.h>

namespace Faunus {

  /**
   * @brief Storage for molecular properties
   *
   * Values can be read from a JSON object with the following format:
   *
   * ~~~~
   * "salt": { "atoms" : "Na,Cl", "init":"RANDOM" },
   * "polymer": { "activity" : 0.05, "atoms" : "MM MM MM MM" }
   * ~~~~
   *
   * The key of type string is the `name` followed, in no particular order,
   * by properties:
   *
   * Key           | Description
   * :------------ | :-------------------------------------------------------
   * `activity`    | Chemical activity for grand canonical MC [mol/l]
   * `atoms`       | List of atoms in molecule (use `AtomData` names)
   * `structure`   | Read conformation and species from AAM file
   *
   */
  template<class Tpvec>
    class MoleculeData  : public PropertyBase {
      private:

        typedef typename Tpvec::value_type Tparticle;
        using PropertyBase::Tjson;
        bool _isAtomic;

      public:

        std::vector<typename Tparticle::Tid> atoms; //!< List of atoms in molecule
        vector<Tpvec> conformations;                //!< Conformations of molecule
        double activity;                            //!< Chemical activity (mol/l)
        double chemPot;                             //!< Chemical potential (kT)
        Average<double> rho;                        //!< Average concentration (1/angstrom^3)

        /** @brief Constructor - by default data is initialized; mass set to unity */
        inline MoleculeData( const Tjson &molecule=Tjson()) : _isAtomic(false) {
          readJSON(molecule);
        }

        /** @brief Get a random conformation */
        p_vec getRandomConformation() const {
          assert(!conformations.empty());
          return conformations[slp_global.range(0, conformations.size()-1) ];
        }

        /**
         * @brief Store a single configuration
         * @param vec of particles
         */
        void pushConfiguration(const Tpvec& vec) {
          conformations.push_back(vec);
        }

        bool isAtomic() const { return _isAtomic; }

        bool operator==(const MoleculeData &d) const { return (*this==d); }

        /** @brief Read data from a picojson object */
        void readJSON(const Tjson &molecule) FOVERRIDE {

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
            Tpvec _vec;
            Tparticle _p;
            for (auto id : atoms) {
              _p = atom[id];
              _vec.push_back( _p );
            }
            conformations.push_back( _vec );
          }
        }

        inline string info() {

          using namespace textio;
          ostringstream o;
          char w=22;

          o << header("Molecule: " + name)
            << pad(SUB, w, "Activity") << activity << " mol/l\n"
            << pad(SUB, w, "Chem. pot.") << chemPot << " kT\n"
            << pad(SUB, w, "Atomic") << std::boolalpha << isAtomic() << "\n"
            << pad(SUB, w, "Conformations") << conformations.size() << "\n"
            << pad(SUB, w, "Number of atoms") << std::max(atoms.size(), conformations.size()) << "\n";

          if (rho.cnt>0) {
            o << pad(SUB, w, "Average conc.") << rho.avg() << " 1/"+angstrom+" = "
              << rho.avg() / 1.0_molar << " mol/l\n"
              << pad(SUB, w, "Activity coefficient") << activity / rho.avg() << "\n";
          }

          return o.str();
        }

    };

  /**
   * @brief Class for loading and storing Molecule properties
   *
   * This will load molecule properties from disk and store them in a
   * vector of `MoleculeData`. The file format is JSON (<http://www.json.org>)
   * and all molecule properties must be inclosed in an object with
   * the keyword `moleculelist`.
   * While not strictly JSON compliant, comments beginning
   * with `//` are allowed.
   *
   * For example:
   *
   * ~~~~
   * {
   *   "moleculelist" : {
   *     "salt" : { "atoms": "Na,Cl" },
   *     "salt2" : { "atoms": "Mg,Cl,Cl" },
   *     "chloride" : { "atoms": "Cl,Cl,Cl,Cl" },
   *     "polymer" : { "activity": 0.05, "atoms": "MM,MM,MM,MM" },
   *     "polymer2" : { "activity": 0.05, "atoms": "MM,MM,MM,MM" }
   *   }
   * }
   * ~~~~
   *
   * Note that faunus currently has a global instance of `MoleculeMap`,
   * simply named `molecule`. This can be accessed from anywhere.
   *
   * @todo More documentation
   */
  template<class Tpvec, class base=PropertyVector<MoleculeData<Tpvec>>>
    class MoleculeMap : public base {
      public:
        using typename base::value_type;

        MoleculeMap() {
          base::jsonsection = "moleculelist";
          base::name = "Molecule Properties";
          base::push_back( typename base::value_type() ); // add default property
        }

        /** @brief Read JSON file given through `InputMap` */
        bool includefile(InputMap &in) {
          string file = in("moleculelist", string("") );
          return base::includefile(file);
        }

        /** @brief Information string */
        string info() {
          char w=15;
          using namespace textio;
          ostringstream o;
          o << header("Molecule list") << "\n"
            << std::left << setw(w+4) << "  Molecule" << setw(w) << "activity"
            << setw(w) << "atoms" << setw(w) << "atomic" << "\n"
            << "  " << string(4*w, '-') << "\n";

          for (auto &i : *this)
            if ( !i.name.empty() )
              o << setw(w+4) << "  "+i.name << setw(w) << i.activity
                << setw(w) << std::max(i.conformations.size(), i.atoms.size())
                << std::boolalpha << i.isAtomic()
                << "\n";
          return o.str();
        }
    };

  /**
   * @brief Molecular combinations for Grand Canonical moves, Widom insertion etc.
   *
   * JSON entry examples:
   *
   * ~~~~
   * "moleculelist" : {
   *    "ion1" : { "atoms" : "Na" },
   *    "ion2" : { "atoms" : "Cl" }
   * },
   *
   * "moleculecombinations" : {
   *    "NaCl" : { "molecules" : "ion1 ion2", "prob" : 0.5 },
   * }
   * ~~~~
   *
   * The key of type string is the `name` followed, in no particular order,
   * by properties:
   *
   * Key           | Description
   * :------------ | :--------------------------------------
   * `molecules`   | list of space separated molecule names
   * `prob`        | insertion probability [0:1] (default: 1=100%)
   *
   * the final probability of combination is: prob_i / (sum_all prob_i)
   */
  template<class Tpvec>
    class MoleculeCombination : public PropertyBase {

      private:

        inline void readJSON(const Tjson &comb) FOVERRIDE {
          molComb.clear();
          name = comb.first;
          probability = json::value<double>( comb.second, "prob", 1.0 );
          string mollist = json::value<string>( comb.second, "molecules", "" );
          molName = textio::words2vec<string>(mollist);
        }

      public:

        MoleculeCombination(const Tjson &comb=Tjson()) { readJSON(comb); }

        vector<PropertyBase::Tid> molComb;/// vector of molecule index in combination

        vector<string> molName;           /// vector of molecule names in combination

        double probability;               /// probability of this combination in GC-MC move

    };

  /**
   * @brief Vector of molecular combinations
   *
   * When examining a JSON file, individual entries must be placed
   * in a section called `moleculecombinations`. When loading
   * from `InputMap` the keyword for the json file is
   * also `moleculecombinations`.
   */
  template<class Tpvec, class base=PropertyVector< MoleculeCombination<Tpvec> >>
    class MoleculeCombinationMap : public base {
      private:

        MoleculeMap<Tpvec> *molmap; // typically points to `Space::molecule`

        void update() { // convert name (string) -> id (int) using MoleculeMap
          for (auto &i : *this) {
            assert( !i.molName.empty() && "Molecule combination cannot be empty");
            for (auto &name : i.molName) {
              auto it = molmap->find(name);
              if ( it!=molmap->end() )
                i.molComb.push_back(it->id);
            }
            assert( i.molName.size() == i.molComb.size() );
            assert( !i.molName.empty() );
          }
        }

      public:

        bool includefile(InputMap &in) {
          bool rc = base::includefile( in(base::jsonsection, string("") ) );
          update();
          return rc;
        }

        bool includefile(const string& file) {
          bool rc = base::includefile(file);
          update();
          return rc;
        }

        MoleculeCombinationMap(MoleculeMap<Tpvec> &moleculemap) : molmap(&moleculemap) {
          assert( molmap != nullptr);
          base::name = "Molecule Combinations";
          base::jsonsection = "moleculecombinations";
        }

        string info() {
          using namespace textio;
          char w=25;
          std::ostringstream o;
          o << header(base::name)
            << pad(SUB,w,"Number of entries:") << base::size() << endl;
          if (!base::jsonfile.empty())
            o << pad(SUB,w,"Input JSON file:") << base::jsonfile << endl;
          o << indent(SUB) << "Element info:\n";
          for (auto &i : *this) {
            o << std::left << setw(w) << "    "+i.name+":";
            for (auto &n : i.molName)
              o << setw(10) << n;
            o << "\n";
          }
          o << endl;
          return o.str();
        }

    };


}//namespace

#endif
