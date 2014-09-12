#ifndef FAUNUS_SPECIES_H
#define FAUNUS_SPECIES_H

#ifndef SWIG
#include <faunus/common.h>
#include <faunus/json.h>
#include <faunus/physconst.h>
#include <faunus/textio.h>
#include <faunus/point.h>
#include <faunus/slump.h>
#endif

namespace Faunus {

  /**
   * @brief Base class for property data
   *
   * Used to collect properties for atoms, molecules etc.
   */
  class PropertyBase {
  public:
    typedef unsigned char Tid;
    typedef picojson::object::value_type Tjson;
    std::string name;  //!< Unique, user-defined name
    Tid id;            //!< Unique id (automatically set by `PropertyVector`)
    virtual void readJSON(const Tjson&)=0; //!< Process single json entry
    virtual ~PropertyBase() {};
    PropertyBase() : name("UNK") {}
  };

  /**
   * @brief Storage for atomic properties
   *
   * This stores properties about single atoms and can be used together
   * with `PropertyVector`. Values can be read from a JSON
   * object with the following format:
   *
   * ~~~~
   * "Na" : { "q": 1.0, "r":1.9, "mw":22.99 }, // sodium ion
   * ~~~~
   *
   * or more advanced,
   *
   * ~~~~
   * "pol" :
   *    {
   *      "r":2.0,
   *      "mu":"1.3 0.1 0",
   *      "alpha":"1.1   0   0  2.3   0    1",    // Matrix input in 
   *      "theta":"2.4 0.3   0  1.8   0  3.2"     // (xx, xy, xz, yy, yz, zz) form
   *    }
   * ~~~~
   *
   * The key of type string is the `name` followed, in no particular order,
   * by properties:
   *
   * Key           | Description
   * :------------ | :-------------------------------------------------------
   * `activity`    | Chemical activity for grand canonical MC [mol/l]
   * `alpha`       | Polarizability in units of [\f$ 4\pi\epsilon_0 \f$]
   * `dp`          | Translational displacement parameter [angstrom] 
   * `dprot`       | Rotational displacement parameter [degrees] (will be converted to radians)
   * `eps`         | Epsilon energy scaling commonly used for Lennard-Jones interactions etc. [kJ/mol] 
   * `hydrophobic` | Is the particle hydrophobic? [true/false]
   * `mu`          | Dipole moment vector [Debye]
   * `theta`       | Quadrupole moment tensor [Debye \f$ \unicode{x212B} \f$]
   * `mw`          | Molecular weight [g/mol]
   * `patchtype`   | Patchtype for sphero-cylinders
   * `q`           | Valency / partial charge number [e]
   * `r`           | Radius = `sigma/2` [angstrom]
   * `sigma`       | `2*r` [angstrom] (overrides radius)
   * `tfe`         | Transfer free energy [J/mol/angstrom^2/M] (default: 0)
   */

  class AtomData : public PropertyBase {
  public:
    using PropertyBase::Tjson;
    double sigma,            //!< LJ diameter [angstrom]
           eps,              //!< LJ epsilon [kJ/mol] (pair potentials should convert to kT)
           radius,           //!< Radius [angstrom]
           muscalar,         //!< Dipole momentscalar [e \f$ \unicode{x212B} \f$]
           mw,               //!< Weight [g/mol]
           charge,           //!< Charge/valency [e]
           activity,         //!< Chemical activity [mol/l]
           chemPot,          //!< Chemical potential, calculated from activity
           dp,               //!< Translational displacement parameter [angstrom]
           dprot,            //!< Rotational displacement parameter [degrees]
           mean,             //!< Mean value... (charge, sasa, etc.)
           variance,         //!< Spread around AtomData::mean
           len,              //!< Spherocylinder length [angstrom]
           half_len,         //!< Spherocylinder half length [angstrom]
           pangl,            //!< Angle of attrative patch on PSC [degrees]
           panglsw,          //!< Angle of angluar switch on sides of patch on PSC [degrees]
           pdis,             //!< Dist. to which attraction is flat on attrative patch on PSC [angstrom]
           pswitch,          //!< Distance on which attraction switches to zero on PSC [angstrom]
           chiral_angle,     //!< Angle of chirality (rotation of patch) on PSC [degrees]
           betaC,            //!< Value of the charge distribution (inverse) width [1/angstrom]
           betaD,            //!< Value of the dipole distribution (inverse) width [1/angstrom] 
           betaQ,            //!< Value of the quadrupole distribution (inverse) width [1/angstrom] 
           tfe;              //!< Transfer free energy (J/mol/angstrom^2/M)
    Point mu;                //!< Dipolemoment vector
    short int patchtype;     //!< If patchy particle, which type of patch
    bool hydrophobic;        //!< Are we hydrophobic?
    Tensor<double>
      alpha,                 //!< Polarizability
      theta;                 //!< Quadrupole moment

    bool operator==(const AtomData &d) const { return (*this==d); }

    /** @brief Constructor - by default data is initialized; mass set to unity */
    inline AtomData(const Tjson &atom=Tjson()) { readJSON(atom); }

    /** @brief Read data from a picojson object */
    inline void readJSON(const Tjson &atom) FOVERRIDE {
      name = atom.first;
      activity = json::value<double>(atom.second, "activity", 0);
      chemPot = log( activity*pc::Nav*1e-27);
      alpha << json::value<std::string>(atom.second, "alpha", "");
      alpha *= 4*pc::pi*pc::e0*(1e-10)*pc::kT()/(pc::e*pc::e);
      theta << json::value<std::string>(atom.second, "theta", "");
      theta *= 0.20819434; // Debye Å -> e Å^2
      dp = json::value<double>(atom.second, "dp", 0);
      dprot = json::value<double>(atom.second, "dprot", 0) * pc::pi / 180.; // deg->rads
      eps = json::value<double>(atom.second, "eps", 0);
      hydrophobic = json::value<bool>(atom.second, "hydrophobic", false);
      mu << json::value<std::string>(atom.second, "mu", "0 0 0");
      muscalar = mu.len()*pc::D2eA();
      if (mu.len()>1e-6)
        mu = mu/mu.len();
      mw = json::value<double>(atom.second, "Mw", 1.);
      charge = json::value<double>(atom.second, "q", 0);
      radius = json::value<double>(atom.second, "r", 0);
      sigma = 2*radius;
      sigma = json::value<double>(atom.second, "sigma", sigma);
      radius = sigma/2;
      tfe = json::value<double>(atom.second, "tfe", 0);
      half_len = 0.5 * json::value<double>(atom.second, "len", 0);
      patchtype = json::value<double>(atom.second, "patchtype", 0);
      pswitch = json::value<double>(atom.second, "patchswitch", 0);
      pdis = json::value<double>(atom.second, "patchdistance", 0);
      pangl = json::value<double>(atom.second, "patchangle", 0)/180.0*pc::pi;
      panglsw = json::value<double>(atom.second, "patchangleswitch", 0)/180.0*pc::pi;
      chiral_angle = json::value<double>(atom.second, "patchchiralangle", 0)/180.0*pc::pi;
      betaC = json::value<double>(atom.second, "betaC", pc::infty);
      betaD = json::value<double>(atom.second, "betaD", pc::infty);
      betaQ = json::value<double>(atom.second, "betaQ", pc::infty);
    }
  };

  /**
   * @brief Container for properties accessible either by index or string
   *
   * @details
   * This is a specialization of `std::vector` where the elements
   * must by derived from `PropertyBase` and can be inserted only
   * using `push_back()`. This enforces that `id` of the data
   * is always equal to the vector index. Elements can be accessed
   * by either index of string matching. Due to compatibility, a
   * default fallback property is added upon construction (index 0);
   *
   * Example:
   *
   * ~~~~
   * struct MyProp : public PropertyBase {
   *   double mass;
   * };
   *
   * PropertyVector<MyProp> v;
   *
   * MyProp d;
   * d.name="oxygen";
   * v.push_back(d);
   * v["oxygen"].mass = 15.999;
   *
   * for (auto i : v)
   *   std::cout << i.name << " " << i.mass << "\n"
   * ~~~~
   *
   * Element access complexity by `string` is linear with size so be
   * careful. For constant speed, use direct access.
   */
  template<class Tproperty, class base=std::vector<Tproperty> >
    class PropertyVector : private base {
      protected:
        string jsonfile; // keep name of json file
        string name;     // name of properties
      public:
        using typename base::value_type;
        using typename base::size_type;
        using typename base::iterator;
        using typename base::reference;
        using typename base::const_iterator;
        using typename base::const_reference;

        /** @brief Add element at the end */
        void push_back(const value_type &d) {
          base::push_back(d);
          base::back().id = base::size()-1;
        }

        iterator begin() { return base::begin(); } //!< Iterator
        iterator end() { return base::end(); }     //!< Iterator
        const_iterator begin() const { return base::begin(); } //!< Iterator
        const_iterator end() const { return base::end(); } //!< Iterator
        size_type size() const { return base::size(); } //!< Size
        bool empty() const { return base::empty(); } //!< Empty test

        /** @brief Find element by name */
        const_iterator find(const string &name) const {
          for (auto &i : *this)
            if (i.name==name)
              return begin()+i.id;
          return end();
        }

        /** @brief Access element */
        const_reference operator[](size_type i) const {
          assert( i==base::operator[](i).id && "Property out of sync");
          return base::operator[](i);
        }

        /** @brief Access element */
        reference operator[](size_type i) {
          assert( i==base::operator[](i).id && "Property out of sync");
          return base::operator[](i);
        }

        /**
         * @brief Access element by string - exception if not found. 
         * @details If not found, return default property at index 0
         */
        reference operator[](const std::string &name) {
          for (auto &i : *this)
            if (i.name==name)
              return i;
          return base::front(); // fallback
        }

        /**
         * @brief Read properties from JSON file
         * @param file Filename
         * @param section Section in JSON file to scan
         * @note All data is reset before loading
         */
        bool includefile(const string& file, const string& section) {
          jsonfile=file;
          base::resize(1); // keep default property
          auto j = json::open(file);
          for (auto &a : json::object(section, j) )
            push_back( value_type(a) );
          return ( empty() ? false : true );
        }

        PropertyVector() {
          push_back( value_type() ); // add default property
        }

        string info() {
          using namespace textio;
          char w=25;
          std::ostringstream o;
          o << header(name)
            << pad(SUB,w,"Number of entries:") << size() << endl;
          if (!jsonfile.empty())
            o << pad(SUB,w,"Input JSON file:") << jsonfile << endl;
          o << indent(SUB) << "Element info:";
          for (auto &i : *this) {
            if (i.id%10==0)
              o << endl << indent(SUBSUB);
            o << setw(SUBSUB+2) << std::left << i.name;
          }
          o << endl;
          return o.str();
        }
    };

  /**
   * @brief Class for loading and storing atomic properties
   * 
   * This will load atom properties from disk and store them in a
   * vector of `AtomData`. The file format is JSON (<http://www.json.org>)
   * and all atom properties must be inclosed in an object with
   * the keyword `atomlist`.
   * While not strictly JSON compliant, comments beginning
   * with `//` are allowed.
   *
   * For example:
   *
   *     {
   *       // All atoms reside inside the "atomlist" object.
   *       "atomlist" :
   *       {
   *         "Na" : { "q": 1.0, "r":1.9, "mw":22.99 }, // sodium ion
   *         "Cl" : { "q":-1.0, "r":1.7, "mw":35.45 }  // chloride ion
   *       }
   *     }
   *
   * Code example:
   *
   *     AtomMap a;
   *     a.includefile("atoms.json");   // load parameters
   *
   *     double s=a["Na"].sigma;        // get LJ sigma for sodium
   *
   *     PointParticle p;               // Copy properties to particle
   *     p = a["Cl"];
   *     std::cout << a[p.id].activity; // Get property via particle id
   *
   *     for (auto &i : a)              // Loop over all atoms
   *       std::cout << i.charge;
   *
   * Note that faunus currently has a global instance of `AtomMap`,
   * simply named `atom`. This can be accessed from anywhere.
   * 
   * @note
   * For simple JSON syntax highlighting in the VIM editor, add
   * the following to `~/.vimrc`:
   *
   *     au! BufRead,BufNewFile *.json set filetype=javascript
   */

  class AtomMap : public PropertyVector<AtomData> {
    public:
      typedef PropertyVector<AtomData> base;

      AtomMap() {
        base::name = "Atom Properties";
      }

      bool includefile(const std::string&); //!< Read JSON file
      bool includefile(InputMap&);     //!< Read JSON file given through `InputMap`

      /**
       * @brief Copy properties into particle vector.
       * @details Positions are left untouched!
       */
      template<typename Tpvec>
        void reset_properties(Tpvec &pvec) const {
          for (auto &i : pvec)
            i = base::at( i.id );
        }
  };

  extern AtomMap atom; //!< Global instance of AtomMap - can be accessed from anywhere


  /**
   * @brief Storage for molecule properties
   *
   * Values can be read from a JSON object with the following format:
   *
   * ~~~~
   * "salt": { "atoms": "Na,Cl", "init":"RANDOM"},
   * "polymer": {"activity": 0.05, "atoms": "MM,MM,MM,MM", "init": "POOL"}
   * ~~~~
   *
   * or more advanced like in atomData type
   *
   * The key of type string is the `name` followed, in no particular order,
   * by properties:
   *
   * Key           | Description
   * :------------ | :-------------------------------------------------------
   * `activity`    | Chemical activity for grand canonical MC [mol/l]
   * `init`        | RANDOM or POOL - using pregenerated configurations
   * `atoms`       | List of atoms in molecule (use atomsTypes names)
   *
   * NOTE: RANDOM for atomic molecules - for example salt
   *       POOL for molecular molecules - for example polymer
   */
  class MoleculeData  : public PropertyBase {
      using PropertyBase::Tjson;
    public:
      // Grand Canonical ensemble - type of initialization of insertion combinations
      enum{RANDOM,POOL};

      /** @brief Constructor - by default data is initialized; mass set to unity */
      inline MoleculeData(const Tjson &molecule=Tjson()) { readJSON(molecule); }

      std::vector<particle::Tid> atoms; ///< \brief List of atoms in molecule
      std::vector<p_vec> conformations; ///< \brief Conformations of molecule

      double activity;
      double chemPot;

      int initType;  ///< \brief sets how inserted combination will be initialized

      bool isAtomic; // set in GCMolecular, used only in GCM, here for conveniency

      bool operator==(const MoleculeData &d) const { return (*this==d); }

      /** @brief Read data from a picojson object */
      inline void readJSON(const Tjson &molecule) FOVERRIDE {
        string::size_type pos = 0;
        string::size_type oldPos = 0;
        string token;
        string line;

        name = molecule.first;
        activity = json::value<double>(molecule.second, "activity", 0);
        chemPot = log( activity*pc::Nav*1e-27);

        line = json::value<string>(molecule.second, "init", "Error");
        initType=-1;
        if(line.compare("POOL") == 0) initType = POOL;
        if(line.compare("RANDOM") == 0) initType = RANDOM;

        line = json::value<string>(molecule.second, "atoms", "Error");

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
      }
    };


  /**
   * @brief Class for loading and storing Molecule properties
   *
   * State topology "filename" in *.input file
   *
   * This will load molecule properties from disk and store them in a
   * vector of `AtomData`. The file format is JSON (<http://www.json.org>)
   * and all molecule properties must be inclosed in an object with
   * the keyword `topology`.
   * While not strictly JSON compliant, comments beginning
   * with `//` are allowed.
   *
   * For example:
   *
   * {
   *   "topology": {
   *     "salt": { "atoms": "Na,Cl", "init":"RANDOM"},
   *     "salt2": { "atoms": "Mg,Cl,Cl", "init":"RANDOM"},
   *     "chloride": { "atoms": "Cl,Cl,Cl,Cl", "init":"RANDOM"},
   *     "polymer": {"activity": 0.05, "atoms": "MM,MM,MM,MM", "init": "POOL"},
   *     "polymer2": {"activity": 0.05, "atoms": "MM,MM,MM,MM", "init": "POOL"}
   *   }
   * }
   *
   * Note that faunus currently has a global instance of `MoleculeMap`,
   * simply named `molecule`. This can be accessed from anywhere.
   *
   * @note
   * For simple JSON syntax highlighting in the VIM editor, add
   * the following to `~/.vimrc`:
   *
   *     au! BufRead,BufNewFile *.json set filetype=javascript
   */
  class MoleculeMap : public PropertyVector<MoleculeData> {
    public:
      typedef PropertyVector<MoleculeData> base;

      MoleculeMap() {
        base::name = "Molecule Properties";
      }

      bool includefile(const string&); //!< Read JSON file
      bool includefile(InputMap&);     //!< Read JSON file given through `InputMap`

      ///
      /// \return count of moleculeTypes stored
      ///
      int molTypeCount() {return size();}

      p_vec& getRandomConformation(PropertyBase::Tid molId) {
        assert(!this->operator [](molId).conformations.empty());
        return this->operator [](molId).conformations[slp_global.rand() % this->operator [](molId).conformations.size() ];
      }

      ///
      /// \brief Store a single configuration for grand canonical POOL insert
      /// \param molName of moleculeType
      /// \param vec of particles
      ///
      void pushConfiguration(string molName, p_vec& vec) {
        for(auto& mol: *this)
          if(molName.compare(mol.name) == 0) {
            pushConfiguration(mol.id, vec);
            break;
          }
      }

      ///
      /// \brief Store a single configuration for grand canonical POOL insert
      /// \param molId of moleculeType
      /// \param vec of particles
      ///
      void pushConfiguration(PropertyBase::Tid molId, p_vec& vec) {
        this->operator [](molId).conformations.push_back(vec);
      }

      /** @brief Information string */
      string info() {
        char s=14;
        using namespace textio;
        ostringstream o;

        o << header("Topology");

        o << setw(4) << "" << std::left
          << setw(s) << "Molecule" << setw(s) << "init"
          << setw(s) << "atoms" << endl;

        for (auto &i : *this) {
          o << setw(4) << "" << setw(s) << i.name;

          if(i.initType == MoleculeData::RANDOM) o << setw(s) << "RANDOM";
          else
            if(i.initType == MoleculeData::POOL) o << setw(s) << "POOL";
            else o << setw(s) << "";

          for ( auto j=i.atoms.begin(); j!=i.atoms.end(); ++j ) {
            o << setw(0) << atom[(*j)].name;

            if(j != i.atoms.end()-1)
              o<<",";
          }
          o<<"\n";
        }
        return o.str();
      }
  };

  extern MoleculeMap molecule;
}//namespace
#endif
