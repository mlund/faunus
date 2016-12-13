#ifndef faunus_molecule_h
#define faunus_molecule_h

#include <faunus/io.h>
#include <faunus/inputfile.h>
#include <faunus/bonded.h>
#include <faunus/geometry.h>

namespace Faunus
{

  /**
   * @brief Random position and orientation - typical for rigid bodies
   *
   * Molecule inserters take care of generating molecules
   * for insertion into space and can be used in Grand Canonical moves,
   * Widom analysis, and for generating initial configurations.
   * Inserters will not actually insert anything, but rather
   * return a particle vector with proposed coordinates.
   *
   * All inserters are function objects, expecting
   * a geometry, particle vector, and molecule data.
   */
  template<typename TMoleculeData>
  struct RandomInserter
  {
      typedef typename TMoleculeData::TParticleVector Tpvec;
      string name;
      Point dir;         //!< Scalars for random mass center position. Default (1,1,1)
      Point offset;      //!< Added to random position. Default (0,0,0)
      bool checkOverlap; //!< Set to true to enable container overlap check
      bool rotate;       //!< Set to true to randomly rotate molecule when inserted. Default: true
      bool keeppos;      //!< Set to true to keep original positions (default: false)
      int maxtrials;     //!< Maximum number of overlap checks if `checkOverlap==true`

      RandomInserter() : dir(1, 1, 1), offset(0, 0, 0), checkOverlap(true), keeppos(false), maxtrials(2e3)
      {
          name = "random";
      }

      Tpvec operator()( Geometry::Geometrybase &geo, const Tpvec &p, TMoleculeData &mol )
      {
          bool _overlap = true;
          Tpvec v;
          int cnt = 0;
          do
          {
              if ( cnt++ > maxtrials )
                  throw std::runtime_error("Max. # of overlap checks reached upon insertion.");

              v = mol.getRandomConformation();

              if ( mol.isAtomic())
              { // insert atomic species
                  for ( auto &i : v )
                  { // for each atom type id
                      Geometry::QuaternionRotate rot;
                      Point u;
                      u.ranunit(slump);
                      if ( rotate )
                      {
                          rot.setAxis(geo, {0, 0, 0}, u, 2 * pc::pi * slump());
                          i.rotate(rot);
                      }
                      geo.randompos(i);
                      i = i.cwiseProduct(dir) + offset;
                      geo.boundary(i);
                  }
              }
              else
              { // insert molecule
                  if ( keeppos )
                  {                   // keep original positions (no rotation/trans)
                      for ( auto &i : v )              // ...but let's make sure it fits
                          if ( geo.collision(i, 0))
                              throw std::runtime_error("Error: Inserted molecule does not fit in container");
                  }
                  else
                  {                         // generate a now position/orientation
                      Point a, b;
                      geo.randompos(a);              // random point in container
                      a = a.cwiseProduct(dir);       // apply user defined directions (default: 1,1,1)
                      Geometry::cm2origo(geo, v);    // translate to origo - obey boundary conditions
                      Geometry::QuaternionRotate rot;
                      b.ranunit(slump);              // random unit vector
                      rot.setAxis(geo, {0, 0, 0}, b, slump() * 2 * pc::pi); // random rot around random vector
                      for ( auto &i : v )
                      {            // apply rotation to all points
                          if ( rotate )
                              i = rot(i) + a + offset;   // ...and translate
                          else
                              i += a + offset;
                          geo.boundary(i);             // ...and obey boundaries
                      }
                  }
              }

              assert(!v.empty());

              _overlap = false;
              if ( checkOverlap )              // check for container overlap
                  for ( auto &i : v )
                      if ( geo.collision(i, i.radius))
                      {
                          _overlap = true;
                          break;
                      }
          }
          while ( _overlap == true );
          return v;
      }
  };

  /**
   * @brief Weight molecule according to deviation from mean charge
   *
   * This is a weight function to be used with Grand Canonical insertion
   * by associating molecular configurations with a weight defined by how
   * much it deviates from the average, bulk charge. It is required that
   * both `MoleculeData::capacitance` and `MoleculeData::meancharge` have
   * been specified.
   * 
   * [Details: J. Am. Chem. Soc. 2010, 132:17337](http://dx.doi.org/10.1021/ja106480a)
   */
  template<class Tpvec, class TMoleculeData>
  double WeightByCharge(
      Geometry::Geometrybase &geo, const Tpvec &p,
      const Group &g, const TMoleculeData &m )
  {
      double dz = netCharge(p, g) - m.meancharge;
      if ( m.capacitance > 1e-9 )
          return exp(-0.5 * dz * dz / m.capacitance);
      return 1;
  }

  /**
   * @brief Storage for molecular properties
   *
   * Values can be read from a JSON object with the following format:
   *
   * ~~~~
   *    "mymolecule" : { "atoms":"ARG ARG GLU", "structure":"peptide.aam" }
   * ~~~~
   *
   * The key of type string is the molecule `name` followed, in no particular order,
   * by properties:
   *
   * Key           | Type    | Description
   * :------------ | :------ | : ------------------------------------------------
   * `activity`    | float   | Chemical activity for grand canonical MC [mol/l]
   * `atoms`       | string  | List of atoms in molecule (use `AtomData` names)
   * `atomic`      | bool    | `true` if molecule a free atomic species (default: false)
   * `bonds`       | string  | List of harmonic bonds - index 0 corresponds to first atom in structure
   * `dihedrals`   | string  | List of dihedrals (under construction!)
   * `fasta`       | string  | Construct bonded chain from fasta sequence (hardcoded k and req)
   * `insdir`      | string  | Directions for generation of random position. Default: "1 1 1" = XYZ
   * `insoffset`   | string  | Translate generated random position. Default: "0 0 0" = no translation
   * `keeppos`     | bool    | Keep original positions (`insdir`, `insoffset` ignored. Default: `false`)
   * `Ninit`       | int     | Initial number of molecules
   * `checkoverlap`| bool    | Check for overlap while inserting. Default: true
   * `rotate`      | bool    | Randomly rotate molecule or anisotropic atom upon insertion. Default: true
   * `structure`   | string  | Read conformation from AAM file
   * `traj`        | string  | Read conformations from PQR trajectory (`structure` will be ignored)
   * `trajweight`  | string  | One column file w. relative weights for each conformations. Must match `traj` file.
   * `Cavg`        | float   | Charge capacitance [e^2]
   * `Zavg`        | float   | Average charge number [e]
   */
  template<class Tpvec>
  class MoleculeData : public PropertyBase
  {
  public:
      typedef Tpvec TParticleVector;
      typedef typename Tpvec::value_type Tparticle;

      /** @brief Signature for inserted function */
      typedef std::function<Tpvec( Geometry::Geometrybase &,
                                   const Tpvec &, MoleculeData<Tpvec> & )> TinserterFunc;

      std::vector<typename Tparticle::Tid> atoms; //!< List of atoms in molecule
      std::vector<Tpvec> conformations;           //!< Conformations of molecule
      std::discrete_distribution<> confDist;      //!< Weight of conformations
      double activity;                            //!< Chemical activity (mol/l)
      double chemPot;                             //!< Chemical potential (kT)
      Average<double> rho;                        //!< Average concentration (1/angstrom^3)
      vector<Bonded::BondData> bonds;             //!< List of harmonic bonds
      vector<Bonded::DihedralData> dihedrals;     //!< List of harmonic bonds
      int Ninit;                                  //!< Initial number of molecules

      double capacitance;                         //!< Charge capacitance
      double meancharge;                          //!< Mean charge

  private:
      bool _isAtomic;

  public:
      TinserterFunc inserterFunctor;              //!< Function for insertion into space

      /** @brief Constructor - by default data is initialized; mass set to unity */
      inline MoleculeData( Tmjson::iterator &molecule )
          : Ninit(0), _isAtomic(false)
      {
          readJSON(molecule);
          auto ins = RandomInserter<MoleculeData<Tpvec> >();
          ins.dir << (molecule.value()["insdir"] | string("1 1 1"));
          ins.offset << (molecule.value()["insoffset"] | string("0 0 0"));
          ins.checkOverlap = molecule.value()["checkoverlap"] | true;
          ins.rotate = molecule.value()["rotate"] | true;
          ins.keeppos = molecule.value()["keeppos"] | false;
          setInserter(ins);
      }

      /** @brief Get list of bonds for molecule */
      std::vector<Bonded::BondData> &getBondList() { return bonds; }

      /**
       * @brief Specify function to be used when inserting into space.
       *
       * By default a random position and orientation is generator and overlap
       * with container is avoided.
       */
      void setInserter( const TinserterFunc &ifunc ) { inserterFunctor = ifunc; };

      /**
       * @brief Get a random conformation
       *
       * This will return the raw coordinates of a random conformation
       * as loaded from a directory file. The propability of a certain
       * conformation is dictated by the weight which, by default,
       * is set to unity. Specify a custom distribution using the
       * `weight` keyword.
       */
      Tpvec getRandomConformation()
      {
          if ( conformations.empty())
              throw std::runtime_error("No configurations for molecule '" + name +
                  "'. Perhaps you forgot to specity the 'atomic' keyword?");

          assert(size_t(confDist.max()) == conformations.size() - 1);
          assert(atoms.size() == conformations.front().size());

          return conformations[confDist(slump.eng)];
      }

      /**
       * @brief Get random conformation that fits in container
       * @param geo Geometry
       * @param otherparticles Typically `spc.p` is insertion depends on other particle
       *
       * By default the molecule is placed at a random position and orientation with
       * no container overlap using the `RandomInserter` class. This behavior can
       * be changed by specifying another inserter using `setInserter()`.
       */
      Tpvec getRandomConformation(
          Geometry::Geometrybase &geo,
          const Tpvec &otherparticles = Tpvec())
      {
          return inserterFunctor(geo, otherparticles, *this);
      }

      /**
       * @brief Store a single conformation
       * @param vec Vector of particles
       * @param weight Relative weight of conformation (default: 1)
       */
      void pushConformation( const Tpvec &vec, double weight = 1 )
      {
          if ( !conformations.empty())
          {     // resize weights
              auto w = confDist.probabilities();// (defaults to 1)
              w.push_back(weight);
              confDist = std::discrete_distribution<>(w.begin(), w.end());
          }
          conformations.push_back(vec);
          assert(confDist.probabilities().size() == conformations.size());
      }

      /** True if molecule holds individual atoms */
      bool isAtomic() const { return _isAtomic; }

      /** True if molecule holds molecular entities */
      bool isMolecular() const { return !_isAtomic; }

      bool operator==( const MoleculeData &d ) const { return (*this == d); }

      /** Read data from json object **/
      void readJSON( Tmjson::iterator &molecule )
      {

          name = molecule.key();
          auto _js = molecule.value();

          _isAtomic = _js["atomic"] | false;
          Ninit = _js["Ninit"] | 0;
          capacitance = _js["Cavg"] | 0.0;
          meancharge = _js["Zavg"] | 0.0;
          bool keeppos = _js["keeppos"] | false;
          activity = _js["activity"] | 0.0;
          chemPot = log(activity * 1.0_molar);

          // create bond list
          auto b = _js["bonds"];
          for ( auto it = b.begin(); it != b.end(); ++it )
              bonds.push_back(Bonded::BondData(it));

          // create dihedral list
          b = _js["dihedrals"];
          for ( auto it = b.begin(); it != b.end(); ++it )
              dihedrals.push_back(Bonded::DihedralData(it));

          // read conformation from disk
          {
              string structure = _js["structure"] | string();
              if ( !structure.empty())
              {
                  Tpvec v;
                  if ( structure.substr(structure.find_last_of(".") + 1) == "aam" )
                      FormatAAM::load(structure, v);
                  if ( structure.substr(structure.find_last_of(".") + 1) == "pqr" )
                      FormatPQR::load(structure, v);
                  if ( structure.substr(structure.find_last_of(".") + 1) == "xyz" )
                      FormatXYZ::load(structure, v);
                  if ( !v.empty())
                  {
                      if ( keeppos == false )
                          Geometry::cm2origo(
                              Geometry::Sphere(1e20), v); // move to origo
                      pushConformation(v);        // add conformation
                      for ( auto &p : v )           // add atoms to atomlist
                          atoms.push_back(p.id);
                  }
                  if ( v.empty())
                      throw std::runtime_error("Structure " + structure + " not loaded. Filetype must be .aam/.pqr");
              }
          }

          // construct flexible peptide from fasta sequence
          string fasta = _js["fasta"] | string();
          if ( !fasta.empty())
          {
              assert(atoms.empty());
              atoms = FastaToAtoms(fasta); // -> atomid
              atoms.insert(atoms.begin(), atom["NTR"].id);
              atoms.push_back(atom["CTR"].id);
              double req = 4.9, k = 0.76;
              Tpvec v;
              Point u;
              typename Tpvec::value_type a;
              v.reserve(atoms.size());
              for ( auto i : atoms )
              {
                  a = atom[i];
                  if ( v.empty())
                      a = Point(0, 0, 0);
                  else
                  {
                      u.ranunit(slump);   // position randomly
                      a = v.back() + req * u; // around previous particle
                  }
                  v.push_back(a);
              }
              pushConformation(v);
              bonds.reserve(v.size() - 1);
              for ( size_t i = 0; i < v.size() - 1; i++ )
                  bonds.push_back(
                      Bonded::BondData(i, i + 1, k, req, Bonded::BondData::Type::HARMONIC));
          }

          // read tracjectory w. conformations from disk
          string traj = _js["traj"] | string();
          if ( !traj.empty())
          {
              conformations.clear();
              FormatPQR::load(traj, conformations);
              if ( !conformations.empty())
              {
                  // create atom list
                  atoms.clear();
                  atoms.reserve(conformations.front().size());
                  for ( auto &p : conformations.front())           // add atoms to atomlist
                      atoms.push_back(p.id);

                  // set default uniform weight
                  vector<float> w(conformations.size(), 1);
                  confDist = std::discrete_distribution<>(w.begin(), w.end());

                  // look for weight file
                  string weightfile = _js["trajweight"] | string();
                  if ( !weightfile.empty())
                  {
                      std::ifstream f(weightfile.c_str());
                      if ( f )
                      {
                          w.clear();
                          w.reserve(conformations.size());
                          double _val;
                          while ( f >> _val )
                              w.push_back(_val);
                          if ( w.size() == conformations.size())
                              confDist = std::discrete_distribution<>(w.begin(), w.end());
                          else
                              throw std::runtime_error("Number of weights does not match conformations.");
                      }
                      else
                          throw std::runtime_error("Weight file " + weightfile + " not found.");
                  }
              }
              else
                  throw std::runtime_error("Error: trajectory " + traj + " not loaded or empty.");
          }

          // add atoms to atom list
          if ( atoms.empty())
          {
              string atomlist = _js["atoms"] | string();
              for ( auto &a : textio::words2vec<string>(atomlist))
                  if ( atom[a].id > 0 )
                      atoms.push_back(atom[a].id);
          }

          // for atomic groups, add particles
          if ( isAtomic())
          {
              Tparticle _p;
              Tpvec _vec;
              _vec.reserve(atoms.size());
              for ( auto id : atoms )
              {
                  _p = atom[id];
                  _vec.push_back(_p);
              }
              pushConformation(_vec);
          }
      }

      inline string info()
      {

          using namespace textio;
          ostringstream o;
          char w = 30;

          o << header("Molecule: " + name)
            << pad(SUB, w, "Activity") << activity << " mol/l\n"
            << pad(SUB, w, "Chem. pot.") << chemPot << " kT\n"
            << pad(SUB, w, "Atomic") << std::boolalpha << isAtomic() << "\n"
            << pad(SUB, w, "Number of atoms") << std::max(atoms.size(), conformations.size()) << "\n"
            << pad(SUB, w, "Number of conformations") << conformations.size() << "\n";
          if ( Ninit > 0 )
              o << pad(SUB, w, "N initial") << Ninit << "\n";
          if ( confDist.probabilities().size() > 1 )
          {
              o << pad(SUB, w, "Number of weigts") << confDist.probabilities().size() << "\n";
              for ( auto w : confDist.probabilities())
                  o << w << " ";
              o << "\n";
          }

          if ( rho.cnt > 0 )
          {
              o << pad(SUB, w, "Average conc.") << rho.avg() << " 1/" + angstrom + " = "
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
   *
   * For example:
   *
   * ~~~~
   * {
   *   "moleculelist" :
   *   {
   *     "salt" : { "atoms":"Na Cl", "atomic":True },
   *     "polymer" :
   *        { "activity":0.05, "atoms":"MM MM MM MM",
   *         "bonds" :
   *         {
   *           "0 1" : { "k":0.5, "req":4.0 },
   *           "1 2" : { "k":0.5, "req":4.0 },
   *           "2 3" : { "k":0.5, "req":4.0 }
   *         },
   *         "dihedrals" :
   *         {
   *           "0 1 2 3" : { "k1":0.01, "k2":2.123 }
   *         }
   *       }
   *   }
   * }
   * ~~~~
   *
   * Note that faunus currently has a global instance of `MoleculeMap`,
   * simply named `molecule`. This can be accessed from anywhere.
   *
   * @todo More documentation. Rename entry for filename to "moleculefile"
   */
  template<class Tpvec, class base=PropertyVector<MoleculeData<Tpvec>>>
  class MoleculeMap : public base
  {
  public:

      MoleculeMap()
      {
          base::jsonsection = "moleculelist";
          base::name = "Molecule Properties";
      }

      /** @brief Information string */
      string info()
      {
          char w = 15;
          using namespace textio;
          ostringstream o;
          o << header("Molecule list") << "\n"
            << std::left << setw(w + 4) << "  Molecule" << setw(w) << "activity"
            << setw(w) << "atoms" << setw(w) << "atomic" << "\n"
            << "  " << string(4 * w, '-') << "\n";

          for ( auto &i : *this )
              if ( !i.name.empty())
                  o << setw(w + 4) << "  " + i.name << setw(w) << i.activity
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
  class MoleculeCombination : public PropertyBase
  {
  public:

      MoleculeCombination( Tmjson::iterator &comb )
      {
          name = comb.key();
          probability = comb.value()["prob"] | 1.0;
          string mollist = comb.value()["molecules"] | string();
          molName = textio::words2vec<string>(mollist);
      }

      vector<PropertyBase::Tid> molComb;/// vector of molecule index in combination

      vector<string> molName;           /// vector of molecule names in combination

      double probability;               /// probability of this combination in GC-MC move

  };

  /**
   * @brief Vector of molecular combinations
   *
   * When examining a JSON file, individual entries must be placed
   * in a section called `moleculecombinations`.
   *
   * @todo Add warning if `moleculecombinations` is not found
   */
  template<class Tpvec, class base=PropertyVector<MoleculeCombination<Tpvec> >>
  class MoleculeCombinationMap : public base
  {
  private:

      MoleculeMap<Tpvec> *molmap; // typically points to `Space::molecule`

      void update()
      { // convert name (string) -> id (int) using MoleculeMap
          for ( auto &i : *this )
          {
              assert(!i.molName.empty() && "Molecule combination cannot be empty");
              for ( auto &name : i.molName )
              {
                  auto it = molmap->find(name);
                  if ( it != molmap->end())
                      i.molComb.push_back(it->id);
              }
              assert(i.molName.size() == i.molComb.size());
              assert(!i.molName.empty());
          }
      }

  public:

      MoleculeCombinationMap( MoleculeMap<Tpvec> &moleculemap ) : molmap(&moleculemap)
      {
          assert(molmap != nullptr);
          base::name = "Molecule Combinations";
          base::jsonsection = "moleculecombinations";
      }

      bool include( Tmjson &j ) override
      {
          bool r = base::include(j);
          update();
          return r;
      }

      string info()
      {
          using namespace textio;
          char w = 25;
          std::ostringstream o;
          o << header(base::name)
            << pad(SUB, w, "Number of entries:") << base::size() << endl;
          o << indent(SUB) << "Element info:\n";
          for ( auto &i : *this )
          {
              o << std::left << setw(w) << "    " + i.name + ":";
              for ( auto &n : i.molName )
                  o << setw(10) << n;
              o << "\n";
          }
          o << endl;
          return o.str();
      }

  };

}//namespace

#endif
