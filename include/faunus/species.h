#ifndef FAUNUS_SPECIES_H
#define FAUNUS_SPECIES_H

#ifndef SWIG
#include <faunus/common.h>
#include <faunus/json.h>
#include <faunus/json.hpp>
#include <faunus/physconst.h>
#include <faunus/textio.h>
#include <faunus/point.h>
#include <faunus/slump.h>
#include <faunus/average.h>

#endif

namespace Faunus
{

  /**
   * @brief Base class for property data
   *
   * Used to collect properties for atoms, molecules etc.
   */
  class PropertyBase
  {
  public:
      typedef unsigned char Tid;
      std::string name;  //!< Unique, user-defined name
      Tid id;            //!< Unique id (automatically set by `PropertyVector`)
      virtual ~PropertyBase() {};

      PropertyBase() : name("UNK") {}
  };

  /**
   * @brief Container for properties accessible either by index or string
   *
   * @details
   * This is a specialization of `std::vector` where the elements
   * must by derived from `PropertyBase` and can be inserted only
   * using `push_back()`. This enforces that `id` of the data
   * is always equal to the vector index. Elements can be accessed
   * by either index of string matching.
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
  class PropertyVector : private base
  {
  protected:
      string jsonsection; // section to look for elements
      string name;        // name of properties
  public:
      using typename base::value_type;
      using typename base::size_type;
      using typename base::iterator;
      using typename base::reference;
      using typename base::const_iterator;
      using typename base::const_reference;

      /** @brief Iterator to random element */
      const_iterator random() const { return slump.element(begin(), end()); }

      /** @brief Iterator to random element */
      iterator random() { return slump.element(begin(), end()); }

      /** @brief Add element at the end */
      void push_back( const value_type &d )
      {
          base::push_back(d);
          base::back().id = base::size() - 1;
      }

      iterator begin() { return base::begin(); } //!< Iterator
      iterator end() { return base::end(); }     //!< Iterator
      const_iterator begin() const { return base::begin(); } //!< Iterator
      const_iterator end() const { return base::end(); } //!< Iterator
      size_type size() const { return base::size(); } //!< Size
      bool empty() const { return base::empty(); } //!< Empty test

      /** @brief Find element by name */
      const_iterator find( const string &name ) const
      {
          return std::find_if(begin(), end(), [&name]( const value_type &i ) { return i.name == name; });
      }

      /** @brief Find element by name */
      iterator find( const string &name )
      {
          return std::find_if(begin(), end(), [&name]( const value_type &i ) { return i.name == name; });
      }

      /** @brief Access element by string */
      const_reference operator[]( const std::string &name ) const { return *find(name); }

      /** @brief Access element */
      const_reference operator[]( size_type i ) const
      {
          assert(i == base::operator[](i).id && "Property out of sync");
#ifdef NDEBUG
          return base::operator[](i);
#else
          return base::at(i);
#endif
      }

      /** @brief Access element */
      reference operator[]( size_type i )
      {
          assert(i == base::operator[](i).id && "Property out of sync");
          assert(i < base::size());
#ifdef NDEBUG
          return base::operator[](i);
#else
          return base::at(i);
#endif
      }

      /** Load data from json object */
      virtual bool include( Tmjson &j )
      {
          assert(!jsonsection.empty());
          auto m = j[jsonsection];
          base::reserve(m.size());
          for ( auto it = m.begin(); it != m.end(); ++it )
              push_back(value_type(it));
          return (empty() ? false : true);
      }

      PropertyVector()
      {
          static_assert(std::is_base_of<PropertyBase, Tproperty>::value,
                        "Elements must be derived from `PropertyBase`");
      }

      virtual ~PropertyVector() {};

      string info()
      {
          using namespace textio;
          char w = 25;
          std::ostringstream o;
          o << header(name)
            << pad(SUB, w, "Number of entries:") << size() << endl;
          o << indent(SUB) << "Element info:";
          for ( auto &i : *this )
          {
              if ( i.id % 10 == 0 )
                  o << endl << indent(SUBSUB);
              o << setw(SUBSUB + 2) << std::left << i.name;
          }
          o << endl;
          return o.str();
      }
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
   * `Ninit`       | Initial number of atoms (used by `MoleculeData` to insert atoms
   * `theta`       | Quadrupole moment tensor [Debye \f$ \unicode{x212B} \f$]
   * `mw`          | Molecular weight [g/mol]
   * `patchtype`   | Patchtype for sphero-cylinders
   * `q`           | Valency / partial charge number [e]
   * `r`           | Radius = `sigma/2` [angstrom]
   * `sigma`       | `2*r` [angstrom] (overrides radius)
   * `tfe`         | Transfer free energy [J/mol/angstrom^2/M] (default: 0)
   * `alphax`      | Excess polarizability in units of [angstrom^3]
   */

  class AtomData : public PropertyBase
  {
  public:
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
          tfe,              //!< Transfer free energy (J/mol/angstrom^2/M)
          alphax;           //!< Excess polarizability [angstrom^3]
      Point mu;                //!< Dipolemoment vector
      int Ninit;               //!<
      short int patchtype;     //!< If patchy particle, which type of patch
      bool hydrophobic;        //!< Are we hydrophobic?
      Tensor<double>
          alpha,                 //!< Polarizability
          theta;                 //!< Quadrupole moment

      bool operator==( const AtomData &d ) const { return (*this == d); }

      /** @brief Constructor - by default data is initialized; mass set to unity */
      inline AtomData( Tmjson::iterator &atom )
      {

          auto _js = atom.value();
          if ( !atom.key().empty())
              name = atom.key();

          activity = _js["activity"] | 0.0;
          chemPot = log(activity * 1.0_molar);
          alpha << (_js["alpha"] | string());
          alpha /= pc::lB(1.0);
          theta << (_js["theta"] | string());
          theta *= 1.0_Debye;
          dp = _js["dp"] | 0.0;
          dprot = (_js["dprot"] | 0.0) * 1._deg; // deg->rads
          eps = (_js["eps"] | 0.0) * 1.0_kJmol;
          hydrophobic = _js["hydrophobic"] | false;
          mu << (_js["mu"] | string("0 0 0"));
          muscalar = mu.len() * 1.0_Debye;
          if ( mu.len() > 1e-6 )
              mu = mu / mu.len();

          mw = _js["mw"] | 1.0;
          Ninit = _js["Ninit"] | 0.0;
          charge = _js["q"] | 0.0;
          radius = (_js["r"] | 0.0) * 1.0_angstrom;
          sigma = (_js["sigma"] | 2 * radius) * 1.0_angstrom;
          radius = 0.5 * sigma;
          tfe = _js["tfe"] | 0.0;
          alphax = (_js["alphax"] | 0.0) * std::pow(radius, 3);
          string unit = _js["alphax_unit"] | string("unitless");
          if ( unit == "angstrom^3" )
              alphax /= std::pow(radius, 3);

          // spherocylindrical properties
          half_len = 0.5 * (_js["len"] | 0.0);
          patchtype = _js["patchtype"] | 0.0;
          pswitch = _js["patchswitch"] | 0.0;
          pdis = _js["patchdistance"] | 0.0;
          pangl = (_js["patchangle"] | 0.0) * 1._deg;
          panglsw = (_js["patchangleswitch"] | 0.0) * 1._deg;
          chiral_angle = (_js["patchchiralangle"] | 0.0) * 1._deg;

          betaC = _js["betaC"] | pc::infty;
          betaD = _js["betaD"] | pc::infty;
          betaQ = _js["betaQ"] | pc::infty;
      }
  };

  /**
   * @brief Class for loading and storing atomic properties
   * 
   * This will load atom properties from disk and store them in a
   * vector of `AtomData`. The file format is JSON (<http://www.json.org>)
   * and all atom properties must be inclosed in an object with
   * the keyword `atomlist`.
   * Due to compatibility, a default fallback property is added upon construction
   * (index 0).
   *
   * For example:
   *
   *     {
   *       "atomlist" :
   *       {
   *         "Na" : { "q": 1.0, "r":1.9, "mw":22.99 },
   *         "Cl" : { "q":-1.0, "r":1.7, "mw":35.45 }
   *       }
   *     }
   *
   * Code example:
   *
   *     AtomMap a;
   *     a.include("atoms.json");   // load parameters
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

  class AtomMap : public PropertyVector<AtomData>
  {
  public:
      typedef PropertyVector<AtomData> base;

      AtomMap()
      {
          base::name = "Atom Properties";
          base::jsonsection = "atomlist";

          // to be removed! This is here only to be compatible with test
          // data that expects an initial dummy atom with id 0.
          Tmjson j;
          j["unk"] = {{"dummy", 0}};
          auto it = j.begin();
          push_back(it); // add default property
      }

  };

  extern AtomMap atom; //!< Global instance of AtomMap - can be accessed from anywhere

}//namespace
#endif
