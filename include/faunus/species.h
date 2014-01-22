#ifndef FAUNUS_SPECIES_H
#define FAUNUS_SPECIES_H

#ifndef SWIG
#include <faunus/common.h>
#include <faunus/point.h>
#endif

namespace Faunus {

  /**
   * @brief Atomic properties
   */
  struct AtomData {
    typedef unsigned char Tid;
    typedef bool Thydrophobic;
    AtomData();
    Tid id;            //!< Identification number
    double sigma,      //!< LJ diameter [angstrom]
           eps,        //!< LJ epsilon [kJ/mol] (pair potentials should convert to kT)
           radius,     //!< Radius [angstrom]
           muscalar,   //!< Dipole momentscalar [eÅ]
           mw,         //!< Weight [g/mol]
           charge,     //!< Charge/valency [e]
           activity,   //!< Chemical activity [mol/l]
           dp,         //!< Translational displacement parameter [angstrom]
           dprot,      //!< Rotational displacement parameter [degrees]
           mean,       //!< Mean value... (charge, sasa, etc.)
           variance,   //!< Spread around AtomData::mean
      
           len,        //!< Spherocylinder length [Å]
           half_len,   //!< Spherocylinder half length [Å]
           pangl,      //!< Angle of attrative patch on PatchySpherocylinder [degrees]
           panglsw,    //!< Angle of angluar switch on sides of patch on PatchySpherocylinder [degrees]
           pdis,       //!< Distance to which attraction is flat (typicaly end of repulsion) on attrative patch on PatchySpherocylinder [Å]
           pswitch,    //!< Distance on which attraction switches to zero on PatchySpherocylinder [Å]
           chiral_angle,//!< Angle of chirality (rotation of patch) on PatchySpherocylinder [degrees]
           betaC,      //!< Value of the charge distribution (inverse) width [Å^-1]
           betaD,      //!< Value of the dipole distribution (inverse) width [Å^-1] 
           betaQ;      //!< Value of the quadrupole distribution (inverse) width [Å^-1] 
    Point mu;
    short int patchtype;     //!< If patchy particle, which type of patch
    Thydrophobic hydrophobic;//!< Are we hydrophobic?
    Tensor<double>
      alpha,           //!< Polarizability
      theta;           //!< Quadrupole moment
    string name;       //!< Name. Avoid spaces.
    bool operator==(const AtomData &d) const { return (*this==d); }
  };

  /**
   * @brief Container for data between pairs
   *
   * This will maintain a symmetric, dynamic NxN matrix for storing data
   * about pairs.
   * Use the `set()` function for setting values and use the function
   * operator for access: 
   *
   *     int i=2,j=3; // particle type, for example
   *     PairMatrix<double> cutoff;
   *     eps2.set(i,j,12.0);
   *     cout << cutoff(i,j); // -> 12.0
   */
  template<class T=double>
    class PairMatrix {
      public:
        vector< vector<T> > m; // symmetric matrix (mem.wasteful - fast access)
        void resize(size_t n) {
          m.resize(n);
          for (auto &i : m)
            i.resize(n,0);
        }
        PairMatrix(size_t n=0) {
          resize(n);
        }
        const T& operator()(size_t i, size_t j) const {
          assert( i<m.size() );
          assert( j<m[i].size() );
          assert( m[i][j]==m[j][i] );
          return m[i][j]; 
        }
        void set(size_t i, size_t j, T val) {
          size_t n=std::max(i,j);
          if (n>=m.size())
            resize(n+1);
          m[i][j]=m[j][i]=val;
        }
    };

  /*!
   * @brief Class for loading and storing atomic properties
   * 
   * This will load atom properties from disk and store them in a
   * map of AtomData. The recommended file format is JSON
   * (<http://www.json.org>)
   * and all atom properties must be inclosed in an object with
   * the keyword "atomlist".
   * While not strictly JSON compliant, comments beginning
   * with `//` are allowed.
   * For example:
   *
   *     {
   *       // All atoms must be inside the "atomlist" object.
   *       "atomlist" :
   *       {
   *         "Na" : { "q": 1.0, "r":1.9, "mw":22.99 }, // sodium ion
   *         "Cl" : { "q":-1.0, "r":1.7, "mw":35.45 }, // chloride ion
   *         "Cigar" :
   *         {
   *           "sigma" : 3.1,
   *           "hydrophobic" : true,
   *           "patchtype" : 1
   *         }
   *       }
   *     }
   *
   * Or: 
   * 
   *     "atomlist" :
   *     {
   *       "sol1" : { "q": 1.0, "r":1.9, "mw":22.99, "mu":"1.3 0.1 0",
   *          "alpha":"1.1   0 0 2.3 0   1",    // Matrix input in 
   *          "theta":"2.4 0.3 0 1.8 0 3.2" }   // (xx, xy, xz, yy, yz, zz) form
   *     }
   * 
   * @note
   * For simple JSON syntax highlighting in the VIM editor, add
   * the following to `.vimrc`:
   *
   *     au! BufRead,BufNewFile *.json set filetype=javascript
   *
   * Loading the input file, the following keywords are read and
   * if not present, a default value of zero, false or unity (mw) is assumed.
   *
   * Key           | Description
   * :------------ | :-------------------------------------------------------
   * `activity`    | Chemical activity for grand canonical MC [mol/l]
   * `alpha`       | Polarizability [\f$ 4\pi\epsilon_0 \f$ cubic angstrom]
   * `dp`          | Translational displacement parameter [angstrom] 
   * `dprot`       | Rotational displacement parameter [degrees] (will be converted to radians)
   * `eps`         | Epsilon energy scaling commonly used for Lennard-Jones interactions etc. [kJ/mol] 
   * `hydrophobic` | Is the particle hydrophobic? [true/false]
   * `mu`          | Dipole moment vector [Debye]
   * `theta`       | Quadrupole moment matrix [Debye Å]
   * `mw`          | Molecular weight [g/mol]
   * `patchtype`   | Patchtype for sphero-cylinders
   * `q`           | Valency / partial charge number [e]
   * `r`           | Radius = `sigma/2` [angstrom]
   * `sigma`       | `2*r` [angstrom] (overrides radius)
   * 
   *
   * Code example:
   *
   *     AtomMap a;
   *     a.include("atoms.json");       // load parameters
   *     double s=a["Na"].sigma;        // get LJ sigma for sodium
   *     particle p = a["Cl"];          // Copy properties to particle
   *     std::cout << a[p.id].activity; // Get property via particle id
   *
   * Note that faunus currently has a global instance of AtomMap,
   * simply named `atom`. This can be accessed from anywhere.
   *
   * If a non-JSON file is parsed, the loader attempts to use a file format
   * from older versions of Faunus - see below. This is discouraged, though,
   * as this is a fixed format that cannot be expanded.
   *
   *     #     name    charge  radius     epsilon      Mw       hydrophobic?
   *     #             (e)     (angstrom) (kJ/mol)     (g/mol)  (yes/no)
   *     Atom  Na      +1      1.665      0.01158968   22       no
   *     Atom  Cl      -1      2.200      0.4184       36       no
   *
   */
  class AtomMap {
    private:
      string filename;
      bool includeJSON(const string&);     //!< Load JSON file
    public:
      AtomMap();                           //!< Constructor
      std::vector<AtomData> list;          //!< List of atoms
      bool includefile(string);            //!< Append atom parameters from file
      bool includefile(InputMap&);         //!< Append atom parameters from file
      int size();                          //!< Number of atom-types
      AtomData& operator[] (string);       //!< Name->data
      AtomData& operator[] (AtomData::Tid);//!< Id->data
      string info();                       //!< Print info

      /**
       * @brief Copy properties into particles vector.
       *
       * Positions are left untouched!
       */
      template<typename Tpvec>
        void reset_properties(Tpvec &pvec) const {
          for (auto &i : pvec)
            i = list.at( i.id );
        }
  };

extern AtomMap atom; //!< Global instance of AtomMap - can be accessed from anywhere
}//namespace
#endif
