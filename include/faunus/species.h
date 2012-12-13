#ifndef FAU_species_h
#define FAU_species_h

#ifndef SWIG
#include <faunus/common.h>
#include <faunus/point.h>
#endif

namespace Faunus {
  //class _inputfile;

  class AtomData {
    public:
      AtomData();
      particle::Tid id;  //!< Identification number
      double sigma,      //!< LJ diameter [A]
             eps,        //!< LJ epsilon [kJ/mol]
             radius,     //!< Radius [A]
             mw,         //!< Weight [g/mol]
             charge,     //!< Charge [valency]
             activity,   //!< Chemical activity "(mol/l)"
             dp,         //!< Translational displacement parameter
             dprot,      //!< Rotational displacement parameter
             mean,       //!< Mean value... (charge, sasa, etc.)
             variance;   //!< ...and the spread around it.
      short int patchtype;  //!< Type of patch
      particle::Thydrophobic hydrophobic;  //!< Are we hydrophobic?
      string name;       //!< Name. Avoid spaces.
      bool operator==(const AtomData &d) const { return (*this==d); }
  };

  /*!
   * \brief Class to load and handle atom types from disk
   * \author Mikael Lund
   * \date Lund 2008
   *
   * Example:\n
   * \code
   * atoms atom;
   * atom.load("atoms.dat");           // load parameters
   * double s=atom["Na"].sigma;        // get LJ sigma for sodium
   * particle p = atom["Cl"];          // Copy properties to particle
   * std::cout << atom[p.id].activity; // Get property via particle id
   * \endcode
   *
   * The file format of the atom list is as follows -- note that atom sizes
   * are given as the radius, i.e. sigma/2:
   \verbatim
       #     name    charge  radius     epsilon      Mw       hydrophobic?
       #             (e)     (angstrom) (kJ/mol)     (g/mol)  (yes/no)
       Atom  Na      +1      1.665      0.01158968   22       no
       Atom  Cl      -1      2.200      0.4184       36       no
   \endverbatim
   */
  class AtomTypes {
    private:
      string filename;
    public:
      AtomTypes();                           //!< Constructor - set UNK atom type (fallback)
      std::vector<AtomData> list;            //!< List of atoms
      bool includefile(string);              //!< Append atom parameters from file
      bool includefile(InputMap&);           //!< Append atom parameters from file
      AtomData& operator[] (string);         //!< Name->data
      AtomData& operator[] (particle::Tid);  //!< Id->data
      string info();                         //!< Print info
      void reset_properties(p_vec&);         //!< Reset particle properties according to particle id
  };

  extern AtomTypes atom; // GLOBAL SCOPE
}//namespace
#endif
