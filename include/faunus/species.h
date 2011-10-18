#ifndef FAU_species_h
#define FAU_species_h

#include <faunus/common.h>

namespace Faunus {
  class inputfile;

  class AtomData {
    public:
      AtomData();
      short  id;         //!< Identification number
      double sigma,      //!< LJ diameter
             eps,        //!< LJ epsilon
             radius,     //!< Radius
             mw,         //!< Weight
             charge,     //!< Charge
             chempot,    //!< Chemical potential
             dp,         //!< Displacement parameter
             mean,       //!< Mean value... (charge, sasa, etc.)
             variance;   //!< ...and the spread around it.
      bool hydrophobic;  //!< Are we hydrophobic?
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
   * atom.load("atoms.dat");          // load parameters
   * double s=atom["Na"].sigma;       // get LJ sigma for sodium
   * particle p = atom["Cl"];         // Copy properties to particle
   * std::cout << atom[p.id].chempot; // Get property via particle id
   * \endcode
   */
  class AtomTypes {
    private:
      void init();                  //!< Recalc eps and sigma vectors
      string filename;
    public:
      AtomTypes();                               //!< Constructor - set UNK atom type (fallback)
      vector<AtomData> list;                 //!< List of atoms
      vector< vector<double> >
        qq,                                  //!< Charge product between atoms i and j
        eps,                                 //!< LJ epsilon between atoms i and j
        sigma;                               //!< LJ sigma between atoms i and j

      bool includefile(string);              //!< Append atom parameters from file
      bool includefile(inputfile&);          //!< Append atom parameters from file
      AtomData& operator[] (string);         //!< Name->data
      AtomData& operator[] (short);          //!< Id->data
      string info();                         //!< Print info
      void reset_properties(vector<particle> &);//!< Reset particle properties according to particle id
  };

  extern AtomTypes atom; // GLOBAL SCOPE
}//namespace
#endif
