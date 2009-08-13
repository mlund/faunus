#ifndef FAU_species_h
#define FAU_species_h
#include "faunus/point.h"
#include "faunus/inputfile.h"

namespace Faunus {
  /*!
   * Class to load and handle atom types from disk
   * Example:\n
   * \code
   *   atoms atom;
   *   atom.load("atoms.dat");         // load parameters
   *   double s=atom["Na"].sigma;      // get LJ sigma for sodium
   *   particle p=atom("Na");          // Initialize a particle
   *   std::cout << atom[p.id].charge; // Get charge via particle id
   * \endcode
   *
   * \author Mikael Lund
   * \date Lund 2008
   */
  class atoms {
    private:
      void init();                  //!< Recalc eps and sigma vectors
      string filename;
    public:
      particle get(char);           //!< Convert n'th atom to a particle
      char find(string);            //!< Find atom id from name
      char find(double, double=0.1);//!< Find atom id from molecular weight
      struct data {
        char id;          //!< id number
        double sigma,     //!< LJ diameter
               eps,       //!< LJ epsilon
               radius,    //!< Radius
               mw,        //!< Weight
               charge,    //!< Charge
               pka;       //!< pKa value
        bool hydrophobic; //!< Are we hydrophobic?
        string name;
        bool operator==(const data &d) const { return (*this==d); }
      };
      atoms();                               //!< Constructor - set UNK atom type (fallback)
      vector<data> list;                     //!< List of atoms
      vector< vector<double> >
        qq,                                  //!< Charge product between atoms i and j
        eps,                                 //!< LJ epsilon between atoms i and j
        sigma;                               //!< LJ sigma between atoms i and j
      bool load(string);                     //!< Load atom parameter from a file
      bool load(inputfile &);                //!< Load atom parameter from a file
      particle set(particle &, char);        //!< Set particle properties
      particle operator() (string);          //!< Name->particle
      particle operator() (char);            //!< Id->particle
      data & operator[] (string);            //!< Name->data
      data & operator[] (char);              //!< Id->data
      string info();                         //!< Print info
      void reset_properties(vector<particle> &);//!< Reset particle properties according to particle id
  };

  extern atoms atom; // GLOBAL SCOPE
}//namespace
#endif
