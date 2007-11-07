/*
 * GROUP (M.Lund 2004)
 *
 * Defines a set of particles in terms of a list of integers
 * that points to the main particle vector in the SPACE class
 *
 * - Maaske klassen skulle hedde "SLICE" i stedet.
 * - der skal vaere en overloaded '+' operator til at tilfoeje partikler
 *   til 'v' vektoren... void +(int)
 * 
 *
 */

#ifndef group_h
#define group_h

#include <string>
#include <iostream>
#include <iomanip>

using std::ostream;
using namespace std;

//! Groups set of particles, such as molecules, salt etc.
/*!
 * A group defines the range of a a protein, the mobile ions,
 * a chain or whatever. If is used extensively in MC movements,
 * energy calculations and many, many other places.
 */
class group {
public:
  friend ostream &operator<<(ostream &, group &);

  short int beg,end;    //define range in particle vector. [beg;end]
  short int cm;         //points to mass-center vector
  short int dipv;       /// Points to the dipole vector. Specified by space::adddipole()
  string name;          /// Informative (arbitrary) name
  
  double radius;        ///< Smallest sphere that can include the group. Set with with space::radius()
  double dip;           ///< Dipole moment scalar. Updated with space::calcdipole()
  double charge;        ///< Total charge. Updated with space::charge()

  bool vdw;             ///< True if we should calculate vdW interactions.
  bool chain;           ///< True if the group is a chain point::len()
  short int graftpoint; ///< Chain is grafted to this point. -1 if free (default)
  
  group(int=0);         ///< Constructor, initialize data.
  
  void set(short int,short int);///< Set particle range, "beg" and "end".
  short int size();             ///< Number of particles in group
  short int random();           ///< Picks a random particle within this group
  bool isingroup(int);          ///< Check if particle i is included in this group
  void operator++(int);         ///< Expand range  [beg;end+1]
  void operator--(int);         ///< Decrease range [beg;end-1]
  
  void operator+=(group);
  group operator+(group);
  

};
#endif

