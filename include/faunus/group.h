#ifndef FAUNUS_GROUP_H
#define FAUNUS_GROUP_H

#ifndef SWIG
#include <faunus/common.h>
#include <faunus/geometry.h>
#include <faunus/range.h>
#endif

namespace Faunus {

  typedef ContinuousRange<int> Range;  //<! Basic, continuous range of integers

  /**
   * @brief Defines a continuous range of particles as described in `ContinuousRange`
   */
  class Group : public Range {
    private:
      virtual Point _massCenter(const Space&) const;
    protected:
      virtual std::ostream& write(std::ostream &) const; //!< Write all Group data to stream
      virtual string _info();
      char w;                           //!< Text padding for info() functions
    public:
      Group(int=-1, int=-1);
      string info();                          //!< Information string
      string name;                            //!< Information time (and short) name
      Point cm_trial;                         //!< mass center vector for trial position
      Point cm;                               //!< mass center vector
      int random() const;                     //!< Pick random particle index in Group
      double charge(const p_vec&) const;      //!< Calculates total charge
      Point massCenter(const Space&) const;   //!< Calculates mass center - does not touch group!
      Point setMassCenter(const Space &);     //!< Calculate AND set mass center (cm and cm_trial)
      Point dipolemoment(const Space&) const; //!< Calculates dipole moment

      virtual void rotate(Space&, const Point&, double);//!< Rotate around a vector
      virtual void translate(Space&, const Point&);     //!< Translate along a vector
      virtual void scale(Space&, double);               //!< Volume scaling for NPT ensemble
      virtual void undo(Space&);                        //!< Undo move operation
      virtual void accept(Space&);                      //!< Accept a trial move

      virtual bool isMolecular() const;                 //!< True if group represents a molecule
      virtual bool isAtomic() const;                    //!< True if group represents atomic species

      virtual int numMolecules() const; //!< Number of molecules in group

      bool operator==(const Group&) const;                   //!< Compare two Groups
      friend std::ostream &operator<<(std::ostream&, Group&);//!< Output Group data to stream
      virtual Group &operator<<(std::istream&);              //!< Get Group data from stream
      virtual ~Group();
  };

  /**
   * @brief Group class for atomic species, typically salt
   */
  class GroupAtomic : public Group {
    public:
      GroupAtomic(int=-1, int=-1);
      GroupAtomic(Space&, InputMap&);        //!< Construct and call add()
      void add(Space&, InputMap&);           //!< Add atomic particles via InputMap parameters
      bool isAtomic() const;                 //!< Always true for GroupAtomic
  };

  /**
   * @brief Class for molecular groups - proteins, polymers etc.
   */
  class GroupMolecular : public Group {
    private:
      string _info();
    public:
      GroupMolecular(int=-1, int=-1);
      void scale(Space&, double);               //!< Mass-center volume scale
      bool isMolecular() const;                 //!< Always true for GroupMolecular
      int numMolecules() const;                 //!< Number of molecules in group
  };

  /**
   * @brief Class for an array of multiatom molecules - solvent, lipids etc.
   */
  class GroupArray : public GroupMolecular {
    private:
      GroupMolecular sel;             //!< A temporary group class
      string _info();                 //!< Show information
    public:
      int N;                          //!< Number of atoms in each molecule
      GroupArray(int);                //!< Constructor - number of atoms per molecule.
      int randomMol() const;          //!< Pick a random molecule
      GroupMolecular& operator[](int);//!< Access i'th molecule
      void add(const GroupMolecular&);//!< Add a molecule to the array - range must be continuous
      void scale(Space&, double);     //!< Mass-center volume scale
      int numMolecules() const;       //!< Number of molecules
  };

}//namespace
#endif

