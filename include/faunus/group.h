#ifndef FAUNUS_GROUP_H
#define FAUNUS_GROUP_H

#include <faunus/common.h>
#include <faunus/average.h>
#include <faunus/geometry.h>

namespace Faunus {
  
  class Group {
    protected:
      virtual std::ostream& write(std::ostream &) const; //!< Write all Group data to stream
      virtual Point _massCenter(const Space&) const;

    public:
      enum type {GROUP,ATOMIC,MOLECULAR,CIGAR,HYPER};
      enum prop {PVEC, TRIALVEC, RECALC};
      type id;
      Group(int=-1, int=-1);
      string info();
      string name;
      Point cm_trial;           //!< mass center vector for trial position
      Point cm;                 //!< mass center vector
      int beg;                  //!< index of first particle
      int end;                  //!< index of last particle
      int size() const;         //!< number of particles in Group.
      int random() const;
      bool find(int) const;     //!< Check if index is part of group
      vector<Move::Movebase*> moves;    //!< pointers to move functions

      virtual double charge(const p_vec&) const;             //!< Calculate total charge

      Point massCenter(const Space&) const;                  //!< Calculate mass center - does not set touch group!
      Point setMassCenter(const Space &);                    //!< Calculate and Set mass center (cm and cm_trial)
      virtual Point dipolemoment(const Space&) const;        //!< Calculate dipole moment

      virtual void rotate(Space&, const Point&, double);     //!< Rotate around a vector
      virtual void translate(Space&, const Point&);          //!< Translate along a vector
      virtual void scale(Space&, double);                    //!< Volume scaling
      virtual void undo(Space&);                             //!< Undo move operation
      virtual void accept(Space&);

      // Operators
      bool operator==(const Group&) const;                     //!< Compare two Groups
      Group& operator+=(const Group&);                         //!< Add two Groups
      const Group operator+(const Group&) const;
      friend std::ostream &operator<<(std::ostream&, Group&);  //!< Output Group data to stream
      virtual Group &operator<<(std::istream &);               //!< Get Group data from stream
      virtual ~Group();
  };

  /*!
   * \brief Group class for atomic species - for example salt particles.
   */
  class GroupAtomic : public Group {
    public:
      GroupAtomic();
      GroupAtomic(Space&, InputMap&);        //!< Construct and call add()
      GroupAtomic &operator<<(std::istream&);
      void add(Space&, InputMap&);      //!< Add atomic particles via InputMap paramters
  };

  class GroupMolecular : public Group {
    private:
      Geometry::VectorRotate vrot;
    protected:
      std::ostream & write(std::ostream&) const;  //!< Write all Group data to stream

    public:
      Average<double> Q;        //!< average net charge
      Average<double> mu;       //!< average dipole moment

      GroupMolecular();
      void translate(Space&, const Point&);
      void rotate(Space&, const Point&, double);     //!< Rotate around a vector
      void scale(Space&, double);

      GroupMolecular &operator<<(std::istream&);                        //!< Get information
  };

}//namespace

#endif

