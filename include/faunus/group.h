#ifndef FAUNUS_GROUP_H
#define FAUNUS_GROUP_H

#include <faunus/common.h>
#include <faunus/average.h>

namespace Faunus {
  class Group {
    protected:
      virtual std::ostream& write(std::ostream &) const; //!< Write all Group data to stream

    public:
      enum type {GROUP,ATOMIC,MOLECULAR,CIGAR,HYPER};
      type id;
      Group(int=-1, int=-1);
      string info();
      string name;
      point cm_trial;           //!< mass center vector for trial position
      point cm;                 //!< mass center vector
      int beg;                  //!< index of first particle
      int end;                  //!< index of last particle
      int size() const;         //!< number of particles in Group.
      int random() const;
      vector<Move::Movebase*> moves;    //!< pointers to move functions

      virtual double charge(const p_vec&) const;      //!< Calculate total charge
      virtual point masscenter(const Space&) const;          //!< Calculate mass center
      virtual point dipolemoment(const Space&) const;        //!< Calculate dipole moment

      virtual void rotate(Space&, const point&, double);     //!< Rotate around a vector
      virtual void translate(Space&, const point&);          //!< Translate along a vector
      virtual void scale(Space&, double);                    //!< Volume scaling
      virtual void undo(Space&);                                 //!< Undo move operation
      virtual void accept(Space&);

      // Operators
      bool operator==(const Group&) const;                     //!< Compare two Groups
      Group& operator+=(const Group&);                         //!< Add two Groups
      const Group operator+(const Group&) const;
      friend std::ostream &operator<<(std::ostream&, Group&);  //!< Output Group data to stream
      virtual Group &operator<<(std::istream &);               //!< Get Group data from stream
      virtual ~Group() {};
  };

  /*!
   * \brief Group class for atomic species - for example salt particles.
   */
  class Atomic : public Group {
    public:
      Atomic();
      Atomic(Space&, InputMap&);        //!< Construct and call add()
      Atomic &operator<<(std::istream&);
      void add(Space&, InputMap&);      //!< Add atomic particles via InputMap paramters
      void scale(Space&, double);       //!< Scale all atomic particles in Group to new volume
  };

  class Molecular : public Group {
    protected:
      std::ostream & write(std::ostream&) const;  //!< Write all Group data to stream

    public:
      average<double> Q;        //!< average net charge
      average<double> mu;       //!< average dipole moment

      Molecular();
      void translate(Space&, const point&);
      void accept(Space&);
      void scale(Space&, double);

      Molecular &operator<<(std::istream&);                        //!< Get information
  };

}//namespace

#endif

