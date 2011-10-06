#ifndef FAUNUS_GROUP_H
#define FAUNUS_GROUP_H

#include <faunus/common.h>
#include <faunus/average.h>

namespace Faunus {
  class group {
    protected:
      virtual std::ostream& write(std::ostream &) const; //!< Write all group data to stream

    public:
      enum type {SALT, MOLECULE, CIGAR, HYPER};
      type id;
      group(int=-1, int=-1);
      string info();
      string name;
      int beg;                  //!< index of first particle
      int end;                  //!< index of last particle
      int size() const;         //!< number of particles in group.
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
      bool operator==(const group&) const;                     //!< Compare two groups
      group& operator+=(const group&);                         //!< Add two groups
      const group operator+(const group&) const;
      friend std::ostream &operator<<(std::ostream&, group&);  //!< Output group data to stream
      virtual group &operator<<(std::istream &);               //!< Get group data from stream
      virtual ~group() {};
  };

  class Atomic : public group {
  };

  class Molecular : public group {
    protected:
      std::ostream & write(std::ostream&) const;  //!< Write all group data to stream
      point cm_trial;                             //!< mass center vector for trial position

    public:
      point cm;                 //!< mass center vector
      average<double> Q;        //!< average net charge
      average<double> mu;       //!< average dipole moment

      void translate(Space&, const point&);
      void accept(Space&);

      Molecular &operator<<(std::istream &);                        //!< Get information
  };

}//namespace

#endif

