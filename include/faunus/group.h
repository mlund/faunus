#ifndef FAUNUS_GROUP_H
#define FAUNUS_GROUP_H

#include <faunus/common.h>
#include <faunus/average.h>
#include <faunus/geometry.h>

namespace Faunus {

  //http://stackoverflow.com/questions/7185437/is-there-a-range-class-in-c0x-aka-c11-for-use-with-range-based-for-loops
  //https://bitbucket.org/AraK/range/wiki/Home
  template<class Tint=int>
  class myrange {
    public:
      class iterator {
        friend class myrange;
        public:
        Tint operator *() const { return i_; }
        const iterator &operator ++() { ++i_; return *this; }
        iterator operator ++(int) { iterator copy(*this); ++i_; return copy; }
        iterator operator +(int i) { iterator copy(*this); copy.i_+=i; return copy; }
        iterator operator -(int i) { iterator copy(*this); copy.i_-=i; return copy; }
        bool operator ==(const iterator &other) const { return i_ == other.i_; }
        bool operator !=(const iterator &other) const { return i_ != other.i_; }
        protected:
        iterator(Tint start) : i_ (start) { }
        private:
        Tint i_;
      };
      iterator begin() const { return begin_; }   //!< Iterator to beginning
      iterator end() const { return end_; }       //!< Iterator to end
      Tint front() const { return begin_.i_; }    //!< Get first value in range
      Tint back() const { return end_.i_-1; }     //!< Get last value in range
      void resize(unsigned int size) { end_.i_=begin_.i_+size; } //!< Resize range, keeping same beginning
      void operator++() { end_.i_++; }
      void operator--() { if (end_.i_>=0) end_.i_--; }
      myrange(Tint begin=0, Tint end=0) : begin_(begin), end_(end) {}
      void set(Tint begin, Tint end) {
        begin_.i_=begin;
        end_.i_=end;
      }
    private:
      iterator begin_;
      iterator end_;
  };

  /*!
   * \brief Defines a continuous range of particles in the Space particle vector.
   * \todo Implement iterator
   *
   * This class defines a range, [beg:end], in the particle vector vector and knows how to
   * perform geometric operations on it - rotate, translate etc.
   *
   * \note http://stackoverflow.com/questions/7185437/is-there-a-range-class-in-c0x-aka-c11-for-use-with-range-based-for-loops
   *
   */
  class Group : public myrange<int> {
    protected:
      virtual std::ostream& write(std::ostream &) const; //!< Write all Group data to stream
      virtual Point _massCenter(const Space&) const;
      vector<Move::Movebase*> moves;    //!< pointers to move functions

    public:
      enum type {GROUP,ATOMIC,MOLECULAR,CIGAR,RIGID,ISOBARIC,GRANDCANONICAL};
      std::set<type> property;
      //enum prop {PVEC, TRIALVEC, RECALC};
      type id;
      Group(int=-1, int=-1);
      string info();
      string name;
      Point cm_trial;           //!< mass center vector for trial position
      Point cm;                 //!< mass center vector
      int beg;                  //!< index of first particle
      int last;                 //!< index of last particle
      int size() const;         //!< number of particles in Group.
      bool empty() const;       //!< True if group is empty.
      int random() const;
      bool find(int) const;     //!< Check if index is part of group

      virtual double charge(const p_vec&) const;             //!< Calculate total charge

      Point massCenter(const Space&) const;                  //!< Calculate mass center - does not set touch group!
      Point setMassCenter(const Space &);                    //!< Calculate and Set mass center (cm and cm_trial)
      virtual Point dipolemoment(const Space&) const;        //!< Calculate dipole moment

      virtual void rotate(Space&, const Point&, double);     //!< Rotate around a vector
      virtual void translate(Space&, const Point&);          //!< Translate along a vector
      virtual void scale(Space&, double);                    //!< Volume scaling
      virtual void undo(Space&);                             //!< Undo move operation
      virtual void accept(Space&);                           //!< Accept a trial move

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

  /*!
   * \brief Class for molecular groups - proteins, polymers etc.
   */
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
      void rotate(Space&, const Point&, double);
      void scale(Space&, double);

      GroupMolecular &operator<<(std::istream&);                        //!< Get information
  };

}//namespace

#endif

