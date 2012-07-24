#ifndef FAUNUS_GROUP_H
#define FAUNUS_GROUP_H

#ifndef SWIG
#include <faunus/common.h>
#include <faunus/average.h>
#include <faunus/geometry.h>
#endif

namespace Faunus {

  /* 
   * \brief Range iterator class for continuous ranges
   * \todo Replace with https://bitbucket.org/AraK/range/wiki/Home ??
   * \comment http://stackoverflow.com/questions/7185437/is-there-a-range-class-in-c0x-aka-c11-for-use-with-range-based-for-loops
   */
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
      iterator end() const { return end_; }       //!< Iterator to end (end points to last element + 1)
      Tint front() const { return begin_.i_; }    //!< Get first value in range
      Tint back() const { return end_.i_-1; }     //!< Get last value in range
      bool empty() const { return (end_.i_<=begin_.i_) ? true : false; } //!< Determines if range is empty.
      void resize(unsigned int size) { end_.i_ = begin_.i_ + size; } //!< Resize range, keeping same beginning
      Tint size() const { return end_.i_-begin_.i_; }  //!< Size of range
      void setfront(Tint front) { begin_.i_=front; }   //!< Set first element
      void setback(Tint back) { end_.i_=back+1; }      //!< Set last element
      void setrange(Tint front, Tint back=-1) {        //!< Set range [front:back]
        setfront(front);
        if (back!=-1)
          setback(back);
        else setback(front-1);
      }
      myrange(Tint first=0, Tint size=0) : begin_(first), end_(first+size) {}
    private:
      iterator begin_;
      iterator end_;
  };

  /*!
   * \brief Defines a continuous range of particles in the Space particle vector.
   * \todo Implement iterator
   * \note http://stackoverflow.com/questions/7185437/is-there-a-range-class-in-c0x-aka-c11-for-use-with-range-based-for-loops
   *
   * This class defines a range, \c [front:back], in the particle vector and behaves much like a standard
   * STL container. Groups know how to perform geometric operations on it - rotate, translate etc. 
   *
   * Example:
   * \code
   *   Group g(2,5);           // first, last particle
   *   for (auto i : g)        // iterator access
   *     cout << i;            // -> 23456
   *   g.front();              // -> 2
   *   g.back();               // -> 6
   *   g.size();               // -> 5
   *   g.resize( g.size()+1 ); // -> size=6, back=7.
   * \endcode
   */
  class Group : public myrange<int> {
    private:
      Geometry::VectorRotate vrot;
      virtual Point _massCenter(const Space&) const;
    protected:
      virtual std::ostream& write(std::ostream &) const; //!< Write all Group data to stream
      vector<Move::Movebase*> moves;    //!< pointers to move functions
      virtual string _info();
      char w;                           //!< Text padding for info() functions

    public:
      //enum type {GROUP,ATOMIC,MOLECULAR,CIGAR,RIGID,ISOBARIC,GRANDCANONICAL};
      //std::set<type> property;
      //type id;
      Group(int=-1, int=-1);
      string info();                                         //!< Information string
      string name;                                           //!< Information time (and short) name
      Point cm_trial;                                        //!< mass center vector for trial position
      Point cm;                                              //!< mass center vector
      int random() const;                                    //!< Pick random particle index in Group
      bool find(int) const;                                  //!< Check if index is part of group

      virtual double charge(const p_vec&) const;             //!< Calculate total charge

      Point massCenter(const Space&) const;                  //!< Calculate mass center - does not set touch group!
      Point setMassCenter(const Space &);                    //!< Calculate AND set mass center (cm and cm_trial)
      virtual Point dipolemoment(const Space&) const;        //!< Calculate dipole moment

      virtual void rotate(Space&, const Point&, double);     //!< Rotate around a vector
      virtual void translate(Space&, const Point&);          //!< Translate along a vector
      virtual void scale(Space&, double);                    //!< Volume scaling for NPT ensemble
      virtual void undo(Space&);                             //!< Undo move operation
      virtual void accept(Space&);                           //!< Accept a trial move

      virtual bool isMolecular() const;                      //!< True if group represents a molecule
      virtual bool isAtomic() const;                         //!< True if group represents atomic species

      bool operator==(const Group&) const;                   //!< Compare two Groups
      friend std::ostream &operator<<(std::ostream&, Group&);//!< Output Group data to stream
      virtual Group &operator<<(std::istream&);              //!< Get Group data from stream
      virtual ~Group();
  };

  /*!
   * \brief Group class for atomic species, typically salt
   */
  class GroupAtomic : public Group {
    public:
      GroupAtomic(int=-1, int=-1);
      GroupAtomic(Space&, InputMap&);        //!< Construct and call add()
      void add(Space&, InputMap&);           //!< Add atomic particles via InputMap parameters
      bool isAtomic() const;                 //!< Always true for GroupAtomic
  };

  /*!
   * \brief Class for molecular groups - proteins, polymers etc.
   */
  class GroupMolecular : public Group {
    private:
      string _info();
    public:
      Average<double> Q;                        //!< average net charge
      Average<double> mu;                       //!< average dipole moment

      GroupMolecular(int=-1, int=-1);
      void scale(Space&, double);               //!< Mass-center volume scale
      bool isMolecular() const;                 //!< Always true for GroupMolecular
  };

}//namespace
#endif

