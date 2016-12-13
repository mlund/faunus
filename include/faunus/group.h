#ifndef FAUNUS_GROUP_H
#define FAUNUS_GROUP_H

#ifndef SWIG
#include <faunus/common.h>
#include <faunus/geometry.h>
#include <faunus/range.h>
#include <faunus/species.h>
#include <faunus/textio.h>
#endif

namespace Faunus
{

  typedef ContinuousRange<int> Range;  //<! Basic, continuous range of integers

  /**
   * @brief Defines a continuous range of particles as described in `ContinuousRange`
   * @todo Is `isRange()` functionality needed?
   */
  class Group : public Range
  {
  private:
      int molsize; // Number of atoms per molecule


  public:
      typedef PropertyBase::Tid Tid;
      int molIndex = -1; /// \brief Index of group in molecular tracker (sum of previous molecules in vector and on lower TID)
      int Index;        /// \brief Index of group in molecular tracker sub vector
      string name;                            //!< Information time (and short) name
      Point cm_trial;                         //!< mass center vector for trial position
      Point cm;                               //!< mass center vector
      Tid molId;                              //!< molecule id

      /**
       * @brief Constructor
       * @param front First index
       * @param back Last index
       */
      Group( int front = -1, int back = -1 ) : Range(front, back - front + 1), name("")
      {
          molId = 0;
          setMolSize(-1);
          if ( front < 0 || back < 0 )
              resize(0);
      }

      Group( string name, int front = -1, int back = -1 ) : Range(front, back - front + 1), name(name)
      {
          molId = 0;
          setMolSize(-1);
          if ( front < 0 || back < 0 )
              resize(0);
      }

      inline int getMolSize() { return molsize; }

      /* @brief Shift front and back index by `n` */
      inline void shift( int n )
      {
          setfront(front() + n);
          setback(back() + n);
          assert(front() >= 0);
      }

      /* @brief Get begin/end iterators to contained particles in particle vector */
      template<class Tpvec>
      std::pair<typename Tpvec::iterator, typename Tpvec::iterator> to_iter( Tpvec &p )
      {
          assert(p.size() >= size());
          assert(size() == std::distance(begin(), end()));
          auto _beg = p.begin() + front();
          auto _end = _beg + std::distance(begin(), end());
          return {_beg, _end};
      }

      /* @brief Get begin/end iterators to contained particles in particle vector */
      template<class Tpvec>
      std::pair<typename Tpvec::const_iterator, typename Tpvec::const_iterator> to_iter( const Tpvec &p )
      {
          assert(p.size() >= size());
          assert(size() == std::distance(begin(), end()));
          auto _beg = p.begin() + front();
          auto _end = _beg + std::distance(begin(), end());
          return {_beg, _end};
      }

      /* @brief Greater than operator w. respect to back element */
      bool operator>( Group &other ) { return this->back() > other.back(); }

      bool operator==( Group &other )
      {
          return (
              (this->front() == other.front()) &&
                  (this->back() == other.back()) &&
                  (this->molId == other.molId) &&
                  (this->molsize == other.molsize));
      }

      /** @brief Convert front/back to vector of ints, i.e. `[front, front+1, ... , back]` */
      vector<int> range() const
      {
          std::vector<int> v(size());
          if ( size() > 0 )
          {
              std::iota(v.begin(), v.end(), front());
              assert(v.front() == front() && v.back() == back());
          }
          return v;
      }

      /** @brief Information string */
      std::string info()
      {
          using namespace textio;
          using namespace std;
          char w = 15;
          ostringstream o;
          o << header("Group: " + name)
            << pad(SUB, w, "Size") << size() << endl
            << pad(SUB, w, "Range")
            << ((empty()) ? "Empty" : "[" + to_string(front()) + "-" + to_string(back()) + "]")
            << endl
            << pad(SUB, w, "Mol size") << molsize << endl
            << pad(SUB, w, "Molecules") << numMolecules() << endl
            << pad(SUB, w, "Mass center") << cm.transpose() << endl;
          return o.str();
      }

      /** @brief Random particle index */
      int random() const
      {
          if ( !empty())
              return *slump.element(begin(), end());
          return -1;
      }

      /** @brief True if group represents a molecule */
      inline bool isMolecular() const { return (molsize == size()); }

      /** @brief True if group represents a molecule */
      inline bool isAtomic() const { return (molsize == 1); }

      /** @brief True if group contains a range of groups */
      inline bool isRange() const { return (molsize > 1 && molsize < size()); }

      /** @brief Number of atoms in contained molecule(s) */
      inline void setMolSize( int N )
      {
          molsize = N;
          assert((size() % molsize) == 0);
      }

      /** @brief Number of molecules in group */
      inline int numMolecules() const
      {
          assert((size() % molsize) == 0);
          return size() / molsize;
      }

      /**
       *  @brief   Calculates mass center - does not touch group!
       *  @warning Intra-molecular distances must not exceed half
       *           the box size for cuboid geometry.
       *  @todo    Implement assertion to catch failure when molecule
       *           is bigger than half the box size.
       */
      template<class Tspace>
      Point massCenter( const Tspace &spc ) const
      {
          return Geometry::massCenter(spc.geo, spc.p, *this);
      }

      /** @brief Calculate AND set mass center (cm and cm_trial) */
      template<class Tspace>
      Point setMassCenter( const Tspace &spc )
      {
          cm = massCenter(spc);
          cm_trial = cm;
          return cm;
      }

      /**
       * @brief Translate along a vector
       * @param spc Simulation Space
       * @param p Vector to translate with
       */
      template<class Tspace>
      void translate( Tspace &spc, const Point &p )
      {
          assert(spc.geo.sqdist(cm, massCenter(spc)) < 1e-6 && "Mass center out of sync.");
          cm_trial.translate(spc.geo, p);
          for ( auto i : *this )
              spc.trial[i].translate(spc.geo, p);
      }

      /**
       * @brief Rotate around a vector
       * @param spc Simulation space
       * @param endpoint End point of rotation axis, starting from the mass center
       * @param angle [rad]
       */
      template<class Tspace>
      void rotate( Tspace &spc, const Point &endpoint, double angle )
      {
          assert(spc.geo.dist(cm, massCenter(spc)) < 1e-6);
          Geometry::QuaternionRotate vrot1;
          cm_trial = cm;
          vrot1.setAxis(spc.geo, cm, endpoint, angle);//rot around CM->point vec
          auto vrot2 = vrot1;
          vrot2.getOrigin() = Point(0, 0, 0);
          for ( auto i : *this )
          {
              spc.trial[i] = vrot1(spc.trial[i]); // rotate coordinates
              spc.trial[i].rotate(vrot2);         // rotate internal coordinates
          }
          assert(spc.geo.dist(cm_trial, massCenter(spc)) < 1e-9
                     && "Rotation messed up mass center. Is the box too small?");
      }

      /**
       * @brief Get the i'th molecule in the group
       * @warning You must manually update the mass center of the returned group
       */
      template<class Tgroup>
      void getMolecule( int i, Tgroup &sel ) const
      {
          sel.setfront(front() + i * molsize);
          sel.setback(sel.front() + molsize - 1);
          sel.setMolSize(sel.size());

          assert(sel.molsize > 0);
          assert((sel.size() % molsize) == 0);
          assert(sel.isMolecular());
          assert(find(sel.front()));
          assert(find(sel.back()));
      }

      /**
       * @brief Get the i'th molecule in the group
       * @warning You must manually update the mass center of the returned group
       */
      Group getMolecule( int i ) const
      {
          Group sel(name, front() + i * molsize, front() + i * molsize + molsize - 1);
          sel.molId = this->molId;
          sel.setMolSize(molsize);

          assert(sel.back() <= this->back());
          assert(sel.molsize > 0);
          assert((sel.size() % molsize) == 0);
          assert(sel.isMolecular());
          assert(find(sel.front()));
          assert(find(sel.back()));
          return sel;
      }

      /** @brief Scaling for isobaric and isochoric moves */
      template<class Tspace>
      void scale( Tspace &spc, Point &s, double xyz = 1, double xy = 1 )
      {
          if ( empty())
              return;

          if ( isAtomic())
          {
              cm_trial = cm;
              cm_trial.scale(spc.geo, s, xyz, xy);
              for ( auto i : *this )
                  spc.trial[i].scale(spc.geo, s, xyz, xy);
              return;
          }

          if ( isMolecular())
          {
              assert(spc.geo.dist(cm, massCenter(spc)) < 1e-6);
              assert(spc.geo.dist(cm, cm_trial) < 1e-7);

              Point newcm = cm;
              newcm.scale(spc.geo, s, xyz, xy);
              translate(spc, -cm);                 // move to origo

              Point oldlen = spc.geo.len; // store original volume
              Point newlen = oldlen;
              newlen.scale(spc.geo, s, xyz, xy);
              spc.geo.setlen(newlen);         // apply trial volume

              for ( auto i : *this )
              {
                  spc.trial[i] += newcm;            // move all particles to new cm
                  spc.geo.boundary(spc.trial[i]);  // respect boundary conditions
              }
              cm_trial = newcm;
              spc.geo.setlen(oldlen);         // restore original volume
              return;
          }

          if ( isRange())
          {
              for ( int i = 0; i != numMolecules(); ++i )
              {
                  Group sel;
                  getMolecule(i, sel);
                  sel.setMassCenter(spc);
                  sel.scale(spc, s, xyz, xy);
              }
              return;
          }
      }

      /** @brief Undo move operation */
      template<class Tspace>
      void undo( Tspace &s )
      {
          for ( auto i : *this )
              s.trial[i] = s.p[i];
          cm_trial = cm;
      }

      /** @brief Accept a trial move */
      template<class Tspace>
      void accept( Tspace &s )
      {
          for ( auto i : *this )
              s.p[i] = s.trial[i];
          cm = cm_trial;
      }

      /** @brief Write group data to stream */
      friend std::ostream &operator<<( std::ostream &o, const Group &g )
      {
          o << g.front() << " " << g.back() << " "
            << int(g.molId) << " " << g.molsize << " " << g.cm.transpose();
          return o;
      }

      /** @brief Read group data from stream */
      Group &operator<<( std::istream &in )
      {
          int front, back, id;
          in >> front >> back >> id >> molsize;
          setrange(front, back);
          molId = PropertyBase::Tid(id);
          assert(size() == back - front + 1 && "Problem with Group range");
          cm.operator<<(in);
          cm_trial = cm;
          return *this;
      }

      /**
       * @brief Add atomic particles via `InputMap` parameters
       *
       * The InputMap is scanned for the following keywords, starting with X=1:
       *
       * Key            | Description
       * :------------- | :---------------------
       * `tionX`        | Name of atom X
       * `nionX`        | Number of type X atoms
       * `overlap`      | Allow overlap of atoms [default: no]
       *
       * @todo Rename to addAtoms; rename 'overlap' keyword - perhaps to 'overlapionX' ?
       */
      template<class Tspace, class Tinputmap>
      void addParticles( Tspace &spc, Tinputmap &in )
      {
          name = "Atomic Species";
          setfront(spc.p.size());
          int size = 0, n = 1, npart;

          auto overlap = Tspace::OVERLAP_CHECK;
          if ( in.get("overlap", false))
              overlap = Tspace::NOOVERLAP_CHECK;

          do
          {
              std::ostringstream nion("nion"), tion("tion");
              nion << "nion" << n;
              tion << "tion" << n++;
              npart = in.get(nion.str(), 0);
              if ( npart > 0 )
              {
                  auto id = atom[in.get(tion.str(), string("UNK"))].id;
                  spc.insert(atom[id].name, npart, overlap);
                  size += npart;
              }
              else
                  break;
          }
          while ( npart > 0 );
          if ( size > 0 )
              resize(size);
          else
              resize(0);
          setMolSize(1);
          setMassCenter(spc);
          spc.enroll(*this);
      }

      /**
       * @brief Add atomic particles, checks overlaps
       * @param spc Space
       * @param name Name of particle type
       * @param count Number of particles
       */
      template<class Tspace>
      void addParticles( Tspace &spc, string &name, int count )
      {
          if ( spc.insert(name, count))
          {
              if ( size() < 0 )
                  resize(count);
              else
                  resize(size() + count);
              setMolSize(1);
              setMassCenter(spc);
          }
          else
          {
              cout << "Error inserting particles, group.h:329";
              cout << " Group::addParticles(Tspace &spc, string& name, int count)" << endl;
          }
      }
  };

  /** @brief Number of hydrophobic sites */
  template<class Tpvec, class Tindex>
  int numHydrophobic( const Tpvec &p, const Tindex &g )
  {
      return std::count_if(g.begin(), g.end(),
                           [&]( int i ) { return p[i].hydrophobic; });
  }

  /**
   * @brief Summed valency of a set of particles
   * @param p Particle vector
   * @param g Range (`Group` or arbitrary container with index)
   * @param Z Starting charge (default: 0)
   */
  template<class Tparticle, class Talloc, class Tindex>
  double netCharge( const std::vector<Tparticle, Talloc> &p, const Tindex &g, double Z = 0 )
  {
      for ( auto i : g )
          Z += p[i].charge;
      return Z;
  }

  template<class Titer>
  double netCharge( const Titer &beg, const Titer &end, double Z = 0 )
  {
      for ( auto i = beg; i != end; i++ )
          Z += i->charge;
      return Z;
  }

}//namespace
#endif

