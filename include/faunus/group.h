#ifndef FAUNUS_GROUP_H
#define FAUNUS_GROUP_H

#ifndef SWIG
#include <faunus/common.h>
#include <faunus/geometry.h>
#include <faunus/range.h>
#include <faunus/species.h>
#include <faunus/textio.h>
#endif

namespace Faunus {

  typedef ContinuousRange<int> Range;  //<! Basic, continuous range of integers

  /**
   * @brief Defines a continuous range of particles as described in `ContinuousRange`
   */
  class Group : public Range {
    private:
      int molsize; // Number of atoms per molecule
    public:

      /**
       * @brief Constructor
       * @param front First index
       * @param back Last index
       */
      Group(int front=-1, int back=-1) : Range(front,back-front+1) {
        setMolSize(-1);
        if (front<0 || back<0)
          resize(0);
      }

      string name;                            //!< Information time (and short) name
      Point cm_trial;                         //!< mass center vector for trial position
      Point cm;                               //!< mass center vector

      /** @brief Information string */
      std::string info() {
        using namespace textio;
        using namespace std;
        char w=15;
        ostringstream o;
        o << header("Group: " + name)
          << pad(SUB,w,"Size") << size() << endl
          << pad(SUB,w,"Range")
          << ((empty()) ? "Empty" : "["+to_string(front())+"-"+to_string(back())+"]")
          << endl
          << pad(SUB,w,"Mol size") << molsize << endl
          << pad(SUB,w,"Molecules") << numMolecules() << endl
          << pad(SUB,w,"Mass center") << cm.transpose() << endl;
        return o.str();
      }

      /** @brief Pick random particle index in Group */
      int random() const {
        if (!empty()) {
          int i = front() + slp_global.rand() % size();
          assert(find(i) && "Generated random element out of range!");
          return i;
        }
        return -1;
      }

      /** @brief True if group represents a molecule */
      inline bool isMolecular() const { return (molsize==size()); }

      /** @brief True if group represents a molecule */
      inline bool isAtomic() const { return (molsize==1); }

      /** @brief True if group contains a range of groups */
      inline bool isRange() const { return (molsize>1 && molsize<size()); }

      /** @brief Number of atoms in contained molecule(s) */
      inline void setMolSize(int N) {
        molsize=N;
        assert( (size()%molsize)==0 );
      }

      /** @brief Number of molecules in group */
      inline int numMolecules() const {
        assert( (size()%molsize)==0 );
        return size()/molsize;
      }

      /**
       *  @brief   Calculates mass center - does not touch group!
       *  @warning Intra-molecular distances must not exceed half
       *           the box size for cubouid geometry.
       *  @todo    Implement assertion to catch failure when molecule
       *           is bigger than half the box size.
       */
      template<class Tspace>
        Point massCenter(const Tspace &spc) const {
          assert(&spc!=nullptr);
          return Geometry::massCenter(spc.geo, spc.p, *this);
        }

      /** @brief Calculate AND set mass center (cm and cm_trial) */
      template<class Tspace>
        Point setMassCenter(const Tspace &spc) {
          cm=massCenter(spc);
          cm_trial=cm;
          return cm;
        }
        
      /** @brief Calculates electric dipole moment */
      template<class Tspace>
        Point dipolemoment(const Tspace &s, Point mu=Point(0,0,0)) const {
          for (auto i : *this) {
            Point t=s.p[i] - cm;
            s.geo.boundary(t);
            mu += t*s.p[i].charge;
          }
          return mu;
        }

      /**
       * @brief Translate along a vector
       * @param spc Simulation Space
       * @param p Vector to translate with
       */
      template<class Tspace>
        void translate(Tspace &spc, const Point &p) {
          assert( spc.geo.sqdist(cm,massCenter(spc))<1e-6
              && "Mass center out of sync.");
          cm_trial.translate(spc.geo, p);
          for (auto i : *this)
            spc.trial[i].translate(spc.geo, p);
        }

      /**
       * @brief Rotate around a vector
       * @param spc Simulation space
       * @param endpoint End point of rotation axis, starting from the mass center
       * @param angle [rad]
       */
      template<class Tspace>
        void rotate(Tspace &spc, const Point &endpoint, double angle) {
          assert( spc.geo.dist(cm,massCenter(spc) )<1e-6 );
          Geometry::QuaternionRotate vrot1;
          cm_trial = cm;
          vrot1.setAxis(spc.geo, cm, endpoint, angle);//rot around CM->point vec
          auto vrot2 = vrot1;
          vrot2.getOrigin() = Point(0,0,0);
          for (auto i : *this) {
            spc.trial[i] = vrot1(spc.trial[i]); // rotate coordinates
            spc.trial[i].rotate(vrot2);         // rotate internal coordinates
          }
          assert( spc.geo.dist(cm_trial, massCenter(spc))<1e-9
              && "Rotation messed up mass center. Is the box too small?");
        }

      /**
       * @brief Get the i'th molecule in the group
       * @warning You must manually update the mass center of the returned group
       */
      template<class Tgroup>
        void getMolecule(int i, Tgroup &sel) const {
          sel.setfront( front()+i*molsize );
          sel.setback( sel.front()+molsize-1 );
          sel.setMolSize(sel.size());

          assert( sel.molsize>0 );
          assert( (sel.size()%molsize)==0 );
          assert( sel.isMolecular() );
          assert( find( sel.front() ) );
          assert( find( sel.back()  ) );
        }

      /** @brief Volume scaling for NPT ensemble */
      template<class Tspace>
        void scale(Tspace &s, double newvol) {
          if (empty()) return;

          if (isAtomic()) {
            cm_trial=cm;
            cm_trial.scale(s.geo, newvol);
            for (auto i : *this)
              s.trial[i].scale(s.geo, newvol);
            return;
          }

          if (isMolecular()) {
            assert( s.geo.dist(cm, massCenter(s))<1e-6);
            assert( s.geo.dist(cm, cm_trial)<1e-7);

            Point newcm=cm;
            newcm.scale(s.geo, newvol);
            translate(s,-cm);                 // move to origo

            double oldvol=s.geo.getVolume(); // store original volume
            s.geo.setVolume(newvol);         // apply trial volume

            for (auto i : *this) {
              s.trial[i] += newcm;            // move all particles to new cm
              s.geo.boundary( s.trial[i] );  // respect boundary conditions
            }
            cm_trial=newcm;
            s.geo.setVolume(oldvol);         // restore original volume
            return;
          }

          if (isRange()) {
            for (int i=0; i!=numMolecules(); ++i) {
              Group sel;
              getMolecule(i,sel);
              sel.setMassCenter(s);
              sel.scale(s,newvol);
            }
            return;
          }
        }

      /** @brief Undo move operation */
      template<class Tspace>
        void undo(Tspace &s) {
          for (auto i : *this)
            s.trial[i]=s.p[i];
          cm_trial=cm;
        }

      /** @brief Accept a trial move */
      template<class Tspace>
        void accept(Tspace &s) {
          for (auto i : *this)
            s.p[i] = s.trial[i];
          cm=cm_trial;
        }

      /** @brief Write group data to stream */
      friend std::ostream& operator<<(std::ostream &o, const Group &g) {
        o << g.front() << " " << g.back() << " " << g.cm;
        return o;
      }

      /** @brief Read group data from stream */
      Group& operator<<(std::istream &in) {
        int front;
        int back;
        in >> front >> back;
        setrange(front,back);
        assert( size()==back-front+1 && "Problem with Group range");
        cm.operator<<(in);
        return *this;
      }

      /** @brief Select random molecule */
      int randomMol() const {
        int i=(random()-front())/molsize;
        assert(molsize>0);
        assert(find(i) && "Out of range!");
        return i;
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
        void addParticles(Tspace &spc, Tinputmap &in) {
          name="Atomic Species";
          setfront( spc.p.size() );
          int size=0, n=1, npart;

          auto overlap=Tspace::OVERLAP_CHECK;
          if (in.get("overlap", false))
            overlap=Tspace::NOOVERLAP_CHECK;

          do {
            std::ostringstream nion("nion"), tion("tion");
            nion << "nion" << n;
            tion << "tion" << n++;
            npart = in.get(nion.str(), 0);
            if (npart>0) {
              auto id = atom[in.get(tion.str(), string("UNK")) ].id;
              spc.insert(atom[id].name, npart, overlap);
              size+=npart;
            } else break;
          } while (npart>0);
          if (size>0)
            resize(size);
          else
            resize(0);
          setMolSize(1);
          setMassCenter(spc);
          spc.enroll(*this);
        }

      /**
       * @todo rename to addGroup or implement operator
       */
      template<class Tgroup>
        void addMolecule(const Tgroup &g) {
          if ((g.size()%molsize)==0) {
            if (empty())
              setrange(g.front(), g.back());
            else if (g.front()==back()+1)
              setback(g.back());
          }
          assert( (size()%molsize)==0 && "GroupArray not a multiple of N");
        }
  };

  /** @brief Number of hydrophobic sites */
  template<class Tpvec, class Tindex>
    int numHydrophobic(const Tpvec &p, const Tindex &g) {
      return std::count_if(g.begin(), g.end(),
          [&](int i) { return p[i].hydrophobic; });
    }

  /**
   * @brief Summed valency of a set of particles
   * @param p Particle vector
   * @param g Range (`Group` or arbitrary container with index)
   * @param Z Starting charge (default: 0)
   */
  template<class Tpvec, class Tindex>
    double netCharge(const Tpvec &p, const Tindex &g, double Z=0) {
      for (auto i : g)
        Z += p[i].charge;
      return Z;
    }

}//namespace
#endif

