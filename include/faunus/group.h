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
   * @todo Is `isRange()` functionality needed?
   */
  class Group : public Range {
    private:
      int molsize; // Number of atoms per molecule
    public:
      typedef PropertyBase::Tid Tid;

      /**
       * @brief Constructor
       * @param front First index
       * @param back Last index
       */
      Group(int front=-1, int back=-1) : Range(front,back-front+1), name("") {
        molId=0;
        setMolSize(-1);
        if (front<0 || back<0)
          resize(0);
      }

      Group(string name, int front=-1, int back=-1) : Range(front,back-front+1), name(name) {
        molId=0;
        setMolSize(-1);
        if (front<0 || back<0)
          resize(0);
      }

      string name;                            //!< Information time (and short) name
      Point cm_trial;                         //!< mass center vector for trial position
      Point cm;                               //!< mass center vector
      Tid molId;                              //!< molecule id

      inline int getMolSize() {return molsize;}

      /* @brief Shift front and back index by `n` */
      inline void shift(int n) {
        setfront( front()+n );
        setback( back()+n );
        assert(front()>=0);
      }

      /* @brief Greater than operator w. respect to back element */
      bool operator> (Group &other) { return this->back() > other.back(); }

      bool operator== (Group& other) {
        return (
            (this->front() == other.front()) &&
            (this->back() == other.back()) &&
            (this->molId == other.molId) &&
            (this->molsize == other.molsize)  );
      }

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

      /** @brief Random particle index */
      int random() const {
        if (!empty())
          return *slump.element(begin(), end());
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
       *           the box size for cuboid geometry.
       *  @todo    Implement assertion to catch failure when molecule
       *           is bigger than half the box size.
       */
      template<class Tgeometry, class Tparticles>
        Point massCenter(const Tgeometry &geo, const Tparticles &p) const {
          return Geometry::massCenter(geo, p, *this);
        }
        
      template<class Tgeometry, class Tparticles>
        Point setTrialMassCenter(const Tgeometry &geo, const Tparticles &p) {
	    cm_trial = massCenter(geo,p);
	    return cm_trial;
        }
        
      template<class Tgeometry, class Tparticles>
        Point setMassCenter(const Tgeometry &geo, const Tparticles &p) {
	    cm = massCenter(geo,p);
	    return cm;
        }
        
      /** @brief Calculate AND set mass center (cm and cm_trial) */
      template<class Tspace>
        void setMassCenters(const Tspace &spc) {
	  setTrialMassCenter(spc.geo_trial,spc.trial);
	  setMassCenter(spc.geo,spc.p);
        }

      /**
       * @brief Translate along a vector
       * @param geo Geometry
       * @param p Particle vector
       * @param vec Vector to translate with
       * @param isTrial Whether 'geo' and 'p' are trial entities
       */
      template<class Tgeometry, class Tparticles>
        void translate(const Tgeometry &geo, Tparticles &p, const Point &vec, bool isTrial = true) {
	  if(isTrial) {
	    cm_trial.translate(geo, vec);
	  } else {
	    cm.translate(geo, vec);
	  }
          for (auto i : *this)
            p[i].translate(geo, vec);
        }
        
      /**
       * @brief Rotate around a vector
       * @param spc Simulation space
       * @param endpoint End point of rotation axis, starting from the mass center
       * @param angle [rad]
       */
      template<class Tgeometry, class Tparticles>
        void rotate(Tgeometry &geo, Tparticles &p, const Point &endpoint, Point cm_in, double angle, bool isTrial = true) {
	  assert( geo.dist(cm_in,massCenter(geo,p) )<1e-9 );
          Geometry::QuaternionRotate vrot1;
          vrot1.setAxis(geo, cm_in, endpoint, angle);//rot around CM->point vec
          auto vrot2 = vrot1;
          vrot2.getOrigin() = Point(0,0,0);
          for (auto i : *this) {
            p[i] = vrot1(p[i]); // rotate coordinates
            p[i].rotate(vrot2);         // rotate internal coordinates
          }
	  assert( geo.dist(cm_in, massCenter(geo,p))<1e-9 && "Rotation messed up mass center. Is the box too small?");
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

      /**
       * @brief Get the i'th molecule in the group
       * @warning You must manually update the mass center of the returned group
       */
      Group getMolecule(int i) const {
        Group sel(name, front()+i*molsize, front()+i*molsize+molsize-1);
        sel.molId = this->molId;
        sel.setMolSize(molsize);

        assert( sel.back() <= this->back());
        assert( sel.molsize>0 );
        assert( (sel.size()%molsize)==0 );
        assert( sel.isMolecular() );
        assert( find( sel.front() ) );
        assert( find( sel.back()  ) );
        return sel;
      }
        
      /** @brief Scaling for isobaric and isochoric moves. Assumes that if 'geo_trial' is used so is the trial particle vector.
       */ 
      template<class Tgeometry, class Tparticles>
        void scale(const Tgeometry &geo, const Tgeometry &geo_old, Tparticles &p, Point &s, double xyz=1, double xy=1, bool isTrial = true) {
          if (empty()) return;

          if (isAtomic()) {
	    if(isTrial) {
	      cm_trial=cm;
	      cm_trial.scale(geo,s,xyz,xy);
	    } else {
	      cm.scale(geo,s,xyz,xy);
	    }
	    for (auto i : *this)
	      p[i].scale(geo,s,xyz,xy);
	    return;
          }

          if (isMolecular()) {
	    if(isTrial) {
	      assert( geo.dist(cm_trial, cm) < 1e-9);
	    } else {
	      assert( geo.dist(cm, Geometry::massCenter(geo, p, *this)) < 1e-9);
	    }

            Point newcm = cm;
            newcm.scale(geo,s,xyz,xy);                
	    
            translate(geo_old,p,-cm,isTrial);                 // move to origo
	    translate(geo,p,newcm,isTrial);               // move to scaled position
	    

            if(isTrial) {
	      cm_trial=newcm;
	      assert( geo.dist(cm_trial, Geometry::massCenter(geo, p, *this)) < 1e-9);
	    } else {
	      cm = newcm;
	    }
            return;
          }

          if (isRange()) {
            for (int i=0; i!=numMolecules(); ++i) {
              Group sel;
              getMolecule(i,sel);
	      if(isTrial) {
		sel.setTrialMassCenter(geo,p);
	      } else {
		sel.setMassCenter(geo,p);
	      }
              sel.scale(geo,geo_old,p,s,xyz,xy,isTrial);
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
        o << g.front() << " " << g.back() << " "
          << int(g.molId) << " " << g.molsize << " " << g.cm.transpose();
        return o;
      }

      /** @brief Read group data from stream */
      Group& operator<<(std::istream &in) {
        int front, back, id;
        in >> front >> back >> id >> molsize;
        setrange(front,back);
        molId=PropertyBase::Tid(id);
        assert( size()==back-front+1 && "Problem with Group range");
        cm.operator<<(in);
        cm_trial=cm;
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
       * @brief Add atomic particles, checks overlaps
       * @param spc Space
       * @param name Name of particle type
       * @param count Number of particles
       */
      template<class Tspace>
        void addParticles(Tspace &spc, string& name, int count) {
          if(spc.insert(name, count) ) {
            if(size() < 0) resize(count);
            else resize(size()+count);
            setMolSize(1);
            setMassCenter(spc);
          } else {
            cout << "Error inserting particles, group.h:329";
            cout << " Group::addParticles(Tspace &spc, string& name, int count)" << endl;
          }
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
  template<class Tparticle, class Talloc, class Tindex>
    double netCharge(const std::vector<Tparticle,Talloc> &p, const Tindex &g, double Z=0) {
      for (auto i : g)
        Z += p[i].charge;
      return Z;
    }

  template<class Titer>
    double netCharge(const Titer &beg, const Titer &end, double Z=0) {
      for (auto i=beg; i!=end; i++)
        Z += i->charge;
      return Z;
    }

}//namespace
#endif

