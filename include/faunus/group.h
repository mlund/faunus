#ifndef FAUNUS_GROUP_H
#define FAUNUS_GROUP_H

#ifndef SWIG
#include <faunus/common.h>
#include <faunus/geometry.h>
#include <faunus/range.h>
#include <faunus/species.h>
#include <faunus/inputfile.h>
#endif

namespace Faunus {

  typedef ContinuousRange<int> Range;  //<! Basic, continuous range of integers

  /**
   * @brief Defines a continuous range of particles as described in `ContinuousRange`
   */
  class Group : public Range {
    private:
      /**
       * @warning Intra-molecular distances must not exceed half
       *          the box size for cubouid geometry.
       *    @todo Implement assertion to catch failure when molecule
       *          is bigger than half the box size.
       */
      template<class Tspace>
        Point _massCenter(const Tspace &spc) const {
          return Geometry::massCenter(*spc.geo, spc.p, *this);
        }
    protected:
      virtual std::ostream& write(std::ostream &) const; //!< Write all Group data to stream
      virtual string _info();
    public:
      Group(int=-1, int=-1);
      string info();                          //!< Information string
      string name;                            //!< Information time (and short) name
      Point cm_trial;                         //!< mass center vector for trial position
      Point cm;                               //!< mass center vector
      int random() const;                     //!< Pick random particle index in Group

      /** @brief Total charge */
      template<class Tpvec>
        double charge(const Tpvec &p, double Z=0) const {
          for (auto i : *this) Z+=p[i].charge;
          return Z;
        }

      /** @brief Calculates mass center - does not touch group! */
      template<class Tspace>
        Point massCenter(const Tspace &spc) const {
          assert(&spc!=nullptr);
          return _massCenter(spc);
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
            s.geo->boundary(t);
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
          assert( spc.geo->sqdist(cm,massCenter(spc))<1e-6
              && "Mass center out of sync.");
          for (auto i : *this)
            spc.trial[i].translate(*spc.geo, p);
          cm_trial.translate(*spc.geo, p);
        }

      /**
       * @brief Rotate around a vector
       * @param spc Simulation space
       * @param endpoint End point of rotation axis, starting from the mass center
       * @param angle [rad]
       */
      template<class Tspace>
        void rotate(Tspace &spc, const Point &endpoint, double angle) {
          assert( spc.geo->dist(cm,massCenter(spc) )<1e-6 );
          Geometry::QuaternionRotate vrot;
          cm_trial = cm;
          vrot.setAxis(*spc.geo, cm, endpoint, angle);//rot around CM->point vec
          for (auto i : *this)
            spc.trial[i].rotate(vrot);
          assert( spc.geo->dist(cm_trial, massCenter(spc))<1e-9
              && "Rotation messed up mass center. Is the box too small?");
        }

      virtual void scale(Space&, double);               //!< Volume scaling for NPT ensemble

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

      virtual bool isMolecular() const; //!< True if group represents a molecule
      virtual bool isAtomic() const;    //!< True if group represents atomic species
      virtual int numMolecules() const; //!< Number of molecules in group

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
      bool isAtomic() const;                 //!< Always true for GroupAtomic

      /**
       * @brief Add atomic particles via InputMap parameters
       * The InputMap is scanned for the following keywords, starting with X=1:
       *
       * Key            | Description
       * :------------- | :---------------------
       * `tionX`        | Name of atom X
       * `nionX`        | Number of type X atoms
       * `dpionX`       | (Displacement parameter of atom type X)
       * `aionX`        | (Activity of atom X (molar scale))
       *
       * If the two latter properties, displacement and activity, are omitted
       * (recommended) the values from AtomTypes is used instead. That is, you
       * should specify these directly in the input JSON file.
       */
      template<class Tspace>
        void add(Tspace &spc, InputMap &in) {
          setfront( spc.p.size() );
          int size=0;
          int n=1, npart;
          do {
            std::ostringstream nion("nion"),
              tion("tion"), dpion("dpion"), aion("aion");
            nion << "nion" << n;
            tion << "tion" << n;
            dpion<< "dpion"<< n;
            aion << "aion" << n++; //activity
            npart = in.get(nion.str(), 0);
            if (npart>0) {
              short id = atom[ in.get<string>(tion.str(), "UNK") ].id;
              atom[id].dp = in.get(dpion.str(), atom[id].dp);
              atom[id].activity = in.get(aion.str(), 0.);
              spc.insert( atom[id].name, npart);
              size+=npart;
            } else break;
          } while (npart>0);
          if (size>0)
            resize(size);
          else
            resize(0);
          setMassCenter(spc);
          spc.enroll(*this);
        }
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

