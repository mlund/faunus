#ifndef FAUNUS_SPACE_H
#define FAUNUS_SPACE_H

#ifndef SWIG
#include <faunus/common.h>
#include <faunus/slump.h>
#endif

namespace Faunus {

  /**
   * @brief Place holder for particles and groups
   *
   * Every simulation must have a `Space` instance as this contains
   * the particles and information about groups (particle ranges).
   * `Space` must also be given a valid reference to a `Geomebrybase`.
   */
  class Space {
    protected:
      std::ifstream fin;
    private:
      slump slp;
      bool overlap_container() const;
      bool overlap() const;
      bool overlap(const particle&) const;   //!< Check hardspheres overlap with particle
      bool checkSanity();                    //!< Check group length and vector sync
      std::vector<Group*> g;                 //!< Pointers to ALL groups in the system

    public:
      enum keys {OVERLAP,NOOVERLAP,RESIZE,NORESIZE};
      Geometry::Geometrybase* geo;               //!< Pointer to a geometry
      p_vec p;                                   //!< Main particle vector
      p_vec trial;                               //!< Trial particle vector. 
      std::vector<Group*>& groupList();          //!< Vector with pointers to all groups

      Space(Geometry::Geometrybase&);
      virtual ~Space();

      virtual bool save(string);                      //!< Save container state to disk
      virtual bool load(string, keys=NORESIZE);       //!< Load container state from disk

      GroupMolecular insert(const p_vec&, int=-1);
      bool insert(const particle&, int=-1);           //!< Insert particle at pos n (old n will be pushed forward).
      bool insert(string, int, keys=NOOVERLAP); 
      bool erase(int);                                //!< Remove n'th particle
      int enroll(Group&);                             //!< Store group pointer
      void reserve(int);                              //!< Reserve space for particles for better memory efficiency

      double charge() const;                          //!< Sum all charges
      string info();                                  //!< Information string
      void displace(const Point&);                    //!< Displace system by a vector
  };
} //namespace
#endif
