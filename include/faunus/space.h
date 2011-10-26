#ifndef FAUNUS_SPACE_H
#define FAUNUS_SPACE_H
#include <faunus/common.h>
#include <faunus/slump.h>
#include <faunus/group.h>
//extern template class std::vector<Faunus::particle>;

namespace Faunus {
  class Space {
    protected:
      std::ifstream fin;
    private:
      slump slp;
    public:
      enum spacekeys {NOOVERLAP};
      Geometry::Geometrybase* geo;
      p_vec p;                                   //!< The main particle vector
      p_vec trial;                               //!< Trial particle vector. 
      vector<Group*> g;                          //!< Pointers to all groups in the system.

      Space(Geometry::Geometrybase&);
      virtual ~Space();

      virtual bool save(string);                      //!< Save container state to disk
      virtual bool load(string, bool=false);          //!< Load container state from disk

      Group insert(const p_vec&, int=-1);
      bool insert(particle, int=-1);                  //!< Insert particle at pos n.
      bool insert(string, int, spacekeys=NOOVERLAP); 
      bool remove(int);                               //!< Remove particle n.
      int enroll(Group&);                             //!< Add group pointer to g vector

      bool overlap(const particle&) const;            //!< Check for hardspheres overlap with particle

      double charge() const;                          //!< Sum all charges in particle vector
      bool check_vector();                            //!< Check if p and trial are equal!
      string info();
  };
} //namespace
#endif
