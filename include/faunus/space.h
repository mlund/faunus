#ifndef FAUNUS_SPACE_H
#define FAUNUS_SPACE_H
#include <faunus/common.h>

namespace Faunus {
  class space {
    protected:
      std::ifstream fin;
    private:
      slump slp;
    public:
      Geometry::geometrybase* geo;
      vector<particle> p;                             //!< The main particle vector
      vector<particle> trial;                         //!< Trial particle vector. 
      vector<group*> g;                               //!< Pointers to all groups in the system.

      space(Geometry::geometrybase&);
      virtual bool save(string);                      //!< Save container state to disk
      virtual bool load(string, bool=false);          //!< Load container state from disk

      bool insert(particle, unsigned int=-1);         //!< Insert particle at pos n.
      bool remove(unsigned int);                      //!< Remove particle n.

      double charge() const;                          //!< Sum all charges in particle vector
      bool check_vector();                            //!< Check if p and trial are equal!
      string info();
  };
} //namespace
#endif
