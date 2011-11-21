#ifndef FAUNUS_SPACE_H
#define FAUNUS_SPACE_H
#include <faunus/common.h>
#include <faunus/slump.h>
#include <faunus/group.h>
//extern template class std::vector<Faunus::particle>;

namespace Faunus {


  /*!
   * \brief Contains the particle vector and takes care of particle insertion and deletion.
   * \todo Take care of bonds too?
   *
   * Every simulation must have a Space instance as this contains the particles. While instantiating
   * Space you must provide a valid Geometry which will be stored as a pointer. Typically the Geometrybase
   * is owned by some Energybase.
   * Space also bookkeeps groups in the system and this information can be requested by classes that needs
   * to know about all possible groups -- two such examples are NmuT and NPT simulations that need to
   * insert and re-scale mass centra, for example.
   */
  class Space {
    protected:
      std::ifstream fin;
    private:
      slump slp;
      bool overlap() const;
      bool overlap(const particle&) const;            //!< Check for hardspheres overlap with particle
      bool checkSanity();                             //!< Do a number of checks to see if eveything is OK

    public:
      enum keys {OVERLAP,NOOVERLAP,RESIZE,NORESIZE};
      Geometry::Geometrybase* geo;               //!< Pointer to a valid Geometry (!=nullptr)
      p_vec p;                                   //!< The main particle vector
      p_vec trial;                               //!< Trial particle vector. 
      vector<Group*> g;                          //!< Pointers to all groups in the system (sum must match particle size!)

      Space(Geometry::Geometrybase&);
      virtual ~Space();

      virtual bool save(string);                      //!< Save container state to disk
      virtual bool load(string, keys=NORESIZE);       //!< Load container state from disk

      GroupMolecular insert(const p_vec&,int=-1, keys=NOOVERLAP);
      bool insert(const particle&, int=-1);           //!< Insert particle at pos n (old n will be pushed forward).
      bool insert(string, int, keys=NOOVERLAP); 
      bool erase(int);                               //!< Remove n'th particle
      int enroll(Group&);                             //!< Add group pointer to g vector

      double charge() const;                          //!< Sum all charges in particle vector
      string info();                                  //!< Print information string
  };
} //namespace
#endif
