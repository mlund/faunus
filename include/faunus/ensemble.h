#ifndef FAU_SLUMP_H
#define FAU_SLUMP_H

#include <faunus/common.h>
#include <faunus/slump.h>

namespace Faunus {
  class container;
  class group;
  class particle;

  /*!
   * Polymorph ensemble class
   * \todo Maybe not really needed...
   */
  class ensemble {
    protected:
      slump slp;
    public:
      ensemble();
      bool metropolis(double);
  };

  /*! \brief NVT ensemble
  */
  class canonical : public ensemble { public: canonical(); };

  /*! \brief NPT ensemble
   *  \todo Silly name due to name collision...
   */
  class isobarical : public ensemble { public: isobarical(); };

  /*!
   * Class for handling Grand Canonical steps. This class can dynamically insert and
   * erase particles while at the same time update the groups in the systems. 
   * \author Mikael Lund, Bjorn Persson
   * \date Lund, 2009
   * \todo A lot.
   */
  class grandcanonical : public ensemble {
    private:
      int           size();  //!< Calculates the number of particles in the specifies groups
      short         findgroup(unsigned int);
      bool          meMbeR(unsigned int&, unsigned int&);
    public:
      vector<group*> gp;                           //!< Contains all groups that are effected by changing size of con.p/trial
      void          addGroupPtr(group &); 
      short         findgroup(string);
      bool          insert(container &, particle &); //!< Insert a particle in group w. name 'particle'
      bool          erase(container &, unsigned int &); //!< Erase --//--     (overload or use virtual?)
      bool          load(container &, string);      //!< Regenerate groups to fit previus configuration
      string        print();                        //!< To print info of last configuration
      string        info();
      int           gcd(int,int);
      
  };
}//namespace
#endif
