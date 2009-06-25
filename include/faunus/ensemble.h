#ifndef FAU_SLUMP_H
#define FAU_SLUMP_H

#include <faunus/common.h>
#include <faunus/slump.h>
#include <faunus/group.h>
#include <faunus/container.h>

namespace Faunus {
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
   * \author Mikael Lund
   * \date Lund, 2009
   * \todo A lot.
   */
  class grandcanonical : public ensemble {
    private:
      int size();  //!< Calculates the number of particles in the specifies groups
      short findgroup(unsigned int);
    public:
      vector<group> g;      //!< Groups with fluctuating species
      void addgroup(group &); 
      short findgroup(string);
      void searchsalt(container &, salt &);
      bool insert(vector<particle> &, unsigned int, particle, unsigned char=1);
      bool erase(vector<particle> &, unsigned int, unsigned char=1);
      string info();
  };
}//namespace
#endif
