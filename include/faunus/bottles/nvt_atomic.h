#ifndef FAU_NVT_ATOMIC_H
#define FAU_NVT_ATOMIC_H

#include "faunus/bottles/base.h"

namespace Faunus {

  /*!
   * \brief NVT simulation system for atomic species
   * \author Mikael Lund
   * \date Lund 2010
   */
  template<typename Tcon, typename Tpot> class nvt_atomic : public simBottle {
    private:
      string dumpfile;  //!< Filename of container dumpfile
    protected:
      canonical nvt;
      io fileio;
      saltmove *sm;
    public:
      Tcon con; //!< Container
      Tpot pot; //!< Interaction scheme (based on energybase)

      nvt_atomic(string pfx) : simBottle(pfx), con(in), pot(in) {
        cPtr=&con;
        pPtr=&pot;
        dumpfile=prefix+"confout.dump";
      }

      void prepare() {}
      double systemEnergy() {}
      string preinfo() {}
      string postinfo() {}
      void microloop() {}
      void macroloop() {}
      void save() {}
  };

} // end of namespace
#endif
