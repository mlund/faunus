#ifndef FAU_SIMBOTTLE_H
#define FAU_SIMBOTTLE_H

#include "faunus/inputfile.h"
#include "faunus/energy/base.h"
#include "faunus/ensemble.h"
#include "faunus/analysis.h"
#include "faunus/mcloop.h"
#include "faunus/io.h"

namespace Faunus {

  /*!
   * \brief Base class for entire simulation systems
   * \author Mikael Lund
   * \date Malmo, 2010
   */
  class simBottle {
    protected:
      iopqr pqr;
      ioaam aam;
      inputfile in;
      mcloop loop;

    public:
      systemenergy usys;
      container* cPtr;                  //!< Pointer to the container for this system
      string name;                      //!< Name of the system (arbitrary)
      string prefix;                    //!< Prefix for input files, output etc.
      virtual void prepare()=0;         //!< Prepare and setup systems
      virtual void microloop()=0;       //!< Things to do in each micro move
      virtual void macroloop()=0;       //!< Things to do in each macro move
      virtual void save()=0;            //!< Save data to disk
      virtual double systemEnergy()=0;  //!< Calculate system energy (of TRIAL particles!)
      virtual string preinfo()=0;       //!< Give initial information string
      virtual string postinfo()=0;      //!< Give final information string
      simBottle(string pfx) : in(pfx+".conf"), loop(in) {
        prefix=pfx;
      }
  };

}//namespace
#endif
