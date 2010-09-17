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
      double T;                         //!< Temperature [K]
      double P;                         //!< Pressure over kT
      systemenergy usys;                //!< Class for tracking total system energy
      container* cPtr;                  //!< Pointer to the container for this system
      energybase* pPtr;                 //!< Pointer to interaction class
      string name;                      //!< Name of the system (arbitrary)
      string prefix;                    //!< Prefix for input files, output etc.
      virtual void prepare()=0;         //!< Prepare and setup systems
      virtual void microloop(int=1)=0;  //!< Things to do (n times) in each micro move
      virtual void macroloop()=0;       //!< Things to do in each macro move
      virtual void save()=0;            //!< Save data to disk
      virtual double systemEnergy()=0;  //!< Calculate system energy (of TRIAL particles!)
      virtual string preinfo()=0;       //!< Give initial information string
      virtual string postinfo()=0;      //!< Give final information string
      simBottle(string pfx, double betaP=0) : in(pfx+".conf"), loop(in) {
        prefix=pfx;
        P=betaP;
      }
  };

}//namespace
#endif
