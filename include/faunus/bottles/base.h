#ifndef FAU_BOTTLE_H
#define FAU_BOTTLE_H

#include "faunus/common.h"
#include "faunus/inputfile.h"
#include "faunus/energy/base.h"
#include "faunus/ensemble.h"
#include "faunus/analysis.h"
#include "faunus/io.h"

namespace Faunus {

  /*!
   * \brief Base class for entire simulation systems
   * \author Mikael Lund
   * \date Malmo, 2010
   *
   * Being entire simulations wrapped in a class, BOTTLES
   * are used to run many simulations simultaneously - replica
   * exchange / parallel tempering, Gibb's ensemble or merely
   * parallellization are some of the uses. The idea is that
   * many instances of the bottle can be created and the calling
   * program will then take care of the micro and macro loops.
   *
   * Things to remember:
   * \li The constructor should set up and prepare the simulation system.
   *     All input and output must be prefixed by the prefix string.
   * \li Output streamed to fout will be saved in "prefix.out"
   * \li It can be useful to design the bottle as a class template
   *     that takes container and interactions types as input.
   */
  class bottle {
    protected:
      iopqr pqr;
      ioaam aam;
      inputfile in;

    public:
      std::ofstream fout;               //!< Output file stream - place all file output here.
      double T;                         //!< Temperature [K] (currently unused))
      double P;                         //!< Pressure over kT
      systemenergy usys;                //!< Class for tracking total system energy
      container* cPtr;                  //!< Pointer to the container for this system
      energybase* pPtr;                 //!< Pointer to interaction class
      string name;                      //!< Name of the system (arbitrary)
      string prefix;                    //!< Prefix for input files, output etc.
      virtual void microloop(int=1)=0;  //!< Things to do (n times) in each micro move
      virtual void macroloop()=0;       //!< Things to do in each macro move
      virtual void save()=0;            //!< Save data to disk (usually called by macroloop())
      virtual void finish()=0;          //!< Clean up and terminate simulation (write output etc.)
      virtual double systemEnergy()=0;  //!< Calculate system energy (of TRIAL particles!)

      bottle(string pfx, double betaP=0) : in(pfx+".conf") {
        prefix=pfx;
        P=betaP;
        fout.open((prefix+".out").c_str());
      }

      ~bottle() {
        fout.close();
      }
  };

}//namespace
#endif
