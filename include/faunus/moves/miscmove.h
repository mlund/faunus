#ifndef FAU_ROTTRANSLATIONAL_H
#define FAU_ROTTRANSLATIONAL_H
#include "faunus/moves/base.h"
namespace Faunus {
  /*! \brief Rotate group around cm plus dr (effectively a combined rotation and translation).
   *  \author Bjorn Persson
   *  \date Lund 2008
   */
  class transrot : public markovmove { 
    public:
      transrot( ensemble&, container&, energybase&);
      double move(macromolecule &);                  //!< Translate/rotate one molecule
      double move(vector<macromolecule> &, int);     //!< Translate/rotate one molecule in vector of molecules
      double dpt;                                    //!< Translational displacement parameter
      double dpr;                                    //!< Rotational displacement parameter
      string info();
      
  };
  /*! \brief Rotate group around cm plus dr (effectively a combined rotation and translation).
   *  \author Bjorn Persson
   *  \date Lund 2008
   */
  class multtr : public markovmove { 
    public:
      multtr( ensemble&, container&, energybase&, int);
      double move(molecules &, vector<int> &);
      double dpt, dpr;      // step parameters
      string info();
      int m;
      
  };
}//namespace
#endif
