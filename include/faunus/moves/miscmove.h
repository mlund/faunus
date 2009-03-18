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
      double move(macromolecule &);
      double dpt, dpr;      // step parameters
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
