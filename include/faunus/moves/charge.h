#ifndef FAU_CHARGEREG_H
#define FAU_CHARGEREG_H

#include <faunus/moves/base.h>
#include <faunus/titrate.h>

namespace Faunus {
  /*! \brief Titrate all titrateable sites
   *  \author Mikael Lund
   *  \note "titrate" used to be private. Changed because iopqr::save()
   */
  class chargereg : public markovmove, public titrate {
    public:
      chargereg( ensemble&, container&, energybase&, group&, float);
      double titrateall();
      string info();
  };

  /*!
   * \brief Grand Canonical titation of all sites
   * \author Bjoern Persson
   * \todo Untested, in principle it must be supplemented with grand canonical salt
   */
  class HAchargereg : public chargereg {
    public: 
      HAchargereg( ensemble&, container&, energybase&, group&, float, float);
      string info();
    private:
      double energy(vector<particle> &, double, titrate::action &); //!< New titrate energy function
      double CatPot;  //!< Chemical potential of coupled cation
  };

  class DHchargereg : public markovmove, public titrate_implicit {
    public:
      DHchargereg( ensemble&, container&, energybase&, float, float);
      double titrateall();
      string info();
   };

  /*!
   * \brief (Implicit titration scheme)
   * \author (Andre Teixeira)
   */
  class ATchargereg : public markovmove, public titrate_implicit {
    public:
      ATchargereg( ensemble&, container&, energybase&, float, float);
      double titrateall();
      string info();
  };

}//namespace
#endif
