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
   * \todo Untested, must be supplemented with grand canonical salt
      This is perhaps not the best forum but a detail regarding grand canonical 
      titration should be mentioned. What ever group that is titrating should
      be declared BEFORE any fluctuating salt group. If not, titrate::sites will
      no longer be in sync with con->p.
   */
  class GCchargereg : public markovmove, public titrate_gc {
    public: 
      GCchargereg( grandcanonical&, container&, energybase&, inputfile&);
      double titrateall();
      string info();
    private:
  };

  class DHchargereg : public markovmove, public titrate_implicit {
    public:
      DHchargereg( ensemble&, container&, energybase&, float, float);
      double titrateall();
      string info();
   };

#ifdef BABEL
   class glu3corechargereg : public chargereg {
     public:
       glu3corechargereg(ensemble &, container &, energybase &, inputfile &, group &);
       string info();
       double move(glu3 &);
     private:
       double energy(vector<particle> &, double, titrate::action &, titrate::action &);
       double porphyrinpKa;
       int porph1, porph2;
       average<double> p1, p2;
       int cntcore;
   };
   class GCglu3corechargereg : public GCchargereg {
     public:
       GCglu3corechargereg(grandcanonical &, container &, energybase &, inputfile &);
       string info();
       double move(glu3 &);
     private:
       double energy(vector<particle> &, double, int i);
       double porphyrinpKa;
       int porph1, porph2;
       average<double> p1, p2;
       int cntcore;
       unsigned int o1,o2;
       particle i1,i2;
   };
#endif

  /*!
   * \brief Implicit ions titration scheme
   * \author Andre Teixeira
   * \date Jul 2009
   * \todo Use templates in order to work with different potential classes.
   */
  class ATchargereg : public markovmove, public titrate_implicit {
    public:
      ATchargereg( ensemble&, container&, energybase&, float, inputfile& , pot_debyehuckel& );
      double titrateall( vector<macromolecule>& );
      string info();
      
    private:
      pot_debyehuckel* pairpot;
      double const_kappa;
      double knew , kold;
	  double ionic_str0 , ionic_str1;
	  double protein_conc;
	  inline double calc_kappa( vector<macromolecule>& , vector<particle>& );
	  inline double prot_ion_u( vector<macromolecule>& , vector<particle>& );
	  inline double prot_ion_u( vector<macromolecule>& , vector<particle>& , double& );
	  inline int who_is_tit( vector<macromolecule>& , int& );
	  inline void set_kappa( double& );
  };

}//namespace
#endif
