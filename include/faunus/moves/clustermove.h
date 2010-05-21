#ifndef FAU_CLUSTERMOVE_H
#define FAU_CLUSTERMOVE_H

#include "faunus/moves/base.h"
#include "faunus/hardsphere.h"
#include "faunus/energy/reactionfield.h"
#include "faunus/potentials/pot_test.h"

namespace Faunus {
  /*!
   * This type of move will attempt to move collective sets of macromolecules that
   * obeys some criteria (here a hardcore overlap) with a symmetric transition
   * matrix (no flow through the clusters).
   * \author Bjoern Persson
   * \brief Quick and dirty translational cluster move
   * \warning This rutine is only compatible for systems containing a set of macromolecules!
   *
   *
  class clustertranslate : public markovmove {
    public: 
      clustertranslate( ensemble&, container&, energybase&, double);
      double move(vector<macromolecule> &); 
      vector<int> cluster;    //! Index vector for macromolecule number in cluster
      vector<int> free;       //! Index vector for 'free' macromolecules
      vector<int> nacc;
      vector<int> ntrial ;
      int np;                 //! Number of particles in config.
      bool decreasing(int, int);
      void sort(vector<macromolecule> &);  //function to sort macromolecules in to two classes, cluster and free
      void flowcheck(vector<macromolecule> &); //ensure detailed balance
      hardsphere coll;
      double sep;             // separation parameter
      int FLOW;               // control variable
      string info();
  };
*/
 /*!
    * This type of move will attempt to move collective sets of macromolecules that
    * obeys some criteria (here a hardcore overlap) with a symmetric transition
    * matrix (no flow through the clusters).
    * \author Bjoern Persson
    * \brief Quick and dirty rotational cluster move
    * \warning This rutine is only compatible for systems containing a set of macromolecules!
    *
    */
    class clusterrotate : public markovmove {
      public: 
        clusterrotate( ensemble&, container&, energybase&);
        double move(vector<macromolecule> &); 
        vector<int> cluster;    //! Index vector for macromolecule number in cluster
        vector<int> free;       //! Index vector for 'free' macromolecules
        vector<int> nacc;       //! Local acceptance counter
        vector<int> ntrial ;    //! Local trial counter
        int np;                                  //! Number of particles in config.
        bool decreasing(int, int);               
        void sort(vector<macromolecule> &);      //! function to sort macromolecules in to two classes, cluster and free, from a random seed
        void flowcheck(vector<macromolecule> &); //! ensure detailed balance, no flux through the cluster
        hardsphere coll;                         //! used to define clusters
        double sep, angle;                       //! separation parameter(defines the clustering distance) and angle to rotate around
        int FLOW;                                //! check parameter
        string info();
        void print();                            //! Debug-info
    };

  // Non-rejective cluster inversion. Minor testing suggest that 
  // the algorithm is ok.

  class clusterinvw : public markovmove {
    public :
      clusterinvw( ensemble&, container&, sphericalimage<pot_test>&);
      double move(molecules &);
      string info();
      average< double > movefrac;
      vector<int> moved;
      vector<int> remaining;
      sphericalimage<pot_test> *ipot;
  };

  // Non-rejective cluster inversion on a partial volume and meteropolis weigth
  // with the rest. Untested and undone.

  class clusterrinvw : public markovmove {
    public :
      clusterrinvw( ensemble&, container&, sphericalimage<pot_test>&,
          double&, double&);
      double r;                 // restricted radius
      double cr,cr2;            // cavity radius
      double move(molecules &);
      string info();
      average< double > movefrac;
      vector<int> moved;
      vector<int> remaining;
      sphericalimage<pot_test> *ipot;
  };

  /*!
   * \brief Non-rejective cluster translation.
   * \author Bjoern Persson
   * \date Lund 2009-2010
   * \brief Quick and dirty translational cluster move
   * \warning This rutine is only compatible for systems containing a set of macromolecules!
   *
   * This type of move will attempt to move collective sets of macromolecules that
   * obeys some criteria (here a hardcore overlap) with a symmetric transition
   * matrix (no flow through the clusters).
   *
   * Setting the boolen "skipEnergyUpdate" to true (default is false) updates of the
   * total energy are skipped to speed up the move.
   * While this has no influence on the Markov chain it will cause an apparent energy
   * drift. It is recommended that this is enabled only for long production runs after
   * having properly checked that no drifts occur with skipEnergyUpdate=false.
   */
  class clustertrans : public markovmove {
    private:
      vector<int> moved;
      vector<int> remaining;
      //vector<macromolecule> *g;
      //distributions dist;
    public:
      clustertrans( ensemble&, container&, energybase&, vector<macromolecule>&);
      double move(vector<macromolecule> &);
      string info();
      average< double > movefrac; //!< Fraction of particles moved
      void check(checkValue &);   //!< Unit testing
      bool skipEnergyUpdate;      //!< True if energy updates should be skipped (faster evaluation!)
  };
}//namespace
#endif

