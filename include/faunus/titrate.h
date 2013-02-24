#ifndef FAU_TITRATE_H
#define FAU_TITRATE_H

#ifndef SWIG
#include <faunus/common.h>
#include <faunus/move.h>
#include <faunus/average.h>
#endif

namespace Faunus {

  namespace Energy {

    class EquilibriumEnergy;

    /**
     * @brief  Class for implicit titration of species with fixed chemical potential.
     *
     * Consider the dissociation process AX=A+X. This class will locate
     * all species of type AX and A and make a MC swap move between them.
     * X is implicit, meaning that it enters only with its chemical potential
     * (activity). The titrating species, their dissociation constants
     * and the chemical potential of the titrant are read from a
     * `processes` JSON object.
     * For example, for proton titration of phosphate one would
     * use the following JSON input (pH 7.0):
     *
     *     {
     *       "processes" :
     *       {
     *         "K1" : { "bound":"H3PO4", "free":"H2PO4", "pKd":2.12,  "pX":7.0 }, 
     *         "K2" : { "bound":"H2PO4", "free":"HPO4",  "pKd":7.21,  "pX":7.0 }, 
     *         "K3" : { "bound":"HPO4",  "free":"PO4",   "pKd":12.67, "pX":7.0 }
     *       }
     *     }
     *
     * All species and their properties must be defined in `AtomMap` before
     * initializing this class.
     *
     * @date Malmo, October 2010
     * @author Mikael Lund and Chris Evers
     */
    class EquilibriumController {
      private:
        bool includeJSON(const string&); //!< Read equilibrium processes 
      public:
        friend class EquilibriumEnergy;
        class processdata {
          friend class EquilibriumController;
          private:
          double mu_AX;                //!< chemical potential of AX
          double mu_A;                 //!< chemical potential of A
          double mu_X;                 //!< chemical potential of X (this is the titrant)
          double ddG;                  //!< ddG = mu_A + mu_X - mu_AX
          int cnt;                     //!< number of sites for this process
          public:
          particle::Tid id_AX, id_A;   //!< Particle id's for AX and A
          bool one_of_us(const particle::Tid&);//!< Does the particle belong to this process?
          double energy(const particle::Tid&); //!< Returns intrinsic energy of particle id
          double swap(particle &);     //!< Swap AX<->A and return intrinsic energy change
          void set(double,double);     //!< Set activity of X and the pKd value
          void set_mu_AX(double);      //!< Set chemical potential of species AX - mu_A then follows.
          void set_mu_A(double);       //!< Set chemical potential of species A  - mu_AX then follows.
        };

        std::map<int, Average<double> > q;       //!< Map of average charges per site
        std::vector<processdata> process;        //!< Vector of processes.

        EquilibriumController(InputMap&, string="eq_");
        bool include(string);                    //!< Read equilibrium processes
        void findSites(const p_vec&);            //!< Locate all titratable sites
        double intrinsicEnergy(const particle::Tid&);    //!< Intrinsic energy of particle id (kT)
        string info(char=25);                    //!< Get information string
        processdata& random(const p_vec&, int&); //!< Random titratable particle and assiciated random process

        std::vector<int> sites;                  //!< List of titratable sites

        void sampleCharge(const p_vec&);         //!< Updates the average charge vector titrate::q
        double applycharges(p_vec &);            //!< Copy average charges to particles in the particle vector
        double avgcharge(const p_vec&, int&);    //!< Print average charges of process i
    };

    /**
     * @brief Energy class for implicit titration of species
     *        used with Move::SwapMove.
     *
     *  This is a Hamiltonian for swapping atomic species according
     *  to their chemical potential and equilibrium constant as
     *  explained in `EquilibriumController`.
     */
    class EquilibriumEnergy : public Energybase {
      private:
        string _info();
      protected:
        std::map<particle::Tid, double> energymap;//!< Intrinsic site energy
      public:
        EquilibriumController eq;                 //!< Process controller
        EquilibriumEnergy(InputMap&);
        double i_internal(const p_vec&, int);
        double g_internal(const p_vec&, Group&);
        int findSites(const p_vec&);
    };

  }//Energy namespace 

  namespace Move {

    /**
     * @brief Move for swapping species types - i.e. implicit titration
     *
     * Upon construction this class will add an instance of
     * Energy::EquilibriumEnergy to the Hamiltonian. For details
     * about the titration procedure see Energy::EquilibriumController.
     */
    class SwapMove : public Movebase {
      private:
        std::map<int, Average<double> > accmap; //!< Site acceptance map
        string _info();
        void _trialMove();
        void _acceptMove();
        void _rejectMove();
      protected:
        double _energyChange();
        int ipart;                              //!< Particle to be swapped
        Energy::EquilibriumEnergy eqpot;
      public:
        SwapMove(InputMap&, Energy::Hamiltonian&, Space&, string="swapmv_"); //!< Constructor
        int findSites(const p_vec&); //!< Search for titratable sites (old ones are discarded)
        double move();
        void applycharges(p_vec &); 
    };

    /**
     * @brief As SwapMove but Minimizes Short Ranged interactions
     *        within a molecule upon swapping
     *
     * Before calculating dU of an attempted swap move, radii on
     * particles within the SAME group are set to minus radius of
     * the swapped particle and hydrophobicity is set to false.
     * This to minimize large interactions in molecules with overlapping
     * particles - i.e LJ will be zero. It can also be used to avoid
     * internal hydrophobic interactions in rigid groups upon swapping
     * between hydrophobic and non-hydrophobic species.
     * Group information is found in Space::g and to avoid energy drifts by
     * ignoring hydrophobic interactions internally in groups the
     * Energy::EnergyRest is used to collect the missing contribution to dU.
     */
    class SwapMoveMSR : public SwapMove {
      private:
        std::map<int, double> radiusbak;    // backup for radii
        std::map<int, bool> hydrophobicbak; // backup for hydrophobic state
        Energy::EnergyRest potrest; // dummy energy class missing contributions to dU in energy drift calc.
        double _energyChange();
        void modify();
        void restore();
      public:
        SwapMoveMSR(InputMap&, Energy::Hamiltonian&, Space&, string="swapmv_");
    };

  }// Move namespace
}//Faunus namespace
#endif
