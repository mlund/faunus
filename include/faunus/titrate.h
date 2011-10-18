#ifndef FAU_TITRATE_H
#define FAU_TITRATE_H

#include <faunus/common.h>
#include <faunus/move.h>
#include <faunus/average.h>
#include <faunus/energy.h>
#include <faunus/move.h>

namespace Faunus {

  namespace Energy {

    class EquilibriumEnergy;

    /*!
     * \brief  Class for implicit titration of species with fixed chemical potential.
     * \author Mikael Lund and Chris Evers
     * \date   Malmo, October 2010
     *
     * Consider the dissociation process AX=A+X. This class will locate
     * all species of type AX and A and make a MC swap move between them.
     * X is implicit, meaning that it enters only with its chemical potential
     * (activity). The titrating species, their dissociation constants
     * and the chemical potential of the titrant are read from an input file with
     * this format:
     * \code 
     * Process type_AX type_A pKd pX
     * \endcode
     * where pKd and pX are the negative logarithms of the dissociation constant
     * and the activity of X, respectively. Make sure the standard states (units) are
     * consistent. For example, for proton titration of the phosphate ion one would
     * use the following input (pH 7):
     * \code
     * Process H3PO4 H2PO4   2.12    7
     * Process H2PO4 HPO4    7.21    7
     * Process HPO4  PO4     12.67   7
     * \endcode
     * All species and their properties must be defined in the faunatoms.dat
     * file before initializing this class.
     */
    class EquilibriumController {
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
          short id_AX, id_A;           //!< Particle id's for AX and A
          bool one_of_us(const short&);//!< Does the particle belong to this process?
          double energy(const short&); //!< Returns intrinsic energy of particle id
          double swap(particle &);     //!< Swap AX<->A and return intrinsic energy change
          void set(double,double);     //!< Set activity of X and the pKd value
          void set_mu_AX(double);      //!< Set chemical potential of species AX - mu_A then follows.
          void set_mu_A(double);       //!< Set chemical potential of species A  - mu_AX then follows.
        };

        std::map<int, Average<double> > q;       //!< Map of average charges per site
        vector<processdata> process;             //!< Vector of processes.

        EquilibriumController(InputMap&, string="eq_");
        bool include(string);                    //!< Read equilibrium processes from file
        void findSites(const p_vec&);            //!< Locate all titratable sites
        double intrinsicEnergy(const short&);    //!< Intrinsic energy of particle id (kT)
        string info(char=25);                    //!< Get information string
        processdata& random(const p_vec&, int&); //!< Random titratable particle and assiciated random process

        vector<int> sites;                       //!< List of titratable sites

        void sampleCharge(const p_vec&);         //!< Updates the average charge vector titrate::q
        double applycharges(p_vec &);            //!< Copy average charges to particles in the particle vector
        double avgcharge(const p_vec&, int&);    //!< Print average charges of process i
    };

    class EquilibriumEnergy : public Energybase {
      protected:
        std::map<short, double> energymap;       //!< Map of intrinsic energy for titratable sites
      public:
        EquilibriumController eq;                //!< Process controller
        EquilibriumEnergy(InputMap&);
        double i_internal(const p_vec&, int);
        double g_internal(const p_vec&, Group&);
        int findSites(const p_vec&);
        string info();
    };

  }//Energy namespace 

  namespace Move {

    class SwapMove : public Movebase {
      private:
        std::map<int, Average<double> > accmap; //!< Site acceptance map
        int ipart;                              //!< Particle to be swapped
      protected:
        void trialMove();
        double energyChange();
        void acceptMove();
        void rejectMove();
        Energy::EquilibriumEnergy eqpot;
      public:
        SwapMove(InputMap&, Energy::Hamiltonian&, Space&, string="swapmv_");
        int findSites(const p_vec&);
        double move();
        double totalEnergy();
        string info();
    };

  }// Move namespace
}//Faunus namespace
#endif
