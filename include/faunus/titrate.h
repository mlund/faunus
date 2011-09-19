#ifndef FAU_EQUILIBRIA_H
#define FAU_EQUILIBRIA_H

#include <faunus/common.h>
#include <faunus/move.h>
#include <faunus/average.h>

namespace Faunus {

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
  class equilibria {
    private:
      class data {
        friend class equilibria;
        private:
        double mu_AX;                //!< chemical potential of AX
        double mu_A;                 //!< chemical potential of A
        double mu_X;                 //!< chemical potential of X (this is the titrant)
        double ddG;                  //!< ddG = mu_A + mu_X - mu_AX
        int cnt;                     //!< number of sites for this process
        public:
        char id_AX, id_A;            //!< Particle id's for AX and A
        bool one_of_us(const short&);  //!< Does the particle belong to this process?
        double energy(const short&);   //!< Returns intrinsic energy of particle
        double swap(particle &);     //!< Swap AX<->A and return intrinsic energy change
        void set(double,double);     //!< Set activity of X and the pKd value
        void set_mu_AX(double);      //!< Set chemical potential of species AX - mu_A then follows.
        void set_mu_A(double);       //!< Set chemical potential of species A  - mu_AX then follows.
      };

      vector<average <double> > q;         //!< List of average charges per site
      vector<data> process;                //!< Vector of eq. processes.

      bool load(string);                   //!< Read equilibrium processes from file
      void findSites(p_vec &);  //!< Locate all titratable sites
      double intrinsicenergy(const short&);   //!< Intrinsic energy of particle

    public:
      vector<int> sites;                           //!< List of titratable sites
      equilibria(space&, inputfile&, string="eqtit_");
      string info();                               //!< Get information string
      string info(p_vec&);              //!< Get extended information string
      double intrinsicenergy(p_vec&);   //!< Intrinsic energy of all titratable sites
      void samplesites(p_vec &);        //!< Updates the average charge vector titrate::q
      double applycharges(p_vec &);     //!< Copy average charges to particles in the particle vector
      double avgcharge(p_vec&, int&);   //!< Print average charges of process i
  };

}//namespace

#endif
