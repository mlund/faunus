#ifndef FAU_POTENTIAL_H
#define FAU_POTENTIAL_H

#include <faunus/common.h>
#include <faunus/point.h>
#include <faunus/geometry.h>
#include <faunus/textio.h>
#include <faunus/inputfile.h>
#include <faunus/physconst.h>

namespace Faunus {

  /*!
   * \brief Namespace for various potentials - pair potentials, external potentials etc.
   */
  namespace Potential {

    class DebyeHuckel;

    using namespace Faunus::textio;

    /*!
     * \brief Base class for pair potential classes
     *
     * This is a base class for pair potentials which must implement the function operator so
     * that the potential can work as a class function.
     */
    class PairPotentialBase {
      private:
        virtual string _brief()=0;
        virtual void _setScale(double);
      protected:
        double _tokT;
      public:  
        PairPotentialBase();
        string name;             //!< Short (preferably one-word) description of the core potential
        string brief();          //!< Brief, one-lined information string
        void setScale(double=1); //!< Set scaling factor
        double tokT();           //!< Convert returned energy to kT.

        /*!
         * \brief Particle-particle energy divided by tokT()
         * \param a First particle
         * \param b Second particle
         * \param r2 Squared distance between them (angstrom squared)
         */
        virtual double operator() (const particle &a, const particle&b, double r2) const=0;
    };

    /*!
     * \brief Harmonic pair potential
     *
     * The harmonic potential has the form \f$ \beta u_{ij} = k(r_{ij}-r_{eq})^2 \f$ where k is the force constant
     * (kT/angstrom^2) and req is the equilibrium distance (angstrom).
     */
    class Harmonic : public PairPotentialBase {
      private:
        string _brief();
        void _setScale(double);
      public:
        double k;   //!< Force constant (kT/angstrom squared)
        double req; //!< Equilibrium distance (angstrom)
        Harmonic(double=0, double=0);
        double operator() (const particle&, const particle&, double) const; //!< Pair interaction energy (kT)
    };

    /*!
     * \brief Hard sphere pair potential
     */
    class HardSphere : public PairPotentialBase {
      private:
        string _brief();
      public:
        HardSphere();
        HardSphere(InputMap&);
        inline double operator() (const particle &a, const particle &b, double r2) const {
          double mindist=a.radius+b.radius;
          if (r2<mindist*mindist)
            return pc::infty;
          return 0;
        }
        string info(char w);
    };

    /*!
     * \brief Lennard-Jones (12-6) pair potential
     *
     * The Lennard-Jones potential has the form:
     *
     * \f$ \beta u = 4\epsilon_{lj} \left (  (\sigma_{ij}/r_{ij})^{12} - (\sigma_{ij}/r_{ij})^6    \right ) \f$
     *
     * where \f$\sigma_{ij} = (\sigma_i+\sigma_j)/2\f$ and \f$\epsilon_{lj}\f$ is fixed for this class.
     */
    class LennardJones : public PairPotentialBase {
      private:
        string _brief();
        void _setScale(double);
      protected:
        double eps;
      public:
        LennardJones();
        LennardJones(InputMap&);
        inline double r6(double sigma, double r2) const {
          double x=sigma*sigma/r2;  // 2
          return x*x*x;             // 6
        }
        inline double r12(double sigma, double r2) const {
          double x=r6(sigma,r2);
          return x*x;               // 12
        }
        inline double energy(double sigma, double r2) const {
          double x=r6(sigma,r2);
          return eps*(x*x - x);
        }
        inline double operator() (const particle &a, const particle &b, double r2) const {
          return energy(a.radius+b.radius, r2);
        }
        string info(char);
    };

    /*!
     * \brief Square well pair potential
     */
    class SquareWell : public PairPotentialBase {
      private:
        string _brief();
        void _setScale(double);
      public:
        double threshold;                           //!< Threshold between particle *surface* [A]
        double depth;                               //!< Energy depth [kT]
        SquareWell(InputMap&, string="squarewell"); //!< Constructor
        inline double operator() (const particle &a, const particle &b, double r2) const {
          if ( sqrt(r2)-a.radius-b.radius<threshold )
            return -depth;
          return 0;
        }
        string info(char);
    };

    class SquareWellHydrophobic : public SquareWell {
      public:
        SquareWellHydrophobic(InputMap&, string="squarewell"); //!< Constructor
        inline double operator() (const particle &a, const particle &b, double r2) const {
          if (a.hydrophobic)
            if (b.hydrophobic)
              return SquareWell::operator()(a,b,r2);
          return 0;
        }
    };
    
    /*!
     * \brief Soft repulsion of the form \f$ \beta u = \sigma^6 / (r_{ij}-r_i-r_j)^6
     * \todo This applies sqrt() and thus may be slow
     */
    class SoftRepulsion : public PairPotentialBase {
      private:
        string _brief();
        void _setScale(double);
        double sigma6;
      public:
        SoftRepulsion(InputMap&);
        inline double operator() (const particle &a, const particle &b, double r2) const {
          double x=sqrt(r2)-a.radius-b.radius;
          x=x*x; //2
          return sigma6 / (x*x*x);
        }
        string info(char);
    };

    /*!
     * \brief Coulomb pair potential
     * 
     * The Coulomb potential has the form:
     * \f[ \beta u_{ij} = \frac{e^2}{4\pi\epsilon_0\epsilon_rkT} \frac{z_i z_j}{r_{ij}} \f]
     * where the first factor is the Bjerrum length, \f$l_B\f$. By default the scaling is set
     * to the Bjerrum length - i.e. returned energies are divided by \f$l_B\f$ and to get the energy
     * in kT simply multiply with tokT(). This to increase performance when summing over many
     * particles.
     */
    class Coulomb : public PairPotentialBase {
      friend class DebyeHuckel;
      private:
      string _brief();
      void _setScale(double);
      double temp, epsilon_r;
      protected:
      double lB;          //!< Bjerrum length (angstrom)
      public:
      Coulomb(InputMap&); //!< Construction from InputMap

      /*!
       * \brief Particle-particle energy
       * \param zz Charge number product i.e. \f$z_az_b = q_aq_b/e^2\f$
       * \param r Distance between charges (angstrom) 
       * \returns \f$\beta u/l_B\f$
       */
      inline double energy(double zz, double r) const {
        return zz/r;
      }

      /*! \returns \f$\beta u/l_B\f$ */
      inline double operator() (const particle &a, const particle &b, double r2) const {
        return energy( a.charge*b.charge, sqrt(r2) );
      }
      string info(char);
    };

    /*!
     * \brief Debye-Huckel pair potential
     *
     * This potential is similar to the plain Coulomb potential but with an extra exponential term to described salt screening:
     * \f[ \beta u_{ij} = \frac{e^2}{4\pi\epsilon_0\epsilon_rkT} \frac{z_i z_j}{r_{ij}} \exp(-\kappa r_{ij}) \f]
     * where \f$\kappa=1/D\f$ is the inverse Debye screening length.
     */
    class DebyeHuckel : public Coulomb {
      private:
        string _brief();
      protected:
        double c,k;
      public:
        DebyeHuckel(InputMap&);       //!< Construction from InputMap
        double ionicStrength() const; //!< Returns the ionic strength (mol/l)
        double debyeLength() const;   //!< Returns the Debye screening length (angstrom)
        inline double energy(double zz, double r) const { return zz/r * exp(-k*r); }
        inline double operator() (const particle &a, const particle &b, double r2) const {
          r2=sqrt(r2);
          return a.charge*b.charge/r2 * exp(-k*r2); //note that k=-kappa
        }
        string info(char);
    };

    /*!
     * \brief Combines two PairPotentialBases
     * \date Lund, 2012
     * \author Mikael Lund
     *
     * This simply combines two PairPotentialBases -- typically short-ranged such as
     * Lennard-Jones and SquareWell, for example. This can then be mixed with electrostatics
     * using the CoulombSR template
     * \code
     *   using namespace Potential;
     *   typedef CombinedPairPotential< LennardJones, SquareWell > srpot;
     * \endcode
     */
    template<class T1, class T2>
      class CombinedPairPotential : public PairPotentialBase {
        protected:
          T1 sr1;
          T2 sr2;
        private:
          string _brief() { return name; };
          void _setScale(double s) {
            sr1.setScale(s);
            sr2.setScale(s);
          }
        public:
          CombinedPairPotential(InputMap &in) : sr1(in), sr2(in) {
            name=sr1.name+"+"+sr2.name;
          }
          inline double operator() (const particle &a, const particle &b, double r2) const {
            return sr1(a,b,r2) + sr2(a,b,r2);
          }
          string info(char w=20) {
            std::ostringstream o;
            o << sr1.info(w) << sr2.info(w);
            return o.str();
          }
      };

    /*!
     * \brief Combined electrostatic/short ranged pair potential
     *
     * PairPotentialBase classes do not need to implement distance calculation functions
     * as the operator() takes the squared distance as input. For nonbonded interactions we
     * typically want to call functions of the type energy(particle1, particle2) and
     * therefore need to calculate distances explicitly. This is handled by templates
     * that explicitly handles the simulation Geometry. Note that the surrounding energy
     * loops in the Energy namespace should point to this Geometry so as to handle volume
     * fluctuations etc.
     */
    template<class Tgeometry, class Tcoulomb=Coulomb, class Tshortranged=LennardJones>
      class CoulombSR : public PairPotentialBase {
        private:
          string _brief() { return name; };
        protected:
          Tshortranged sr;
          Tcoulomb el;
        public:
          Tgeometry geo;
          CoulombSR(InputMap &in) : sr(in), el(in), geo(in) {
            setScale( el.tokT() );
            sr.setScale( el.tokT() );
            name=sr.name+"+"+el.name;
          }
          inline double energy(const particle &a, const particle &b) const {
            double r2=geo.sqdist(a,b);
            return el.energy( a.charge*b.charge, sqrt(r2) ) + sr(a,b,r2);
          }
          inline double operator() (const particle &a, const particle &b, double r2) const {
            return el.energy( a.charge*b.charge, sqrt(r2) ) + sr(a,b,r2);
          }
          string info(char w=20) {
            std::ostringstream o;
            o << el.info(w) << sr.info(w);
            return o.str();
          }
      };

    class MultipoleEnergy {
      public:
        double lB;
        double ionion(double, double, double);
        double iondip(double, const Point&, double);
        double dipdip(const Point&, const Point&, double);
    };

  } //end of potential namespace

} //end of Faunus namespace
#endif
