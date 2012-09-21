#ifndef FAU_POTENTIAL_H
#define FAU_POTENTIAL_H

#ifndef SWIG
#include <faunus/common.h>
#include <faunus/point.h>
#include <faunus/textio.h>
#include <faunus/physconst.h>
#include <faunus/auxiliary.h>
#include <faunus/inputfile.h>
#include <faunus/species.h>
#endif

namespace Faunus {

  /*!
   * \brief Namespace for various potentials - pair potentials, external potentials etc.
   */
  namespace Potential {

    class DebyeHuckel;

    /*!
     * \brief Base class for pair potential classes
     *
     * This is a base class for all pair potentials which must implement the function operator so
     * that the potential can work as a class function. To make a new pair potential you must
     * implement 1) a function that takes two particles as arguments as well as the squared distance
     * between them (i.e. the function operator), and 2)
     * a brief information string (both are pure virtual). The unit of the
     * returned energy is arbitrary but you *must* ensure that when multiplied by tokT() that
     * it is converted to kT units (thermal energy). By default _tokT=1 and it is a good policy
     * to return energies in kT. Several pair potentials can be combined by the class template
     * Potential::CombinedPairPotential.
     */
    class PairPotentialBase {
      private:
        virtual string _brief()=0;
        virtual void _setScale(double);
      protected:
        double _tokT;
      public:  
        PairPotentialBase();
        virtual ~PairPotentialBase();
        string name;             //!< Short (preferably one-word) description of the core potential
        string brief();          //!< Brief, one-lined information string
        void setScale(double=1); //!< Set scaling factor
        /*! \brief Convert returned energy to kT.*/
        inline double tokT() {
          return _tokT;
        }
        virtual void setTemperature(double); //!< Set temperature [K]

        /*!
         * \brief Particle-particle energy divided by tokT()
         * \param a First particle
         * \param b Second particle
         * \param r2 Squared distance between them (angstrom squared)
         */
        virtual double operator() (const particle &a, const particle &b, double r2) const=0;

        bool save(string, particle::Tid, particle::Tid); //!< Save table of pair potential to disk
    };

    /*!
     * \brief Harmonic pair potential
     * \note We do not multiply with 1/2 which must be included in the supplied force constant, k
     *
     * The harmonic potential has the form \f$ \beta u_{ij} = k(r_{ij}-r_{eq})^2 \f$ where k is the force constant
     * (kT/angstrom^2) and req is the equilibrium distance (angstrom).
     */
    class Harmonic : public PairPotentialBase {
      private:
        string _brief();
        void _setScale(double);
      public:
        double k;   //!< Force constant (kT/A^2) - Did you rember to divide by two? See note.
        double req; //!< Equilibrium distance (angstrom)
        Harmonic(double=0, double=0);
        Harmonic(InputMap&, string="harmonic_");
        inline double operator() (const particle &a, const particle &b, double r2) const FOVERRIDE {
          double d=sqrt(r2)-req;
          return k*d*d;
        }
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
        inline double operator() (const particle &a, const particle &b, double r2) const FOVERRIDE {
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
        inline double r6(double sigma, double r2) const {
          double x=sigma*sigma/r2;  // 2
          return x*x*x;             // 6
        }
        double eps;
      public:
        LennardJones();
        LennardJones(InputMap&, string="lj_");
        inline double operator() (const particle &a, const particle &b, double r2) const FOVERRIDE {
          double x=r6(a.radius+b.radius,r2);
          return eps*(x*x - x);
        }
        string info(char);
    };

    /*!
     * \brief As LennardJones but only the repulsive R12 part.
     */
    class LennardJonesR12 : public LennardJones {
      public:
        LennardJonesR12(InputMap&, string="r12rep_");
        inline double operator() (const particle &a, const particle &b, double r2) const FOVERRIDE {
          double x=r6(a.radius+b.radius,r2);
          return eps*x*x;
        }
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
        double depth;                               //!< Energy depth [kT] (positive number)
        SquareWell(InputMap&, string="squarewell"); //!< Constructor
        inline double operator() (const particle &a, const particle &b, double r2) const FOVERRIDE {
          double d=a.radius+b.radius+threshold;
          if ( r2 < d*d )
            return -depth;
          return 0;
        }
        string info(char);
    };
    /*!
     * \brief Hydrophobic pair potential based on SASA and surface tension
     *
     * The potential is not zero iff the distance between hydrophobic particles is smaller than size of solvent molecule (2*Rs)  
     *
     * Potential has the form:
     *
     * \f$ u = Surface tension * (\Delta SASA_i + \Delta SASA_j) \f$
     *
     * Surface area which is not accesible for solvent \f$ \Delta SASA_i = (SASA_i(r_{ij})-SASA_i(\inf)) \f$ is calculated based on surface of a sphere cap
     *
     * \f$ SA_{cap,i}=2\pi(R_i+R_s)h_i \f$ where h is dependent on distance between the particles as 
     *
     * \f$ h_i=(R_i+R_s)*(\frac{(R_i+R_s)^2-(R_j+R_s)^2)+r_{ij}^2}{2r_{ij}(R_i+R_s)}\f$
     *
     */
    class SquareWellHydrophobic : public SquareWell {
      public:
        SquareWellHydrophobic(InputMap&, string="squarewell"); //!< Constructor
        inline double operator() (const particle &a, const particle &b, double r2) const FOVERRIDE {
          if (a.hydrophobic)
            if (b.hydrophobic)
              return SquareWell::operator()(a,b,r2);
          return 0;
        }
    };

    /*!
     * \brief Soft repulsion of the form \f$ \beta u = \sigma^6 / (r_{ij}-r_i-r_j)^6 \f$
     * \todo This applies sqrt() and thus may be slow. Also remove floating point comparison.
     */
    class SoftRepulsion : public PairPotentialBase {
      private:
        string _brief();
        void _setScale(double);
        double sigma6;
      public:
        SoftRepulsion(InputMap&);
        inline double operator() (const particle &a, const particle &b, double r2) const FOVERRIDE {
          double d=a.radius+b.radius;
          if (r2<=d*d) //fp comparison!
            return pc::infty;
          d = sqrt(r2)-d;
          d=d*d; //d^2
          return sigma6 / (d*d*d); // s^6/d^6
        }
        string info(char);
    };

    /*!
     * \brief r12-Repulsion of the form
     *
     * \f$ \beta u = 4\epsilon_{lj} \left (  (\sigma_{ij}/r_{ij})^{12}  \right ) \f$
     *
     * where \f$\sigma_{ij} = (\sigma_i+\sigma_j)/2\f$ and \f$\epsilon_{lj}\f$ is fixed for this class.
     * \todo Inherit from Potential::LennardJones
     */
    class R12Repulsion : public PairPotentialBase {
      private:
        string _brief();
        void _setScale(double);
        inline double r6(double sigma, double r2) const {
          double x=sigma*sigma/r2;  // 2
          return x*x*x;             // 6
        }
        inline double energy(double sigma, double r2) const {
          double x=r6(sigma,r2);
          return eps*(x*x);
        }
      protected:
        double eps;
      public:
        R12Repulsion();
        R12Repulsion(InputMap&, string="r12rep_");
        inline double operator() (const particle &a, const particle &b, double r2) const FOVERRIDE {
          return energy(a.radius+b.radius, r2);
        }
        string info(char);
    };

    /*!
     * \brief Coulomb pair potential
     * 
     * The Coulomb potential has the form:
     * \f[ \beta u_{ij} = \frac{e^2}{4\pi\epsilon_0\epsilon_rk_BT} \frac{z_i z_j}{r_{ij}} \f]
     * where the first factor is the Bjerrum length, \f$l_B\f$. By default the scaling is set
     * to the Bjerrum length - i.e. returned energies are divided by \f$l_B\f$ and to get the energy
     * in kT simply multiply with tokT(). This to increase performance when summing over many
     * particles.
     */
    class Coulomb : public PairPotentialBase {
      friend class Potential::DebyeHuckel;
      friend class Energy::GouyChapman;
      private:
      string _brief();
      void _setScale(double);
      double epsilon_r;
      protected:
      double depsdt;      //!< \f$ T\partial \epsilon_r / \epsilon_r \partial T = -1.37 \f$
      double lB;          //!< Bjerrum length (angstrom)

      public:
      Coulomb(InputMap&); //!< Construction from InputMap
      double bjerrumLength() const;  //!< Returns Bjerrum length [AA]

      /*! \returns \f$\beta u/l_B\f$ */
      inline double operator() (const particle &a, const particle &b, double r2) const FOVERRIDE {
#ifdef FAU_APPROXMATH
        return lB*a.charge*b.charge * invsqrtQuake(r2);
#else
        return lB*a.charge*b.charge / sqrt(r2);
#endif
      }
      string info(char);
    };

    /*!
     * \brief Debye-Huckel pair potential
     *
     * This potential is similar to the plain Coulomb potential but with an extra exponential term to described salt screening:
     * \f[ \beta w_{ij} = \frac{e^2}{4\pi\epsilon_0\epsilon_rk_BT} \frac{z_i z_j}{r_{ij}} \exp(-\kappa r_{ij}) \f]
     * where \f$\kappa=1/D\f$ is the inverse Debye screening length.
     */
    class DebyeHuckel : public Coulomb {
      private:
        string _brief();
      protected:
        double c,k;
      public:
        DebyeHuckel(InputMap&);                       //!< Construction from InputMap
        /*! \returns \f$\beta w/l_B\f$ */
        inline double operator() (const particle &a, const particle &b, double r2) const FOVERRIDE {
#ifdef FAU_APPROXMATH
          double rinv = invsqrtQuake(r2);
          return lB * a.charge * b.charge * rinv * exp_cawley(-k/rinv);
#else
          double r=sqrt(r2);
          return lB * a.charge * b.charge / r * exp(-k*r);
#endif
        }
        double entropy(double, double) const;         //!< Returns the interaction entropy 
        double ionicStrength() const;                 //!< Returns the ionic strength (mol/l)
        double debyeLength() const;                   //!< Returns the Debye screening length (angstrom)
        double excessChemPot(double, double=0) const; //!< Single ion excess chemical potential (kT)
        double activityCoeff(double, double=0) const; //!< Single ion activity coefficient (molar scale) 
        string info(char);
    };

    /*!
     * \brief as DebyeHuckel but shifted to reaach zero at a user specified cut-off distance
     *
     * The cut-off distance is read from the InputMap with the following keyword:
     * \li \c pairpot_cutoff Spherical cut-off in angstroms
     */
    class DebyeHuckelShift : public DebyeHuckel {
      private:
        double shift;
      public:
        DebyeHuckelShift(InputMap&);                       //!< Construction from InputMap
        inline double operator() (const particle &a, const particle &b, double r2) const {
#ifdef FAU_APPROXMATH
          double rinv = invsqrtQuake(r2);
          return lB * a.charge * b.charge * ( exp_cawley(-k/rinv)*rinv - shift );
#else
          double r=sqrt(r2);
          return lB * a.charge * b.charge * ( exp(-k*r)/r - shift );
#endif
        }
    };

    /*!
     * \brief Combines two pair potentials
     * \date Lund, 2012
     * \author Mikael Lund
     * \note tokT() functions are *ignored* and the two pair potentials must return
     *       energies in kT directly. The plan is to remove tokT() permanently.
     *
     * This simply combines two PairPotentialBases. The combined potential can subsequently
     * be used as a normal pair potential and even be combined with a third potential and
     * so forth. A number of typedefs such as Potential::CoulombHS is aliasing this.
     *
     * \code
     *   // mix two and three pair potentials
     *   using namespace Potential;
     *   typedef CombinedPairPotential< LennardJones, SquareWell > Tpairpot1;
     *   typedef CombinedPairPotential< Tpairpot1, Coulomb > Tpairpot2;
     *   Tpairpot2 mypairpot;
     *   std::cout << mypairpot.info();
     * \endcode
     */
    template<class T1, class T2>
      class CombinedPairPotential : public PairPotentialBase {
        protected:
          T1 sr1;
          T2 sr2;
        private:
          string _brief() { return sr1.brief() + " " + sr2.brief(); }

          void _setScale(double s) {
            //sr2.setScale( s*sr1.tokT() );
            //sr1.setScale( s );
          }
        public:
          CombinedPairPotential(InputMap &in) : sr1(in), sr2(in) {
            name=sr1.name+"+"+sr2.name;
          }
          CombinedPairPotential(InputMap &in, string pfx1, string pfx2) : sr1(in,pfx1), sr2(in,pfx2) {
            name=sr1.name+"+"+sr2.name;
          }
          inline double operator() (const particle &a, const particle &b, double r2) const FOVERRIDE {
            return sr1(a,b,r2) + sr2(a,b,r2);
          }
          string info(char w=20) {
            std::ostringstream o;
            o << sr1.info(w) << sr2.info(w);
            return o.str();
          }
      };

    /*!
     * \brief As CombinedPairPotential but with a spherical cut-off
     *
     * The following keywords are read from the InputMap:
     * \li \c pairpot_cutoff Spherical cut-off distance (Default: infinity)
     */
    template<class T1, class T2>
      class CombinedPairPotentialCutoff : public CombinedPairPotential<T1,T2> {
        private:
          typedef CombinedPairPotential<T1,T2> Tbase;
          double cutoff2; // squared, spherical cutoff
          void getin(InputMap &in) {
            cutoff2 = pow(in.get<double>("pairpot_cutoff",pc::infty),2);
          }
        public:
          CombinedPairPotentialCutoff(InputMap &in) : Tbase(in) { getin(in); }

          CombinedPairPotentialCutoff(InputMap &in, string pfx1, string pfx2) : Tbase(in,pfx1,pfx2) {
            getin(in);
          }

          inline double operator() (const particle &a, const particle &b, double r2) const {
            return (r2>cutoff2) ? 0 : Tbase::operator()(a,b,r2);
          }

          string info(char w=20) {
            std::ostringstream o;
            o << Tbase::info(w)
              << textio::pad(textio::SUB,w,"Pair potential cutoff") << sqrt(cutoff2) << textio::_angstrom << endl;
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

    /*!
     * \brief Combined Coulomb / HardSphere potential
     */
    typedef CombinedPairPotential<Coulomb, HardSphere> CoulombHS;

    /*!
     * \brief Combined Coulomb / LennardJones potential
     */
    typedef CombinedPairPotential<Coulomb, LennardJones> CoulombLJ;

    /*!
     * \brief Combined DebyeHuckel / HardSphere potential
     */
    typedef CombinedPairPotential<DebyeHuckel, HardSphere> DebyeHuckelHS;

    /*!
     * \brief Combined DebyeHuckel / LennardJones potential
     */
    typedef CombinedPairPotential<DebyeHuckel, LennardJones> DebyeHuckelLJ;

    /*!
     * \brief Combined DebyeHuckel / R12Repulsion potential
     */
    typedef CombinedPairPotential<DebyeHuckel, R12Repulsion> DebyeHuckelr12;

  } //end of Potential namespace

} //end of Faunus namespace
#endif
