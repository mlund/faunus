#ifndef FAUNUS_POTENTIAL_H
#define FAUNUS_POTENTIAL_H

#ifndef SWIG
#include <faunus/common.h>
#include <faunus/point.h>
#include <faunus/textio.h>
#include <faunus/physconst.h>
#include <faunus/auxiliary.h>
#include <faunus/inputfile.h>
#include <faunus/species.h>
#endif

using namespace Eigen;

namespace Faunus {

  /**
   * @brief Namespace for pair potentials
   *
   * This namespace contains classes and templates that calculates the
   * pair potential between particles as well as particles with external
   * potentials. The majority of these classes/templates are derived
   * from
   *
   * 1. `Potential::PairPotentialBase`
   * 2. `Potential::ExternalPotentialBase`
   *
   * and thus have common interfaces.
   * Several pair potentials can be combined into
   * one by the template Potential::CombinedPairPotential and a number of
   * common combinations are already defined as typedefs.
   */
  namespace Potential {

    class DebyeHuckel;

    /**
     * @brief Base class for pair potential classes
     *
     * This is a base class for all pair potentials. All derived classes
     * are class functions, i.e. they the function operator returns the
     * energy. For example:
     *
     *     PointParticle a,b;
     *     Coulomb pot(in);
     *     pot(a,b,16.); // -> energy in kT
     * 
     */
    class PairPotentialBase {
      private:
        virtual string _brief();
      protected:
        Eigen::MatrixXf rcut2;                //!< Squared cut-off distance (angstrom)
        void initCutoff(size_t, float);       //!< Initialize all cut-off distances to single value
        void setCutoff(size_t, size_t, float);//!< Specialized cut-off for a pair
      public:  
        PairPotentialBase(std::string="");
        virtual ~PairPotentialBase();
        string name;       //!< Name of potential
        string prefix;     //!< Prefix used for user input specific to the potential
        string brief();    //!< Brief, one-lined information string
        virtual void setTemperature(double); //!< Set temperature [K]

        /**
         * @brief Particle-particle force in units of \c kT/AA
         * @param a First particle
         * @param b Second particle
         * @param r2 Squared distance between them (angstrom squared)
         * @param p Vector from: p=b-a
         */
        template<typename Tparticle>
          Point force(const Tparticle &a, const Tparticle &b, double r2, const Point &p) {
            return Point(0,0,0);
          }

        /** @brief Electric field at spatial position */
        template<typename Tparticle>
          Point field(const Tparticle &a, const Point &r) const {
            return Point(0,0,0);
          }

        template<class Tparticle>
          inline bool save(string, typename Tparticle::Tid, typename Tparticle::Tid); //!< Save table of pair potential to disk
        virtual void test(UnitTest&);                    //!< Perform unit test

        virtual std::string info(char=20);
        
    };

    template<class Tparticle>
      bool PairPotentialBase::save(string filename, typename Tparticle::Tid ida, typename Tparticle::Tid idb) {
        std::ofstream f(filename.c_str());
        if (f) {
          double min=0.9 * (atom[ida].radius+atom[idb].radius);
          particle a,b;
          a = atom[ida];
          b = atom[idb];
          f << "# Pair potential: " << brief() << endl
            << "# Atoms: " << atom[ida].name << "<->" << atom[idb].name << endl;
          for (double r=min; r<=150; r+=0.5)
            f << std::left << std::setw(10) << r << " " ;//<< operator()(a,b,r*r) << endl; 
          return true;
        }
        return false;
      }

    /**
     * @brief Harmonic pair potential
     * @details The harmonic potential has the form
     * \f$ \beta u_{ij} = k(r_{ij}-r_{eq})^2 \f$ where k is the force constant
     * (kT/angstrom^2) and req is the equilibrium distance (angstrom).
     * @note We do not multiply with 1/2 which must be included in the supplied force constant, k
     */
    class Harmonic : public PairPotentialBase {
      private:
        string _brief();
      public:
        double k;   //!< Force constant (kT/A^2) - Did you rember to divide by two? See note.
        double req; //!< Equilibrium distance (angstrom)
        Harmonic(double=0, double=0);
        Harmonic(InputMap&, string="harmonic");
        template<class Tparticle>
          double operator()(const Tparticle &a, const Tparticle &b, double r2) {
            double d=sqrt(r2)-req;
            return k*d*d;
          }
    };

    /**
     * @brief Cosine attraction
     * @details This is an attractive potential used for coarse grained lipids
     * and has the form:
     * @f[
     *     \beta u(r) = -\epsilon \cos^2 [ \pi(r-r_c)/2w_c ]
     * @f]
     * for \f$r_c\leq r \leq r_c+w_c\f$. For \f$r<r_c\f$, \f$\beta u=-\epsilon\f$,
     * while zero for \f$r>r_c+w_c\f$.
     *
     * The InputMap parameters are:
     *
     * Key                | Description
     * :----------------- | :---------------------------
     * `cosattract_eps`   | Depth, \f$\epsilon\f$ [kT]
     * `cosattract_rc`    | Width, r_c [angstrom]
     * `cosattract_wc`    | Decay range, w_c [angstrom] 
     *
     * @warning Untested!
     */
    class CosAttract : public PairPotentialBase {
      private:
        double eps, wc, rc, rc2, c, rcwc2;
        string _brief();
      public:
        CosAttract(InputMap&, string="cosattract"); // Constructor from InputMap
        template<class Tparticle>
          double operator() (const Tparticle &a, const Tparticle &b, double r2) const {
            if (r2<rc2)
              return -eps;
            if (r2>rcwc2)
              return 0;
            double x=cos( c*( sqrt(r2)-rc ) );
            return -eps*x*x;
          }
        string info(char); // More verbose information
    };

    /**
     * @brief Finite Extensible nonlinear elastic (FENE) potential
     * @details This is an anharmonic bonding potential with the form:
     * \f[
     *     \beta u(r) = -\frac{k r_0^2}{2}\ln \left [ 1-(r/r_0)^2 \right ]
     * \f]
     * for \f$r<r_0\f$, otherwise infinity. The input parameters read by InputMap
     * are as follows:
     * - `fene_stiffness` Bond stiffness, `k` [kT]
     * - `fene_maxsep` Maximum separation, `r_0` [angstrom]
     *
     * More info: doi:10.1103/PhysRevE.59.4248
     * 
     */
    class FENE : public PairPotentialBase {
      private:
        double k,r02,r02inv;
        string _brief();
      public:
        FENE(double,double); // Constructor
        FENE(InputMap&, string="fene"); // Constructor
        template<class Tparticle>
          inline double operator() (const Tparticle &a, const Tparticle &b, double r2) {
            return (r2>r02) ? pc::infty : -0.5*k*r02*std::log(1-r2*r02inv);
          }
    };

    /**
     * @brief Hertz pair potential
     */
    class Hertz : public PairPotentialBase {
      private:
        double E;
        string _brief();
      public:
        Hertz(InputMap&, string="hertz");
        template<class Tparticle>
          double operator()(const Tparticle &a, const Tparticle &b, double r2) {
            double m = a.radius+b.radius;
            double diameter = 2.*a.radius;
            if(r2 <= m*m) {
              return E*pow((1-(sqrt(r2)/diameter)),(5./2.));
            }
            return 0;
          }
        template<class Tparticle>
          double operator()(const Tparticle &a, const Tparticle &b,const Point &r2) {
            return operator()(a,b, r2.squaredNorm());

          }
        string info(char); 
    };


    /**
     * @brief Hard sphere pair potential
     */
    class HardSphere : public PairPotentialBase {
      public:
        HardSphere();

        /** @brief Compatibility constructor - no data is read from the InputMap. */
        inline HardSphere(InputMap& in) { name="Hardsphere"; }

        template<class Tparticle>
          double operator() (const Tparticle &a, const Tparticle &b, double r2) {
            double m=a.radius+b.radius;
            return (r2<m*m) ? pc::infty : 0;
          }

        template<class Tparticle>
          double operator() (const Tparticle &a, const Tparticle &b, const Point &r) {
            return operator()(a,b,r.squaredNorm());
          }

        string info(char w);
    };

    /**
     * @brief Lennard-Jones (12-6) pair potential
     * @details The Lennard-Jones potential has the form:
     * @f$
     * \beta u=4\epsilon_{lj} \left ((\sigma_{ij}/r_{ij})^{12}-(\sigma_{ij}/r_{ij})^6\right )
     * @f$
     * where \f$\sigma_{ij} = (\sigma_i+\sigma_j)/2\f$ and \f$\epsilon_{lj}=\epsilon\f$
     * is the same for all pairs in this class.
     */
    class LennardJones : public PairPotentialBase {
      private:
        string _brief();
      protected:
        inline double r6(double sigma, double r2) const {
          double x(sigma*sigma/r2);  // 2
          return x*x*x;             // 6
        }
        double eps;
      public:
        LennardJones();
        LennardJones(InputMap&, string="lj");

        template<class Tparticle>
          double operator() (const Tparticle &a, const Tparticle &b, double r2) const {
            double x(r6(a.radius+b.radius,r2));
            return eps*(x*x - x);
          }
        template<class Tparticle>
          double operator() (const Tparticle &a, const Tparticle &b, const Point &r) {
            return operator()(a,b,r.squaredNorm());
          }

        string info(char);
    };

    /** @brief Lorentz-Berthelot Mixing Rule for sigma and epsilon */
    struct LorentzBerthelot {
      string name;
      LorentzBerthelot();
      double mixSigma(double,double) const;
      double mixEpsilon(double,double) const;
    };

    /**
     * @brief Lennard-Jones with arbitrary mixing rule
     *
     * @details This is a template for Lennard-Jones pair interactions where the template parameter
     * must be a class for the epsilon and sigma mixed rules. The atomic values for 
     * sigma and epsilon are taken from `AtomMap` via the global instance
     * `atom`. In your InputMap configuration file you would typically set the atom list file using
     * the keyword `atomlist`. Note that sigma for each atom is set to two times the radius found in
     * `AtomMap`.
     *
     * For example:
     * 
     *     InputMap mcp("myconfig");
     *     LennardJonesMixed<LorentzBerthelot> lj(mcp);
     *
     */
    template<class Tmixingrule>
      class LennardJonesMixed : public PairPotentialBase {
        private:
          Tmixingrule mixer; // mixing rule class for sigma and epsilon
          string _brief() { return name + " w. " + mixer.name; }
        protected:
          PairMatrix<double> s2,eps; // matrix of sigma_ij^2 and eps_ij
        public:
          LennardJonesMixed(InputMap &in) {
            name="Lennard-Jones";
            size_t n=atom.list.size(); // number of atom types
            s2.resize(n); // not required...
            eps.resize(n);// ...but possible reduced mem. fragmentation
            for (size_t i=0; i<n; i++)
              for (size_t j=0; j<n; j++) {
                s2.set(i,j,
                    pow( mixer.mixSigma( atom.list[i].sigma, atom.list[j].sigma), 2));
                eps.set(i,j,
                    4*mixer.mixEpsilon( atom.list[i].eps, atom.list[j].eps ));
                eps.set(i,j,
                    pc::kJ2kT( eps(i,j) ) ); // convert to kT
              }
          }

          template<class Tparticle>
            double operator()(const Tparticle &a, const Tparticle &b, double r2) const {
              double x=s2(a.id,b.id)/r2; //s2/r2
              x=x*x*x; // s6/r6
              return eps(a.id,b.id) * (x*x - x);
            }

          /**
           * @brief This will set a custom epsilon for a pair of particles
           * @param i Particle id of first particle
           * @param j Particle id of second particle
           * @param eps_kT epsilon in units of kT
           */
          void customEpsilon(particle::Tid i, particle::Tid j, double eps_kT) {
            eps.set(i,j,4*eps_kT);
          }

          void customSigma(particle::Tid i, particle::Tid j, double sigma) {
            s2.set(i,j,sigma*sigma);
          }

          string info(char w=0) {
            using namespace Faunus::textio;
            std::ostringstream o;
            o << indent(SUB) << name+" pair parameters:\n";
            o.precision(4);
            int n=(int)atom.list.size();
            for (int i=0; i<n; i++)
              for (int j=0; j<n; j++)
                if (i>=j)
                  if (i!=0 && j!=0) // ignure first "UNK" particle type
                    o << indent(SUBSUB) << setw(12) << atom[i].name+"<->"+atom[j].name
                      << indent(SUB) << sigma+" = " << setw(7) << sqrt( s2(i,j) ) << _angstrom
                      << indent(SUB) << epsilon+" = " << setw(7) << eps(i,j)/4 << kT+" = "
                      << pc::kT2kJ(eps(i,j)/4) << " kJ/mol"
                      << endl;
            return o.str();
          }
      };

    /**
     * @brief Weeks-Chandler-Andersen pair potential
     * @details This is a Lennard-Jones type potential, cut and shifted to zero
     * at @f$r_c=2^{1/6}\sigma@f$. More info can be found in at
     * <http://doi.org/ct4kh9> and the functional form is:
     * @f[
     * \beta u = 4 \epsilon \left ( (b/r)^{12} - (b/r)^6 + \frac{1}{4} \right )
     * @f]
     * where sigma, epsilon per default are set
     * using Lorentz-Berthelot mixing rules.
     */
    class WeeksChandlerAndersen : public LennardJonesMixed<LorentzBerthelot> {
      protected:
        double onefourth, twototwosixth;
      public:
        typedef LennardJonesMixed<LorentzBerthelot> Tbase;
        WeeksChandlerAndersen(InputMap&);

        template<class Tparticle>
          inline double operator() (const Tparticle &a, const Tparticle &b, double r2) {
            double x=s2(a.id,b.id); // s^2
            if (r2>x*twototwosixth)
              return 0;
            x=x/r2;  // (s/r)^2
            x=x*x*x;// (s/r)^6
            return eps(a.id,b.id)*(x*x - x + onefourth);
          }

        template<class Tparticle>
          double operator() (const Tparticle &a, const Tparticle &b, const Point &r) const {
            return operator()(a,b,r.squaredNorm());
          }
    };


    /*!
     * \brief Square well pair potential
     */
    class SquareWell : public PairPotentialBase {
      private:
        string _brief();
      public:
        double threshold;                           //!< Threshold between particle *surface* [A]
        double depth;                               //!< Energy depth [kT] (positive number)
        SquareWell(InputMap&, string="squarewell"); //!< Constructor

        template<class Tparticle>
          inline double operator() (const Tparticle &a, const Tparticle &b, double r2) {
            double d=a.radius+b.radius+threshold;
            if ( r2 < d*d )
              return -depth;
            return 0;
          }
        string info(char);
    };

    /*!
     * \brief Square well pair potential shifted
     * \author Anil Kurut
     */
    class SquareWellShifted : public SquareWell {
      private:
        string _brief();
      public:
        double threshold_lower;
        SquareWellShifted(InputMap&, string="squarewell"); //!< Constructor

        template<class Tparticle>
          double operator() (const Tparticle &a, const Tparticle &b, double r2) const {
            double d=a.radius+b.radius+threshold_lower;
            if ( r2 > d*d )
              return SquareWell::operator()(a,b,r2);
            return 0;
          }
        string info(char);
    };

    /**
     * @brief Hydrophobic pair potential based on SASA and surface tension
     * @todo Documentation is incorrect.
     * @details The potential is not zero if the distance between hydrophobic particles
     * is smaller than size of solvent molecule (2*Rs)  
     * Potential has the form:
     *
     * \f$ u = Surface tension * (\Delta SASA_i + \Delta SASA_j) \f$
     *
     * Surface area which is not accesible for solvent
     * \f$ \Delta SASA_i = (SASA_i(r_{ij})-SASA_i(\inf))
     * \f$ is calculated based on surface of a sphere cap
     *
     * \f$ SA_{cap,i}=2\pi(R_i+R_s)h_i \f$
     * where h is dependent on distance between the particles as 
     *
     * \f$ h_i=(R_i+R_s)*(\frac{(R_i+R_s)^2-(R_j+R_s)^2)+r_{ij}^2}{2r_{ij}(R_i+R_s)}\f$
     *
     */
    class SquareWellHydrophobic : public SquareWell {
      public:
        SquareWellHydrophobic(InputMap&, string="squarewell"); //!< Constructor

        template<class Tparticle>
          double operator() (const Tparticle &a, const Tparticle &b, double r2) {
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
        double sigma6;
      public:
        SoftRepulsion(InputMap&);

        template<class Tparticle>
          double operator() (const Tparticle &a, const Tparticle &b, double r2) {
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
     * \details \f$ \beta u = 4\epsilon_{lj} \left (  (\sigma_{ij}/r_{ij})^{12}  \right ) \f$
     * where
     * \f$\sigma_{ij} = (\sigma_i+\sigma_j)/2\f$ and \f$\epsilon_{lj}\f$
     * is fixed for this class.
     * \todo Same as LennardJonesR12. Remove?
     */
    class R12Repulsion : public PairPotentialBase {
      private:
        string _brief();
      protected:
        double eps;
      public:
        R12Repulsion(); // __attribute__ ((deprecated));
        R12Repulsion(InputMap&, string="r12rep"); // __attribute__ ((deprecated));
        template<class Tparticle>
          double operator() (const Tparticle &a, const Tparticle &b, double r2) {
            double x=(a.radius+b.radius);
            x=x*x/r2; // r2
            x=x*x*x; // r6
            return eps*x*x;
          }
        string info(char);
    };

    /**
     * @brief Repulsive part of LennardJones
     */
    class LennardJonesR12 : public LennardJones {
      public:
        LennardJonesR12(InputMap&, string="r12rep");
        template<class Tparticle>
          double operator() (const Tparticle &a, const Tparticle &b, double r2) {
            double x=r6(a.radius+b.radius,r2);
            return eps*x*x;
          }
    };

    /**
     * @brief Lennard-Jones truncated and shifted to sigma.
     */
    class LennardJonesTrunkShift : public LennardJones {
      public:
        LennardJonesTrunkShift(InputMap&, string="ljts");
        template<class Tparticle>
          double operator() (const Tparticle &a, const Tparticle &b, double r2) {
            double sigma = a.radius+b.radius;
            if (r2 > sigma*sigma)
              return 0;
            double x=r6(sigma,r2)*0.5;
            return eps*(x*x - x + 0.25);
          }
        template<class Tparticle>
          Point force(const Tparticle &a, const Tparticle &b, double r2, const Point &p) {
            double sigma = a.radius+b.radius;
            if (r2 > sigma*sigma)
              return Point(0,0,0);
            double x=r6(sigma,r2)*0.5;
            return eps*6*(2*x*x - x) / r2 * p;
          }
    };


    /**
     * @brief Coulomb pair potential between charges in a dielectric medium.
     * @details The Coulomb potential has the form
     * @f[
     * \beta u_{ij} = \frac{e^2}{4\pi\epsilon_0\epsilon_rk_BT} \frac{z_i z_j}{r_{ij}}
     *              = \lambda_B \frac{z_i z_j}{r_{ij}}
     * @f]
     * where \f$\lambda_B\f$ is the Bjerrum length and \c z are the valencies.
     */
    class Coulomb : public PairPotentialBase {
      friend class Potential::DebyeHuckel;
      //friend class Energy::GouyChapman;
      private:
      string _brief();

      double epsilon_r;
      protected:
      double depsdt;      //!< \f$ T\partial \epsilon_r / \epsilon_r \partial T = -1.37 \f$
      double lB;          //!< Bjerrum length (angstrom)

      public:
      Coulomb(InputMap&); //!< Construction from InputMap
      double bjerrumLength() const;  //!< Returns Bjerrum length [AA]

      template<class Tparticle>
        double operator() (const Tparticle &a, const Tparticle &b, double r2) const {
#ifdef FAU_APPROXMATH
          return lB*a.charge*b.charge * invsqrtQuake(r2);
#else
          return lB*a.charge*b.charge / sqrt(r2);
#endif
        }

      template<class Tparticle>
        double operator() (const Tparticle &a, const Tparticle &b, const Point &r) {
          return operator()(a,b,r.squaredNorm());
        }

      template<class T>
        Point force(const T &a, const T &b, double r2, const Point &p) {
#ifdef FAU_APPROXMATH
          return lB*a.charge*b.charge * invsqrtQuake(r2) / r2 * p;
#else
          return lB*a.charge*b.charge * p / (sqrt(r2)*r2);
#endif
        }

      /** @brief Electric field at `r` due to charge `p` 
       * Gets returned in [e/Å] (\f$\beta eE \f$)
       */
      template<class Tparticle>
        Point field (const Tparticle &p, const Point &r) const {
          double R2 = 1.0/r.squaredNorm();
          return p.charge*R2*r*sqrt(R2)*lB;
        }

      string info(char);
      void test(UnitTest&); //!< Perform unit test
    };

    /**
     * @brief Coulomb pair potential shifted according to
     *        Wolf/Yonezawaa (doi:10/j97)
     * @details The Coulomb potential has the form:
     * @f[
     * \beta u_{ij} = \frac{e^2}{4\pi\epsilon_0\epsilon_rk_BT}
     * z_i z_j \left (
     * \frac{1}{r} - \frac{1}{R_c} + \frac{r-R_c}{R_c^2}
     * \right )
     * @f]
     */
    class CoulombWolf : public Coulomb {
      private:
        double Rc2, Rcinv;
      public:
        CoulombWolf(InputMap&); //!< Construction from InputMap

        template<class Tparticle>
          double operator() (const Tparticle &a, const Tparticle &b, double r2) {
            if (r2>Rc2)
              return 0;
#ifdef FAU_APPROXMATH
            r2=invsqrtQuake(r2);  // 1/r
            return lB * a.charge * b.charge * (r2 - Rcinv + (Rcinv/r2-1)*Rcinv );
#else
            r2=sqrt(r2); // r
            return lB * a.charge * b.charge * (1/r2 - Rcinv + (r2*Rcinv-1)*Rcinv );
#endif
          }
        string info(char);
    };

    /**
     * @brief Charge-nonpolar pair interaction
     * @details This accounts for polarization of
     * \f[
     * \beta u_{ij} = -\frac{\lambda_B z_i^2 \delta a_j^3}{2r_{ij}^4}
     * \f]
     * where a is the radius of the nonpolar particle.
     * Note that this version requires that one of the particles
     * is charged, while the other is neutral.
     * Delta is a unitless scaling parameter of the excess
     * polarizability.
     * For non-polar particles in a polar medium, this is a negative number.
     * For more information, see Israelachvili, Chapter 5.
     *
     * The InputMap is scanned for
     *
     * - The parameters from `Potential::Coulomb`
     * - `excess_polarization` for the delta value
     *
     */
    class ChargeNonpolar : public Coulomb {
      private:
        double c;
      public:
        ChargeNonpolar(InputMap&); //!< Construction from InputMap
        template<class Tparticle>
          double operator() (const Tparticle &a, const Tparticle &b, double r2) {
            double qq=a.charge * a.charge;
            if (qq>1e-6)
              return -c * qq / (r2*r2) * (b.radius*b.radius*b.radius);
            qq=b.charge * b.charge;
            if (qq>1e-6)
              return -c * qq / (r2*r2) * (a.radius*a.radius*a.radius);
            return 0;
          }
        string info(char);
    };


    class YukawaGel : public Coulomb {
      private:
        string _brief();
        double Z, nc,ns, v, k, Z2e2overER,d,kd,k2d2,ekd, braket7;

      public:

        YukawaGel(InputMap&);
        template<class Tparticle>
          double operator()(const Tparticle &a, const Tparticle &b, double r2) {
            double m = a.radius+b.radius;
            double r = sqrt(r2);
            double kr = k*r;


            double ekr=exp(-kr);

            if(r2 <= m*m) {
              double roverd = r/d;

              double ekdsinhkr = ekd*sinh(kr);

              double A = (2./d)*Z2e2overER;

              double braket = (6./5.)-(2.*pow(roverd,2))+((3./2.)*pow(roverd,3))-((1./5.)*pow(roverd,5));

              double B = (72./((k2d2*k2d2)*r))*Z2e2overER;

              double braket2 = (((1.-ekr+(0.5*kr*kr)+((1./24.)*pow(kr,4)))*(1.-(4./k2d2)))+((4.*ekdsinhkr)/kd));

              double braket3 = ( ekdsinhkr + (k*k*d*r) + ( ((k*k*k*k)/6.) * ((d*d*d*r)+(r*r*r*d)) ) )* (1.+(4./k2d2));

              double braket4 = ((4.*r)/d)*(1.+(0.5*k2d2)+((1./30.)*(k2d2*k2d2)));
              double braket5 = ((8.*r*r*r)/(3.*d*d*d))*((k2d2/4.)+((k2d2*k2d2)/12.));
              double braket6 = (((1./180.)*((k*k*k*k)/(d*d)))*(r*r*r*r*r*r));

              double pot = (A*braket)-(B*(braket2+braket3-braket4-braket5-braket6));
              return pot;
            }
            else{

              double pot2 = ((144./(k2d2*k2d2))*(Z2e2overER)*(braket7*braket7)*(ekr/r));

              return pot2;
            }
          }
        template<class Tparticle>
          double operator() (const Tparticle &a, const Tparticle &b, const Point &r) {
            return operator()(a,b,r.squaredNorm());
          }
        string info(char); 
    };

    /**
     * @brief Debye-Huckel/Yukawa potential
     *1     * @details Similar to the plain Coulomb potential
     *          but with an exponential term to described salt screening:
     * \f[ \beta w_{ij} = \frac{e^2}{4\pi\epsilon_0\epsilon_rk_BT}
     * \frac{z_i z_j}{r_{ij}} \exp(-\kappa r_{ij}) \f]
     * where \f$\kappa=1/D\f$ is the inverse Debye screening length.
     */
    class DebyeHuckel : public Coulomb {
      private:
        string _brief();
      protected:
        double c,k;
      public:
        DebyeHuckel(InputMap&);                       //!< Construction from InputMap
        template<class Tparticle>
          double operator()(const Tparticle &a, const Tparticle &b, double r2) const {
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

    /**
     * @brief DebyeHuckel shifted to reaach zero at given cut-off
     * @details The cut-off distance is read from the InputMap with the following keyword:
     * - `pairpot_cutoff` Spherical cut-off in angstroms
     */
    class DebyeHuckelShift : public DebyeHuckel {
      private:
        double shift;    // offset at cutoff distance
        double sqcutoff; // squared cutoff distance
      public:
        DebyeHuckelShift(InputMap&); //!< Construction from InputMap
        template<class Tparticle>
          inline double operator() (const Tparticle &a, const Tparticle &b, double r2) const {
            if (r2>sqcutoff)
              return 0;
#ifdef FAU_APPROXMATH
            double rinv = invsqrtQuake(r2);
            return lB * a.charge * b.charge * ( exp_cawley(-k/rinv)*rinv - shift );
#else
            double r=sqrt(r2);
            return lB * a.charge * b.charge * ( exp(-k*r)/r - shift );
#endif
          }
    };

    /**
     * @brief Custom potentials between specific particle types
     *
     * If the pair is not recognized, i.e. not added with the
     * `add()` function, the `Tdefault` pair potential is used.
     *
     * Example:
     *
     *     PotentialMap<CoulombLJ> pot(...);
     *     pot.add( atom["Na"].id ,atom["CH4"].id, ChargeNonpolar(...) );
     *     pot.add( atom["Cl"].id ,atom["CH4"].id, ChargeNonpolar(...) );
     *
     */
    template<typename Tdefault, typename Tparticle=PointParticle, typename Tdist=double>
      class PotentialMap : public Tdefault {
        protected:
          typedef opair<int> Tpair;
          typedef std::function<double(const Tparticle&,const Tparticle&,Tdist)> Tfunc;
          std::map<Tpair,Tfunc> m;
          std::string _info; // info for the added potentials (before turning into functors)

        public:
          PotentialMap(InputMap &in) : Tdefault(in) {
            Tdefault::name += " (default)";
          } 

          template<class Tpairpot>
            void add(AtomData::Tid id1, AtomData::Tid id2, Tpairpot pot) {
              pot.name=atom[id1].name + "<->" + atom[id2].name + ": " + pot.name;
              _info+="\n  " + pot.name + ":\n" + pot.info(20);
              m[Tpair(id1,id2)] = pot;
            }

          double operator()(const Tparticle &a, const Tparticle &b, const Tdist &r2) {
            auto i=m.find( Tpair(a.id,b.id) );
            if (i!=m.end())
              return i->second(a,b,r2);
            return Tdefault::operator()(a,b,r2);
          }

          std::string info(char w=20) {
            return Tdefault::info(w) + _info;
          }
      };



    /**
     * @brief Combines two pair potentials
     * @details This combines two PairPotentialBases. The combined potential
     * can subsequently be used as a normal pair potential and even be
     * combined with a third potential and so forth.
     *
     *     // mix two and three pair potentials
     *     using namespace Potential;
     *     typedef CombinedPairPotential< LennardJones, SquareWell > Tpairpot1;
     *     typedef CombinedPairPotential< Tpairpot1, Coulomb > Tpairpot2;
     *     Tpairpot2 mypairpot;
     *     std::cout << mypairpot.info();
     *
     * @date Lund, 2012
     */
    template<class T1, class T2>
      class CombinedPairPotential : public PairPotentialBase {
        private:
          string _brief() {
            return first.brief() + " " + second.brief();
          }
        public:
          T1 first;  //!< First pair potential of type T1
          T2 second; //!< Second pair potential of type T2

          CombinedPairPotential(T1 a, T2 b) : first(a), second(b) {
            name=first.name+"+"+second.name;
          }

          CombinedPairPotential(InputMap &in) : first(in), second(in) {
            name=first.name+"+"+second.name;
          }

          CombinedPairPotential(InputMap &in, string pfx1, string pfx2) :
            PairPotentialBase(pfx1+pfx2), first(in,pfx1), second(in,pfx2) {
              name=first.name+"+"+second.name;
            }

          template<class Tparticle, class Tdist>
            double operator()(const Tparticle &a, const Tparticle &b, const Tdist &r2) {
              return first(a,b,r2) + second(a,b,r2);
            }

          template<typename Tparticle>
            Point field(const Tparticle &a, const Point &r) const {
              return first.field(a,r) + second.field(a,r);
            }

          string info(char w=20) {
            return first.info(w) + second.info(w);
          }

          void test(UnitTest &t) {
            first.test(t);
            second.test(t);
          }
      };

    /**
     * @brief Adds pair potentials together
     *
     * Example:
     *
     *     auto primitiveModel = Potential::Coulomb() + Potential::HardSphere();
     */
    template<class T1, class T2,
      class = typename std::enable_if<std::is_base_of<PairPotentialBase,T1>::value>::type,
      class = typename std::enable_if<std::is_base_of<PairPotentialBase,T2>::value>::type>
        Potential::CombinedPairPotential<T1,T2>& operator+(const T1 &pot1, const T2 &pot2) {
          return *(new Potential::CombinedPairPotential<T1,T2>(pot1,pot2));
        }

    class MultipoleEnergy {
      public:
        double lB;
        double ionion(double, double, double);
        double iondip(double, const Point&, double);
        double dipdip(const Point&, const Point&, double);
    };

    /**
     * @brief Lennard-Jones potential with Lorentz-Berthelot mixing rule
     */
    typedef LennardJonesMixed<LorentzBerthelot> LennardJonesLB;

    /**
     * @brief Combined Coulomb / HardSphere potential
     */
    typedef CombinedPairPotential<Coulomb, HardSphere> CoulombHS;

    /**
     * @brief Combined Coulomb / LennardJones potential
     */
    typedef CombinedPairPotential<Coulomb, LennardJones> CoulombLJ;

    /**
     * @brief Combined Coulomb / WeeksChandlerAndersen potential
     */
    typedef CombinedPairPotential<Coulomb, WeeksChandlerAndersen> CoulombWCA;

    /**
     * @brief Combined Coulomb / LennardJonesTrunkShift potential
     */
    typedef CombinedPairPotential<Coulomb, LennardJonesTrunkShift> CoulombLJTS;

    /**
     * @brief Combined Coulomb / LennardJones potential
     */
    typedef CombinedPairPotential<CoulombWolf, LennardJones> CoulombWolfLJ;

    /**
     * @brief Combined DebyeHuckel / HardSphere potential
     */
    typedef CombinedPairPotential<DebyeHuckel, HardSphere> DebyeHuckelHS;

    /**
     * @brief Combined DebyeHuckel / LennardJones potential
     */
    typedef CombinedPairPotential<DebyeHuckel, LennardJones> DebyeHuckelLJ;

    /**
     * @brief Combined DebyeHuckel / R12Repulsion potential
     */
    typedef CombinedPairPotential<DebyeHuckel, R12Repulsion> DebyeHuckelr12;

    typedef CombinedPairPotential<Hertz,YukawaGel> HertzYukawa;

  } //end of Potential namespace

} //end of Faunus namespace
#endif
