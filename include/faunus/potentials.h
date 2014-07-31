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
#include <faunus/average.h>
#endif

namespace Faunus {

  /**
   * @brief Namespace for pair potentials
   *
   * This namespace contains classes and templates that calculates the
   * pair interaction energy and force between particles.
   * All pair potentials are based
   * on `PairPotentialBase` and thus have common interfaces.
   * Several pair potentials can be combined into
   * a new pair potential and a number of
   * common combinations are pre-defined as typedefs.
   *
   * ~~~
   * CoulombHS primitiveModel();                // Note that constructor
   * auto nonbond = Coulomb() + LennardJones(); // arguments are here
   * auto bond    = Harmonic() - nonbond;       // omitted for clarity
   *
   * PointParticle a,b;                         // two particles
   * double r2 = 100;                           // a<->b squared distance
   * double u  = nonbond(a,b,r2);               // a<->b energy in kT
   * ~~~
   *
   * As shown in the last example, pair potentials can also be subtracted
   * which can be used to for example exclude non-bonded interactions
   * between bonded pairs.
   *
   * If the pair interaction depends on particle types, use
   * `PotentialMap` pair interaction.
   *
   * *Behind the scene:*
   *
   * The above combination of pair potentials is done at compile time
   * using templates. This means that there is a good chance that
   * the mixing overhead can be optimized out by the compiler.
   * For example, when adding two potentials
   * we construct a new `CombinedPairPotential`. Likewise when
   * subtracting, a new `Minus` template is created and then combined.
   *
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
        //void initCutoff(size_t, float);       //!< Initialize all cut-off distances to single value
        //void setCutoff(size_t, size_t, float);//!< Specialized cut-off for a pair
      public:
        typedef PairMatrix<double> Tcutoff;
        Tcutoff rcut2;                        //!< Squared cut-off distance (angstrom^2)
        PairPotentialBase(std::string="");
        virtual ~PairPotentialBase();
        string name;       //!< Name of potential
        string prefix;     //!< Prefix used for user input specific to the potential
        string brief();    //!< Brief, one-lined information string
        virtual void setTemperature(double); //!< Set temperature [K]

        /**
         * @brief Particle-particle force in units of `kT/AA`
         * @param a First particle
         * @param b Second particle
         * @param r2 Squared distance between them (angstrom squared)
         * @param p Vector from: p=b-a
         */
        template<typename Tparticle>
          Point force(const Tparticle &a, const Tparticle &b, double r2, const Point &p) {
            std::cerr << "Using unimplemented force!" << std::endl;
            return Point(0,0,0);
          }

        /** @brief Electric field at spatial position */
        template<typename Tparticle>
          Point field(const Tparticle &a, const Point &r) const {
            return Point(0,0,0);
          }

        /**
         * @brief Set space dependent features such as density dependent potentials
         *
         * The base-class version does nothing.
         */
        template<class Tspace>
          void setSpace(Tspace&) {} 

        virtual void test(UnitTest&);                    //!< Perform unit test

        virtual std::string info(char=20);

    };

    /**
     * @brief Save pair potential and force table to disk
     *
     * This will save a pair potential to disk as three columns:
     * `distance`, `energy`, and `force`. The distance interval
     * is hard coded to `dmin=a.radius+b.radius` to `5*dmin`.
     *
     * Example:
     * ~~~~
     * using namespace Potential;
     * CoulombLJ pot(...);
     * save(pot, atom["Na"].id, atom["Cl"].id, "mytable.dat");
     * ~~~~
     */
    template<class Tpairpot, class Tid>
      bool save(Tpairpot pot, Tid ida, Tid idb, string file) {
        std::ofstream f(file.c_str());
        if (f) {
          double min=atom[ida].radius+atom[idb].radius;
          particle a,b;
          a = atom[ida];
          b = atom[idb];
          f << "# Pair potential: " << pot.brief() << endl
            << "# Atoms: " << atom[ida].name << "<->" << atom[idb].name
            << endl;
          f.precision(12);
          for (double r=min; r<=5*min; r+=0.01)
            f << std::left << std::setw(12) << r << " "
              << pot(a,b,r*r) << " " << pot.force(a,b,r*r,Point(r,0,0)).x() << endl; 
          return true;
        }
        return false;
      }

    /**
     * @brief Harmonic pair potential
     * @details The harmonic potential has the form
     * \f$ \beta u_{ij} = k(r_{ij}-r_{eq})^2 \f$ where k is the force constant
     * (kT/angstrom^2) and req is the equilibrium distance (angstrom).
     * @note We do not multiply with 1/2 which must be included in the
     * supplied force constant, k
     */
    class Harmonic : public PairPotentialBase {
      private:
        string _brief();
      public:
        double k;   //!< Force constant (kT/A^2) - Remember to divide by two!
        double req; //!< Equilibrium distance (angstrom)

        Harmonic(double=0, double=0);

        Harmonic(InputMap&, string="harmonic");

        template<class Tparticle>
          double operator()(const Tparticle &a, const Tparticle &b, double r2) {
            double d=sqrt(r2)-req;
            return k*d*d;
          }

        template<class Tparticle>
          Point force(const Tparticle &a, const Tparticle &b, double r2, const Point &p) {
            double r=sqrt(r2), d=r-req;
            return -2 * k * d / r * p;
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
     */
    class CosAttract : public PairPotentialBase {
      protected:
        double eps, wc, rc, rc2, c, rcwc2;
        string _brief();
      public:
        CosAttract(InputMap&, string="cosattract"); // Constructor from InputMap

        /**
         * @todo
         * The function `x(c,r2,rc)` defined herein could be approximated
         * by a series expansion for `r2=rcwc2`. In this way one can
         * avoid `cos()` and `sqrt()`. C code for this could be generated
         * in Matlab:
         *
         * ~~~
         * with(CodeGeneration)
         * x := series(cos(c*(sqrt(r2)-rc)), r2 = rcwc2, 2)
         * convert(x, polynom)
         * C(%, resultname = "x")
         * ~~~
         */
        template<class Tparticle>
          double operator() (const Tparticle &a, const Tparticle &b, double r2) const {
            if (r2<rc2)
              return -eps;
            if (r2>rcwc2)
              return 0;
#ifdef FAU_APPROXMATH
            double x=cosApprox( c*( sqrt(r2)-rc ) );
#else
            double x=cos( c*( sqrt(r2)-rc ) );
#endif
            return -eps*x*x;
          }

        string info(char); // More verbose information

        template<class Tparticle>
          Point force(const Tparticle &a, const Tparticle &b, double r2, const Point &p) {
            if (r2<rc2 || r2>rcwc2)
              return Point(0,0,0);
            double r=sqrt(r2);
#ifdef FAU_APPROXMATH
            double x1=cosApprox( c*( r-rc ) );
            double x2=sinApprox( c*( r-rc ) );
#else
            double x1=cos( c*( r-rc ) );
            double x2=sin( c*( r-rc ) );
#endif
            return -2*c*eps*x1*x2/r*p;
          }
    };

    /**
     * @brief Finite Extensible Nonlinear Elastic (FENE) potential
     * @details This is an anharmonic bonding potential with the form:
     * @f[
     *     \beta u(r) = -\frac{k r_0^2}{2}\ln \left [ 1-(r/r_0)^2 \right ]
     * @f]
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

        template<class Tparticle>
          Point force(const Tparticle &a, const Tparticle &b, double r2, const Point &p) {
            return (r2>r02) ? -pc::infty*p : -k * r02 / (r02-r2) * p;
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
     * \beta u=4\epsilon_{lj}
     * \left ((\sigma_{ij}/r_{ij})^{12}-(\sigma_{ij}/r_{ij})^6\right )
     * @f$
     * where
     * \f$\sigma_{ij} = (\sigma_i+\sigma_j)/2\f$ and \f$\epsilon_{lj}=\epsilon\f$
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
        template<typename Tparticle>
          Point force(const Tparticle &a, const Tparticle &b, double r2, const Point &p) {
            double s6=pow(a.radius+b.radius,6); //pow() is slow
            double r6=r2*r2*r2;
            double r14=r6*r6*r2;
            return 6.*eps*s6*(2*s6-r6)/r14*p;
          }

        string info(char);
    };

    /**
     * @brief Cuts a pair-potential and shift to zero at cutoff
     *
     * This will cut any pair potential at `prefix_cutoff` and shift to
     * zero at that distance. Slow but general.
     *
     * Example:
     *
     * ~~~~
     * using namespace Faunus::Potential;
     * typedef CutShift<LennardJones> Tpairpot;
     * ~~~~
     *
     * Upon construction with an `InputMap`, the `lj_cutoff`
     * keyword will be used to set the cut-off.
     *
     *
     * @todo Implement force calculation
     */
    template<class Tpairpot>
      class CutShift : public Tpairpot {
        private:
          double rc2; // squared cut-off distance
        public:
          CutShift(InputMap &in) : Tpairpot(in) {
            rc2=pow(in.get<double>(Tpairpot::prefix+"_cutoff", pc::infty),2);
            Tpairpot::name+=" (rcut="
              + std::to_string(sqrt(rc2)) + textio::_angstrom + ")";
          }
          template<class Tparticle>
            double operator() (const Tparticle &a, const Tparticle &b, double r2) const {
              return (r2>rc2) ? 0 : Tpairpot::operator()(a,b,r2) - Tpairpot::operator()(a,b,rc2);
            }
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
     * `AtomMap`. Epsilon is stored internally in `kT*4`.
     *
     * For example:
     * 
     *     InputMap mcp("myconfig");
     *     LennardJonesMixed<LorentzBerthelot> lj(mcp);
     *
     */
    template<class Tmixingrule>
      class LennardJonesMixed : public PairPotentialBase {
        protected:
          Tmixingrule mixer; // mixing rule class for sigma and epsilon
          string _brief() { return name + " w. " + mixer.name; }
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

          template<typename Tparticle>
            Point force(const Tparticle &a, const Tparticle &b, double r2, const Point &p) {
              double s6=pow(s2(a.id,b.id),3);
              double r6=r2*r2*r2;
              double r14=r6*r6*r2;
              return 6.*eps(a.id,b.id) * s6 * (2*s6-r6) / r14 * p;
            }

          template<class Tparticle>
            double operator()(const Tparticle &a, const Tparticle &b, const Point &r) const {
              return operator()(a,b,r.squaredNorm()); 
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
            o << indent(SUB) << name+" pair parameters:\n\n";
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

    template<class Tmixingrule=LorentzBerthelot>
      class CosAttractMixed : public LennardJonesMixed<Tmixingrule> {
        protected:
          typedef LennardJonesMixed<Tmixingrule> base;
          PairMatrix<double> rc2,rc,c; // matrix of sigma_ij^2 and eps_ij
        public:
          CosAttractMixed(InputMap &in) : base(in) {
            base::name="Cos" + textio::squared + " attraction (mixed)";
            size_t n=atom.list.size(); // number of atom types
            c.resize(n);
            rc.resize(n);
            rc2.resize(n);
            for (size_t i=0; i<n; i++)
              for (size_t j=i; j<n; j++) {
                rc.set(i,j,base::mixer.mixSigma( atom.list[i].pdis, atom.list[j].pdis));
                base::rcut2.set(i,j,base::mixer.mixSigma( atom.list[i].pswitch, atom.list[j].pswitch));
                c.set(i,j, 0.5*pc::pi/base::rcut2(i,j) );
                base::rcut2.set(i,j, base::rcut2(i,j) + rc(i,j) );
                base::rcut2.set(i,j, base::rcut2(i,j) * base::rcut2(i,j) );
                rc2.set(i,j, rc(i,j)*rc(i,j) );
              }
          }

          template<class Tparticle>
            double operator()(const Tparticle &a, const Tparticle &b, double r2) const {
              if (r2<rc2(a.id,b.id)) {
                //epsilon is from LJmixed stored as multiplied by 4
                return -0.25*base::eps(a.id,b.id);
              }
              if (r2>base::rcut2(a.id,b.id))
                return 0;
              double x=cos( c(a.id,b.id)*( sqrt(r2)-rc(a.id,b.id) ) );
              //epsilon is from LJmixed stored as multiplied by 4
              return -0.25*base::eps(a.id,b.id)*x*x;
            }
          template<class Tparticle>
            double operator() (const Tparticle &a, const Tparticle &b, const Point &r) const {
              return operator()(a,b,r.squaredNorm());
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

        template<class Tparticle>
          Point force(const Tparticle &a, const Tparticle &b, double r2, const Point &p) {
            double x=s2(a.id,b.id); // s^2
            if (r2>x*twototwosixth)
              return Point(0,0,0);
            x=x/r2;  // (s/r)^2
            x=x*x*x;// (s/r)^6
            return eps(a.id,b.id)*6*(2*x*x - x)/r2*p;
          }
    };

    /**
     * @brief Square well pair potential
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
     *
     * Upon construction, the following parameters are read from the
     * `InputMap`:
     *
     * Keyword        |  Description
     * :------------- | :---------------------
     *  `temperature` | Absolute temperature (K) [default: 298.15 K]
     *  `epsilon_r`   | Relative dielectric constant [default: 80]
     *  `depsdt`      | Dielectric constant temperature derivative [default: -0.368]
     */
    class Coulomb : public PairPotentialBase {
      friend class Potential::DebyeHuckel;
      private:
      string _brief();

      double epsilon_r;
      protected:
      double depsdt;      //!< \f$ T\partial \epsilon_r / \epsilon_r \partial T = -1.37 \f$
      double lB;          //!< Bjerrum length (angstrom)

      public:
      Coulomb(InputMap&, string="coulomb"); //!< Construction from InputMap
      double bjerrumLength() const;  //!< Returns Bjerrum length [AA]

      template<class Tparticle>
        double operator() (const Tparticle &a, const Tparticle &b, double r2) {
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

      template<class Tparticle>
        Point force(const Tparticle &a, const Tparticle &b, double r2, const Point &p) {
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
          double r2i = 1.0/r.squaredNorm();
          return p.charge*r2i*r*sqrt(r2i)*lB;
        }

      string info(char);
      void test(UnitTest&); //!< Perform unit test
    };

    /**
     * @brief Coulomb pair potential shifted according to Wolf/Yonezawa
     * @details The potential has the form:
     * @f[
     * \beta u_{ij} = \frac{e^2}{4\pi\epsilon_0\epsilon_rk_BT}
     * z_i z_j \left (
     * \frac{1}{r} - \frac{1}{R_c} + \frac{r-R_c}{R_c^2}
     * \right )
     * @f]
     *
     * and is hence a particular simple form of the original Wolf
     * formulation. This potential is expected to work reasonably well
     * for dense liquids, see [here](http://dx.doi.org/10/j97).
     *
     * Upon construction using an `InputMap` the keywords from `Potential::Coulomb`
     * are used in addition to `coulomb_cut` to specifify the cut-off.
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

        template<class Tparticle>
          double operator() (const Tparticle &a, const Tparticle &b, const Point &r) {
            return operator()(a,b,r.squaredNorm());
          }

        template<typename Tparticle>
          Point force(const Tparticle &a, const Tparticle &b, double r2, const Point &p) {
            if (r2>Rc2) return Point(0,0,0);
            return lB*a.charge*b.charge*(1-Rcinv*Rcinv*r2)/(r2*sqrt(r2))*p;
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

        template<class Tparticle>
          Point force(const Tparticle &a, const Tparticle &b, double r2, const Point &p) {
            double qq=a.charge * a.charge;
            if (qq>1e-6)
              return -4*c * qq / (r2*r2) * (b.radius*b.radius*b.radius) / r2 * p;
            qq=b.charge * b.charge;
            if (qq>1e-6)
              return -4*c * qq / (r2*r2) * (a.radius*a.radius*a.radius) / r2 * p;
            return Point(0,0,0);
          }
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
     * @details Coulomb with an exponential term to describe salt screening:
     * \f[ \beta w_{ij} = \frac{e^2}{4\pi\epsilon_0\epsilon_rk_BT}
     * \frac{z_i z_j}{r_{ij}} \exp(-\kappa r_{ij}) \f]
     * where \f$\kappa=1/D\f$ is the inverse Debye screening length.
     */
    class DebyeHuckel : public Coulomb {
      private:
        string _brief();
      protected:
        double c,k,k2_count, z_count;
        Average<double> k2_count_avg;
      public:
        DebyeHuckel(InputMap&);                       //!< Construction from InputMap
        template<class Tparticle>
          double operator()(const Tparticle &a, const Tparticle &b, double r2) {
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

        template<class Tparticle>
          Point force(const Tparticle &a, const Tparticle &b, double r2, const Point &p) {
#ifdef FAU_APPROXMATH
            double rinv = invsqrtQuake(r2);
            return lB * a.charge * b.charge * rinv * exp_cawley(-k/rinv) * ( 1/r2 + k/rinv ) * p;
#else
            double r=sqrt(r2);
            return lB * a.charge * b.charge / (r*r2) * exp(-k*r) * ( 1 + k*r ) * p;
#endif
          }

        /**
         * @brief Adds counter ions to kappa
         */
        template<class Tspace>
          void setSpace(Tspace &s) {
            if (std::fabs(z_count)>1e-6) {
              double N=s.charge() / std::fabs(z_count);
              double V=s.geo.getVolume();
              double k2=k*k - k2_count; // salt contribution
              k2_count = 4*pc::pi*lB*N/V*std::pow(z_count,2); // counter ion contrib
              k=sqrt( k2+k2_count );    // total 
              k2_count_avg+=k2_count;   // sample average
            }
          } 
    };

    /**
     * DebyeHuckel ala Chung and Denton, <http://dx.doi.org/10/nkc>
     */
    class DebyeHuckelDenton : public DebyeHuckel {
      private:
        double lB_org; // original Bjerrum length without prefactor

        // Eq. 37 prefactor
        void setBjerrum(double a_m, double a_n) {
          double ka_m=k*a_m, ka_n=k*a_n;
          lB=lB_org * exp(ka_m+ka_n) / ( (1+ka_m)*(1+ka_n) );
        }

        // Eq. A7, first term
        double fmn(double m, double n) {
          return k*(m*m + n*n + k*(m+n)*m*n) / ( (1+k*m)*(1+k*n) );
        }

        // Eq. 11
        template<class Tpvec>
          double eta(Tpvec &p, double V) {
            double v=0;
            for (auto &i : p)
              v += pow(i.radius,3);
            return 4*pc::pi/(3*V) * v;
          }

        // Eq. 41 - volume energy contrib. to pressure (microions)
        // (currently salt free case, only!)
        template<class Tpvec>
          double p0(Tpvec &p, double V) {
            double Z=0,sum=0;
            for (auto &i : p) {
              Z+=i.charge; // total charge
              sum += pow(i.charge / (1+k*i.radius),2) / V; // Eq.41 sum
            }
            double p_id = std::fabs(Z)/V; // assume microion valency is unity
            double p_ex = -k*lB_org/4/(1-eta(p,V))*sum; // Eq.41, complete
            cout << "eta = " << eta(p,V) << endl;
            cout << "id = " << p_id*1660*1e3 << " ex = " << p_ex*1660*1e3 << "\n";
            return p_id + p_ex;
          }
 
      public:
        DebyeHuckelDenton(InputMap &in) : DebyeHuckel(in) {
          name+="-Denton";
          lB_org=lB;
        }

        // Effective macroion-macroion interaction (Eq. 37)
        template<class Tparticle>
          double operator()(const Tparticle &m, const Tparticle &n, double r2) {
            setBjerrum(m.radius, n.radius);
            return DebyeHuckel::operator()(m,n,r2);
          }

        template<class Tparticle>
          Point force(const Tparticle &m, const Tparticle &n, double r2, const Point &p) {
            setBjerrum(m.radius, n.radius);
            return DebyeHuckel::force(m,n,r2,p);
          }

        string info(char w) {
          lB=lB_org;
          return DebyeHuckel::info(w);
        }

        template<class Tpvec>
          string info(Tpvec &p, double V) {
            cout << "p0 = " << p0(p,V)*1660*1e3 << " mM\n";
            return info(20);
          }

        /**
         * @brief Contribution to pressure due to (dU/dk)(dk/dV), Eq.A3 / Eq.42
         */
        template<class Tpvec, class Tgeo>
          double virial(Tpvec &p, Tgeo &geo) {
            int n=int(p.size());
            double P=0, V=geo.getVolume();
            for (int i=0; i<n-1; i++)
              for (int j=i+1; j<n; j++) {
                double r2=geo.sqdist(p[i],p[j]);
                double vmn=operator()(p[i],p[j],r2);
                P += vmn * (fmn(p[i].radius, p[j].radius) - sqrt(r2)); // Eq. A5
              }
            return P * -k/(2*V*( 1-eta(p,V) )); // Eq. A3 = A5*A4
          }

        template<class Tspace>
          void setSpace(Tspace &s) {}
    };

    /**
     * @brief DebyeHuckel shifted to reach zero at given cut-off
     *
     * Shifted and truncated Yukawa potential of the form,
     * @f[
     * u(r) = u^{dh}(r) - u^{dh}(r_c)
     *         - \frac{\partial u^{dh}(r_c)}{\partial r_c}
     *         \left ( r-r_c \right )
     * @f]
     * where @f$r_c@f$ is a spherical cut-off beyond which the
     * energy is zero.
     * See more in for example <http://dx.doi.org/10/fm7qm5>.
     * The cut-off distance is read from the InputMap with the following
     * keyword:
     * - `dh_cutoff` Spherical cut-off in angstroms
     */
    class DebyeHuckelShift : public DebyeHuckel {
      private:
        double rc2,rc; // (squared) cutoff distance
        double u_rc,dudrc;
      public:
        DebyeHuckelShift(InputMap &in) : DebyeHuckel(in) {
          rc=in.get<double>("dh_cutoff",sqrt(std::numeric_limits<double>::max()));
          rc2=rc*rc;
#ifdef FAU_APPROXMATH
          u_rc = exp_cawley(-k*rc)/rc; // use approx. func even here!
#else
          u_rc = exp(-k*rc)/rc;
#endif
          dudrc = -k*u_rc - u_rc/rc; // 1st derivative of u(r) at r_c
          std::ostringstream o;
          o << " (shifted, rcut=" << rc << textio::_angstrom << ")";
          name+=o.str();
        }

        template<class Tparticle>
          inline double operator() (const Tparticle &a, const Tparticle &b, double r2) const {
            if (r2>rc2)
              return 0;
#ifdef FAU_APPROXMATH
            double r = 1./invsqrtQuake(r2);
            return lB * a.charge * b.charge
              * ( exp_cawley(-k*r)/r - u_rc - (dudrc*(r-rc)) );
#else
            double r=sqrt(r2);
            return lB * a.charge * b.charge
              * ( exp(-k*r)/r - u_rc - (dudrc*(r-rc))  );
#endif
          }
        template<class Tparticle>
          Point force(const Tparticle &a, const Tparticle &b, double r2, const Point &p) {
            if (r2>rc2)
              return Point(0,0,0);
#ifdef FAU_APPROXMATH
            double rinv = invsqrtQuake(r2);
            return lB * a.charge * b.charge * ( exp_cawley(-k/rinv) / r2 * (k+rinv) + dudrc*rinv) * p;
#else
            double r=sqrt(r2);
            return lB * a.charge * b.charge * ( exp(-k*r) / r2 * (k + 1/r) + dudrc / r) * p;
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
          typedef std::function<Point(const Tparticle&,const Tparticle&,double,const Point&)> Tforce;
          std::map<Tpair,Tfunc> m;
          std::map<Tpair,Tforce> mforce;
          std::string _info; // info for the added potentials (before turning into functors)

          // Force function object wrapper class
          template<class Tpairpot>
            struct ForceFunctionObject {
              Tpairpot pot;
              ForceFunctionObject(const Tpairpot &p) : pot(p) {}
              Point operator()(const particle &a, const particle &b, double r2, const Point &r) {
                return pot.force(a,b,r2,r);
              }
            };

        public:
          PotentialMap(InputMap &in) : Tdefault(in) {
            Tdefault::name += " (default)";
          } 

          /**
           * @warning Templating `Tpairpot` may cause code-bloat if used on
           * many different potentials
           */
          template<class Tpairpot>
            void add(AtomData::Tid id1, AtomData::Tid id2, Tpairpot pot) {
              pot.name=atom[id1].name + "<->" + atom[id2].name + ": " + pot.name;
              _info+="\n  " + pot.name + ":\n" + pot.info(20);
              m[Tpair(id1,id2)] = pot;
              mforce[Tpair(id1,id2)] = ForceFunctionObject<decltype(pot)>(pot);
            }

          double operator()(const Tparticle &a, const Tparticle &b, const Tdist &r2) {
            auto i=m.find( Tpair(a.id,b.id) );
            if (i!=m.end())
              return i->second(a,b,r2);
            return Tdefault::operator()(a,b,r2);
          }

          Point force(const Tparticle &a, const Tparticle &b, double r2, const Point &p) {
            auto i=mforce.find( Tpair(a.id,b.id) );
            if (i!=mforce.end())
              return i->second(a,b,r2,p);
            return Tdefault::force(a,b,r2,p);
          }

          std::string info(char w=20) {
            return Tdefault::info(w) + _info;
          }
      };

    /**
     * @brief Base class for manybody potentials (dihedrals etc.)
     *
     * This is an inter-particle potential for interactions involving
     * many particles such as angles and dihedrals.
     * It works much in the same way as pair potentials, but instead
     * of taking particles as input, it maintains a list of particle
     * index that interact.
     * Derived classes are expected to implement a function operator
     * that takes as arguments a geometry and particle vector.
     */
    class ManybodyBase {
      protected:
        virtual string _brief() const=0;
        typedef vector<int> Tivec;
        Tivec v; // particle index
        string ndxStr() const {
          std::ostringstream o;
          for (auto i : v)
            o << i << " ";
          return "[ " + o.str() + "]";
        }
      public:
        string brief() const { return _brief(); }
        void setIndex(const Tivec &index) { v = index; }
        Tivec& getIndex() { return v; }
        Tivec getIndex() const { return v; }
        virtual ~ManybodyBase() {}
    };

    /**
     * @brief Angular potential
     * @warning unfinished!
     */
    class Angular : public ManybodyBase {
      private:
        static const string name;
        typedef ManybodyBase base;
        double c1,c2; // pot parameters
        string _brief() const FOVERRIDE {
          std::ostringstream o;
          o << ndxStr() << " - " << name << ": " << c1 << " " << c2;
          return o.str();
        }
      public:
        Angular() : c1(0), c2(0) {}
        Angular(const Tivec& ndx, double angle, double energy) {
          assert(ndx.size()==3);
          setIndex(ndx);
        }

        /** @brief Energy function (kT) */
        template<class Tgeometry, class Tpvec>
          double operator()(Tgeometry &g, const Tpvec &p) {
            assert(v.size()==3);
            assert((int)p.size()>*std::max_element(v.begin(), v.end()));
            auto ab = g.vdist(p[v[0]], p[v[1]]);
            auto ac = g.vdist(p[v[0]], p[v[2]]);
            auto bc = g.vdist(p[v[1]], p[v[2]]);
            return c1 * ab.dot(bc) / ac.norm(); // for example...
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
          void setCutoff() {
            for (size_t i=0; i<atom.list.size(); i++)
              for (size_t j=0; j<atom.list.size(); j++) {
                if (first.rcut2(i,j) > second.rcut2(i,j))
                  PairPotentialBase::rcut2.set(i,j,first.rcut2(i,j));
                else
                  PairPotentialBase::rcut2.set(i,j,second.rcut2(i,j));
              }
          }

        public:
          T1 first;  //!< First pair potential of type T1
          T2 second; //!< Second pair potential of type T2

          CombinedPairPotential(T1 a, T2 b) : first(a), second(b) {
            name=first.name+"+"+second.name;
            setCutoff();
          }

          CombinedPairPotential(InputMap &in) : first(in), second(in) {
            name=first.name+"+"+second.name;
            setCutoff();
          }

          CombinedPairPotential(InputMap &in, string pfx1, string pfx2) :
            PairPotentialBase(pfx1+pfx2), first(in,pfx1), second(in,pfx2) {
              name=first.name+"+"+second.name;
              setCutoff();
            }

          template<class Tparticle, class Tdist>
            double operator()(const Tparticle &a, const Tparticle &b, const Tdist &r2) {
              return first(a,b,r2) + second(a,b,r2);
            }

          template<typename Tparticle>
            Point force(const Tparticle &a, const Tparticle &b, double r2, const Point &p) {
              return first.force(a,b,r2,p) + second.force(a,b,r2,p);
            }

          template<typename Tparticle>
            Point field(const Tparticle &a, const Point &r) {
              return first.field(a,r) + second.field(a,r);
            }

          template<class Tspace>
            void setSpace(Tspace &s) {
              first.setSpace(s);
              second.setSpace(s);
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
     * @brief Creates a new pair potential with opposite sign
     */
    template<class T>
      struct Minus : public T {
        Minus(const T &pot) : T(pot) { T::name = "minus " + T::name; }
        template<class Tparticle> // isotropic energy
          double operator()(const Tparticle &a, const Tparticle &b, double r2) {
            return -T::operator()(a,b,r2);
          }
        template<class Tparticle> // anisotropic energy
          double operator()(const Tparticle &a, const Tparticle &b, const Point &r2) {
            return -T::operator()(a,b,r2);
          }
        template<class Tparticle> // force
          Point force(const Tparticle &a, const Tparticle &b, double r2, const Point &p) {
            return -T::force(a,b,r2,p);
          }
      };

    /**
     * @brief Adds two pair potentials
     *
     * Example:
     *
     *     auto PrimitiveModel = Potential::Coulomb(...) + Potential::HardSphere(...);
     */
    template<class T1, class T2,
      class = typename std::enable_if<std::is_base_of<PairPotentialBase,T1>::value>::type,
      class = typename std::enable_if<std::is_base_of<PairPotentialBase,T2>::value>::type>
        CombinedPairPotential<T1,T2>& operator+(const T1 &pot1, const T2 &pot2) {
          return *(new CombinedPairPotential<T1,T2>(pot1,pot2));
        }

    /**
     * @brief Subtracts two pair potentials
     *
     * This can be useful for excluding non-bonded interactions between bonded pairs.
     * Beware, though, that in the case of strongly repulsive interactions
     * for example due to particle overlap,
     * first adding then subtracting
     * may lead to numerical issues, often manifested in a system energy
     * drift.
     *
     * Example:
     *
     *     auto mypot = Potential::Harmonic(...) - Potential::LennardJones(...);
     */
    template<class T1, class T2,
      class = typename std::enable_if<std::is_base_of<PairPotentialBase,T1>::value>::type,
      class = typename std::enable_if<std::is_base_of<PairPotentialBase,T2>::value>::type>
        CombinedPairPotential<T1,Minus<T2>>& operator-(const T1 &pot1, const T2 &pot2) {
          return *(new CombinedPairPotential<T1,Minus<T2>>(pot1,Minus<T2>(pot2)));
        }

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
