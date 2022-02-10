#pragma once
#include "potentials_base.h"
#include "tabulate.h"
#include "functionparser.h"
#include "multipole.h"
#include "spherocylinder.h"
#include <coulombgalore.h>

namespace Faunus::Potential {

class CosAttract;
class WeeksChandlerAndersen;

using CigarCigarCosAttractWCA = PatchyCigarCigar<CosAttract, WeeksChandlerAndersen>;

struct Dummy : public PairPotentialBase {
    Dummy();
    inline double operator()(const Particle &, const Particle &, double, const Point &) const override { return 0.0; }
    void from_json(const json &) override;
    void to_json(json &) const override;
}; //!< A dummy pair potential that always returns zero

/**
 * @brief Lennard-Jones potential with an arbitrary combination rule.
 * @note Mixing data is _shared_ upon copying
 */
class LennardJones : public MixerPairPotentialBase {
    TExtractorFunc extract_sigma, extract_epsilon;

  protected:
    TPairMatrixPtr sigma_squared;     // sigma_ij * sigma_ij
    TPairMatrixPtr epsilon_quadruple; // 4 * epsilon_ij
    void initPairMatrices() override;
    void extractorsFromJson(const json &j) override;

  public:
    LennardJones(const std::string& name = "lennardjones", const std::string& cite = std::string(),
                 CombinationRuleType combination_rule = CombinationRuleType::LORENTZ_BERTHELOT);

    /**
     * @brief Calculates force on particle a due to another particle, b
     * @param particle_a Particle a ("target")
     * @param particle_b Particle b
     * @param squared_distance Squared norm |𝐚-𝐛|²
     * @param b_towards_a Distance vector 𝐛 -> 𝐚 = 𝐚 - 𝐛
     * @return Force on particle a due to particle b
     */
    inline Point force(const Particle& a, const Particle& b, double squared_distance,
                       const Point& b_towards_a) const override {
        const auto s6 = powi((*sigma_squared)(a.id, b.id), 3);
        const auto r6 = squared_distance * squared_distance * squared_distance;
        const auto r14 = r6 * r6 * squared_distance;
        return 6.0 * (*epsilon_quadruple)(a.id, b.id) * s6 * (2.0 * s6 - r6) / r14 * b_towards_a; // force in a
    }

    inline double operator()(const Particle &a, const Particle &b, double r2, const Point &) const override {
        double x = (*sigma_squared)(a.id, b.id) / r2;              // s2/r2
        x = x * x * x;                                             // s6/r6
        return (*epsilon_quadruple)(a.id, b.id) * (x * x - x);
    }
};

/**
 * @brief Weeks-Chandler-Andersen pair potential
 * @details This is a Lennard-Jones type potential, cut and shifted to zero
 * at @f$r_c=2^{1/6}\sigma@f$. More info can be found in at
 * <http://doi.org/ct4kh9> and the functional form is:
 * @f[
 * \beta u = 4 \beta \epsilon \left ( (b/r)^{12} - (b/r)^6 + \frac{1}{4} \right )
 * @f]
 * where sigma, epsilon per default are set
 * using Lorentz-Berthelot mixing rules.
 *
 * @note Mixing data is _shared_ upon copying
 */
class WeeksChandlerAndersen : public LennardJones {
    static constexpr double onefourth = 0.25, twototwosixth = 1.2599210498948732;

    inline double operator()(const Particle &a, const Particle &b, double r2) const {
        double x = (*sigma_squared)(a.id, b.id); // s^2
        if (r2 > x * twototwosixth)
            return 0;
        x = x / r2;    // (s/r)^2
        x = x * x * x; // (s/r)^6
        return (*epsilon_quadruple)(a.id, b.id) * (x * x - x + onefourth);
    }

  public:
    WeeksChandlerAndersen(const std::string& name = "wca", const std::string& cite = "doi:ct4kh9",
                          CombinationRuleType combination_rule = CombinationRuleType::LORENTZ_BERTHELOT);

    inline double operator()(const Particle &a, const Particle &b, double r2, const Point &) const override {
        return operator()(a, b, r2);
    }

    inline Point force(const Particle &a, const Particle &b, const double r2, const Point &p) const override {
        auto x = (*sigma_squared)(a.id, b.id); // s^2
        if (r2 > x * twototwosixth) {
            return {0.0, 0.0, 0.0};
        }
        x = x / r2;    // (s/r)^2
        x = x * x * x; // (s/r)^6
        return (*epsilon_quadruple)(a.id, b.id) * 6.0 * (2.0 * x * x - x) / r2 * p;
    }
}; // Weeks-Chandler-Andersen potential

/**
 * @brief Hardsphere potential
 *
 * Uses arithmetic mean for sigma as a default combination rule.
 * @note `PairMatrix` is _shared_ upon copying
 */
class HardSphere : public MixerPairPotentialBase {
    TExtractorFunc extract_sigma;

  protected:
    TPairMatrixPtr sigma_squared; // sigma_ij * sigma_ij
    void initPairMatrices() override;
    void extractorsFromJson(const json &j) override;

  public:
    HardSphere(const std::string& name = "hardsphere");

    inline double operator()(const Particle &a, const Particle &b, double r2, const Point &) const override {
        return r2 < (*sigma_squared)(a.id, b.id) ? pc::infty : 0.0;
    }
};

/**
 * @brief Hertz potential
 * @details This is a repulsive potential, that for example, describes the change in elastic energy
 * of two deformable objects when subjected to an axial compression.
 * @f[
 *     u(r) = \epsilon \left(1 - \frac{r}{\sigma}\right)^{5/2}
 * @f]
 * where r_H corresponds to the particle's radius.
 *
 * More info: doi:10.1063/1.3186742
 */
class Hertz : public MixerPairPotentialBase {
    TExtractorFunc extract_sigma, extract_epsilon;

  protected:
    TPairMatrixPtr sigma_squared; // sigma_ij * sigma_ij
    TPairMatrixPtr epsilon;       // epsilon_ij
    void initPairMatrices() override;
    void extractorsFromJson(const json &j) override;

  public:
    Hertz(const std::string &name = "hertz");
    inline double operator()(const Particle &a, const Particle &b, double r2, const Point &) const override {
        if (r2 <= (*sigma_squared)(a.id, b.id)) {
            return (*epsilon)(a.id, b.id) * pow((1 - (sqrt(r2 / (*sigma_squared)(a.id, b.id)))), 2.5);
        }
        return 0.0;
    }
};

/**
 * @brief Square-well potential
 * @details This is an attractive potential described by
 * @f[
 *     u(r) = -\epsilon
 * @f]
 * when r < \sigma, and zero otherwise.
 */
class SquareWell : public MixerPairPotentialBase {
    TExtractorFunc extract_sigma, extract_epsilon;

  protected:
    TPairMatrixPtr sigma_squared; // sigma_ij * sigma_ij
    TPairMatrixPtr epsilon;       // epsilon_ij
    void extractorsFromJson(const json &j) override;
    void initPairMatrices() override;

  public:
    SquareWell(const std::string &name = "squarewell");
    inline double operator()(const Particle &a, const Particle &b, double r2, const Point &) const override {
        return (r2 < (*sigma_squared)(a.id, b.id)) ? -(*epsilon)(a.id, b.id) : 0.0;
    }
};

struct RepulsionR3 : public PairPotentialBase {
    double f = 0, s = 0, e = 0;

    RepulsionR3(const std::string &name = "repulsionr3") : PairPotentialBase(name) {};
    void from_json(const json &j) override;
    void to_json(json &j) const override;

    inline double operator()(const Particle &, const Particle &, double r2, const Point &) const override {
        const auto r = sqrt(r2);
        return f / (r * r2) + e * std::pow(s / r, 12);
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
 * JSON keywords:
 *
 * Key     | Description
 * :-------| :---------------------------
 * `eps`   | Depth, \f$\epsilon\f$ [kJ/mol]
 * `rc`    | Width, r_c [angstrom]
 * `wc`    | Decay range, w_c [angstrom]
 *
 */
class CosAttract : public PairPotentialBase {
    double eps, wc, rc, rc2, c, rcwc2;

  public:
    CosAttract(const std::string &name = "cos2");
    double cutOffSquared() const; //!< Squared cutoff distance where potential goes to zero

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
    inline double operator()(const Particle &, const Particle &, double r2, const Point &) const override {
        if (r2 < rc2) {
            return -eps;
        }
        if (r2 > rcwc2) {
            return 0;
        }
        const auto x = std::cos(c * (sqrt(r2) - rc));
        return -eps * x * x;
    }

    inline Point force(const Particle &, const Particle &, double r2, const Point &p) const override {
        if (r2 > rcwc2 || r2 < rc2) {
            return {0.0, 0.0, 0.0};
        }
        const auto r = sqrt(r2);
        const auto x1 = std::cos(c * (r - rc));
        const auto x2 = std::sin(c * (r - rc));
        return -2.0 * c * eps * x1 * x2 / r * p;
    }

    void to_json(json &j) const override;
    void from_json(const json &j) override;
};

/**
 * @brief Pairwise SASA potential calculating the surface area of inter-secting spheres
 */
class SASApotential : public PairPotentialBase {
  private:
    bool shift = true; // shift potential to zero at large separations?
    double proberadius = 0, conc = 0;

    double area(double R, double r, double d_squared)
    const; //!< Total surface area of two intersecting spheres or radii R and r as a function of separation

  public:
    inline double operator()(const Particle &a, const Particle &b, const double r2, const Point &) const override {
        const auto tfe = 0.5 * (atoms[a.id].tfe + atoms[b.id].tfe);
        const auto tension = 0.5 * (atoms[a.id].tension + atoms[b.id].tension);
        if (fabs(tfe) > 1e-6 or fabs(tension) > 1e-6)
            return (tension + conc * tfe) * area(0.5 * atoms[a.id].sigma, 0.5 * atoms[b.id].sigma, r2);
        return 0.0;
    }
    SASApotential(const std::string &name = "sasa", const std::string &cite = std::string());
    void to_json(json &j) const override;
    void from_json(const json &j) override;
};

/**
 * @brief Plain Coulomb potential
 */
struct Coulomb : public PairPotentialBase {
    Coulomb(const std::string &name = "coulomb");
    double bjerrum_length = 0.0; //!< Bjerrum length
    inline double operator()(const Particle &a, const Particle &b, const double r2, const Point &) const override {
        return bjerrum_length * a.charge * b.charge / std::sqrt(r2);
    }
    void to_json(json &j) const override;
    void from_json(const json &j) override;
};

struct DipoleDipole : public PairPotentialBase {
    DipoleDipole(const std::string &name = "dipoledipole", const std::string &cite = std::string());
    double lB; //!< Bjerrum length
    inline double operator()(const Particle &a, const Particle &b, double, const Point &r) const override {
        return lB*mu2mu(a.getExt().mu, b.getExt().mu, a.getExt().mulen*b.getExt().mulen, r,1.0,0.0);
    }
    void to_json(json &j) const override;
    void from_json(const json &j) override;
};

/**
 * @brief Charge-nonpolar pair interaction
 * @note Pair data is _shared_ upon copying
 */
class Polarizability : public Coulomb {
  private:
    double epsr;
    std::shared_ptr<PairMatrix<double>> m_neutral, m_charged;

  public:
    Polarizability(const std::string &name = "polar");
    void from_json(const json &j) override;
    void to_json(json &j) const override;

    inline double operator()(const Particle &a, const Particle &b, double r2, const Point &) const override {
        double r4inv = 1 / (r2 * r2);
        if (fabs(a.charge) > 1e-9 or fabs(b.charge) > 1e-9)
            return (*m_charged)(a.id, b.id) * r4inv;
        else
            return (*m_neutral)(a.id, b.id) / r2 * r4inv;
    }

    inline Point force(const Particle &a, const Particle &b, double r2, const Point &p) const override {
        double r6inv = 1 / (r2 * r2 * r2);
        if (fabs(a.charge) > 1e-9 or fabs(b.charge) > 1e-9)
            return 4 * m_charged->operator()(a.id, b.id) * r6inv * p;
        else
            return 6 * m_neutral->operator()(a.id, b.id) / r2 * r6inv * p;
    }
};

/**
 * @brief Finite Extensible Nonlinear Elastic (FENE) potential
 *
 * This is an anharmonic bonding potential with the form:
 * @f[
 *     \beta u(r) = -\frac{k r_0^2}{2}\ln \left [ 1-(r/r_0)^2 \right ]
 * @f]
 * for \f$r<r_0\f$, otherwise infinity. JSON keywords:
 *
 * - `stiffness` Bond stiffness, `k` [kT]
 * - `maxsep` Maximum separation, `r_0` [angstrom]
 *
 * More info: doi:10.1103/PhysRevE.59.4248
 */
class FENE : public PairPotentialBase {
    double k, r02, r02inv;

  public:
    FENE(const std::string &name = "fene");
    void from_json(const json &j) override;
    void to_json(json &j) const override;

    /** @brief Energy in kT between two particles, r2 = squared distance */
    inline double operator()(const Particle &, const Particle &, double r2, const Point &) const override {
        return (r2 > r02) ? pc::infty : -0.5 * k * r02 * std::log(1 - r2 * r02inv);
    }

    inline Point force(const Particle &, const Particle &, double r2, const Point &p) const override {
        return (r2 > r02) ? -pc::infty * p : -k * r02 / (r02 - r2) * p;
    }
};

/**
 * @brief Wrapper for external CoulombGalore library
 */
class NewCoulombGalore : public PairPotentialBase {
  protected:
    ::CoulombGalore::Splined pot;
    virtual void setSelfEnergy();

  public:
    NewCoulombGalore(const std::string & = "coulomb");
    inline double operator()(const Particle &a, const Particle &b, const double r2, const Point &) const override {
        return bjerrum_length *
               pot.ion_ion_energy(a.charge, b.charge, sqrt(r2) + std::numeric_limits<double>::epsilon());
    }
    Point force(const Particle& particle_a, const Particle& particle_b, double squared_distance,
                const Point& b_towards_a) const override;
    void from_json(const json &) override;
    void to_json(json &) const override;
    double dielectric_constant(double M2V);
    double bjerrum_length; //!< Bjerrum length (angstrom)
    const CoulombGalore::Splined& getCoulombGalore() const; //!< Access to full coulomb galore class
};

/**
 * @brief Multipole interactions
 * @note Only dipole-dipole interactions are currently enabled
 */
class Multipole : public NewCoulombGalore {
  private:
    void setSelfEnergy() override;

  public:
    Multipole(const std::string & = "multipole");
    inline double operator()(const Particle &a, const Particle &b, double, const Point &r) const override {
        // Only dipole-dipole for now!
        Point mua = a.getExt().mu * a.getExt().mulen;
        Point mub = b.getExt().mu * b.getExt().mulen;
        double dipdip = pot.dipole_dipole_energy(mua, mub, r);
        return bjerrum_length * (dipdip);
    }

    Point force(const Particle &, const Particle &, double, const Point &) const override;
};

/**
 * @brief Custom pair-potential taking math. expressions at runtime
 * @note `symbols` is a shared_ptr as this allows it to be modified by `operator() const`.
 *       A hack, but would otherwise require const-removal in all pair-potentials.
 */
class CustomPairPotential : public PairPotentialBase {
  private:
    // Only ExprFunction<double> is explicitly instantiated in functionparser.cpp. Other types as well as
    // the implicit template instantiation is disabled to save reasources during the compilation/build.
    ExprFunction<double> expression;
    struct Symbols {
        double distance = 0.0; // available as "r"
        double charge1 = 0.0;  // available as "charge1"
        double charge2 = 0.0;  // available as "charge2"
        double sigma1 = 0.0;   // available as "s1"
        double sigma2 = 0.0;   // available as "s2"
    };
    std::shared_ptr<Symbols> symbols;
    double squared_cutoff_distance;
    json original_input;

    inline void setSymbols(const Particle& particle1, const Particle& particle2, double squared_distance) const {
        symbols->distance = std::sqrt(squared_distance);
        symbols->charge1 = particle1.charge;
        symbols->charge2 = particle2.charge;
        symbols->sigma1 = Faunus::atoms[particle1.id].sigma;
        symbols->sigma2 = Faunus::atoms[particle2.id].sigma;
    }

  public:
    inline double operator()(const Particle& particle1, const Particle& particle2, double squared_distance,
                             [[maybe_unused]] const Point& r) const override {
        if (squared_distance < squared_cutoff_distance) {
            setSymbols(particle1, particle2, squared_distance);
            return expression();
        }
        return 0.0;
    }

    inline Point force(const Particle& particle1, const Particle& particle2, double squared_distance,
                       const Point& r) const override {
        if (squared_distance < squared_cutoff_distance) {
            setSymbols(particle1, particle2, squared_distance);
            return -expression.derivative(symbols->distance) / symbols->distance * r;
        }
        return Point::Zero();
    }

    CustomPairPotential(const std::string &name = "custom");

    void from_json(const json &j) override;
    void to_json(json &j) const override;
};


/**
 * @brief Arbitrary potentials for specific atom types
 *
 * This maintains a species x species matrix with function pointers (`std::function`)
 * that wraps pair potentials. Flexibility over performance.
 *
 * @todo `to_json` should retrieve info from potentials instead of merely passing input
 * @warning Each atom pair will be assigned an instance of a pair-potential. This *could* be
 *          problematic if these have large memory requirements.
 */
class FunctorPotential : public PairPotentialBase {
    using uFunc = std::function<double (const Particle &, const Particle &, double, const Point &)>;
    using PrimitiveModel = CombinedPairPotential<Coulomb, HardSphere>;
    using PrimitiveModelWCA = CombinedPairPotential<Coulomb, WeeksChandlerAndersen>;
    json _j; // storage for input json
    bool have_monopole_self_energy = false;
    bool have_dipole_self_energy = false;
    void registerSelfEnergy(PairPotentialBase *); //!< helper func to add to selv_energy_vector

    // List of pair-potential instances used when constructing functors.
    // Note that potentials w. large memory requirements (LJ, WCA etc.)
    // typically use `shared_ptr` so that the created functors _share_
    // the data. That is *only* put the pair-potential here if you can
    // share internal (shared) pointers.
    std::tuple<NewCoulombGalore,      // 0
               CosAttract,            // 1
               Polarizability,        // 2
               HardSphere,            // 3
               LennardJones,          // 4
               RepulsionR3,           // 5
               SASApotential,         // 6
               WeeksChandlerAndersen, // 7
               PrimitiveModel,        // 8
               PrimitiveModelWCA,     // 9
               Hertz,                 // 10
               SquareWell,            // 11
               Multipole,             // 12
               HardSpheroCylinder,     // 13
               CigarCigarCosAttractWCA // 14
               >
        potlist;

    uFunc combineFunc(json &j); // parse json array of potentials to a single potential function object

  protected:
    PairMatrix<uFunc, true> umatrix; // matrix with potential for each atom pair; cannot be Eigen matrix

  public:
    FunctorPotential(const std::string &name = "functor potential");
    void to_json(json &j) const override;
    void from_json(const json &j) override;

    /**
     * @brief Potential energy between two particles
     * @param a First particle
     * @param b Second particle
     * @param r2 Squared distance
     * @param r Distance vector
     * @return Energy in kT
     */
    inline double operator()(const Particle &a, const Particle &b, double r2,
                             const Point &r = {0, 0, 0}) const override {
        return umatrix(a.id, b.id)(a, b, r2, r);
    }
};

/**
 * @brief Splined pair potentials
 *
 * This maintains a species x species matrix as in `FunctorPotential`
 * but with splined pair potentials. This avoids the functor lookup
 * and renders all potentials roughly the same speed.
 *
 * The spline range is automatically detected based on user-defined
 * energy thresholds. If below the range, the default behavior is to return
 * the EXACT energy, while if above ZERO is returned.
 */
class SplinedPotential : public FunctorPotential {
    /** @brief Expand spline data class to hold information about the sign of values for r<rmin */
    class KnotData : public Tabulate::TabulatorBase<double>::data {
      public:
        using base = Tabulate::TabulatorBase<double>::data;
        bool hardsphere_repulsion = false; //!< Use hardsphere repulsion for r smaller than rmin
        KnotData() = default;
        KnotData(const base &);
    };

    PairMatrix<KnotData> matrix_of_knots;                 //!< Matrix with tabulated potential for each atom pair
    Tabulate::Andrea<double> spline;                      //!< Spline method
    bool hardsphere_repulsion = false;                    //!< Use hardsphere repulsion for r smaller than rmin
    const int max_iterations = 1e6;                       //!< Max number of iterations when determining spline interval
    void stream_pair_potential(std::ostream&, size_t id1, size_t id2); //!< Stream pair potential to output stream
    void save_potentials();                               //!< Save splined and exact pair potentials to disk
    double findLowerDistance(int, int, double, double);   //!< Find lower distance for splining (rmin)
    double findUpperDistance(int, int, double, double);   //!< Find upper distance for splining (rmax)
    double dr = 1e-2;                                     //!< Distance interval when searching for rmin and rmax
    void createKnots(int, int, double, double);           //!< Create spline knots for pair of particles in [rmin:rmax]

  public:
    explicit SplinedPotential(const std::string &name = "splined");

    /**
     * Policies:
     *
     * 1. return splined potential if rmin>r<rmax
     * 2. return zero if r>=rmax
     * 3. return infinity if r<=rmin AND `hardsphere_repulsion` has been set to true
     * 4. return exact energy if r<=rmin
     */
    inline double operator()(const Particle &p1, const Particle &p2, double r2, const Point &) const override {
        const auto &knots = matrix_of_knots(p1.id, p2.id);
        if (r2 >= knots.rmax2) {
            return 0.0;
        }
        if (r2 > knots.rmin2) {
            return spline.eval(knots, r2); // spline energy
        }
        if (knots.hardsphere_repulsion) {
            return pc::infty;
        }
        return FunctorPotential::operator()(p1, p2, r2, {0, 0, 0}); // exact energy
    }

    void from_json(const json &) override;
};

} // end of namespace