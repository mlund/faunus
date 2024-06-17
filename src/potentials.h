#pragma once
#include "potentials_base.h"
#include "tabulate.h"
#include "functionparser.h"
#include "multipole.h"
#include "spherocylinder.h"
#include <coulombgalore.h>

namespace Faunus::pairpotential {

/**
 * @brief Lennard-Jones potential with an arbitrary combination rule.
 * @note Mixing data is _shared_ upon copying
 */
class LennardJones : public MixerPairPotentialBase
{
  private:
    TExtractorFunc extract_sigma;
    TExtractorFunc extract_epsilon;

  protected:
    TPairMatrixPtr sigma_squared;     // sigma_ij * sigma_ij
    TPairMatrixPtr epsilon_quadruple; // 4 * epsilon_ij
    void initPairMatrices() override;
    void extractorsFromJson(const json& j) override;

  public:
    explicit LennardJones(
        const std::string& name = "lennardjones", const std::string& cite = std::string(),
        CombinationRuleType combination_rule = CombinationRuleType::LORENTZ_BERTHELOT);

    inline Point force(const Particle& particle_a, const Particle& particle_b,
                       double squared_distance, const Point& b_towards_a) const override
    {
        const auto s6 = powi((*sigma_squared)(particle_a.id, particle_b.id), 3);
        const auto r6 = squared_distance * squared_distance * squared_distance;
        const auto r14 = r6 * r6 * squared_distance;
        return 6.0 * (*epsilon_quadruple)(particle_a.id, particle_b.id) * s6 * (2.0 * s6 - r6) /
               r14 * b_towards_a; // force in a
    }

    inline double operator()(const Particle& particle_a, const Particle& particle_b,
                             double squared_distance,
                             [[maybe_unused]] const Point& b_towards_a) const override
    {
        auto x = (*sigma_squared)(particle_a.id, particle_b.id) / squared_distance; // s2/r2
        x = x * x * x;                                                              // s6/r6
        return (*epsilon_quadruple)(particle_a.id, particle_b.id) * (x * x - x);
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
class WeeksChandlerAndersen : public LennardJones
{
    static constexpr double onefourth = 0.25, twototwosixth = 1.2599210498948732;

    inline double operator()(const Particle& a, const Particle& b, double squared_distance) const
    {
        auto x = (*sigma_squared)(a.id, b.id); // s^2
        if (squared_distance > x * twototwosixth) {
            return 0;
        }
        x = x / squared_distance; // (s/r)^2
        x = x * x * x;            // (s/r)^6
        return (*epsilon_quadruple)(a.id, b.id) * (x * x - x + onefourth);
    }

  public:
    explicit WeeksChandlerAndersen(
        const std::string& name = "wca", const std::string& cite = "doi:ct4kh9",
        CombinationRuleType combination_rule = CombinationRuleType::LORENTZ_BERTHELOT);

    inline double operator()(const Particle& a, const Particle& b, double squared_distance,
                             [[maybe_unused]] const Point& b_towards_a) const override
    {
        return operator()(a, b, squared_distance);
    }

    inline Point force(const Particle& a, const Particle& b, const double squared_distance,
                       const Point& b_towards_a) const override
    {
        auto x = (*sigma_squared)(a.id, b.id); // s^2
        if (squared_distance > x * twototwosixth) {
            return {0.0, 0.0, 0.0};
        }
        x = x / squared_distance; // (s/r)^2
        x = x * x * x;            // (s/r)^6
        return (*epsilon_quadruple)(a.id, b.id) * 6.0 * (2.0 * x * x - x) / squared_distance *
               b_towards_a;
    }
}; // Weeks-Chandler-Andersen potential

/**
 * @brief Hardsphere potential
 *
 * Uses arithmetic mean for sigma as a default combination rule.
 * @note `PairMatrix` is _shared_ upon copying
 */
class HardSphere : public MixerPairPotentialBase
{
    TExtractorFunc extract_sigma;
    TPairMatrixPtr sigma_squared; // sigma_ij * sigma_ij
    void initPairMatrices() override;
    void extractorsFromJson(const json& j) override;

  public:
    explicit HardSphere(const std::string& name = "hardsphere");

    inline double operator()(const Particle& particle_a, const Particle& particle_b,
                             double squared_distance, const Point&) const override
    {
        return squared_distance < (*sigma_squared)(particle_a.id, particle_b.id) ? pc::infty : 0.0;
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
class Hertz : public MixerPairPotentialBase
{
    TExtractorFunc extract_sigma, extract_epsilon;

  protected:
    TPairMatrixPtr sigma_squared; // sigma_ij * sigma_ij
    TPairMatrixPtr epsilon;       // epsilon_ij
    void initPairMatrices() override;
    void extractorsFromJson(const json& j) override;

  public:
    explicit Hertz(const std::string& name = "hertz");

    inline double operator()(const Particle& a, const Particle& b, double squared_distance,
                             [[maybe_unused]] const Point& b_towards_a) const override
    {
        if (squared_distance <= (*sigma_squared)(a.id, b.id)) {
            return (*epsilon)(a.id, b.id) *
                   pow((1 - (sqrt(squared_distance / (*sigma_squared)(a.id, b.id)))), 2.5);
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
class SquareWell : public MixerPairPotentialBase
{
    TExtractorFunc extract_sigma, extract_epsilon;

  protected:
    TPairMatrixPtr sigma_squared; // sigma_ij * sigma_ij
    TPairMatrixPtr epsilon;       // epsilon_ij
    void extractorsFromJson(const json& j) override;
    void initPairMatrices() override;

  public:
    explicit SquareWell(const std::string& name = "squarewell");

    inline double operator()(const Particle& a, const Particle& b, double squared_distance,
                             const Point&) const override
    {
        return (squared_distance < (*sigma_squared)(a.id, b.id)) ? -(*epsilon)(a.id, b.id) : 0.0;
    }
};

class RepulsionR3 : public PairPotential
{
  private:
    double f = 0, s = 0, e = 0;
    void from_json(const json& j) override;
    void to_json(json& j) const override;

  public:
    explicit RepulsionR3(const std::string& name = "repulsionr3")
        : PairPotential(name) {};

    inline double operator()(const Particle&, const Particle&, double squared_distance,
                             const Point&) const override
    {
        const auto r = sqrt(squared_distance);
        return f / (r * squared_distance) + e * std::pow(s / r, 12);
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
class CosAttract : public PairPotential
{
    double eps = 0.0;
    double wc = 0.0;
    double rc = 0.0;
    double rc2 = 0.0;
    double c = 0.0;
    double rcwc2 = 0.0; // (rc + wc)^2 ~ "rcut2" in faunus v1
    void from_json(const json& j) override;

  public:
    explicit CosAttract(const std::string& name = "cos2");
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
    inline double operator()(const Particle&, const Particle&, double squared_distance,
                             const Point&) const override
    {
        if (squared_distance < rc2) {
            return -eps;
        }
        if (squared_distance > rcwc2) {
            return 0;
        }
        const auto x = std::cos(c * (sqrt(squared_distance) - rc));
        return -eps * x * x;
    }

    inline Point force(const Particle&, const Particle&, double squared_distance,
                       const Point& b_towards_a) const override
    {
        if (squared_distance > rcwc2 || squared_distance < rc2) {
            return {0.0, 0.0, 0.0};
        }
        const auto r = sqrt(squared_distance);
        const auto x1 = std::cos(c * (r - rc));
        const auto x2 = std::sin(c * (r - rc));
        return -2.0 * c * eps * x1 * x2 / r * b_towards_a;
    }

    void to_json(json& j) const override;
};

/**
 * @brief Cosine attraction using combination rules (Eq. 4 in doi:10/chqzjk)
 *
 * This will collect `rc`, `wc`, and `eps` from the atom topology and mix using arbitrary
 * combination rules. `EPSILON` is used to mix the energy, `eps`. `SIGMA` is used to mix
 * distances, `rc` and `wc`.
 */
class CosAttractMixed : public MixerPairPotentialBase
{
  private:
    TExtractorFunc extract_rc;
    TExtractorFunc extract_wc;
    TExtractorFunc extract_eps;

    TPairMatrixPtr switching_distance; //!< Switching region begins here (r_c)
    TPairMatrixPtr switching_width;    //!< Width of switching region (w_c)
    TPairMatrixPtr epsilon;            //!< Energy depth in kT

    void initPairMatrices() override;
    void extractorsFromJson(const json& j) override;

  public:
    explicit CosAttractMixed(
        const std::string& name = "cos2", const std::string& cite = "doi:10/chqzjk"s,
        CombinationRuleType combination_rule = CombinationRuleType::LORENTZ_BERTHELOT);

    double cutOffSquared(AtomData::index_type id1, AtomData::index_type id2) const; //!< (r_c+w_c)^2

    inline Point force(const Particle& a, const Particle& b, double squared_distance,
                       const Point& b_towards_a) const override
    {
        const auto rc = (*switching_distance)(a.id, b.id);
        const auto wc = (*switching_width)(a.id, b.id);
        const auto cutoff_squared = (rc + wc) * (rc + wc);
        if (squared_distance > cutoff_squared || squared_distance < (rc * rc)) {
            return {0.0, 0.0, 0.0};
        }
        const auto c = pc::pi / (2.0 * wc);
        const auto r = sqrt(squared_distance);
        const auto x1 = std::cos(c * (r - rc));
        const auto x2 = std::sin(c * (r - rc));
        return -2.0 * c * (*epsilon)(a.id, b.id) * x1 * x2 / r * b_towards_a;
    }

    inline double operator()(const Particle& a, const Particle& b, const double squared_distance,
                             [[maybe_unused]] const Point& b_towards_a) const override
    {
        const auto rc = (*switching_distance)(a.id, b.id);
        if (squared_distance < (rc * rc)) {
            return -(*epsilon)(a.id, b.id);
        }
        const auto wc = (*switching_width)(a.id, b.id);
        const auto cutoff_squared = (rc + wc) * (rc + wc);
        if (squared_distance > cutoff_squared) {
            return 0.0;
        }
        const auto c = pc::pi / (2.0 * wc);
        const auto x = std::cos(c * (sqrt(squared_distance) - rc));
        return -(*epsilon)(a.id, b.id) * x * x;
    }
};

/**
 * @brief Pairwise SASA potential calculating the surface area of inter-secting spheres
 */
class SASApotential : public PairPotential
{
  private:
    bool shift = true; // shift potential to zero at large separations?
    double proberadius = 0, conc = 0;
    void from_json(const json& j) override;

    double area(double R, double r, double center_center_distance_squared)
        const; //!< Total surface area of two intersecting spheres or radii R and r as a function of
               //!< separation

  public:
    inline double operator()(const Particle& a, const Particle& b, const double squared_distance,
                             [[maybe_unused]] const Point& b_towards_a) const override
    {
        const auto tfe = 0.5 * (atoms[a.id].tfe + atoms[b.id].tfe);
        const auto tension = 0.5 * (atoms[a.id].tension + atoms[b.id].tension);
        if (fabs(tfe) > 1e-6 or fabs(tension) > 1e-6)
            return (tension + conc * tfe) *
                   area(0.5 * atoms[a.id].sigma, 0.5 * atoms[b.id].sigma, squared_distance);
        return 0.0;
    }

    explicit SASApotential(const std::string& name = "sasa",
                           const std::string& cite = std::string());
    void to_json(json& j) const override;
};

/**
 * @brief Plain Coulomb potential
 */
class Coulomb : public PairPotential
{
  private:
    void from_json(const json& j) override;

  public:
    explicit Coulomb(const std::string& name = "coulomb");
    double bjerrum_length = 0.0; //!< Bjerrum length

    inline double operator()(const Particle& a, const Particle& b, const double squared_distance,
                             const Point&) const override
    {
        return bjerrum_length * a.charge * b.charge / std::sqrt(squared_distance);
    }

    void to_json(json& j) const override;
};

class DipoleDipole : public PairPotential
{
  private:
    void from_json(const json& j) override;

  public:
    explicit DipoleDipole(const std::string& name = "dipoledipole",
                          const std::string& cite = std::string());
    double bjerrum_length{};

    inline double operator()(const Particle& a, const Particle& b, double,
                             const Point& b_towards_a) const override
    {
        return bjerrum_length * mu2mu(a.getExt().mu, b.getExt().mu,
                                      a.getExt().mulen * b.getExt().mulen, b_towards_a, 1.0, 0.0);
    }

    void to_json(json& j) const override;
};

/**
 * @brief Charge-nonpolar pair interaction
 * @note Pair data is _shared_ upon copying
 */
class Polarizability : public Coulomb
{
  private:
    double epsr;
    std::shared_ptr<PairMatrix<double>> m_neutral, m_charged;
    void from_json(const json& j) override;

  public:
    explicit Polarizability(const std::string& name = "polar");
    void to_json(json& j) const override;

    inline double operator()(const Particle& a, const Particle& b, double squared_distance,
                             [[maybe_unused]] const Point& b_towards_a) const override
    {
        double r4inv = 1 / (squared_distance * squared_distance);
        if (fabs(a.charge) > 1e-9 or fabs(b.charge) > 1e-9) {
            return (*m_charged)(a.id, b.id) * r4inv;
        }
        return (*m_neutral)(a.id, b.id) / squared_distance * r4inv;
    }

    inline Point force(const Particle& a, const Particle& b, double squared_distance,
                       const Point& b_towards_a) const override
    {
        double r6inv = 1 / (squared_distance * squared_distance * squared_distance);
        if (fabs(a.charge) > 1e-9 or fabs(b.charge) > 1e-9) {
            return 4 * m_charged->operator()(a.id, b.id) * r6inv * b_towards_a;
        }
        return 6 * m_neutral->operator()(a.id, b.id) / squared_distance * r6inv * b_towards_a;
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
class FENE : public PairPotential
{
    double k = 0;
    double r02 = 0;
    double r02inv = 0;
    void from_json(const json& j) override;

  public:
    explicit FENE(const std::string& name = "fene");
    void to_json(json& j) const override;

    inline double operator()([[maybe_unused]] const Particle& particle1,
                             [[maybe_unused]] const Particle& particle2, double r2,
                             [[maybe_unused]] const Point& b_towards_a) const override
    {
        return (r2 > r02) ? pc::infty : -0.5 * k * r02 * std::log(1 - r2 * r02inv);
    }

    inline Point force([[maybe_unused]] const Particle& particle_a,
                       [[maybe_unused]] const Particle& particle_b, double squared_distance,
                       const Point& b_towards_a) const override
    {
        return (squared_distance > r02) ? -pc::infty * b_towards_a
                                        : -k * r02 / (r02 - squared_distance) * b_towards_a;
    }
};

/**
 * @brief Wrapper for external CoulombGalore library
 */
class NewCoulombGalore : public PairPotential
{
  protected:
    ::CoulombGalore::Splined pot;
    virtual void setSelfEnergy();
    void from_json(const json& j) override;

  public:
    explicit NewCoulombGalore(const std::string& = "coulomb");

    inline double operator()(const Particle& particle_a, const Particle& particle_b,
                             const double squared_distance,
                             [[maybe_unused]] const Point& b_towards_a) const override
    {
        return bjerrum_length *
               pot.ion_ion_energy(particle_a.charge, particle_b.charge,
                                  sqrt(squared_distance) + std::numeric_limits<double>::epsilon());
    }

    inline Point force(const Particle& particle_a, const Particle& particle_b,
                       [[maybe_unused]] double squared_distance,
                       const Point& b_towards_a) const override
    {
        return bjerrum_length *
               pot.ion_ion_force(particle_a.charge, particle_b.charge, b_towards_a); // force on "a"
    }

    void to_json(json& j) const override;
    [[maybe_unused]] double dielectric_constant(double M2V);
    double bjerrum_length;
    const CoulombGalore::Splined& getCoulombGalore() const; //!< Access to full coulomb galore class
};

/**
 * @brief Multipole interactions
 * @todo Only dipole-dipole interactions are currently enabled
 */
class Multipole : public NewCoulombGalore
{
  private:
    void setSelfEnergy() override;

  public:
    explicit Multipole(const std::string& = "multipole");

    inline double operator()(const Particle& particle_a, const Particle& particle_b,
                             [[maybe_unused]] double squared_distance,
                             const Point& b_towards_a) const override
    {
        // Only dipole-dipole for now!
        Point mua = particle_a.getExt().mu * particle_a.getExt().mulen;
        Point mub = particle_b.getExt().mu * particle_b.getExt().mulen;
        return bjerrum_length * pot.dipole_dipole_energy(mua, mub, b_towards_a);
    }

    inline Point force(const Faunus::Particle& particle1, const Faunus::Particle& particle2,
                       [[maybe_unused]] double squared_distance,
                       const Faunus::Point& b_towards_a) const override
    {
        Point mua = particle1.getExt().mu * particle1.getExt().mulen;
        Point mub = particle2.getExt().mu * particle2.getExt().mulen;
        Point ionion = pot.ion_ion_force(particle1.charge, particle2.charge, b_towards_a);
        Point iondip = pot.ion_dipole_force(particle1.charge, mub, b_towards_a) +
                       pot.ion_dipole_force(particle2.charge, mua, b_towards_a);
        Point dipdip = pot.dipole_dipole_force(mua, mub, b_towards_a);
        return bjerrum_length * (ionion + iondip + dipdip);
    }
};

/**
 * @brief Custom pair-potential taking math. expressions at runtime
 * @note `symbols` is a shared_ptr as this allows it to be modified by `operator() const`.
 *       A hack, but would otherwise require const-removal in all pair-potentials.
 */
class CustomPairPotential : public PairPotential
{
  private:
    // Only ExprFunction<double> is explicitly instantiated in functionparser.cpp. Other types as
    // well as the implicit template instantiation is disabled to save reasources during the
    // compilation/build.
    ExprFunction<double> expression;

    struct Symbols
    {
        double distance = 0.0; // available as "r"
        double charge1 = 0.0;  // available as "charge1"
        double charge2 = 0.0;  // available as "charge2"
        double sigma1 = 0.0;   // available as "s1"
        double sigma2 = 0.0;   // available as "s2"
    };

    std::shared_ptr<Symbols> symbols;
    double squared_cutoff_distance;
    json original_input;
    void from_json(const json& j) override;

    inline void setSymbols(const Particle& particle1, const Particle& particle2,
                           double squared_distance) const
    {
        symbols->distance = std::sqrt(squared_distance);
        symbols->charge1 = particle1.charge;
        symbols->charge2 = particle2.charge;
        symbols->sigma1 = Faunus::atoms[particle1.id].sigma;
        symbols->sigma2 = Faunus::atoms[particle2.id].sigma;
    }

  public:
    inline double operator()(const Particle& particle1, const Particle& particle2,
                             double squared_distance,
                             [[maybe_unused]] const Point& b_towards_a) const override
    {
        if (squared_distance < squared_cutoff_distance) {
            setSymbols(particle1, particle2, squared_distance);
            return expression();
        }
        return 0.0;
    }

    inline Point force(const Particle& particle_a, const Particle& particle_b,
                       double squared_distance, const Point& b_towards_b) const override
    {
        if (squared_distance < squared_cutoff_distance) {
            setSymbols(particle_a, particle_b, squared_distance);
            return -expression.derivative(symbols->distance) / symbols->distance * b_towards_b;
        }
        return Point::Zero();
    }

    explicit CustomPairPotential(const std::string& name = "custom");
    void to_json(json& j) const override;
};

using PrimitiveModel =
    CombinedPairPotential<Coulomb, HardSphere>; //!< Primitive model of electrolytes
using PrimitiveModelWCA = CombinedPairPotential<Coulomb, WeeksChandlerAndersen>;
using CigarCosAttractWCA = CompleteCigarPotential<CosAttractMixed, WeeksChandlerAndersen>;

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
class FunctorPotential : public PairPotential
{
    using EnergyFunctor =
        std::function<double(const Particle&, const Particle&, double, const Point&)>;
    json backed_up_json_input; // storage for input json
    bool have_monopole_self_energy = false;
    bool have_dipole_self_energy = false;
    void registerSelfEnergy(PairPotential*); //!< helper func to add to selv_energy_vector
    EnergyFunctor combinePairPotentials(
        json& potential_array); // parse json array of potentials to a single pair-energy functor

  protected:
    PairMatrix<EnergyFunctor, true>
        umatrix; // matrix with potential for each atom pair; cannot be Eigen matrix
    void from_json(const json& j) override;

  public:
    explicit FunctorPotential(const std::string& name = "functor potential");
    void to_json(json& j) const override;

    inline double operator()(const Particle& particle_a, const Particle& particle_b,
                             const double squared_distance,
                             const Point& b_towards_a = {0, 0, 0}) const override
    {
        return umatrix(particle_a.id, particle_b.id)(particle_a, particle_b, squared_distance,
                                                     b_towards_a);
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
 *
 * @todo Add force
 */
class SplinedPotential : public FunctorPotential
{
    /** @brief Expand spline data class to hold information about the sign of values for r<rmin */
    class KnotData : public Tabulate::TabulatorBase<double>::data
    {
      public:
        using base = Tabulate::TabulatorBase<double>::data;
        bool hardsphere_repulsion = false; //!< Use hardsphere repulsion for r smaller than rmin
        KnotData() = default;
        KnotData(const base&);
    };

    PairMatrix<KnotData> matrix_of_knots; //!< Matrix with tabulated potential for each atom pair
    Tabulate::Andrea<double> spline;      //!< Spline method
    bool hardsphere_repulsion = false;    //!< Use hardsphere repulsion for r smaller than rmin
    const int max_iterations = 1e6; //!< Max number of iterations when determining spline interval
    void streamPairPotential(std::ostream& stream, const size_t id1,
                             const size_t id2); //!< Stream pair potential to output stream
    void savePotentials();                      //!< Save splined and exact pair potentials to disk
    double findLowerDistance(int, int, double, double); //!< Find lower distance for splining (rmin)
    double findUpperDistance(int, int, double, double); //!< Find upper distance for splining (rmax)
    double dr = 1e-2; //!< Distance interval when searching for rmin and rmax
    void createKnots(int, int, double,
                     double); //!< Create spline knots for pair of particles in [rmin:rmax]
    void from_json(const json& j) override;

  public:
    explicit SplinedPotential(const std::string& name = "splined");

    /**
     * Policies:
     *
     * 1. return splined potential if rmin>r<rmax
     * 2. return zero if r>=rmax
     * 3. return infinity if r<=rmin AND `hardsphere_repulsion` has been set to true
     * 4. return exact energy if r<=rmin
     */
    inline double operator()(const Particle& particle_a, const Particle& particle_b,
                             double squared_distance,
                             [[maybe_unused]] const Point& b_towards_a) const override
    {
        const auto& knots = matrix_of_knots(particle_a.id, particle_b.id);
        if (squared_distance >= knots.rmax2) {
            return 0.0;
        }
        if (squared_distance > knots.rmin2) {
            return spline.eval(knots, squared_distance); // spline energy
        }
        if (knots.hardsphere_repulsion) {
            return pc::infty;
        }
        return FunctorPotential::operator()(particle_a, particle_b, squared_distance,
                                            {0, 0, 0}); // exact energy
    }
};

} // namespace Faunus::pairpotential