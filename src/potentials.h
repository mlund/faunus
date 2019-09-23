#pragma once

#include "core.h"
#include "units.h"
#include "particle.h"
#include "auxiliary.h"
#include "tabulate.h"
#include "functionparser.h"
#include "multipole.h"
#include <coulombgalore.h>
#include <array>
#include <functional>

/*
namespace CoulombGalore {
    class SchemeBase;
    class Splined;
}*/

namespace Faunus {
namespace Potential {

using namespace std::string_literals;

//! type of a matrix containing pair potential coefficients
typedef Eigen::MatrixXd TPairMatrix;
typedef std::shared_ptr<TPairMatrix> TPairMatrixPtr;
//! type of a function extracting a potential coefficient from the AtomData, e.g., sigma or eps
typedef std::function<double(const AtomData &)> TExtractorFunc;
//! type of a function defining a combination rule of a heterogeneous pair interaction
typedef std::function<double(double, double)> TCombinatorFunc;
//! type of a function modifying combinator's output
typedef std::function<double(double)> TModifierFunc;

/**
 * @brief Data for a custom (heterogeneous) interaction between two given atom types.
 *
 * The very same format is used as for a homogeneous interaction specified directly on an atom type.
 */
struct InteractionData {
    std::array<AtomData::Tid, 2> atom_id;
    AtomData interaction;
};

void from_json(const json &j, std::vector<InteractionData> &interactions);
void to_json(json &j, const std::vector<InteractionData> &interactions);

/**
 * @brief Known named combination rules for parameters of pair potential interaction.
 *
 * When adding a new one, add a json mapping. Also consider appending the PairMixer::getCombinator()
 * method to recognize the new rule.
 */
enum CombinationRuleType { COMB_UNDEFINED, COMB_ARITHMETIC, COMB_GEOMETRIC, COMB_LORENTZ_BERTHELOT };
NLOHMANN_JSON_SERIALIZE_ENUM(CombinationRuleType, {
                                                      {COMB_UNDEFINED, "undefined"},
                                                      {COMB_ARITHMETIC, "arithmetic"},
                                                      {COMB_GEOMETRIC, "geometric"},
                                                      {COMB_LORENTZ_BERTHELOT, "lorentz_berthelot"},
                                                      {COMB_LORENTZ_BERTHELOT, "LB"}, // alternative non-canonical name
                                                  })

/**
 * @brief Exception for handling pair potential initialization.
 */
struct PairPotentialException : public std::runtime_error {
    PairPotentialException(const std::string msg) : std::runtime_error(msg){};
};

/**
 * @brief PairMixer creates a matrix of pair potential coefficients based on the atom properties
 * and/or custom values using an arbitrary combination rule.
 *
 * PairMixer holds three functions that are applied in order extractor → combinator → modifier to create
 * a coefficient matrix for all possible interactions. The function createPairMatrix applies the functions
 * on all atom type pairs, and optionally also on the list of custom pair parameters (not the combinator
 * function).
 */
class PairMixer {
    TExtractorFunc extractor;   //!< Function extracting the coefficient from the AtomData structure
    TCombinatorFunc combinator; //!< Function combining two values
    TModifierFunc modifier;     //!< Function modifying the result for fast computations, e.g., a square of

  public:
    PairMixer(TExtractorFunc extractor, TCombinatorFunc combinator, TModifierFunc modifier = &modIdentity)
        : extractor(extractor), combinator(combinator), modifier(modifier){};

    //! @return a square matrix of atoms.size()
    TPairMatrixPtr createPairMatrix(const std::vector<AtomData> &atoms);
    //! @return a square matrix of atoms.size()
    TPairMatrixPtr createPairMatrix(const std::vector<AtomData> &atoms,
                                    const std::vector<InteractionData> &interactions);

    enum CoefficientType {COEF_ANY, COEF_SIGMA, COEF_EPSILON};
    static TCombinatorFunc getCombinator(CombinationRuleType combination_rule, CoefficientType coefficient = COEF_ANY);

    // when explicit custom pairs are the only option
    inline static constexpr double combUndefined(double = 0.0, double = 0.0) {
        return std::numeric_limits<double>::signaling_NaN();
    };
    inline static double combArithmetic(double a, double b) { return 0.5 * (a + b); }
    inline static double combGeometric(double a, double b) { return std::sqrt(a * b); }
    inline static double modIdentity(double x) { return x; }
    inline static double modSquared(double x) { return x * x; }
};

/**
 * @brief Base for all pair-potentials
 */
struct PairPotentialBase {
    std::string name; //!< unique name per polymorphic call; used in FunctorPotential::combineFunc
    std::string cite;
    bool isotropic; //!< true if pair-potential is independent of particle orientation
    std::function<double(const Particle &)> selfEnergy; //!< self energy of particle (kT)
    virtual void to_json(json &) const = 0;
    virtual void from_json(const json &) = 0;
    virtual ~PairPotentialBase() = default;
    virtual Point force(const Particle &, const Particle &, double, const Point &) const;
    virtual double operator()(const Particle &a, const Particle &b, const Point &r) const = 0;

  protected:
    PairPotentialBase(const std::string &name = std::string(), const std::string &cite = std::string(),
                      bool isotropic = true);
};

void to_json(json &j, const PairPotentialBase &base);   //!< Serialize any pair potential to json
void from_json(const json &j, PairPotentialBase &base); //!< Serialize any pair potential from json

/**
 * @brief A common ancestor for potentials that use parameter matrices computed from atomic
 * properties and/or custom atom pair properties.
 *
 * The class and their descendants have now also a responsibility to create themselves from a json object
 * and store back. This is gradually becoming a complex task which shall be moved into other class.
 */
class MixerPairPotentialBase : public PairPotentialBase {
  protected:
    CombinationRuleType combination_rule;
    std::shared_ptr<std::vector<InteractionData>> custom_pairs = std::make_shared<std::vector<InteractionData>>();
    json json_extra_params;              //!< pickled extra parameters like a coefficient names mapping
    void init();                         //!< initialize the potential when data, e.g., atom parameters, are available
    virtual void initPairMatrices() = 0; //!< potential-specific initialization of parameter matrices
    virtual void extractorsFromJson(const json &) {}; //!< potential-specific assignment of coefficient extracting functions
  public:
    MixerPairPotentialBase(const std::string &name = std::string(), const std::string &cite = std::string(),
                           CombinationRuleType combination_rule = COMB_UNDEFINED, bool isotropic = true)
        : PairPotentialBase(name, cite, isotropic), combination_rule(combination_rule) {};
    virtual ~MixerPairPotentialBase() = default;
    void from_json(const json &) override;
    void to_json(json &) const override;
};

/**
 * @brief Statically combines two pair potentials at compile-time
 *
 * This is the most efficient way to combining pair-potentials due
 * to the possibility for compile-time optimisation.
 */
template <class T1, class T2> struct CombinedPairPotential : public PairPotentialBase {
    T1 first;  //!< First pair potential of type T1
    T2 second; //!< Second pair potential of type T2
    CombinedPairPotential(const std::string &name = "") : PairPotentialBase(name) {};
    inline double operator()(const Particle &a, const Particle &b, const Point &r) const override {
        return first(a, b, r) + second(a, b, r);
    } //!< Combine pair energy

    inline Point force(const Particle &a, const Particle &b, double r2, const Point &p) const override {
        return first.force(a, b, r2, p) + second.force(a, b, r2, p);
    } //!< Combine force

    void from_json(const json &j) override {
        first = j;
        second = j;
        // combine self-energies
        if (first.selfEnergy or second.selfEnergy) {
            selfEnergy = [u1 = first.selfEnergy, u2 = second.selfEnergy](const Particle &p) {
                if (u1 and u2)
                    return u1(p) + u2(p);
                if (u1)
                    return u1(p);
                return u2(p);
            };
        } else
            selfEnergy = nullptr;
    }

    void to_json(json &j) const override { j = {first, second}; }
};

template <class T1, class T2, class = typename std::enable_if<std::is_base_of<PairPotentialBase, T1>::value>::type,
          class = typename std::enable_if<std::is_base_of<PairPotentialBase, T2>::value>::type>
CombinedPairPotential<T1, T2> &operator+(const T1 &pot1, const T2 &) {
    return *(new CombinedPairPotential<T1, T2>(pot1.name));
} //!< Statically add two pair potentials at compile-time


struct Dummy : public PairPotentialBase {
    Dummy();
    inline double operator()(const Particle &, const Particle &, const Point &) const override { return 0.0; }
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
    LennardJones(const std::string &name = "lennardjones", const std::string &cite = std::string(),
                 CombinationRuleType combination_rule = COMB_LORENTZ_BERTHELOT)
        : MixerPairPotentialBase(name, cite, combination_rule){};

    inline Point force(const Particle &a, const Particle &b, double r2, const Point &p) const override {
        double s6 = powi((*sigma_squared)(a.id, b.id), 3);
        double r6 = r2 * r2 * r2;
        double r14 = r6 * r6 * r2;
        return 6. * (*epsilon_quadruple)(a.id, b.id) * s6 * (2 * s6 - r6) / r14 * p;
    }

    inline double operator()(const Particle &a, const Particle &b, const Point &r) const override {
        double x = (*sigma_squared)(a.id, b.id) / r.squaredNorm(); // s2/r2
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
    WeeksChandlerAndersen(const std::string &name = "wca", const std::string &cite = "doi:ct4kh9",
                          CombinationRuleType combination_rule = COMB_LORENTZ_BERTHELOT)
        : LennardJones(name, cite, combination_rule) {};

    inline double operator()(const Particle &a, const Particle &b, const Point &r) const override {
        return operator()(a, b, r.squaredNorm());
    }

    inline Point force(const Particle &a, const Particle &b, double r2, const Point &p) const override {
        double x = (*sigma_squared)(a.id, b.id); // s^2
        if (r2 > x * twototwosixth)
            return Point(0, 0, 0);
        x = x / r2;    // (s/r)^2
        x = x * x * x; // (s/r)^6
        return (*epsilon_quadruple)(a.id, b.id) * 6 * (2 * x * x - x) / r2 * p;
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
    HardSphere(const std::string &name = "hardsphere")
        : MixerPairPotentialBase(name, std::string(), COMB_ARITHMETIC) {};

    inline double operator()(const Particle &a, const Particle &b, const Point &r) const override {
        return r.squaredNorm() < (*sigma_squared)(a.id, b.id) ? pc::infty : 0.0;
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
    Hertz(const std::string &name = "hertz")
        : MixerPairPotentialBase(name) {};
    inline double operator()(const Particle &a, const Particle &b, const Point &r) const override {
        double r2 = r.squaredNorm();
        if (r2 <= (*sigma_squared)(a.id, b.id))
            return (*epsilon)(a.id, b.id) * pow((1 - (sqrt(r2 / (*sigma_squared)(a.id, b.id)))), 2.5);
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
    SquareWell(const std::string &name = "squarewell")
        : MixerPairPotentialBase(name) {};
    inline double operator()(const Particle &a, const Particle &b, const Point &r) const override {
        return (r.squaredNorm() < (*sigma_squared)(a.id, b.id)) ? -(*epsilon)(a.id, b.id) : 0.0;
    }
};

struct RepulsionR3 : public PairPotentialBase {
    double f = 0, s = 0, e = 0;

    RepulsionR3(const std::string &name = "repulsionr3") : PairPotentialBase(name) {};
    void from_json(const json &j) override;
    void to_json(json &j) const override;

    inline double operator()(const Particle &, const Particle &, const Point &_r) const override {
        double r2 = _r.squaredNorm(), r = sqrt(r2);
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
    CosAttract(const std::string &name = "cos2") : PairPotentialBase(name) {};

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
    inline double operator()(const Particle &, const Particle &, const Point &r) const override {
        double r2 = r.squaredNorm();
        if (r2 < rc2)
            return -eps;
        if (r2 > rcwc2)
            return 0;
        double x = std::cos(c * (sqrt(r2) - rc));
        return -eps * x * x;
    }

    inline Point force(const Particle &, const Particle &, double r2, const Point &p) const override {
        if (r2 < rc2 || r2 > rcwc2)
            return Point(0, 0, 0);
        double r = sqrt(r2);
        double x1 = std::cos(c * (r - rc));
        double x2 = std::sin(c * (r - rc));
        return -2 * c * eps * x1 * x2 / r * p;
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
    inline double operator()(const Particle &a, const Particle &b, const Point &r_ab) const override {
        double tfe = 0.5 * (atoms[a.id].tfe + atoms[b.id].tfe);
        double tension = 0.5 * (atoms[a.id].tension + atoms[b.id].tension);
        if (fabs(tfe) > 1e-6 or fabs(tension) > 1e-6)
            return (tension + conc * tfe) * area(0.5 * atoms[a.id].sigma, 0.5 * atoms[b.id].sigma, r_ab.squaredNorm());
        return 0;
    }
    SASApotential(const std::string &name = "sasa", const std::string &cite = std::string()) :
        PairPotentialBase(name, cite) {};
    void to_json(json &j) const override;
    void from_json(const json &j) override;
};

struct Coulomb : public PairPotentialBase {
    Coulomb(const std::string &name = "coulomb") : PairPotentialBase(name) {};
    double lB; //!< Bjerrum length
    inline double operator()(const Particle &a, const Particle &b, const Point &r) const override {
        return lB * a.charge * b.charge / r.norm();
    }
    void to_json(json &j) const override;
    void from_json(const json &j) override;
};

struct DipoleDipole : public PairPotentialBase {
    DipoleDipole(const std::string &name = "dipoledipole", const std::string &cite = std::string()) :
        PairPotentialBase(name, cite, false) {};
    double lB; //!< Bjerrum length
    inline double operator()(const Particle &a, const Particle &b, const Point &r) const override {
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
    Polarizability(const std::string &name = "polar") : Coulomb(name) {
        m_neutral = std::make_shared<PairMatrix<double>>();
        m_charged = std::make_shared<PairMatrix<double>>();
    };

    void from_json(const json &j) override;

    void to_json(json &j) const override { j = {{"epsr", epsr}}; }

    inline double operator()(const Particle &a, const Particle &b, const Point &r) const override {
        double r2 = r.squaredNorm();
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
    FENE(const std::string &name = "fene") : PairPotentialBase(name) {};
    void from_json(const json &j) override;
    void to_json(json &j) const override;

    /** @brief Energy in kT between two particles, r2 = squared distance */
    inline double operator()(const Particle &, const Particle &, const Point &r) const override {
        double r2 = r.squaredNorm();
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
    double lB;

  public:
    NewCoulombGalore(const std::string & = "coulomb");
    inline double operator()(const Particle &a, const Particle &b, const Point &r) const override {
        return lB * pot.ion_ion_energy(a.charge, b.charge, r.squaredNorm());
    }
    Point force(const Particle &, const Particle &, double, const Point &) const override;
    void from_json(const json &) override;
    void to_json(json &) const override;
    template <class Tpvec, class Tgroup> double internal(const Tgroup &g) const {
        double Eq = 0;
        for (auto i : g)
            Eq += i.charge * i.charge;
        return pot.self_energy({Eq, 0.0}) * lB;
    }
    double dielectric_constant(double M2V) { return pot.calc_dielectric(M2V); }
};

/**
 * @brief Multipole interactions
 * @note Only dipole-dipole interactions are currently enabled
 */
class Multipole : public NewCoulombGalore {
  public:
    Multipole(const std::string & = "coulomb");
    inline double operator()(const Particle &a, const Particle &b, const Point &r) const override {
        // for this to reach reasonable speed, the energy calculations should
        // be aggregated to avoid SRF lookup multiple times
        Point mua = a.getExt().mu * a.getExt().mulen;
        Point mub = b.getExt().mu * b.getExt().mulen;
        // double ionion = pot.ion_ion_energy(a.charge, b.charge, r.squaredNorm());
        // double iondip = pot.ion_dipole_energy(a.charge, mub, r) + pot.ion_dipole_energy(b.charge, mua, r);
        double dipdip = pot.dipole_dipole_energy(mua, mub, r);
        return lB * (dipdip);
        // return lB*(ionion + iondip + dipdip);
    }

    Point force(const Particle &, const Particle &, double, const Point &) const override;
};

/** @brief Coulomb type potentials with spherical cutoff */
class CoulombGalore : public PairPotentialBase {
    std::shared_ptr<PairMatrix<double>> ecs;      // effective charge-scaling
    Tabulate::Andrea<double> sf;                  // splitting function
    Tabulate::TabulatorBase<double>::data table;  // data for splitting function
    std::function<double(double)> calcDielectric; // function for dielectric const. calc.
    std::string type;
    double selfenergy_prefactor = 0;
    double lB, depsdt, rc, rc2, rc1i, epsr, epsrf, alpha, kappa, I;
    int order;
    unsigned int C, D;

    void sfYukawa(const json &j);
    void sfReactionField(const json &j);
    void sfQpotential(const json &j);
    void sfYonezawa(const json &j);
    void sfFanourgakis(const json &j);
    void sfYukawaPoisson(const json &j);
    void sfPoisson(const json &j);
    void sfFennel(const json &j);
    void sfEwald(const json &j);
    void sfWolf(const json &j);
    void sfPlain(const json &j, double val = 1);

    inline double operator()(const Particle &a, const Particle &b, double r2) const {
        if (r2 < rc2) {
            double r = std::sqrt(r2);
            return lB * ecs->operator()(a.id, b.id) * a.charge * b.charge / r * sf.eval(table, r * rc1i);
        }
        return 0;
    }

  public:
    CoulombGalore(const std::string &name = "coulomb") : PairPotentialBase(name) {};
    void from_json(const json &j) override;

    inline double operator()(const Particle &a, const Particle &b, const Point &r) const override {
        return operator()(a, b, r.squaredNorm());
    }

    inline Point force(const Particle &a, const Particle &b, double r2, const Point &p) const override {
        if (r2 < rc2) {
            double r = sqrt(r2);
            return lB * a.charge * b.charge * (-sf.eval(table, r * rc1i) / r2 + sf.evalDer(table, r * rc1i) / r) * p;
        }
        return Point(0, 0, 0);
    }

    /**
     * @brief Self-energy of the potential
     */
    template <class Tpvec, class Tgroup> double internal(const Tgroup &g) const {
        double Eq = 0;
        for (auto i : g)
            Eq += i.charge * i.charge;
        return selfenergy_prefactor * Eq * lB / rc;
    }

    double dielectric_constant(double M2V);

    void to_json(json &j) const override;
};

/** @brief Dipole-dipole type potentials with spherical cutoff */
class DipoleDipoleGalore : public PairPotentialBase {
    Tabulate::Andrea<double> sfA, sfB;            // splitting functions
    Tabulate::TabulatorBase<double>::data tableA, tableB;  // data for splitting function
    std::function<double(double)> calcDielectric; // function for dielectric const. calc.
    std::string type;
    double selfenergy_prefactor;
    double lB, depsdt, rc, rc2, rc1i, epsr, epsrf, alpha, kappa;
    int order;
    // unsigned int C, D;

    void sfEwald(const json &j);
    void sfReactionField(const json &j);
    void sfQ0potential(const json &j);
    void sfQ2potential(const json &j);
    void sfFanourgakis(const json &j);
    void sfFennell(const json &j);
    void sfWolf(const json &j);
    void sfPlain(const json &j, double val = 1);

  public:
    DipoleDipoleGalore(const std::string &name = "dipoledipole", const std::string &cite = std::string()) :
            PairPotentialBase(name, cite, false) {};
    void from_json(const json &j) override;

    inline double operator()(const Particle &a, const Particle &b, const Point &r) const override {
        double r1 = r.norm();
        if (r1 < rc) {
            double af = sfA.eval(tableA,r1*rc1i);
            double bf = sfB.eval(tableB,r1*rc1i);
            return lB*mu2mu(a.getExt().mu, b.getExt().mu, a.getExt().mulen*b.getExt().mulen, r,af,bf);
        }
        return 0.0;
    }

    inline Point force(const Particle &, const Particle &, double, const Point &) const override { return {0, 0, 0}; }

    /**
     * @brief Self-energy of the potential
     */
    template <class Tpvec, class Tgroup> double internal(const Tgroup &g) const {
        double Emu = 0;
        for (auto i : g)
            Emu += i.getExt().mulen * i.getExt().mulen;
        return selfenergy_prefactor * Emu * lB / pow(rc,3.0);
    }

    double dielectric_constant(double M2V);

    void to_json(json &j) const override;
};

/**
 * @brief Custom pair-potential taking math. expressions at runtime
 */
class CustomPairPotential : public PairPotentialBase {
  private:
    // Only ExprFunction<double> is explicitly instantiated in functionparser.cpp. Other types as well as
    // the implicit template instantiation is disabled to save reasources during the compilation/build.
    ExprFunction<double> expr;
    struct Data {
        double r = 0, q1 = 0, q2 = 0, s1 = 0, s2 = 0;
    };
    double Rc2;
    std::shared_ptr<Data> d;
    json jin; // initial json input
  public:
    inline double operator()(const Particle &a, const Particle &b, const Point &r) const override {
        double r2 = r.squaredNorm();
        if (r2 > Rc2)
            return 0;
        d->r = sqrt(r2);
        d->q1 = a.charge;
        d->q2 = b.charge;
        d->s1 = atoms[a.id].sigma;
        d->s2 = atoms[b.id].sigma;
        return expr();
    }
    CustomPairPotential(const std::string &name = "custom") : PairPotentialBase(name), d(std::make_shared<Data>()) {};

    void from_json(const json &) override;
    void to_json(json &) const override;
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
    typedef std::function<double(const Particle &, const Particle &, const Point &)> uFunc;
    json _j; // storage for input json
    typedef CombinedPairPotential<Coulomb, HardSphere> PrimitiveModel;
    typedef CombinedPairPotential<Coulomb, WeeksChandlerAndersen> PrimitiveModelWCA;
    typedef CombinedPairPotential<DipoleDipole, LennardJones> Stockmayer;

    std::vector<std::function<double(const Particle &)>> self_energy_vector;
    bool have_monopole_self_energy = false;
    bool have_dipole_self_energy = false;
    void registerSelfEnergy(PairPotentialBase *); //!< helper func to add to selv_energy_vector

    // List of pair-potential instances used when constructing functors.
    // Note that potentials w. large memory requirements (LJ, WCA etc.)
    // typically use `shared_ptr` so that the created functors _share_
    // the data. That is *only* put the pair-potential here if you can
    // share internal (shared) pointers.
    std::tuple<CoulombGalore,         // 0
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
               DipoleDipoleGalore,    // 12
               Stockmayer,            // 13
               NewCoulombGalore,      // 14
               Multipole              // 15
               >
        potlist;

    uFunc combineFunc(json &j); // parse json array of potentials to a single potential function object

  protected:
    PairMatrix<uFunc, true> umatrix; // matrix with potential for each atom pair; cannot be Eigen matrix

  public:
    FunctorPotential(const std::string &name = "functor") : PairPotentialBase(name) {};
    void to_json(json &j) const override;
    void from_json(const json &j) override;

    inline double operator()(const Particle &a, const Particle &b, const Point &r) const override {
        return umatrix(a.id, b.id)(a, b, r); // pc::infty;
    }
};

/**
 * @brief Tabulated arbitrary potentials for specific atom types
 *
 * This maintains a species x species matrix as in FunctorPotential
 * but with tabulated pair potentials to improve performance.
 *
 */
class TabulatedPotential : public FunctorPotential {

    // expand spline data class to hold information about
    // the sign of values for r<rmin
    struct Ttable : public Tabulate::TabulatorBase<double>::data {
        typedef Tabulate::TabulatorBase<double>::data base;
        bool isNegativeBelowRmin = false;
        Ttable(){};
        Ttable(const base &b) : base(b) {}
    };
    PairMatrix<Ttable, true> tmatrix; // matrix with tabulated potential for each atom pair; cannot be Eigen Matrix
    Tabulate::Andrea<double> tblt;    // spline class
    bool hardsphere = false;          // use hardsphere for r<rmin?

  public:
    TabulatedPotential(const std::string &name = "splined") : FunctorPotential(name) {};

    inline double operator()(const Particle &a, const Particle &b, const Point &r) const override {
        double r2 = r.squaredNorm();
        const Ttable &knots = tmatrix(a.id, b.id);
        if (r2 >= knots.rmax2)
            return 0.0;
        else if (r2 <= knots.rmin2) {
            if (knots.isNegativeBelowRmin or (not hardsphere))
                return this->umatrix(a.id, b.id)(a, b, r); // exact energy
            else
                return pc::infty; // assume extreme repulsion
        }
        return tblt.eval(knots, r2); // we are in splined interval
    }

    void from_json(const json &j) override;
};

} // end of namespace Potential
} // end of namespace Faunus
