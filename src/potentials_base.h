#pragma once
#include "core.h"
#include "particle.h"
#include "aux/pairmatrix.h"
#include <array>
#include <functional>

/** Namespace for particle pair-potentials */
namespace Faunus::Potential {

//! type of a matrix containing pair potential coefficients
using TPairMatrix = Eigen::MatrixXd;
using TPairMatrixPtr = std::shared_ptr<TPairMatrix>;
//! type of a function extracting a potential coefficient from the InteractionData, e.g., sigma or eps
using TExtractorFunc = std::function<double(const InteractionData&)>;
//! type of a function defining a combination rule of a heterogeneous pair interaction
using TCombinatorFunc = std::function<double(double, double)>;
//! type of a function modifying combinator's output
using TModifierFunc = std::function<double(double)>;

/**
 * @brief Data for a custom (heterogeneous) interaction between two given atom types.
 *
 * The very same format is used as for a homogeneous interaction specified directly on an atom type.
 */
struct CustomInteractionData {
    std::array<AtomData::index_type, 2> atom_id;
    InteractionData interaction;
};

void from_json(const json& j, CustomInteractionData& interaction);
void to_json(json& j, const CustomInteractionData& interaction);
void from_json(const json& j, std::vector<CustomInteractionData>& interactions);

/**
 * @brief Known named combination rules for parameters of pair potential interaction.
 *
 * When adding a new one, add a json mapping. Also consider appending the PairMixer::getCombinator()
 * method to recognize the new rule.
 */
enum class CombinationRuleType { UNDEFINED, ARITHMETIC, GEOMETRIC, LORENTZ_BERTHELOT };
NLOHMANN_JSON_SERIALIZE_ENUM(CombinationRuleType, {{CombinationRuleType::UNDEFINED, "undefined"},
                                                   {CombinationRuleType::ARITHMETIC, "arithmetic"},
                                                   {CombinationRuleType::GEOMETRIC, "geometric"},
                                                   {CombinationRuleType::LORENTZ_BERTHELOT, "lorentz_berthelot"},
                                                   {CombinationRuleType::LORENTZ_BERTHELOT, "LB"}})

/**
 * @brief Exception for handling pair potential initialization.
 */
struct PairPotentialException : public std::runtime_error {
    PairPotentialException(const std::string msg);
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
    TExtractorFunc extractor;   //!< Function extracting the coefficient from the InteractionData structure
    TCombinatorFunc combinator; //!< Function combining two values
    TModifierFunc modifier;     //!< Function modifying the result for fast computations, e.g., a square of

  public:
    PairMixer(TExtractorFunc extractor, TCombinatorFunc combinator, TModifierFunc modifier = &modIdentity);

    //! @return a square matrix of atoms.size()
    TPairMatrixPtr createPairMatrix(const std::vector<AtomData>& atoms);
    //! @return a square matrix of atoms.size()
    TPairMatrixPtr createPairMatrix(const std::vector<AtomData>& atoms,
                                    const std::vector<CustomInteractionData>& interactions);

    enum class CoefficientType { ANY, SIGMA, EPSILON };
    static TCombinatorFunc getCombinator(CombinationRuleType combination_rule,
                                         CoefficientType coefficient = CoefficientType::ANY);

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
 *
 * The `selfEnergy` functor is by default a nullptr. Coulombic potentials
 * may change this to add a self-energy stemming from charges and dipoles.
 * This term is important for for example Widom insertion and grand canonical
 * Monte Carlo schemes.
 */
struct PairPotentialBase {
    std::string name;      //!< unique name per polymorphic call; used in FunctorPotential::combineFunc
    std::string cite;      //!< Typically a short-doi litterature reference
    bool isotropic = true; //!< true if pair-potential is independent of particle orientation
    std::function<double(const Particle&)> selfEnergy = nullptr; //!< self energy of particle (kT)
    virtual void to_json(json&) const = 0;
    virtual void from_json(const json&) = 0;
    virtual ~PairPotentialBase() = default;
    virtual Point force(const Particle& a, const Particle& b, double squared_distance, const Point& b_towards_a) const;
    virtual double operator()(const Particle&, const Particle&, double, const Point&) const = 0;

  protected:
    PairPotentialBase(const std::string& name = std::string(), const std::string& cite = std::string(),
                      bool isotropic = true);
};

void to_json(json& j, const PairPotentialBase& base);   //!< Serialize any pair potential to json
void from_json(const json& j, PairPotentialBase& base); //!< Serialize any pair potential from json

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
    std::shared_ptr<std::vector<CustomInteractionData>> custom_pairs =
        std::make_shared<std::vector<CustomInteractionData>>();
    json json_extra_params;              //!< pickled extra parameters like a coefficient names mapping
    void init();                         //!< initialize the potential when data, e.g., atom parameters, are available
    virtual void initPairMatrices() = 0; //!< potential-specific initialization of parameter matrices
    virtual void extractorsFromJson(const json&); //!< potential-specific assignment of coefficient extracting functions
  public:
    MixerPairPotentialBase(const std::string& name = std::string(), const std::string& cite = std::string(),
                           CombinationRuleType combination_rule = CombinationRuleType::UNDEFINED,
                           bool isotropic = true);
    virtual ~MixerPairPotentialBase() = default;
    void from_json(const json&) override;
    void to_json(json&) const override;
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
    CombinedPairPotential(const std::string& name = "")
        : PairPotentialBase(name){};
    inline double operator()(const Particle& a, const Particle& b, double r2,
                             const Point& r = {0, 0, 0}) const override {
        return first(a, b, r2, r) + second(a, b, r2, r);
    } //!< Combine pair energy

    /**
     * @brief Calculates force on particle a due to another particle, b
     * @param particle_a Particle a
     * @param particle_b Particle b
     * @param squared_distance Squared norm |𝐚-𝐛|²
     * @param b_towards_a Distance vector 𝐛 -> 𝐚 = 𝐚 - 𝐛
     * @return Force on particle a due to particle b
     */
    inline Point force(const Particle& a, const Particle& b, double squared_distance,
                       const Point& b_towards_a) const override {
        return first.force(a, b, squared_distance, b_towards_a) + second.force(a, b, squared_distance, b_towards_a);
    } //!< Combine force

    void from_json(const json& j) override {
        Faunus::Potential::from_json(j, first);
        Faunus::Potential::from_json(j, second);
        name = first.name + "/" + second.name;
        if (first.selfEnergy or second.selfEnergy) { // combine self-energies
            selfEnergy = [u1 = first.selfEnergy, u2 = second.selfEnergy](const Particle& p) {
                if (u1 and u2) {
                    return u1(p) + u2(p);
                }
                if (u1) {
                    return u1(p);
                }
                return u2(p);
            };
        } else {
            selfEnergy = nullptr;
        }
    }
    void to_json(json& j) const override {
        assert(j.is_object());
        auto& _j = j["default"] = json::array();
        _j.push_back(first);
        _j.push_back(second);
    }
};

template <class T1, class T2, class = typename std::enable_if<std::is_base_of<PairPotentialBase, T1>::value>::type,
          class = typename std::enable_if<std::is_base_of<PairPotentialBase, T2>::value>::type>
CombinedPairPotential<T1, T2>& operator+(const T1& pot1, const T2&) {
    return *(new CombinedPairPotential<T1, T2>(pot1.name));
} //!< Statically add two pair potentials at compile-time

} // namespace Faunus::Potential