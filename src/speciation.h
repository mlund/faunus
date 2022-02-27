#pragma once

#include "move.h"

namespace Faunus {
namespace Move {

/**
 * Helper class to check if a reaction is possible, i.e.
 * that there's sufficient reactant and product capacity
 */
class ValidateReaction {
  private:
    const Space& spc;
    bool enoughImplicitMolecules(const ReactionData& reaction) const;
    bool canSwapAtoms(const ReactionData& reaction) const;
    bool canReduceMolecularGrups(const ReactionData& reaction) const;
    bool canProduceMolecularGroups(const ReactionData& reaction) const;
    bool canReduceAtomicGrups(const ReactionData& reaction) const;
    bool canProduceAtomicGroups(const ReactionData& reaction) const;

  public:
    explicit ValidateReaction(const Space& spc);
    bool reactionIsPossible(const ReactionData& reaction) const;
};

/**
 * Helper base class for (de)activating groups in speciation move
 */
class GroupDeActivator {
  public:
    using ChangeAndBias = std::pair<Change::GroupChange, double>; //!< Group change and possible bias energy
    using OptionalInt = std::optional<int>;
    virtual ChangeAndBias activate(Group& group, OptionalInt num_particles = std::nullopt) = 0;
    virtual ChangeAndBias deactivate(Group& group, OptionalInt num_particles = std::nullopt) = 0;
    virtual ~GroupDeActivator() = default;
};

/**
 * Helper class for contracting and expanding atomic groups
 */
class AtomicGroupDeActivator : public GroupDeActivator {
  private:
    Space& spc;
    Space& old_spc;
    Random& slump;

  public:
    AtomicGroupDeActivator(Space& spc, Space& old_spc, Random& random);
    ChangeAndBias activate(Group& group, OptionalInt number_to_insert) override;
    ChangeAndBias deactivate(Group& group, OptionalInt number_to_delete) override;
};

/**
 * Helper class to (de)activate a single molecular group
 *
 * Activation policy:
 * - Set random mass center position and orientation
 * - Apply PBC wrapping
 * - Bias is set to *negative* internal bond energy
 *
 * Deactivation policy:
 * - Remove PBC before deactivation
 * - Bias is set to *positive* internal bond energy
 */
class MolecularGroupDeActivator : public GroupDeActivator {
  private:
    Space& spc;
    std::function<void(Group&)> setPositionAndOrientation; //!< Sets position and rotation of inserted group
    double getBondEnergy(const Group& group) const;

  public:
    MolecularGroupDeActivator(Space& spc, Random& random);
    ChangeAndBias activate(Group& group, OptionalInt num_particles = std::nullopt) override;
    ChangeAndBias deactivate(Group& group, OptionalInt num_particles = std::nullopt) override;
};

/**
 * @brief Generalised Grand Canonical Monte Carlo Move
 *
 * This move handles insertion and deletion of atomic and
 * molecular groups, including swap moves and implicit atomic
 * species. Flow of events:
 *
 * 1. pick random `Faunus::ReactionData` object
 * 2. pick random direction (left or right)
 * 3. perform appropriate action:
 *    - atomic swap
 *    - deactivate reactants
 *    - activate products
 */
class SpeciationMove : public MoveBase {
  private:
    using reaction_iterator = decltype(Faunus::reactions)::iterator;
    reaction_iterator reaction;         //!< Randomly selected reaction
    Space* old_spc = nullptr;           //!< Old space (particles, groups)
    double bias_energy = 0.0;           //!< Accumulated bond energy if inserted/deleted molecule
    ValidateReaction validate_reaction; //!< Helper to check if reaction is doable
    std::unique_ptr<GroupDeActivator> molecularGroupDeActivator; //!< (de)activator for molecular groups
    std::unique_ptr<GroupDeActivator> atomicGroupDeActivator;    //!< (de)activator for atomic groups

    struct AcceptanceData {
        Average<double> right; //!< Original left side of reaction
        Average<double> left;  //!< Original right side of reaction
        void update(ReactionData::Direction direction, bool accept);
    };
    std::map<reaction_iterator, AcceptanceData> acceptance;
    std::map<int, Average<double>> average_reservoir_size; //!< Average number of implicit molecules

    void _to_json(json& j) const override;
    void _from_json(const json&) override;
    void _move(Change& change) override;   //!< Perform move
    void _accept(Change& change) override; //!< Called when accepted
    void _reject(Change& change) override; //!< Called when rejected

    void setReaction();                                      //!< Set random reaction and direction
    bool enoughImplicitMolecules() const;                    //!< enough implicit matter for reaction?
    void atomicSwap(Change& change);                         //!< Swap atom type
    void deactivateAllReactants(Change& change);             //!< Delete all reactants
    void activateAllProducts(Change& change);                //!< Insert all products
    void updateGroupMassCenters(const Change& change) const; //!< Update affected molecular mass centers

    SpeciationMove(Space& spc, Space& old_spc, std::string_view name, std::string_view cite);

  public:
    SpeciationMove(Space& spc, Space& old_spc);
    double bias(Change& change, double old_energy, double new_energy) override;
};

} // namespace Move
} // namespace Faunus
