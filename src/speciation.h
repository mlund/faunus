#pragma once

#include "move.h"

namespace Faunus::Speciation {

/**
 * Helper class to check if a reaction is possible, i.e.
 * that there's sufficient reactant and product capacity
 */
class ReactionValidator {
  private:
    const Space& spc;
    bool canSwapAtoms(const ReactionData& reaction) const;
    bool canReduceImplicitGroups(const ReactionData& reaction) const;
    bool canReduceMolecularGroups(const ReactionData& reaction) const;
    bool canProduceMolecularGroups(const ReactionData& reaction) const;
    bool canReduceAtomicGrups(const ReactionData& reaction) const;
    bool canProduceAtomicGroups(const ReactionData& reaction) const;

  public:
    explicit ReactionValidator(const Space& spc);
    bool isPossible(const ReactionData& reaction) const; //!< Enough reactants and product capacity?
};

/**
 * Helper base class for (de)activating groups in speciation move
 */
class GroupDeActivator {
  public:
    using ChangeAndBias = std::pair<Change::GroupChange, double>; //!< Group change and bias energy
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
    Space& spc;     //!< Trial space
    Space& old_spc; //!< Old (accepted) space
    Random& random;

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
    Random& random;
    const bool apply_bond_bias; //!< Set to true to use internal bond energy as bias
    double getBondEnergy(const Group& group) const;
    virtual void setPositionAndOrientation(Group& group) const; //!< Applied to newly activated groups

  public:
    MolecularGroupDeActivator(Space& spc, Random& random, bool apply_bond_bias);
    ChangeAndBias activate(Group& group, OptionalInt num_particles = std::nullopt) override;
    ChangeAndBias deactivate(Group& group, OptionalInt num_particles = std::nullopt) override;
};

/**
 * Helper class to keep track of acceptance in left or right direction
 */
class ReactionDirectionRatio {
  private:
    using reaction_iterator = decltype(Faunus::reactions)::iterator;
    struct AcceptanceData {
        Average<double> right; //!< Acceptance ratio left -> right
        Average<double> left;  //!< Acceptance ratio right -> left
        void update(ReactionData::Direction direction, bool accept);
    };
    std::map<reaction_iterator, AcceptanceData> acceptance;

  public:
    AcceptanceData& operator[](reaction_iterator iter);
    void to_json(json& j) const;
};

} // end of namespace Faunus::Speciation

namespace Faunus::Move {

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
 *
 * @todo Split atom-swap functionality to separate helper class
 */
class SpeciationMove : public MoveBase {
  private:
    using reaction_iterator = decltype(Faunus::reactions)::iterator;
    reaction_iterator reaction;                                            //!< Randomly selected reaction
    double bias_energy = 0.0;                                              //!< Group (de)activators may add bias
    Speciation::ReactionValidator reaction_validator;                      //!< Helper to check if reaction is doable
    Speciation::ReactionDirectionRatio direction_ratio;                    //!< Track acceptance in each direction
    std::unique_ptr<Speciation::GroupDeActivator> molecular_group_bouncer; //!< (de)activator for molecular groups
    std::unique_ptr<Speciation::GroupDeActivator> atomic_group_bouncer;    //!< (de)activator for atomic groups

    std::map<MoleculeData::index_type, Average<double>>
        average_reservoir_size; //!< Average number of implicit molecules

    void _to_json(json& j) const override;
    void _from_json(const json& j) override;
    void _move(Change& change) override;
    void _accept(Change& change) override;
    void _reject(Change& change) override;

    void setRandomReactionAndDirection();                       //!< Set random reaction and direction
    void atomicSwap(Change& change);          //!< Swap atom type
    void deactivateReactants(Change& change); //!< Delete all reactants
    void activateProducts(Change& change);    //!< Insert all products
    void deactivateAtomicGroups(Change& change);
    void activateAtomicGroups(Change& change);
    void deactivateMolecularGroups(Change& change);
    void activateMolecularGroups(Change& change);
    void updateGroupMassCenters(const Change& change) const; //!< Update affected molecular mass centers
    void swapParticleProperties(Particle& particle, int new_atomid) const;
    SpeciationMove(Space& spc, Space& old_spc, std::string_view name, std::string_view cite);

  public:
    SpeciationMove(Space& spc, Space& old_spc);
    double bias(Change& change, double old_energy, double new_energy) override;
};

} // namespace Faunus::Move
