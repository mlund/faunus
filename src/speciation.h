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
 * Base class for (de)activating groups
 *
 * @todo Not yet implemented
 */
class GroupDeActivator {
  public:
    using ChangeAndBias = std::pair<Change::GroupChange, double>; //!< Group change and possible bias energy
    using OptionalInt = std::optional<int>;
    virtual ChangeAndBias activate(Group& group, OptionalInt num_particles = std::nullopt) = 0;
    virtual ChangeAndBias deactivate(Group& group, OptionalInt num_particles = std::nullopt) = 0;
    virtual ~GroupDeActivator() = default;
};

class AtomicGroupDeActivator : public GroupDeActivator {};

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
    std::function<void(Group&)> setGroupCoordinates; //!< Sets position and rotation of inserted group
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
 *
 * @todo This class is still overly messy and needs general refactoring
 */
class SpeciationMove : public MoveBase {
  private:
    using reaction_iterator = decltype(Faunus::reactions)::iterator;
    Space* old_spc = nullptr;           //!< Old space (particles, groups)
    double bond_energy = 0;             //!< Accumulated bond energy if inserted/deleted molecule
    reaction_iterator reaction;         //!< Randomly selected reaction
    ValidateReaction validate_reaction; //!< Helper to check if reaction is doable
    std::unique_ptr<GroupDeActivator> molecularGroupDeActivator; //!< (de)activator for molecular groups
    std::unique_ptr<GroupDeActivator> atomicGroupDeActivator;    //!< (de)activator for atomic groups

    std::function<void(Group&)> setActivatedGroupCoordinates; //!< Sets position and rotation of inserted group

    class AcceptanceData {
      public:
        Average<double> right, left;
        void update(ReactionData::Direction direction, bool accept);
    };
    std::map<reaction_iterator, AcceptanceData> acceptance;
    std::map<int, Average<double>> average_reservoir_size; //!< Average number of implicit molecules

    void _to_json(json&) const override;
    void _from_json(const json&) override;
    void _move(Change&) override;                            //!< Perform move
    void _accept(Change&) override;                          //!< Called when accepted
    void _reject(Change&) override;                          //!< Called when rejected
    bool enoughImplicitMolecules() const;                    //!< Check if we have enough implicit matter for reaction
    void atomicSwap(Change&);                                //!< Swap atom type
    void deactivateAllReactants(Change&);                    //!< Delete reactant species
    void activateAllProducts(Change&);                       //!< Insert product species
    void updateGroupMassCenters(const Change& change) const; //!< Update affected molecular mass centers
    double getBondEnergy(const Group& group) const;          //!< Summed bond energy of group

    Change::GroupChange contractAtomicGroup(Space::GroupType&, Space::GroupType&, int); //!< Contract atomic group
    Change::GroupChange expandAtomicGroup(Space::GroupType&, int);                      //!< Expand atomic group
    Change::GroupChange activateMolecularGroup(Space::GroupType&);                      //!< Activate molecular group
    Change::GroupChange deactivateMolecularGroup(Space::GroupType&);                    //!< Deactivate molecular group

    SpeciationMove(Space& spc, Space &old_spc, std::string_view name, std::string_view cite);

  public:
    SpeciationMove(Space& spc, Space &old_spc);
    double bias(Change& change, double old_energy,
                double new_energy) override; //!< adds extra energy change not captured by the Hamiltonian

}; // End of class SpeciationMove

} // namespace Move
} // namespace Faunus
