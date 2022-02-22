#pragma once

#include "move.h"

namespace Faunus {
namespace Move {

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
    using MoveBase::spc;        //!< Trial space (particles, groups)
    Space *other_spc = nullptr; //!< Old space (particles, groups)
    double bond_energy = 0;     //!< Accumulated bond energy if inserted/deleted molecule
    reaction_iterator reaction; //!< Randomly selected reaction

    class AcceptanceData {
      public:
        Average<double> right, left;
        inline void update(ReactionData::Direction direction, bool accept) {
            if (direction == ReactionData::Direction::RIGHT) {
                right += double(accept);
            } else {
                left += double(accept);
            }
        }
    };
    std::map<reaction_iterator, AcceptanceData> acceptance;
    std::map<int, Average<double>> average_reservoir_size; //!< Average number of implicit molecules

    void _to_json(json &) const override;
    void _from_json(const json &) override;
    void _move(Change &) override;         //!< Perform move
    void _accept(Change &) override;       //!< Called when accepted
    void _reject(Change &) override;       //!< Called when rejected
    bool enoughImplicitMolecules() const;  //!< Check if we have enough implicit matter for reaction
    void atomicSwap(Change &);             //!< Swap atom type
    void deactivateAllReactants(Change &); //!< Delete reactant species
    void activateAllProducts(Change &);    //!< Insert product species
    void updateGroupMassCenters(const Change& change) const; //!< Update affected molecular mass centers
    void swapParticleProperties(Particle& target_particle, int new_atomid) const; //!< Swap particle properties

    Change::GroupChange contractAtomicGroup(Space::GroupType&, Space::GroupType&, int); //!< Contract atomic group
    Change::GroupChange expandAtomicGroup(Space::GroupType&, int);                      //!< Expand atomic group
    Change::GroupChange activateMolecularGroup(Space::GroupType&);                      //!< Activate molecular group
    Change::GroupChange deactivateMolecularGroup(Space::GroupType&);                    //!< Deactivate molecular group

    SpeciationMove(Space& spc, std::string_view name, std::string_view cite);

  public:
    SpeciationMove(Space& spc);
    void setOther(Space& spc);
    double bias(Change& change, double old_energy,
                double new_energy) override; //!< adds extra energy change not captured by the Hamiltonian

}; // End of class SpeciationMove

} // namespace Move
} // namespace Faunus
