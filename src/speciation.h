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
 */
class SpeciationMove : public Movebase {
  private:
    typedef decltype(Faunus::reactions)::iterator reaction_iterator;
    Space &spc;                 //!< Trial space (particles, groups)
    Space *other_spc;           //!< Old space (particles, groups)
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
    bool checkBeforeInsert();              //!< Performs checks before move
    bool atomicSwap(Change &);             //!< Swap atom type
    bool deactivateAllReactants(Change &); //!< Delete reactant species
    bool activateAllProducts(Change &);    //!< Insert product species

    Change::data contractAtomicGroup(Space::Tgroup &, Space::Tgroup &, int);  //!< Contract atomic group
    Change::data expandAtomicGroup(Space::Tgroup &, int);                     //!< Expand atomic group
    Change::data activateMolecularGroup(Space::Tgroup &);                     //!< Activate molecular group
    Change::data deactivateMolecularGroup(Space::Tgroup &);                   //!< Deactivate molecular group

  public:
    SpeciationMove(Space &);
    void setOther(Space &);
    double bias(Change &, double, double) override; //!< adds extra energy change not captured by the Hamiltonian

}; // End of class SpeciationMove

} // namespace Move
} // namespace Faunus
