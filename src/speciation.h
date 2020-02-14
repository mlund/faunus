#pragma once

#include "move.h"

namespace Faunus {
namespace Move {

/*
 * @brief Establishes equilibrium of matter
 * Establishes equilibrium of matter between all species
 *
 * Consider the dissociation process AX=A+X. This class will locate
 * all species of type AX and A and make a MC swap move between them.
 * X can be implicit, meaning that it enters only with its chemical potential
 * (activity). The reacting species, the equilibrium constant,
 * and the activities are read from the JSON input file.
 *
 * @todo Messy code requires substantial cleanup and documentation
 */
class SpeciationMove : public Movebase {
  private:
    Space &spc;
    Space *otherspc;
    ReactionData *trialprocess; // todo: eliminate raw pointer
    std::map<std::string, Average<double>> acceptance_map;

    double lnK;
    double bond_energy = 0;
    bool forward = true;  // reaction moving forwards or backwards? (or left or right)
    bool neutral = false; // true if only neutral molecules are involved in the reaction

    void _to_json(json &) const override;
    void _from_json(const json &) override;

    bool checkInsertProducts(ReactionData &);
    bool swapReaction(Change &, ReactionData &);

    Change::data contractAtomicGroup(Space::Tgroup &, Tspace::Tgroup &, int); //!< Contract atomic group
    Change::data expandAtomicGroup(Space::Tgroup &, int);                     //!< Expand atomic group
    Change::data activateMolecularGroup(Space::Tgroup &);                     //!< Activate molecular group
    Change::data deactivateMolecularGroup(Space::Tgroup &);                   //!< Deactivate molecular group

    void activateProducts(Change &, std::vector<Faunus::ReactionData>::iterator);
    void deactivateReactants(Change &, std::vector<ReactionData>::iterator);

  public:
    SpeciationMove(Space &);
    void setOther(Space &);
    void _move(Change &) override;
    double bias(Change &, double, double) override; //!< adds extra energy change not captured by the Hamiltonian
    void _accept(Change &) override;
    void _reject(Change &) override;

}; // End of class SpeciationMove

} // namespace Move
} // namespace Faunus
