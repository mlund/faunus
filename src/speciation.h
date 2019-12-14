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
 */
class SpeciationMove : public Movebase {
  private:
    Space &spc;
    Space *otherspc;
    ReactionData *trialprocess;
    std::map<std::string, Average<double>> accmap;

    double lnK;
    double bondenergy = 0;
    bool forward = true;  // reaction moving forwards or backwards? (or left or right)
    bool neutral = false; // true if only neutral molecules are involved in the reaction

    void _to_json(json &) const override;
    void _from_json(const json &) override;

    bool insertProducts(std::vector<ReactionData>::iterator);
    bool swapReaction(Change &, std::vector<ReactionData>::iterator);

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
