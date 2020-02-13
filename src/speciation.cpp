#include "bonds.h"
#include "speciation.h"
#include "aux/iteratorsupport.h"

namespace Faunus {
namespace Move {

void SpeciationMove::_to_json(json &j) const {
    json &_j = j["reactions"];
    _j = json::object();
    for (auto &m : acceptance_map)
        _j[m.first] = {{"attempts", m.second.cnt}, {"acceptance", m.second.avg()}};
    Faunus::_roundjson(_j, 3);
}

void SpeciationMove::setOther(Tspace &ospc) { otherspc = &ospc; }

/*
 * This function is only performing checks
 */
bool SpeciationMove::checkInsertProducts(std::vector<ReactionData>::iterator reaction) {
    // Check whether it is possible to insert products (are there any inactive ones?)
    for (auto [molid, number_to_insert] : reaction->moleculesToAdd(forward)) { // Additional checks
        auto mollist = spc.findMolecules(molid, Tspace::ALL);
        if (molecules[molid].atomic) {
            if (rng_size(mollist) != 1) // There can be only one
                throw std::runtime_error("Bad definition: One group per atomic molecule!");
            auto git = mollist.begin();
            if ((git->size() + number_to_insert) > git->capacity()) { // Assure that there are atoms enough in the group
                faunus_logger->warn("molecule {} has reached its maximum capacity", Faunus::molecules[molid].name);
                return false;                                 // Slip out the back door
            }
        } else {
            if (neutral) // Only neutral molecules react
                mollist = spc.findMolecules(molid, Tspace::INACTIVE_NEUTRAL);
            else
                mollist = spc.findMolecules(molid, Tspace::INACTIVE);
            if (rng_size(mollist) < number_to_insert) {
                faunus_logger->warn("molecule {} has reached its maximum capacity", Faunus::molecules[molid].name);
                return false; // Not possible to perform change, escape through the back door
            }
        }
    }
    return true;
}

bool SpeciationMove::swapReaction(Change &change, std::vector<ReactionData>::iterator reaction) {
    // Perform a swap reaction, e.g., the (de)protonation of an atomic species
    if (reaction->swap == true) {
        auto m1 = reaction->atomsToAdd(not forward); // Swap checks
        auto m2 = reaction->atomsToAdd(forward);
        assert((m1.size() == 1) and (m2.size() == 1) &&
               "Bad definition: Only 2 explicit atoms per reaction!"); // Swap A = B
        auto atomlist = spc.findAtoms(m1.begin()->first);
        if (ranges::cpp20::empty(atomlist))                        // Make sure that there are any active atoms to swap
            return false;                                          // Slip out the back door
        auto ait = slump.sample(atomlist.begin(), atomlist.end()); // Random particle iterator
        auto git = spc.findGroupContaining(*ait);

        Change::data d;
        d.atoms.push_back(Faunus::distance(git->begin(), ait)); // Index of particle rel. to group
        d.index = Faunus::distance(spc.groups.begin(), git);
        d.internal = true;
        d.dNswap = true;
        change.groups.push_back(d); // Add to list of moved groups

        // AtomData --> json --> Particle (performance warning!)
        Particle p = atoms.at(m2.begin()->first);

        p.pos = ait->pos;
        *ait = p;
        assert(ait->id == m2.begin()->first);
    }
    return true;
}

void SpeciationMove::deactivateReactants(Change &change, std::vector<ReactionData>::iterator reaction) {
    // Deactivate reactants
    for (auto [molid, N_delete] : reaction->moleculesToAdd(not forward)) { // Delete
        auto mollist = spc.findMolecules(molid, Tspace::ALL);

        // The reagents is an atom

        if (molecules[molid].atomic) {
            auto git = mollist.begin();

            // bail out in case the group is empty and nothing can be removed
            if (git->empty()) {
                faunus_logger->trace("attempting to remove atom(s) from empty group ({})",
                                     Faunus::atoms.at(molid).name);
                return;
            }

            auto othermollist = otherspc->findMolecules(molid, Tspace::ALL); // implies that new and old are in sync
            auto othergit = othermollist.begin();
            Change::data d;
            d.index = Faunus::distance(spc.groups.begin(), git); // integer *index* of moved group
            d.internal = true;
            d.dNatomic = true;
            // m.second is the stoichiometric coefficient
            for (int N = 0; N < N_delete; N++) {                              // deactivate m.second m.first atoms
                auto it_random_atom = slump.sample(git->begin(), git->end()); // iterator to random atom
                // Shuffle back to end, both in trial and new
                auto it_last_atom = git->end() - 1;                      // iterator to last atom
                int dist = Faunus::distance(it_random_atom, git->end()); // distance to random atom from end
                if (Faunus::distance(it_random_atom, it_last_atom) > 1) {
                    std::iter_swap(it_random_atom, it_last_atom);
                    std::iter_swap(othergit->end() - dist - N, othergit->end() - (1 + N));
                }
                d.atoms.push_back(Faunus::distance(git->begin(), it_last_atom));
                git->deactivate(it_last_atom, git->end());
            }
            std::sort(d.atoms.begin(), d.atoms.end());
            change.groups.push_back(d); // add to list of moved groups
        } else
        // The reactant is a molecule
        {
            if (neutral) // Only neutral molecules react
                mollist = spc.findMolecules(molid, Tspace::ACTIVE_NEUTRAL);
            else
                mollist = spc.findMolecules(molid, Tspace::ACTIVE);
            // m.second is the stoichiometric coefficient
            for (int N = 0; N < N_delete; N++) {
                auto git = slump.sample(mollist.begin(), mollist.end());
                git->unwrap(spc.geo.getDistanceFunc());
                // We store the bonded energy of the deactivated molecule
                // The change in bonded energy should not affect the acceptance/rejection of the move
                for (auto &bond : molecules.at(molid).bonds) {
                    auto bondclone = bond->clone();
                    bondclone->shift(std::distance(spc.p.begin(), git->begin()));
                    Potential::setBondEnergyFunction(bondclone, spc.p);
                    bond_energy += bondclone->energy(spc.geo.getDistanceFunc());
                }
                git->deactivate(git->begin(), git->end());
                Change::data d;
                d.index = Faunus::distance(spc.groups.begin(), git); // integer *index* of moved group
                d.all = true;                                        // *all* atoms in group were moved
                d.internal = true;
                for (int i = 0; i < git->capacity(); i++)
                    d.atoms.push_back(i);
                change.groups.push_back(d); // add to list of moved groups
            }
        }
    }
}

void SpeciationMove::activateProducts(Change &change, std::vector<ReactionData>::iterator reaction) {
    // Activate products
    for (auto [molid, number_to_insert] : reaction->moleculesToAdd(forward)) { // Add
        const MoleculeData &molecule_data = molecules.at(molid);
        auto mollist = spc.findMolecules(molid, Tspace::ALL);

        // The product is an atom
        if (molecule_data.atomic) {
            auto group_it = mollist.begin();
            Change::data d;
            d.index = Faunus::distance(spc.groups.begin(), group_it);
            d.internal = true;
            d.dNatomic = true;
            // We insert the atomic product at a random position
            // m.second is the stoichiometric coefficient
            for (int N = 0; N < number_to_insert; N++) { // Activate m.second m.first atoms
                group_it->activate(group_it->end(), group_it->end() + 1);
                auto it_last_atom = group_it->end() - 1;
                spc.geo.randompos(it_last_atom->pos, slump);
                spc.geo.getBoundaryFunc()(it_last_atom->pos);
                d.atoms.push_back(Faunus::distance(group_it->begin(), it_last_atom)); // Index of particle rel. to group
            }
            std::sort(d.atoms.begin(), d.atoms.end());
            change.groups.push_back(d); // Add to list of moved groups
        } else
        // The product is a molecule
        {
            if (neutral) // Only neutral molecules react
                mollist = spc.findMolecules(molid, Tspace::INACTIVE_NEUTRAL);
            else
                mollist = spc.findMolecules(molid, Tspace::INACTIVE);

            // m.second is the stoichiometric coefficient
            for (int N = 0; N < number_to_insert; N++) {
                auto git = slump.sample(mollist.begin(), mollist.end());
                git->activate(git->inactive().begin(), git->inactive().end());
                Point cm = git->cm;
                spc.geo.randompos(cm, slump);
                git->translate(cm, spc.geo.getBoundaryFunc());
                Point u = ranunit(slump);
                Eigen::Quaterniond Q(Eigen::AngleAxisd(2 * pc::pi * (slump() - 0.5), u));
                git->rotate(Q, spc.geo.getBoundaryFunc());

                // We store the bonded energy of the activated molecule
                // The change in bonded energy should not affect the acceptance/rejection of the move
                for (auto &bond : molecule_data.bonds) {
                    auto bondclone = bond->clone();
                    bondclone->shift(std::distance(spc.p.begin(), git->begin()));
                    Potential::setBondEnergyFunction(bondclone, spc.p);
                    bond_energy -= bondclone->energy(spc.geo.getDistanceFunc());
                }

                Change::data d;
                d.index = Faunus::distance(spc.groups.begin(), git); // Integer *index* of moved group
                d.all = true;                                        // All atoms in group were moved
                d.internal = true;

                for (int i = 0; i < git->capacity(); i++)
                    d.atoms.push_back(i);

                change.groups.push_back(d); // Add to list of moved groups
                assert(spc.geo.sqdist(git->cm, Geometry::massCenter(git->begin(), git->end(), spc.geo.getBoundaryFunc(),
                                                                    -git->cm)) < 1e-9);
            }
        }
    }
}

void SpeciationMove::_move(Change &change) {
    if (not reactions.empty()) {
        auto reaction = slump.sample(reactions.begin(), reactions.end()); // select random reaction
        lnK = reaction->lnK;
        neutral = reaction->neutral;       // If true, only neutral molecules participate in the reaction
        forward = (bool)slump.range(0, 1); // Random boolean
        trialprocess = &(*reaction);

        if (reaction->empty(forward)) // Enforce canonic constraint if invoked
            return;              // Out of material, slip out the back door
        if (checkInsertProducts(reaction) == false)
            return;
        if (swapReaction(change, reaction) == false)
            return;

        bond_energy = 0;
        change.dN = true; // Attempting to change the number of atoms / molecules

        deactivateReactants(change, reaction);
        activateProducts(change, reaction);

        std::sort(change.groups.begin(), change.groups.end()); // why?
    } else
        throw std::runtime_error("No reactions in list, disable rcmc or add reactions");
}

double SpeciationMove::bias(Change &, double, double) {
    // The acceptance/rejection of the move is affected by the equilibrium constant
    // but unaffected by the change in bonded energy
    if (forward)
        return -lnK + bond_energy;
    return lnK + bond_energy;
}

void SpeciationMove::_accept(Change &) {
    acceptance_map[trialprocess->reaction] += 1;
    trialprocess->N_reservoir += (forward == true) ? -1 : 1;
    if (trialprocess->N_reservoir < 0 && trialprocess->canonic)
        throw std::runtime_error("There are no negative number of molecules");
}

void SpeciationMove::_reject(Change &) { acceptance_map[trialprocess->reaction] += 0; }

SpeciationMove::SpeciationMove(Tspace &spc) : spc(spc) {
    name = "rcmc";
    cite = "doi:10/fqcpg3";
}
void SpeciationMove::_from_json(const json &) {}

} // end of namespace Move
} // end of namespace Faunus
