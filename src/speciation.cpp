#include "bonds.h"
#include "speciation.h"
#include "aux/iteratorsupport.h"
#include <algorithm>

namespace Faunus {
namespace Move {

void SpeciationMove::_to_json(json &j) const {
    json &_j = j["reactions"];
    _j = json::object();
    for (auto [reaction, data] : acceptance) {
        _j[reaction->reaction_str] = {{"attempts", data.left.cnt + data.right.cnt},
                                      {"acceptance -->", data.right.avg()},
                                      {"acceptance <--", data.left.avg()}};
    }
    for (auto [molid, size] : average_reservoir_size) {
        j["implicit_reservoir"][molecules[molid].name] = size.avg();
    }
}

void SpeciationMove::setOther(Tspace &ospc) { other_spc = &ospc; }

/**
 * Convert from one atom type to another in any group (atomic/molecular).
 * The reaction requires that `swap` is true and there must be *exactly*
 * one atomic reactant and one atomic product in the reaction.
 *
 * The function checks if there are sufficient atomic and molecular
 * reactants and products to perform the move. If not, false is returned
 * and the system is left untouched.
 *
 * @todo If particle has extended properties, make sure to copy the state of those
 */
bool SpeciationMove::atomicSwap(Change &change) {
    if (reaction->swap) {
        [[maybe_unused]] auto [atomic_products, molecular_products] = reaction->getProducts();
        [[maybe_unused]] auto [atomic_reactants, molecular_reactants] = reaction->getReactants();

        assert(atomic_products.size() == 1 and atomic_reactants.size() == 1);

        auto atomlist = spc.findAtoms(atomic_reactants.begin()->first); // search all active molecules
        if (ranges::cpp20::empty(atomlist)) { // Make sure that there are any active atoms to swap
            return false;                     // Slip out the back door
        }

        if (!molecular_reactants.empty()) {          // enough molecular reactants?
            assert(molecular_reactants.size() == 1); // only one allowed this far
            auto [molid, N] = *molecular_reactants.begin();
            if (Faunus::molecules[molid].atomic) {
                auto mollist = spc.findMolecules(molid, Tspace::ALL); // look for a single atomic group
                assert(range_size(mollist) == 1);
                if (mollist.begin()->empty()) {
                    return false;
                }
            } else { // reactant is a molecular group
                auto mollist = spc.findMolecules(molid, Tspace::ACTIVE);
                if (range_size(mollist) < N) {
                    return false;
                }
            }
        }

        if (!molecular_products.empty()) {          // enough inactive molecular products?
            assert(molecular_products.size() == 1); // only one allowed this far
            auto [molid, N] = *molecular_products.begin();
            if (Faunus::molecules[molid].atomic) { // product is an atomic group
                auto mollist = spc.findMolecules(molid, Tspace::ALL);
                assert(range_size(mollist) == 1);
                if (mollist.begin()->capacity() - mollist.begin()->size() < N) {
                    return false;
                }
            } else { // we're producing a molecular group
                auto mollist = spc.findMolecules(molid, Tspace::INACTIVE);
                if (range_size(mollist) < N) {
                    return false;
                }
            }
        }
        auto random_particle = slump.sample(atomlist.begin(), atomlist.end()); // target particle to swap
        auto group = spc.findGroupContaining(*random_particle);                // find enclosing group

        Change::data d; // describe what has change - used for energy cal.
        d.atoms.push_back(Faunus::distance(group->begin(), random_particle)); // Index of particle rel. to group
        d.index = Faunus::distance(spc.groups.begin(), group); // index of particle in group (starting from zero)
        d.internal = true;
        d.dNswap = true;
        change.groups.push_back(d); // Add to list of moved groups

        int atomid = atomic_products.begin()->first; // atomid of new atom type
        Particle p = Faunus::atoms[atomid];          // temporary particle of new type
        p.pos = random_particle->pos;                // get position from old particle
        // todo: extended properties, dipole etc?
        *random_particle = p; // copy new particle onto old particle
        assert(random_particle->id == atomid);
    }
    return true;
}

/**
 * Reduce an atomic group by `number_to_delete` particles. The deleted particles are
 * picked by random; moved the the end of the group; then deactivated. In order for
 * the Hamiltonian to pick up the energy change, particles in the reference Space
 * (`old_target`) are swapped to the same index, albeit not deactivated.
 *
 * @warning Directly modifying the groups in spc and otherspc might interfere with
 *          a future neighbour list implementation.
 */
Change::data SpeciationMove::contractAtomicGroup(Space::Tgroup &target, Space::Tgroup &old_target,
                                                 int number_to_delete) {
    assert(target.atomic);
    Change::data change_data; // describes what has changed

    if ((int)target.size() - number_to_delete >= 0) {
        change_data.index = &target - &spc.groups.front(); // index of moved group
        change_data.internal = true;
        change_data.dNatomic = true;
        for (int i = 0; i < number_to_delete; i++) {
            auto atom_to_delete = slump.sample(target.begin(), target.end()); // iterator to atom to delete
            auto last_atom = target.end() - 1;                                // iterator to last atom
            int dist = std::distance(atom_to_delete, target.end());           // distance to atom from end

            if (std::distance(atom_to_delete, last_atom) > 1) { // Shuffle back to end, both in trial and old target
                std::iter_swap(atom_to_delete, last_atom);
                std::iter_swap(old_target.end() - dist - i, old_target.end() - (1 + i));
            }

            change_data.atoms.push_back(std::distance(target.begin(), last_atom));
            target.deactivate(last_atom, target.end()); // deactivate a single atom at the time
        }
        std::sort(change_data.atoms.begin(), change_data.atoms.end());
    } else {
        faunus_logger->warn("atomic group {} is depleted; increase simulation volume?",
                            Faunus::molecules[target.id].name);
    }
    return change_data;
}

/**
 * Deactivate a single, active molecular group. If there are internal bonds, the total
 * bond energy is stored and used to avoid that the bond energy affect acceptance.
 * When deactivated, the molecule is made whole, i.e. periodic boundary conditions
 * are removed as this cannot be achieved later if the system volume changes.
 */
Change::data SpeciationMove::deactivateMolecularGroup(Space::Tgroup &target) {
    assert(target.atomic == false); // group must be molecular
    assert(not target.empty());
    assert(target.size() == target.capacity()); // group must be active

    target.unwrap(spc.geo.getDistanceFunc()); // when in storage, remove PBC

    // Store internal bond energy of the deactivated molecule
    for (auto &bond : Faunus::molecules.at(target.id).bonds) {
        auto bond_clone = bond->clone();
        bond_clone->shift(std::distance(spc.p.begin(), target.begin()));
        Potential::setBondEnergyFunction(bond_clone, spc.p);
        bond_energy += bond_clone->energyFunc(spc.geo.getDistanceFunc());
    }

    target.deactivate(target.begin(), target.end()); // deactivate whole group
    assert(target.empty());

    Change::data change_data; // describes the change
    change_data.internal = true;
    change_data.index = &target - &spc.groups.front(); // index of moved group
    change_data.all = true;                            // all atoms in group were moved
    change_data.atoms.resize(target.capacity());       // list of changed atom index
    std::iota(change_data.atoms.begin(), change_data.atoms.end(), 0);

    return change_data;
}

/**
 * Expand a single molecular group by `number_to_insert` particles by activating inactive
 * particles at the end of the group. The activated particles are assigned new
 * random positions, guaranteed to fall within the simulation box.
 *
 * If the capacity of the group will be exceeded, a warning is issued and
 * the returned Change::data object will be empty.
 */
Change::data SpeciationMove::expandAtomicGroup(Space::Tgroup &target, int number_to_insert) {
    assert(target.atomic);

    Change::data change_data;
    if (target.size() + number_to_insert <= target.capacity()) {
        change_data.index = &target - &spc.groups.front();
        change_data.internal = true;
        change_data.dNatomic = true;
        for (int i = 0; i < number_to_insert; i++) {
            target.activate(target.end(), target.end() + 1); // activate one particle
            auto last_atom = target.end() - 1;
            spc.geo.randompos(last_atom->pos, slump);                              // give it a random position
            spc.geo.getBoundaryFunc()(last_atom->pos);                             // apply PBC if needed
            change_data.atoms.push_back(std::distance(target.begin(), last_atom)); // index relative to group
        }
    } else {
        faunus_logger->warn("atomic group {} is full; increase capacity?", Faunus::molecules[target.id].name);
    }
    return change_data;
}

/**
 * Activate a single inactive molecule and assign a new random position and orientation.
 * If the molecule has internal bonds, the bond-energy is calculated to ensure that
 * the bond-energy does not affect the insertion acceptance.
 */
Change::data SpeciationMove::activateMolecularGroup(Space::Tgroup &target) {
    assert(not target.atomic);  // must be a molecule group
    assert(target.empty());     // must be inactive
    target.activate(target.inactive().begin(), target.inactive().end()); // activate all particles
    assert(not target.empty());

    Point cm = target.cm;
    spc.geo.randompos(cm, slump);                    // generate random position
    target.translate(cm, spc.geo.getBoundaryFunc()); // assign random position to mass-center
    Point u = ranunit(slump);                        // random unit vector
    Eigen::Quaterniond Q(Eigen::AngleAxisd(2 * pc::pi * (slump() - 0.5), u));
    target.rotate(Q, spc.geo.getBoundaryFunc()); // assign random orientation

    assert(spc.geo.sqdist(target.cm, Geometry::massCenter(target.begin(), target.end(), spc.geo.getBoundaryFunc(),
                                                          -target.cm)) < 1e-9);

    // Store internal bond energy of activated molecule
    for (auto &bond : Faunus::molecules[target.id].bonds) {
        auto bondclone = bond->clone();
        bondclone->shift(std::distance(spc.p.begin(), target.begin()));
        Potential::setBondEnergyFunction(bondclone, spc.p);
        bond_energy -= bondclone->energyFunc(spc.geo.getDistanceFunc());
    }

    Change::data d;                          // describes the changed - used for energy evaluation
    d.index = &target - &spc.groups.front(); // index* of moved group
    d.all = true;                            // all atoms in group were moved
    d.internal = true;
    d.atoms.resize(target.capacity()); // list of changed atom index
    std::iota(d.atoms.begin(), d.atoms.end(), 0);

    return d;
}

bool SpeciationMove::activateAllProducts(Change &change) {
    auto [atomic_products, molecular_products] = reaction->getProducts();

    for (auto [molid, number_to_insert] : molecular_products) {
        if (number_to_insert == 0 or Faunus::molecules[molid].isImplicit()) {
            continue;                                 // implicit molecules are added/removed after the move
        } else if (Faunus::molecules[molid].atomic) { // The product is an atom
            auto mollist = spc.findMolecules(molid, Tspace::ALL);
            if (range_size(mollist) > 0) {
                Change::data change_data = expandAtomicGroup(*mollist.begin(), number_to_insert);
                if (not change_data.atoms.empty()) {
                    change.groups.push_back(change_data);
                } else {
                    return false;
                }
            } else {
                assert(false); // we should never reach here
            }
        } else { // The product is a molecule
            auto selection = (reaction->only_neutral_molecules) ? Tspace::INACTIVE_NEUTRAL : Tspace::INACTIVE;
            auto inactive = spc.findMolecules(molid, selection); // all inactive molecules
            std::vector<std::reference_wrapper<Tspace::Tgroup>> molecules_to_activate;
            std::sample(inactive.begin(), inactive.end(), std::back_inserter(molecules_to_activate), number_to_insert,
                        slump.engine);
            if (molecules_to_activate.size() == number_to_insert) {
                for (auto &target : molecules_to_activate) {
                    if (auto change_data = activateMolecularGroup(target); not change_data.atoms.empty()) {
                        change.groups.push_back(change_data); // Add to list of moved groups
                    } else {
                        return false;
                    }
                }
            } else {
                faunus_logger->warn("maximum number of {} molecules reached; increase capacity?",
                                    Faunus::molecules[molid].name);
                return false;
            }
        }
    }
    return true;
}

/**
 * If possible, deactivate all atomic and molecular reactants.
 *
 * If impossible due to lack of reactant - either
 * implicit or explicit - the function returns `false`.
 */
bool SpeciationMove::deactivateAllReactants(Change &change) {

    if (not enoughImplicitMolecules()) {
        return false;
    }

    auto [atomic_reactants, molecular_reactants] = reaction->getReactants();

    for (auto [molid, N_delete] : molecular_reactants) { // Delete
        if (N_delete <= 0 or Faunus::molecules[molid].isImplicit()) {
            continue;                         // implicit molecules are added/deleted after move
        } else if (molecules[molid].atomic) { // reactant is an atomic group
            auto mollist = spc.findMolecules(molid, Tspace::ALL);
            assert(range_size(mollist) == 1);
            auto target = spc.findMolecules(molid, Tspace::ALL).begin();
            auto other_target = other_spc->findMolecules(molid, Tspace::ALL).begin();
            auto change_data = contractAtomicGroup(*target, *other_target, N_delete);
            if (not change_data.atoms.empty()) {
                change.groups.push_back(change_data);
            } else {
                return false;
            }
        } else { // molecular reactant (non-atomic)
            auto selection = (reaction->only_neutral_molecules) ? Tspace::ACTIVE_NEUTRAL : Tspace::ACTIVE;
            auto active = spc.findMolecules(molid, selection);
            std::vector<std::reference_wrapper<Tspace::Tgroup>> molecules_to_deactivate;
            std::sample(active.begin(), active.end(), std::back_inserter(molecules_to_deactivate), N_delete,
                        slump.engine); // pick random molecules to delete
            if (molecules_to_deactivate.size() == N_delete) {
                for (auto &target : molecules_to_deactivate) {
                    if (auto change_data = deactivateMolecularGroup(target); not change_data.atoms.empty()) {
                        change.groups.push_back(change_data); // add to list of moved groups
                    } else {
                        return false;
                    }
                }
            } else {
                if (Faunus::molecules[molid].activity > 0) { // if molecule has an activity, issue a warning
                    faunus_logger->warn("all {} molecules have been deleted; increase system volume?",
                                        Faunus::molecules[molid].name);
                }
                return false;
            }
        }
    }
    return true;
}

/**
 * Checks if there is enough implicit molecules to carry out the reaction
 */
bool SpeciationMove::enoughImplicitMolecules() const {
    auto tooSmall = [&](auto &i) {                                // check if species has enough implicit
        auto [molid, nu] = i;                                     // molid and stoichiometric coefficient
        if (Faunus::molecules[molid].isImplicit()) {              // matter to perform process
            assert(spc.getImplicitReservoir().count(molid) == 1); // must be registered in space!
            if (spc.getImplicitReservoir()[molid] - nu < 0) {
                return true; // nope, no more material
            }
        }
        return false;
    };

    [[maybe_unused]] auto [atomic_reactants, molecular_reactants] = reaction->getReactants();

    if (std::count_if(molecular_reactants.begin(), molecular_reactants.end(), tooSmall) > 0) {
        return false;
    } else {
        return true;
    }
}

void SpeciationMove::_move(Change &change) {
    assert(other_spc != nullptr);        // knowledge of other space should be provided by now
    if (not Faunus::reactions.empty()) { // global list of reactions
        reaction = slump.sample(Faunus::reactions.begin(), Faunus::reactions.end()); // random reaction
        assert(reaction != Faunus::reactions.end());
        auto direction = static_cast<ReactionData::Direction>((char)slump.range(0, 1)); // random direction
        reaction->setDirection(direction);

        bond_energy = 0;

        if (atomicSwap(change)) {
            if (deactivateAllReactants(change)) {
                if (activateAllProducts(change)) {
                    assert(not change.empty());
                    change.dN = true; // Attempting to change the number of atoms / molecules
                    std::sort(change.groups.begin(), change.groups.end()); // why?
                    return;
                }
            }
        }
    }
    change.clear();
}

double SpeciationMove::bias(Change &change, double, double) {
    // The acceptance/rejection of the move is affected by the equilibrium constant
    // but unaffected by the change in bonded energy
    for (auto change_data : change.groups) {         // loop over changed groups/molecules
        auto &group = spc.groups[change_data.index]; // current group
        for (int i : change_data.atoms) {            // loop over changed particles in group
            double z = group[i].pos.z();             // i'th particle z-position (Å)
            if (z > -15.0 and z < 15.0) {            // is it inside membrane? (hard-coded!)
                return pc::infty;                    // if yes, reject using infinite energy
            }
        }
    }
    return -reaction->lnK + bond_energy;
}

void SpeciationMove::_accept(Change &) {
    acceptance[reaction].update(reaction->getDirection(), true);

    [[maybe_unused]] auto [atomic_products, molecular_products] = reaction->getProducts();
    [[maybe_unused]] auto [atomic_reactants, molecular_reactants] = reaction->getReactants();

    // adjust amount of implicit matter
    for (auto [molid, nu] : molecular_reactants) {
        if (Faunus::molecules[molid].isImplicit()) {
            spc.getImplicitReservoir()[molid] -= nu;
            other_spc->getImplicitReservoir()[molid] -= nu;
            average_reservoir_size[molid] += spc.getImplicitReservoir()[molid];
            assert(spc.getImplicitReservoir()[molid] >= 0);
            assert(spc.getImplicitReservoir()[molid] == other_spc->getImplicitReservoir()[molid]);
        }
    }
    for (auto [molid, nu] : molecular_products) {
        if (Faunus::molecules[molid].isImplicit()) {
            spc.getImplicitReservoir()[molid] += nu;
            other_spc->getImplicitReservoir()[molid] += nu;
            average_reservoir_size[molid] += spc.getImplicitReservoir()[molid];
            assert(spc.getImplicitReservoir()[molid] == other_spc->getImplicitReservoir()[molid]);
        }
    }
}

void SpeciationMove::_reject(Change &) {
    acceptance[reaction].update(reaction->getDirection(), false);

    [[maybe_unused]] auto [atomic_products, molecular_products] = reaction->getProducts();
    [[maybe_unused]] auto [atomic_reactants, molecular_reactants] = reaction->getReactants();

    // average number of implicit molecules
    for (auto [molid, nu] : molecular_reactants) {
        if (Faunus::molecules[molid].isImplicit()) {
            average_reservoir_size[molid] += spc.getImplicitReservoir()[molid];
        }
    }
    for (auto [molid, nu] : molecular_products) {
        if (Faunus::molecules[molid].isImplicit()) {
            average_reservoir_size[molid] += spc.getImplicitReservoir()[molid];
        }
    }
}

SpeciationMove::SpeciationMove(Tspace &spc) : spc(spc) {
    name = "rcmc";
    cite = "doi:10/fqcpg3";
}
void SpeciationMove::_from_json(const json &) {}

} // end of namespace Move
} // end of namespace Faunus
