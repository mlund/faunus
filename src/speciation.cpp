#include "bonds.h"
#include "speciation.h"
#include "aux/iteratorsupport.h"

namespace Faunus {
namespace Move {

void SpeciationMove::_to_json(json &j) const {
    json &_j = j["reactions"];
    _j = json::object();
    for (auto [reaction, acceptance] : acceptance_map)
        _j[reaction] = {{"attempts", acceptance.cnt}, {"acceptance", acceptance.avg()}};
    Faunus::_roundjson(_j, 3);
}

void SpeciationMove::setOther(Tspace &ospc) { other_spc = &ospc; }

/**
 * This function is only performing checks
 * Check whether it is possible to insert products (are there any inactive ones?)
 * @todo Redundant as most (all?) of these checks are checkes in later functions
 */
bool SpeciationMove::checkInsertProducts() {
    [[maybe_unused]] auto [atomic_products, molecular_products] = reaction->getProducts();

    for (auto [molid, number_to_insert] : molecular_products) {
        auto molecule_list = spc.findMolecules(molid, Tspace::ALL);
        if (molecules[molid].atomic) {
            if (range_size(molecule_list) != 1) // There can be only one
                throw std::runtime_error("Bad definition: One group per atomic molecule!");
            auto group = molecule_list.begin();
            if ((group->size() + number_to_insert) > group->capacity()) { // insufficient atoms?
                faunus_logger->warn("molecule {} has reached its maximum capacity", Faunus::molecules[molid].name);
                return false; // Slip out the back door
            }
        } else {
            auto selection = (reaction->only_neutral_molecules) ? Tspace::INACTIVE_NEUTRAL : Tspace::INACTIVE;
            molecule_list = spc.findMolecules(molid, selection);
            if (range_size(molecule_list) < number_to_insert) {
                faunus_logger->warn("molecule {} has reached its maximum capacity", Faunus::molecules[molid].name);
                return false; // Not possible to perform change, escape through the back door
            }
        }
    }
    return true;
}

/**
 * Convert from one atom type to another in any group (atomic/molecular).
 * The reaction require that the `swap` is true and there must be *exactly*
 * one atomic reactant and one atomic product in the reaction.
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
 *          a future neightbour list implementation.
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
        assert(std::is_sorted(change_data.atoms.begin(), change_data.atoms.end())); // redundant
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
    assert(target.atomic == false);             // group must be molecular
    assert(target.size() == target.capacity()); // group must be active

    target.unwrap(spc.geo.getDistanceFunc()); // when in storage, remove PBC

    // Store internal bond energy of the deactivated molecule
    for (auto &bond : Faunus::molecules.at(target.id).bonds) {
        auto bond_clone = bond->clone();
        bond_clone->shift(std::distance(spc.p.begin(), target.begin()));
        Potential::setBondEnergyFunction(bond_clone, spc.p);
        bond_energy += bond_clone->energy(spc.geo.getDistanceFunc());
    }

    target.deactivate(target.begin(), target.end()); // deactivate whole group

    Change::data change_data; // describes the change
    change_data.internal = true;
    change_data.index = &target - &spc.groups.front(); // index of moved group
    change_data.all = true;                            // all atoms in group were moved
    for (int i = 0; i < target.capacity(); i++) {      // list of all changed atom index
        change_data.atoms.push_back(i);
    }
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
    assert(target.size() == 0); // must be inactive

    target.activate(target.inactive().begin(), target.inactive().end()); // activate all particles
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
        bond_energy -= bondclone->energy(spc.geo.getDistanceFunc());
    }

    Change::data d;                          // describes the changed - used for energy evaluation
    d.index = &target - &spc.groups.front(); // index* of moved group
    d.all = true;                            // all atoms in group were moved
    d.internal = true;
    d.atoms.reserve(target.capacity());

    for (size_t i = 0; i < target.capacity(); i++) { // list of changed atom index
        d.atoms.push_back(i);
    }
    return d;
}

void SpeciationMove::activateAllProducts(Change &change) {
    auto [atomic_products, molecular_products] = reaction->getProducts();

    for (auto [molid, number_to_insert] : molecular_products) {
        if (Faunus::molecules[molid].atomic) { // The product is an atom
            auto mollist = spc.findMolecules(molid, Tspace::ALL);
            if (range_size(mollist) > 0) {
                Change::data change_data = expandAtomicGroup(*mollist.begin(), number_to_insert);
                if (not change_data.atoms.empty()) {
                    change.groups.push_back(change_data);
                }
            } else {
                assert(false); // we should never reach here
            }
        } else { // The product is a molecule
            auto selection = (reaction->only_neutral_molecules) ? Tspace::INACTIVE_NEUTRAL : Tspace::INACTIVE;
            auto mollist = spc.findMolecules(molid, selection);
            if (range_size(mollist) >= number_to_insert) {
                for (int i = 0; i < number_to_insert; i++) {
                    auto target = slump.sample(mollist.begin(), mollist.end());
                    auto change_data = activateMolecularGroup(*target);
                    change.groups.push_back(change_data); // Add to list of moved groups
                }
            } else {
                faunus_logger->warn("maximum number of {} molecules reached; increase capacity?",
                                    Faunus::molecules[molid].name);
            }
        }
    }
}

/**
 * Deactivate all atomic and molecular reactants
 */
void SpeciationMove::deactivateAllReactants(Change &change) {

    auto [atomic_reactants, molecular_reactants] = reaction->getReactants();

    for (auto [molid, number_to_delete] : molecular_reactants) { // Delete
        if (molecules[molid].atomic) {                           // The reagents is an atomic group
            auto target = spc.findMolecules(molid, Tspace::ALL).begin();
            auto other_target = other_spc->findMolecules(molid, Tspace::ALL).begin();
            auto change_data = contractAtomicGroup(*target, *other_target, number_to_delete);
            if (not change_data.atoms.empty()) {
                change.groups.push_back(change_data);
            }
        } else { // The reactant is a molecule
            auto selection = (reaction->only_neutral_molecules) ? Tspace::ACTIVE_NEUTRAL : Tspace::ACTIVE;
            auto mol_list = spc.findMolecules(molid, selection);
            if (range_size(mol_list) >= number_to_delete) {
                for (int i = 0; i < number_to_delete; i++) {
                    auto target = slump.sample(mol_list.begin(), mol_list.end());
                    auto change_data = deactivateMolecularGroup(*target);
                    change.groups.push_back(change_data); // add to list of moved groups
                }
            } else {
                faunus_logger->warn("all {} molecules have been deleted; increase system volume?",
                                    Faunus::molecules[molid].name);
            }
        }
    }
}

void SpeciationMove::_move(Change &change) {
    assert(other_spc != nullptr);

    if (not Faunus::reactions.empty()) {                                                 // global list of reactions
        reaction = &(*slump.sample(Faunus::reactions.begin(), Faunus::reactions.end())); // random reaction
        auto direction = static_cast<ReactionData::Direction>((char)slump.range(0, 1));  // random direction
        reaction->setDirection(direction);

        if (reaction->empty()) { // Enforce canonic constraint if invoked
            return;              // Out of material, slip out the back door
        }
        if (checkInsertProducts() == false) {
            return;
        }
        if (atomicSwap(change) == false) {
            return;
        }

        bond_energy = 0;
        change.dN = true; // Attempting to change the number of atoms / molecules

        deactivateAllReactants(change);
        activateAllProducts(change);
        std::sort(change.groups.begin(), change.groups.end()); // why?
    } else {
        throw std::runtime_error("no reactions in list; add some or disable the `rcmc` move");
    }
}

double SpeciationMove::bias(Change &, double, double) {
    // The acceptance/rejection of the move is affected by the equilibrium constant
    // but unaffected by the change in bonded energy
    return -reaction->lnK + bond_energy;
}

void SpeciationMove::_accept(Change &) {
    acceptance_map[reaction->reaction_str] += 1;
    if (reaction->getDirection() == ReactionData::Direction::RIGHT) {
        reaction->N_reservoir -= 1;
    } else {
        reaction->N_reservoir += 1;
    }
    if (reaction->N_reservoir < 0 && reaction->canonic) {
        throw std::runtime_error("negative number of molecules");
    }
}

void SpeciationMove::_reject(Change &) { acceptance_map[reaction->reaction_str] += 0; }

SpeciationMove::SpeciationMove(Tspace &spc) : spc(spc) {
    name = "rcmc";
    cite = "doi:10/fqcpg3";
}
void SpeciationMove::_from_json(const json &) {}

} // end of namespace Move
} // end of namespace Faunus
