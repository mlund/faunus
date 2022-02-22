#include "bonds.h"
#include "speciation.h"
#include "aux/iteratorsupport.h"
#include <algorithm>

namespace Faunus {
namespace Move {

/**
 * @brief Internal exception to carry errors preventing the speciation move.
 */
struct SpeciationMoveException : public std::exception {};

void SpeciationMove::_to_json(json &j) const {
    json &_j = j["reactions"];
    _j = json::object();
    for (auto [reaction, data] : acceptance) {
        _j[reaction->reaction_str] = {{"attempts", data.left.size() + data.right.size()},
                                      {"acceptance -->", data.right.avg()},
                                      {"acceptance <--", data.left.avg()}};
    }
    for (auto [molid, size] : average_reservoir_size) {
        j["implicit_reservoir"][molecules.at(molid).name] = size.avg();
    }
}

/**
 * @param other Space representing the "old" state in a MC move
 */
void SpeciationMove::setOther(Space& other) { other_spc = &other; }

/**
 * Convert from one atom type to another in any group (atomic/molecular).
 * The reaction requires that `swap` is true and there must be *exactly*
 * one atomic reactant and one atomic product in the reaction.
 *
 * The function checks if there are sufficient atomic and molecular
 * reactants and products to perform the move. If not, an exception is
 * thrown and the system is left untouched.
 *
 * @throw SpeciationMoveException  when atomic swap cannot be performed
 *
 * @todo If particle has extended properties, make sure to copy the state of those
 */
void SpeciationMove::atomicSwap(Change& change) {
    if (!reaction->swap) {
        return;
    }
    [[maybe_unused]] const auto [atomic_products, molecular_products] = reaction->getProducts();
    [[maybe_unused]] const auto [atomic_reactants, molecular_reactants] = reaction->getReactants();

    assert(atomic_products.size() == 1 and atomic_reactants.size() == 1);

    auto atomlist = spc.findAtoms(atomic_reactants.begin()->first); // search all active molecules
    if (ranges::cpp20::empty(atomlist)) {                           // Make sure that there are any active atoms to swap
        throw SpeciationMoveException();
    }
    auto& target_particle = *slump.sample(atomlist.begin(), atomlist.end()); // target particle to swap
    auto& group = *spc.findGroupContaining(target_particle);                 // find enclosing group

    if (!molecular_reactants.empty()) {          // enough molecular reactants?
        assert(molecular_reactants.size() == 1); // only one allowed this far
        const auto [molid, N] = *molecular_reactants.begin();
        if (Faunus::molecules.at(molid).isAtomic()) {
            auto mollist = spc.findMolecules(molid, Space::Selection::ALL); // look for a single atomic group
            assert(range_size(mollist) == 1);
            if (mollist.begin()->empty()) {
                throw SpeciationMoveException();
            }
        } else { // reactant is a molecular group
            auto mollist = spc.findMolecules(molid, Space::Selection::ACTIVE);
            if (range_size(mollist) < N) {
                throw SpeciationMoveException();
            }
        }
    }

    if (!molecular_products.empty()) {          // enough inactive molecular products?
        assert(molecular_products.size() == 1); // only one allowed this far
        const auto [molid, N] = *molecular_products.begin();
        if (Faunus::molecules.at(molid).isAtomic()) { // product is an atomic group
            auto mollist = spc.findMolecules(molid, Space::Selection::ALL);
            assert(range_size(mollist) == 1);
            if (mollist.begin()->capacity() - mollist.begin()->size() < N) {
                throw SpeciationMoveException();
            }
        } else { // we're producing a molecular group
            auto mollist = spc.findMolecules(molid, Space::Selection::INACTIVE);
            if (range_size(mollist) < N) {
                throw SpeciationMoveException();
            }
        }
    }

    auto& group_change = change.groups.emplace_back();
    group_change.relative_atom_indices.push_back(group.getParticleIndex(target_particle));
    group_change.group_index = spc.getGroupIndex(group);
    group_change.internal = true;
    group_change.dNswap = true;

    const auto new_atomid = atomic_products.begin()->first;  // atomid of new atom type
    swapParticleProperties(target_particle, new_atomid);
}

/**
 * @brief Swap particle to another type
 * @param target_particle The particle to update
 * @param new_atomid The new atomtype to apply
 *
 * This will keep original positions and internal orientation of particles, while swapping
 * the atomid and other properties related to the new atom type
 *
 * @todo This could make use of policies to customize how particles should be updated
 */
void SpeciationMove::swapParticleProperties(Particle& target_particle, const int new_atomid) const {
    Particle source_particle = atoms.at(new_atomid); // temporary particle of new type
    source_particle.pos = target_particle.pos;       // keep original positio
    if (source_particle.hasExtension() || target_particle.hasExtension()) { // keep other orientational data
        auto& source_ext = source_particle.getExt();
        auto& target_ext = target_particle.getExt();
        if (target_ext.isCylindrical() || source_ext.isCylindrical()) {
            source_ext.setDirections(source_particle.traits().sphero_cylinder, target_ext.scdir, target_ext.patchdir);
        }
        if (target_ext.isDipolar() || source_ext.isQuadrupolar()) {
            throw std::runtime_error("dipolar/quadrupolar properties not yet implemented");
        }
    }
    target_particle = source_particle; // copy new particle onto old particle
    assert(target_particle->id == new_atomid);
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
Change::GroupChange SpeciationMove::contractAtomicGroup(Space::GroupType& target, Space::GroupType& old_target,
                                                        int number_to_delete) {
    assert(target.isAtomic());
    Change::GroupChange change_data; // describes what has changed

    if ((int)target.size() - number_to_delete >= 0) {
        change_data.group_index = &target - &spc.groups.front(); // index of moved group
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

            change_data.relative_atom_indices.push_back(std::distance(target.begin(), last_atom));
            target.deactivate(last_atom, target.end()); // deactivate a single atom at the time
        }
        std::sort(change_data.relative_atom_indices.begin(), change_data.relative_atom_indices.end());
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
Change::GroupChange SpeciationMove::deactivateMolecularGroup(Space::GroupType& target) {
    assert(target.isAtomic() == false); // group must be molecular
    assert(not target.empty());
    assert(target.size() == target.capacity()); // group must be active

    target.unwrap(spc.geometry.getDistanceFunc()); // when in storage, remove PBC

    // Store internal bond energy of the deactivated molecule
    for (auto &bond : Faunus::molecules.at(target.id).bonds) {
        auto bond_clone = bond->clone();
        bond_clone->shiftIndices(std::distance(spc.particles.begin(), target.begin()));
        bond_clone->setEnergyFunction(spc.particles);
        bond_energy += bond_clone->energyFunc(spc.geometry.getDistanceFunc());
    }

    target.deactivate(target.begin(), target.end()); // deactivate whole group
    assert(target.empty());

    Change::GroupChange change_data; // describes the change
    change_data.internal = true;
    change_data.group_index = &target - &spc.groups.front(); // index of moved group
    change_data.all = true;                            // all atoms in group were moved
    change_data.relative_atom_indices.resize(target.capacity()); // list of changed atom index
    std::iota(change_data.relative_atom_indices.begin(), change_data.relative_atom_indices.end(), 0);

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
Change::GroupChange SpeciationMove::expandAtomicGroup(Space::GroupType& target, int number_to_insert) {
    assert(target.isAtomic());

    Change::GroupChange change_data;
    if (target.size() + number_to_insert <= target.capacity()) {
        change_data.group_index = &target - &spc.groups.front();
        change_data.internal = true;
        change_data.dNatomic = true;
        for (int i = 0; i < number_to_insert; i++) {
            target.activate(target.end(), target.end() + 1); // activate one particle
            auto last_atom = target.end() - 1;
            spc.geometry.randompos(last_atom->pos, slump);  // give it a random position
            spc.geometry.getBoundaryFunc()(last_atom->pos); // apply PBC if needed
            change_data.relative_atom_indices.push_back(
                std::distance(target.begin(), last_atom)); // index relative to group
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
Change::GroupChange SpeciationMove::activateMolecularGroup(Space::GroupType& target) {
    assert(not target.isAtomic()); // must be a molecule group
    assert(target.empty());     // must be inactive
    target.activate(target.inactive().begin(), target.inactive().end()); // activate all particles
    assert(not target.empty());

    Point cm = target.mass_center;
    spc.geometry.randompos(cm, slump);                    // generate random position
    target.translate(cm, spc.geometry.getBoundaryFunc()); // assign random position to mass-center
    Point u = randomUnitVector(slump);               // random unit vector
    Eigen::Quaterniond Q(Eigen::AngleAxisd(2 * pc::pi * (slump() - 0.5), u));
    target.rotate(Q, spc.geometry.getBoundaryFunc()); // assign random orientation

    assert(spc.geometry.sqdist(target.mass_center,
                               Geometry::massCenter(target.begin(), target.end(), spc.geometry.getBoundaryFunc(),
                                                    -target.mass_center)) < 1e-9);

    // Store internal bond energy of activated molecule
    for (auto &bond : Faunus::molecules[target.id].bonds) {
        auto bondclone = bond->clone();
        bondclone->shiftIndices(std::distance(spc.particles.begin(), target.begin()));
        bondclone->setEnergyFunction(spc.particles);
        bond_energy -= bondclone->energyFunc(spc.geometry.getDistanceFunc());
    }

    Change::GroupChange d;                         // describes the changed - used for energy evaluation
    d.group_index = &target - &spc.groups.front(); // index* of moved group
    d.all = true;                            // all atoms in group were moved
    d.internal = true;
    d.relative_atom_indices.resize(target.capacity()); // list of changed atom index
    std::iota(d.relative_atom_indices.begin(), d.relative_atom_indices.end(), 0);

    return d;
}

/**
 * @brief Activate all atomic and molecular products.
 *
 * @throw SpeciationMoveException  when maximal capacity for products is exceeded
 */
void SpeciationMove::activateAllProducts(Change &change) {
    auto [atomic_products, molecular_products] = reaction->getProducts();

    // First we need to check if there are enough inactive products
    // for *all* participating molecules. Check only, *no* actual activation.
    for (auto [molid, number_to_insert] : molecular_products) {
        if (number_to_insert == 0 or Faunus::molecules[molid].isImplicit()) {
            continue;                                 // implicit molecules are added/removed after the move
        } else if (Faunus::molecules[molid].atomic) { // The product is an atom
            auto mollist = spc.findMolecules(molid, Space::Selection::ALL);
            if (range_size(mollist) > 0) {
                if (mollist.begin()->size() + number_to_insert > mollist.begin()->capacity()) {
                    faunus_logger->warn("atomic molecule {} is full; increase capacity?",
                                        Faunus::molecules[molid].name);
                    throw SpeciationMoveException();
                }
            } else {
                assert(false); // we should never reach here
            }
        } else { // The product is a molecule
            auto selection =
                (reaction->only_neutral_molecules) ? Space::Selection::INACTIVE_NEUTRAL : Space::Selection::INACTIVE;
            auto inactive = spc.findMolecules(molid, selection); // all inactive molecules
            std::vector<std::reference_wrapper<Space::GroupType>> molecules_to_activate;
            std::sample(inactive.begin(), inactive.end(), std::back_inserter(molecules_to_activate), number_to_insert,
                        slump.engine);
            if (molecules_to_activate.size() != number_to_insert) {
                faunus_logger->warn("maximum number of {} molecules reached; increase capacity?",
                                    Faunus::molecules[molid].name);
                throw SpeciationMoveException();
            }
        }
    }

    // Actually activate the products. This could be optimized as it also checks
    // the availability of products which is also done above
    for (auto [molid, number_to_insert] : molecular_products) {
        if (number_to_insert == 0 or Faunus::molecules[molid].isImplicit()) {
            continue;                                 // implicit molecules are added/removed after the move
        } else if (Faunus::molecules[molid].atomic) { // The product is an atom
            auto mollist = spc.findMolecules(molid, Space::Selection::ALL);
            if (range_size(mollist) > 0) {
                Change::GroupChange change_data = expandAtomicGroup(*mollist.begin(), number_to_insert);
                if (not change_data.relative_atom_indices.empty()) {
                    change.groups.push_back(change_data);
                } else {
                    assert(false); // we should never reach here
                }
            } else {
                assert(false); // we should never reach here
            }
        } else { // The product is a molecule
            auto selection =
                (reaction->only_neutral_molecules) ? Space::Selection::INACTIVE_NEUTRAL : Space::Selection::INACTIVE;
            auto inactive = spc.findMolecules(molid, selection); // all inactive molecules
            std::vector<std::reference_wrapper<Space::GroupType>> molecules_to_activate;
            std::sample(inactive.begin(), inactive.end(), std::back_inserter(molecules_to_activate), number_to_insert,
                        slump.engine);
            if (molecules_to_activate.size() == number_to_insert) {
                for (auto &target : molecules_to_activate) {
                    if (auto change_data = activateMolecularGroup(target);
                        not change_data.relative_atom_indices.empty()) {
                        change.groups.push_back(change_data); // Add to list of moved groups
                    } else {
                        assert(false); // we should never reach here
                    }
                }
            } else {
                assert(false); // we should never reach here
            }
        }
    }
}

/**
 * @brief Deactivate all atomic and molecular reactants.
 *
 * @throw SpeciationMoveException  when impossible due to the lack of reactant – either implicit or explicit
 */
void SpeciationMove::deactivateAllReactants(Change &change) {
    if (not enoughImplicitMolecules()) {
        throw SpeciationMoveException();
    }

    auto [atomic_reactants, molecular_reactants] = reaction->getReactants();
    // check (only!) if there are anything to deactivate
    for (auto [molid, N_delete] : molecular_reactants) { // Delete
        if (N_delete <= 0 or Faunus::molecules[molid].isImplicit()) {
            continue;                         // implicit molecules are added/deleted after move
        } else if (Faunus::molecules[molid].atomic) { // reactant is an atomic group
            auto mollist = spc.findMolecules(molid, Space::Selection::ALL);
            assert(range_size(mollist) == 1);
            auto target = spc.findMolecules(molid, Space::Selection::ALL).begin();
            if ((int)target->size() - N_delete < 0) {
                faunus_logger->warn("atomic group {} is depleted; increase simulation volume?",
                                    Faunus::molecules[molid].name);
                throw SpeciationMoveException();
            }
        } else { // molecular reactant (non-atomic)
            auto selection =
                (reaction->only_neutral_molecules) ? Space::Selection::ACTIVE_NEUTRAL : Space::Selection::ACTIVE;
            auto active = spc.findMolecules(molid, selection);
            std::vector<std::reference_wrapper<Space::GroupType>> molecules_to_deactivate;
            std::sample(active.begin(), active.end(), std::back_inserter(molecules_to_deactivate), N_delete,
                        slump.engine); // pick random molecules to delete
            if (molecules_to_deactivate.size() != N_delete) {
                throw SpeciationMoveException();
            }
        }
    }

    // perform actual deactivation
    for (auto [molid, N_delete] : molecular_reactants) { // Delete
        if (N_delete <= 0 or Faunus::molecules[molid].isImplicit()) {
            continue;                         // implicit molecules are added/deleted after move
        } else if (molecules[molid].atomic) { // reactant is an atomic group
            auto mollist = spc.findMolecules(molid, Space::Selection::ALL);
            assert(range_size(mollist) == 1);
            auto target = spc.findMolecules(molid, Space::Selection::ALL).begin();
            auto other_target = other_spc->findMolecules(molid, Space::Selection::ALL).begin();
            auto change_data = contractAtomicGroup(*target, *other_target, N_delete);
            if (not change_data.relative_atom_indices.empty()) {
                change.groups.push_back(change_data);
            } else {
                assert(false); // we should never reach here
            }
        } else { // molecular reactant (non-atomic)
            auto selection =
                (reaction->only_neutral_molecules) ? Space::Selection::ACTIVE_NEUTRAL : Space::Selection::ACTIVE;
            auto active = spc.findMolecules(molid, selection);
            std::vector<std::reference_wrapper<Space::GroupType>> molecules_to_deactivate;
            std::sample(active.begin(), active.end(), std::back_inserter(molecules_to_deactivate), N_delete,
                        slump.engine); // pick random molecules to delete
            if (molecules_to_deactivate.size() == N_delete) {
                for (auto &target : molecules_to_deactivate) {
                    if (auto change_data = deactivateMolecularGroup(target);
                        not change_data.relative_atom_indices.empty()) {
                        change.groups.push_back(change_data); // add to list of moved groups
                    } else {
                        assert(false); // we should never reach here
                    }
                }
            } else {
                if (Faunus::molecules[molid].activity > 0) { // if molecule has an activity, issue a warning
                    faunus_logger->warn("all {} molecules have been deleted; increase system volume?",
                                        Faunus::molecules[molid].name);
                }
                throw SpeciationMoveException();
            }
        }
    }
}

/**
 * Checks if there is enough implicit molecules to carry out the reaction
 */
bool SpeciationMove::enoughImplicitMolecules() const {
    auto isExhausted = [&](auto &i) {                             // check if species has enough implicit
        auto [molid, nu] = i;                                     // molid and stoichiometric coefficient
        if (Faunus::molecules[molid].isImplicit()) {              // matter to perform process
            assert(spc.getImplicitReservoir().count(molid) == 1); // must be registered in space!
            return nu > spc.getImplicitReservoir()[molid];        // not enough material?
        } else {
            return false;                                         // non-implicit molecules cannot be exhausted
        }
    };

    [[maybe_unused]] auto [atomic_reactants, molecular_reactants] = reaction->getReactants();

    return std::find_if(molecular_reactants.begin(), molecular_reactants.end(), isExhausted) ==
           molecular_reactants.end(); // no molecule exhausted?
}

void SpeciationMove::_move(Change &change) {
    assert(other_spc != nullptr);        // knowledge of other space should be provided by now
    try {
        if (Faunus::reactions.empty()) {     // global list of reactions
            throw SpeciationMoveException();
        }
        reaction = slump.sample(Faunus::reactions.begin(), Faunus::reactions.end());    // random reaction
        auto direction = static_cast<ReactionData::Direction>((char)slump.range(0, 1)); // random direction
        reaction->setDirection(direction);
        bond_energy = 0.0;
        atomicSwap(change);
        deactivateAllReactants(change);
        activateAllProducts(change);
        if (!change.empty()) {
            change.matter_change = true; // Attempting to change the number of atoms / molecules
            std::sort(change.groups.begin(), change.groups.end()); // change groups *must* be sorted!
            updateGroupMassCenters(change);
        }
    } catch (SpeciationMoveException&) { change.clear(); }
}

/**
 * Speciation move may induce a change in molecular mass centers
 */
void SpeciationMove::updateGroupMassCenters(const Change& change) const {
    for (const auto& change_data : change.groups) {
        if (change_data.dNatomic || change_data.dNswap) {
            auto& group = spc.groups.at(change_data.group_index);
            if (group.massCenter()) { // update only if group has a well-defined mass center
                group.updateMassCenter(spc.geometry.getBoundaryFunc(), group.mass_center);
            }
        }
    }
}

double SpeciationMove::bias([[maybe_unused]] Change& change, [[maybe_unused]] double old_energy,
                            [[maybe_unused]] double new_energy) {
    // The acceptance/rejection of the move is affected by the equilibrium constant
    // but unaffected by the change in bonded energy
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

SpeciationMove::SpeciationMove(Space& spc, std::string_view name, std::string_view cite) : MoveBase(spc, name, cite) {}

SpeciationMove::SpeciationMove(Space &spc) : SpeciationMove(spc, "rcmc", "doi:10/fqcpg3") {}

void SpeciationMove::_from_json(const json &) {}

} // end of namespace Move
} // end of namespace Faunus
