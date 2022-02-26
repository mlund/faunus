#include "bonds.h"
#include "speciation.h"
#include "aux/iteratorsupport.h"
#include <algorithm>
#include <range/v3/algorithm/all_of.hpp>
#include <range/v3/algorithm/any_of.hpp>
#include <range/v3/view/take.hpp>
#include <range/v3/view/sample.hpp>

namespace Faunus {
namespace Move {

/**
 * @brief Internal exception to carry errors preventing the speciation move.
 */
struct SpeciationMoveException : public std::exception {};

void SpeciationMove::_to_json(json& j) const {
    json& _j = j["reactions"];
    _j = json::object();
    for (auto [reaction, data] : acceptance) {
        _j[reaction->getReactionString()] = {{"attempts", data.left.size() + data.right.size()},
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
    if (!reaction->containsAtomicSwap()) {
        return;
    }
    const auto& atomic_products = reaction->getProducts().first;
    const auto& atomic_reactants = reaction->getReactants().first;

    assert(atomic_products.size() == 1 and atomic_reactants.size() == 1);

    auto atomlist = spc.findAtoms(atomic_reactants.begin()->first);        // search all active molecules
    auto random_particle = slump.sample(atomlist.begin(), atomlist.end()); // target particle to swap
    auto group = spc.findGroupContaining(*random_particle);                // find enclosing group

    Change::GroupChange d; // describe what has change - used for energy cal.
    d.relative_atom_indices.push_back(
        Faunus::distance(group->begin(), random_particle));      // Index of particle rel. to group
    d.group_index = Faunus::distance(spc.groups.begin(), group); // index of particle in group (starting from zero)
    d.internal = true;
    d.dNswap = true;
    change.groups.push_back(d); // Add to list of moved groups

    int atomid = atomic_products.begin()->first; // atomid of new atom type
    Particle p = Faunus::atoms[atomid];          // temporary particle of new type
    p.pos = random_particle->pos;                // get position from old particle
    // todo: extended properties, dipole etc?
    assert(!p.hasExtension() && "extended properties not yet implemented");
    *random_particle = p; // copy new particle onto old particle
    assert(random_particle->id == atomid);
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
    for (auto& bond : Faunus::molecules.at(target.id).bonds) {
        auto bond_clone = bond->clone();
        bond_clone->shiftIndices(std::distance(spc.particles.begin(), target.begin()));
        bond_clone->setEnergyFunction(spc.particles);
        bond_energy += bond_clone->energyFunc(spc.geometry.getDistanceFunc());
    }

    target.deactivate(target.begin(), target.end()); // deactivate whole group
    assert(target.empty());

    Change::GroupChange change_data; // describes the change
    change_data.internal = true;
    change_data.group_index = &target - &spc.groups.front();     // index of moved group
    change_data.all = true;                                      // all atoms in group were moved
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
        change_data.group_index = spc.getGroupIndex(target);
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
    assert(not target.isAtomic());                                       // must be a molecule group
    assert(target.empty());                                              // must be inactive
    target.activate(target.inactive().begin(), target.inactive().end()); // activate all particles
    assert(not target.empty());

    Point cm = target.mass_center;
    spc.geometry.randompos(cm, slump);                    // generate random position
    target.translate(cm, spc.geometry.getBoundaryFunc()); // assign random position to mass-center
    Point u = randomUnitVector(slump);                    // random unit vector
    Eigen::Quaterniond Q(Eigen::AngleAxisd(2 * pc::pi * (slump() - 0.5), u));
    target.rotate(Q, spc.geometry.getBoundaryFunc()); // assign random orientation

    assert(spc.geometry.sqdist(target.mass_center,
                               Geometry::massCenter(target.begin(), target.end(), spc.geometry.getBoundaryFunc(),
                                                    -target.mass_center)) < 1e-9);

    // Store internal bond energy of activated molecule
    for (auto& bond : Faunus::molecules[target.id].bonds) {
        auto bondclone = bond->clone();
        bondclone->shiftIndices(std::distance(spc.particles.begin(), target.begin()));
        bondclone->setEnergyFunction(spc.particles);
        bond_energy -= bondclone->energyFunc(spc.geometry.getDistanceFunc());
    }

    Change::GroupChange d;                     // describes the changed - used for energy evaluation
    d.group_index = spc.getGroupIndex(target); // index* of moved group
    d.all = true;                              // all atoms in group were moved
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
void SpeciationMove::activateAllProducts(Change& change) {
    namespace rv = ranges::cpp20::views;
    auto atomic_groups = reaction->getProducts().second | rv::filter(ReactionData::not_implicit_group) |
                         rv::filter(ReactionData::is_atomic_group);
    auto molecular_groups = reaction->getProducts().second | rv::filter(ReactionData::not_implicit_group) |
                         rv::filter(ReactionData::is_molecular_group);

    for (auto [molid, number_to_insert] : atomic_groups) {
        auto mollist = spc.findMolecules(molid, Space::Selection::ALL);
        change.groups.emplace_back(expandAtomicGroup(*mollist.begin(), number_to_insert));
    }

    for (auto [molid, number_to_insert] : molecular_groups) {
        auto selection =
            (reaction->only_neutral_molecules) ? Space::Selection::INACTIVE_NEUTRAL : Space::Selection::INACTIVE;
        auto inactive = spc.findMolecules(molid, selection) | ranges::views::sample(number_to_insert, slump.engine) |
                        ranges::to<std::vector<std::reference_wrapper<Group>>>;

        assert(inactive.size() == number_to_insert); // should be ok if passed by reaction validator
        ranges::cpp20::for_each(inactive,
                                [&](auto& group) { change.groups.emplace_back(activateMolecularGroup(group)); });
    }
}

/**
 * @brief Deactivate all atomic and molecular reactants.
 *
 * @throw SpeciationMoveException  when impossible due to the lack of reactant – either implicit or explicit
 */
void SpeciationMove::deactivateAllReactants(Change& change) {
    namespace rv = ranges::cpp20::views;
    if (not enoughImplicitMolecules()) {
        throw SpeciationMoveException();
    }

    auto molecular_reactants = reaction->getReactants().second | rv::filter(ReactionData::not_implicit_group) |
                               rv::filter([](auto& i) { return i.second > 0; });

    // perform actual deactivation
    for (auto [molid, N_delete] : molecular_reactants) { // Delete
        if (molecules[molid].isAtomic()) {               // reactant is an atomic group
            auto mollist = spc.findMolecules(molid, Space::Selection::ALL);
            assert(range_size(mollist) == 1);
            auto target = mollist.begin();
            auto other_target = other_spc->findMolecules(molid, Space::Selection::ALL).begin();
            auto change_data = contractAtomicGroup(*target, *other_target, N_delete);
            assert(!change_data.relative_atom_indices.empty());
            change.groups.push_back(change_data);
        } else { // molecular reactant (non-atomic)
            const auto selection =
                (reaction->only_neutral_molecules) ? Space::Selection::ACTIVE_NEUTRAL : Space::Selection::ACTIVE;
            auto active = spc.findMolecules(molid, selection) | ranges::views::sample(N_delete, slump.engine) |
                          ranges::to<std::vector<std::reference_wrapper<Group>>>;
            assert(active.size() == N_delete);
            std::for_each(active.begin(), active.end(),
                          [&](auto& group) { change.groups.emplace_back(deactivateMolecularGroup(group)); });
        }
    }
}

TEST_CASE("[Faunus] Speciation - Ranges::sample") {
    std::vector<int> vec = {1, 2, 3, 4};
    auto take_nothing = vec | ranges::views::sample(0);
    auto take_less = vec | ranges::views::sample(2);
    auto take_all = vec | ranges::views::sample(4);
    auto take_too_much = vec | ranges::views::sample(10);

    CHECK(range_size(take_nothing) == 0);
    CHECK(range_size(take_less) == 2);
    CHECK(range_size(take_all) == 4);
    CHECK(range_size(take_too_much) == 4);
}

/**
 * Checks if there is enough implicit molecules to carry out the reaction
 */
bool SpeciationMove::enoughImplicitMolecules() const {
    auto isExhausted = [&](auto& i) {                             // check if species has enough implicit
        auto [molid, nu] = i;                                     // molid and stoichiometric coefficient
        if (Faunus::molecules[molid].isImplicit()) {              // matter to perform process
            assert(spc.getImplicitReservoir().count(molid) == 1); // must be registered in space!
            return nu > spc.getImplicitReservoir()[molid];        // not enough material?
        } else {
            return false; // non-implicit molecules cannot be exhausted
        }
    };

    [[maybe_unused]] auto [atomic_reactants, molecular_reactants] = reaction->getReactants();

    return std::find_if(molecular_reactants.begin(), molecular_reactants.end(), isExhausted) ==
           molecular_reactants.end(); // no molecule exhausted?
}

void SpeciationMove::_move(Change& change) {
    assert(other_spc != nullptr); // knowledge of other space should be provided by now
    try {
        if (Faunus::reactions.empty()) { // global list of reactions
            throw SpeciationMoveException();
        }
        reaction = slump.sample(Faunus::reactions.begin(), Faunus::reactions.end());    // random reaction
        auto direction = static_cast<ReactionData::Direction>((char)slump.range(0, 1)); // random direction
        reaction->setDirection(direction);
        bond_energy = 0.0;
        if (not validate_reaction.reactionIsPossible(*reaction)) {
            throw SpeciationMoveException();
        }
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
    return -reaction->freeEnergy() + bond_energy;
}

void SpeciationMove::_accept(Change&) {
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

void SpeciationMove::_reject(Change&) {
    acceptance[reaction].update(reaction->getDirection(), false);

    const auto& molecular_products = reaction->getProducts().second;
    const auto& molecular_reactants = reaction->getReactants().second;

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

SpeciationMove::SpeciationMove(Space& spc, std::string_view name, std::string_view cite)
    : MoveBase(spc, name, cite)
    , validate_reaction(spc) {}

SpeciationMove::SpeciationMove(Space& spc)
    : SpeciationMove(spc, "rcmc", "doi:10/fqcpg3") {}

void SpeciationMove::_from_json(const json&) {}

void SpeciationMove::AcceptanceData::update(const ReactionData::Direction direction, const bool accept) {
    if (direction == ReactionData::Direction::RIGHT) {
        right += static_cast<double>(accept);
    } else {
        left += static_cast<double>(accept);
    }
}

// -----------------------------------------

ValidateReaction::ValidateReaction(const Space& spc)
    : spc(spc) {}

bool ValidateReaction::reactionIsPossible(const ReactionData& reaction) const {
    return enoughImplicitMolecules(reaction) && canSwapAtoms(reaction) && canReduceMolecularGrups(reaction) &&
           canProduceMolecularGroups(reaction) && canReduceAtomicGrups(reaction) && canProduceAtomicGroups(reaction);
}

bool ValidateReaction::enoughImplicitMolecules(const ReactionData& reaction) const {
    auto isExhausted = [&](auto& _pair) {                     // check if species has enough implicit
        const auto [molid, nu] = _pair;                       // molid and stoichiometric coefficient
        if (Faunus::molecules.at(molid).isImplicit()) {       // matter to perform process
            return nu > spc.getImplicitReservoir().at(molid); // not enough material?
        }
        return false; // non-implicit molecules cannot be exhausted
    };

    [[maybe_unused]] auto [atomic_reactants, molecular_reactants] = reaction.getReactants();

    return std::find_if(molecular_reactants.begin(), molecular_reactants.end(), isExhausted) ==
           molecular_reactants.end(); // no molecule exhausted?
}

bool ValidateReaction::canSwapAtoms(const ReactionData& reaction) const {
    namespace rv = ranges::cpp20::views;
    if (reaction.containsAtomicSwap()) {
        auto reactive_atoms = reaction.getReactants().first | rv::filter(ReactionData::not_implicit_atom);
        for (const auto [atomid, number_to_swap] : reactive_atoms) {
            auto particles = spc.findAtoms(atomid) | rv::take(number_to_swap);
            if (range_size(particles) != number_to_swap) {
                return false;
            }
        }
    }
    return true;
}

bool ValidateReaction::canReduceMolecularGrups(const ReactionData& reaction) const {
    namespace rv = ranges::cpp20::views;
    auto molecular_groups = reaction.getReactants().second | rv::filter(ReactionData::not_implicit_group) |
                            rv::filter(ReactionData::is_molecular_group);

    auto cannot_reduce_molecular_group = [&](const auto& _pair) { // molecular reactant (non-atomic)
        const auto [molid, number_to_delete] = _pair;
        if (number_to_delete > 0) {
            const auto selection =
                reaction.only_neutral_molecules ? Space::Selection::ACTIVE_NEUTRAL : Space::Selection::ACTIVE;
            auto active = spc.findMolecules(molid, selection) | rv::take(number_to_delete);
            if (range_size(active) != number_to_delete) {
                if (Faunus::molecules.at(molid).activity > 0.0) {
                    faunus_logger->warn("all grand canonical {} molecules have been deleted; increase system volume?",
                                        Faunus::molecules.at(molid).name);
                }
                return true;
            }
        }
        return false;
    };
    return !ranges::cpp20::any_of(molecular_groups, cannot_reduce_molecular_group);
}

bool ValidateReaction::canReduceAtomicGrups(const ReactionData& reaction) const {
    namespace rv = ranges::cpp20::views;
    auto atomic_groups = reaction.getReactants().second | rv::filter(ReactionData::not_implicit_group) |
                         rv::filter(ReactionData::is_atomic_group);

    auto cannot_reduce_atomic_group = [&](const auto& _pair) {
        const auto [molid, number_to_delete] = _pair;
        if (number_to_delete > 0) {
            auto mollist = spc.findMolecules(molid, Space::Selection::ALL);
            if (range_size(mollist) != 1) {
                return true;
            }
            auto target = mollist.begin();
            if ((int)target->size() < number_to_delete) {
                faunus_logger->warn("atomic group {} is depleted; increase simulation volume?",
                                    Faunus::molecules.at(molid).name);
                return true;
            }
        }
        return false;
    };
    return !ranges::cpp20::any_of(atomic_groups, cannot_reduce_atomic_group);
}

/** enough molecular products? */
bool ValidateReaction::canProduceMolecularGroups(const ReactionData& reaction) const {
    namespace rv = ranges::cpp20::views;
    auto molecular_groups = reaction.getProducts().second | rv::filter(ReactionData::not_implicit_group) |
                            rv::filter(ReactionData::is_molecular_group);

    auto cannot_create_molecular_group = [&](const auto& _pair) {
        const auto [molid, number_to_insert] = _pair;
        if (number_to_insert > 0) {
            const auto selection =
                (reaction.only_neutral_molecules) ? Space::Selection::INACTIVE_NEUTRAL : Space::Selection::INACTIVE;
            auto inactive = spc.findMolecules(molid, selection) | rv::take(number_to_insert);
            if (range_size(inactive) != number_to_insert) {
                faunus_logger->warn("maximum number of {} molecules reached; increase capacity?",
                                    Faunus::molecules.at(molid).name);
                return true;
            }
        }
        return false;
    };
    return !ranges::cpp20::any_of(molecular_groups, cannot_create_molecular_group);
}

bool ValidateReaction::canProduceAtomicGroups(const ReactionData& reaction) const {
    namespace rv = ranges::cpp20::views;
    auto atomic_groups = reaction.getProducts().second | rv::filter(ReactionData::not_implicit_group) |
                         rv::filter(ReactionData::is_atomic_group);

    auto cannot_create_atomic_group = [&](const auto& _pair) {
        const auto [molid, number_to_insert] = _pair;
        auto mollist = spc.findMolecules(molid, Space::Selection::ALL);
        if (number_to_insert > 0 && range_size(mollist) > 0) {
            if (mollist.begin()->size() + number_to_insert > mollist.begin()->capacity()) {
                faunus_logger->warn("atomic molecule {} is full; increase capacity?", Faunus::molecules.at(molid).name);
                return true;
            }
        }
        return false;
    };
    return !ranges::cpp20::any_of(atomic_groups, cannot_create_atomic_group);
}

} // end of namespace Move
} // end of namespace Faunus
