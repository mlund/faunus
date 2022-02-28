#include "bonds.h"
#include "speciation.h"
#include "aux/iteratorsupport.h"
#include <algorithm>
#include <range/v3/algorithm/all_of.hpp>
#include <range/v3/view/take.hpp>
#include <range/v3/view/sample.hpp>

namespace Faunus::Move {

/**
 * @brief Internal exception to carry errors preventing the speciation move.
 */
struct SpeciationMoveException : public std::exception {};

// ----------------------------------

void ReactionDirectionRatio::AcceptanceData::update(ReactionData::Direction direction, bool accept) {
    if (direction == ReactionData::Direction::RIGHT) {
        right += static_cast<double>(accept);
    } else {
        left += static_cast<double>(accept);
    }
}
void ReactionDirectionRatio::to_json(json& j) const {
    auto& _j = j["reactions"] = json::object();
    for (const auto& [reaction, data] : acceptance) {
        _j[reaction->getReactionString()] = {{"attempts", data.left.size() + data.right.size()},
                                             {"acceptance -->", data.right.avg()},
                                             {"acceptance <--", data.left.avg()}};
    }
}

ReactionDirectionRatio::AcceptanceData&
ReactionDirectionRatio::operator[](ReactionDirectionRatio::reaction_iterator iter) {
    return acceptance[iter];
}

// ----------------------------------

void SpeciationMove::_to_json(json& j) const {
    direction_ratio.to_json(j);
    for (auto [molid, size] : average_reservoir_size) {
        j["implicit_reservoir"][molecules.at(molid).name] = size.avg();
    }
}

/**
 * Convert from one atom type to another in any group (atomic/molecular).
 * The reaction requires that `swap` is true and there must be *exactly*
 * one atomic reactant and one atomic product in the reaction.
 *
 * The function checks if there are sufficient atomic and molecular
 * reactants and products to perform the move. If not, an exception is
 * thrown and the system is left untouched.
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
        const auto& [change_data, bias] = atomic_group_bouncer->activate(*mollist.begin(), number_to_insert);
        change.groups.emplace_back(change_data);
        bias_energy += bias;
    }

    for (auto [molid, number_to_insert] : molecular_groups) {
        auto selection =
            (reaction->only_neutral_molecules) ? Space::Selection::INACTIVE_NEUTRAL : Space::Selection::INACTIVE;
        auto inactive = spc.findMolecules(molid, selection) | ranges::views::sample(number_to_insert, slump.engine) |
                        ranges::to<std::vector<std::reference_wrapper<Group>>>;

        assert(inactive.size() == number_to_insert); // should be ok if passed by reaction validator
        ranges::cpp20::for_each(inactive, [&](auto& group) {
            const auto& [change_data, bias] = molecular_group_bouncer->activate(group);
            change.groups.emplace_back(change_data);
            bias_energy += bias;
        });
    }
}

/**
 * @brief Deactivate all atomic and molecular reactants.
 *
 * @throw SpeciationMoveException  when impossible due to the lack of reactant – either implicit or explicit
 */
void SpeciationMove::deactivateAllReactants(Change& change) {
    namespace rv = ranges::cpp20::views;
    auto molecular_reactants = reaction->getReactants().second | rv::filter(ReactionData::not_implicit_group) |
                               rv::filter([](auto& i) { return i.second > 0; });

    for (auto [molid, N_delete] : molecular_reactants) { // Delete
        if (molecules.at(molid).isAtomic()) {            // reactant is an atomic group
            auto mollist = spc.findMolecules(molid, Space::Selection::ALL);
            assert(range_size(mollist) == 1);
            auto group = mollist.begin();
            auto const& [change_data, bias] = atomic_group_bouncer->deactivate(*group, N_delete);
            assert(!change_data.relative_atom_indices.empty());
            change.groups.push_back(change_data);
            bias_energy += bias;
        } else { // molecular reactant (non-atomic)
            const auto selection =
                (reaction->only_neutral_molecules) ? Space::Selection::ACTIVE_NEUTRAL : Space::Selection::ACTIVE;
            auto active = spc.findMolecules(molid, selection) | ranges::views::sample(N_delete, slump.engine) |
                          ranges::to<std::vector<std::reference_wrapper<Group>>>;
            assert(active.size() == N_delete);
            std::for_each(active.begin(), active.end(), [&](auto& group) {
                const auto& [change_data, bias] = molecular_group_bouncer->deactivate(group);
                change.groups.emplace_back(change_data);
                bias_energy += bias;
            });
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

void SpeciationMove::_move(Change& change) {
    if (Faunus::reactions.empty()) {
        return;
    }
    bias_energy = 0.0;
    try {
        setReaction();
        if (reaction_validator.isPossible(*reaction)) {
            atomicSwap(change);
            deactivateAllReactants(change);
            activateAllProducts(change);
            if (change) {
                change.matter_change = true;
                std::sort(change.groups.begin(), change.groups.end()); // change groups *must* be sorted!
                updateGroupMassCenters(change);
            }
        }
    } catch (SpeciationMoveException&) { change.clear(); }
}

void SpeciationMove::setReaction() {
    reaction = slump.sample(reactions.begin(), reactions.end());
    reaction->setRandomDirection(slump);
}

/**
 * Speciation move may induce a change in molecular mass centers.
 * Curently updated for the following moves:
 *
 * - Swap moves (as the swapped atom may have a different mass)
 */
void SpeciationMove::updateGroupMassCenters(const Change& change) const {
    namespace rv = ranges::cpp20::views;

    auto atomic_or_swap = [](const Change::GroupChange& c) { return c.dNatomic || c.dNswap; };
    auto to_group = [&](const Change::GroupChange& c) -> Group& { return spc.groups.at(c.group_index); };
    auto has_mass_center = [](Group& group) { return group.massCenter().has_value(); };

    auto groups = change.groups | rv::filter(atomic_or_swap) | rv::transform(to_group) | rv::filter(has_mass_center);
    ranges::cpp20::for_each(groups, [&](Group& group) {
        group.updateMassCenter(spc.geometry.getBoundaryFunc(), group.massCenter().value());
    });
}

/**
 * The acceptance/rejection of the move is affected by the equilibrium constant,
 * but unaffected by the change in internal bond energy
 */
double SpeciationMove::bias([[maybe_unused]] Change& change, [[maybe_unused]] double old_energy,
                            [[maybe_unused]] double new_energy) {
    return -reaction->freeEnergy() + bias_energy;
}

void SpeciationMove::_accept(Change&) {
    direction_ratio[reaction].update(reaction->getDirection(), true);
    const auto& molecular_products = reaction->getProducts().second;
    const auto& molecular_reactants = reaction->getReactants().second;

    // adjust amount of implicit matter
    for (const auto [molid, nu] : molecular_reactants) {
        if (Faunus::molecules.at(molid).isImplicit()) {
            spc.getImplicitReservoir()[molid] -= nu;
            average_reservoir_size[molid] += spc.getImplicitReservoir().at(molid);
        }
    }
    for (const auto [molid, nu] : molecular_products) {
        if (Faunus::molecules.at(molid).isImplicit()) {
            spc.getImplicitReservoir()[molid] += nu;
            average_reservoir_size[molid] += spc.getImplicitReservoir().at(molid);
        }
    }
}

void SpeciationMove::_reject([[maybe_unused]] Change& change) {
    direction_ratio[reaction].update(reaction->getDirection(), false);

    const auto& molecular_products = reaction->getProducts().second;
    const auto& molecular_reactants = reaction->getReactants().second;

    // average number of implicit molecules
    for (auto [molid, nu] : molecular_reactants) {
        if (Faunus::molecules[molid].isImplicit()) {
            average_reservoir_size[molid] += spc.getImplicitReservoir().at(molid);
        }
    }
    for (auto [molid, nu] : molecular_products) {
        if (Faunus::molecules[molid].isImplicit()) {
            average_reservoir_size[molid] += spc.getImplicitReservoir().at(molid);
        }
    }
}

SpeciationMove::SpeciationMove(Space& spc, Space& old_spc, std::string_view name, std::string_view cite)
    : MoveBase(spc, name, cite)
    , reaction_validator(spc) {
    molecular_group_bouncer = std::make_unique<MolecularGroupDeActivator>(spc, slump, true);
    atomic_group_bouncer = std::make_unique<AtomicGroupDeActivator>(spc, old_spc, slump);
}

SpeciationMove::SpeciationMove(Space& spc, Space& old_spc)
    : SpeciationMove(spc, old_spc, "rcmc", "doi:10/fqcpg3") {}

void SpeciationMove::_from_json(const json&) {}

// -----------------------------------------

ReactionValidator::ReactionValidator(const Space& spc)
    : spc(spc) {}

bool ReactionValidator::isPossible(const ReactionData& reaction) const {
    return canReduceImplicitGroups(reaction) && canSwapAtoms(reaction) && canReduceMolecularGroups(reaction) &&
           canProduceMolecularGroups(reaction) && canReduceAtomicGrups(reaction) && canProduceAtomicGroups(reaction);
}

bool ReactionValidator::canReduceImplicitGroups(const ReactionData& reaction) const {
    auto has_enough = [&](const auto& key_value) {
        const auto& [molid, number_to_delete] = key_value;
        return spc.getImplicitReservoir().at(molid) >= number_to_delete;
    };
    auto implicit_reactants =
        reaction.getReactants().second | ranges::cpp20::views::filter(ReactionData::is_implicit_group);
    return ranges::cpp20::all_of(implicit_reactants, has_enough);
}

bool ReactionValidator::canSwapAtoms(const ReactionData& reaction) const {
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

bool ReactionValidator::canReduceMolecularGroups(const ReactionData& reaction) const {
    namespace rv = ranges::cpp20::views;
    auto can_reduce = [&](auto key_value) { // molecular reactant (non-atomic)
        const auto [molid, number_to_delete] = key_value;
        if (number_to_delete > 0) {
            const auto selection =
                reaction.only_neutral_molecules ? Space::Selection::ACTIVE_NEUTRAL : Space::Selection::ACTIVE;
            auto active = spc.findMolecules(molid, selection) | rv::take(number_to_delete);
            if (range_size(active) != number_to_delete) {
                if (Faunus::molecules.at(molid).activity > 0.0) {
                    faunus_logger->warn("all grand canonical {} molecules have been deleted; increase system volume?",
                                        Faunus::molecules.at(molid).name);
                }
                return false;
            }
        }
        return true;
    };

    auto molecular_groups = reaction.getReactants().second | rv::filter(ReactionData::not_implicit_group) |
                            rv::filter(ReactionData::is_molecular_group);

    return ranges::cpp20::all_of(molecular_groups, can_reduce);
}

bool ReactionValidator::canReduceAtomicGrups(const ReactionData& reaction) const {
    namespace rv = ranges::cpp20::views;
    auto can_reduce = [&](auto key_value) {
        const auto [molid, number_to_delete] = key_value;
        if (number_to_delete > 0) {
            auto mollist = spc.findMolecules(molid, Space::Selection::ALL);
            if (range_size(mollist) != 1) {
                return false;
            }
            auto target = mollist.begin();
            if ((int)target->size() < number_to_delete) {
                faunus_logger->warn("atomic group {} is depleted; increase simulation volume?",
                                    Faunus::molecules.at(molid).name);
                return false;
            }
        }
        return true;
    };

    auto atomic_groups = reaction.getReactants().second | rv::filter(ReactionData::not_implicit_group) |
                         rv::filter(ReactionData::is_atomic_group);
    return ranges::cpp20::all_of(atomic_groups, can_reduce);
}

/** enough molecular products? */
bool ReactionValidator::canProduceMolecularGroups(const ReactionData& reaction) const {
    namespace rv = ranges::cpp20::views;
    auto can_create = [&](auto key_value) {
        const auto [molid, number_to_insert] = key_value;
        if (number_to_insert > 0) {
            const auto selection =
                (reaction.only_neutral_molecules) ? Space::Selection::INACTIVE_NEUTRAL : Space::Selection::INACTIVE;
            auto inactive = spc.findMolecules(molid, selection) | rv::take(number_to_insert);
            if (range_size(inactive) != number_to_insert) {
                faunus_logger->warn("maximum number of {} molecules reached; increase capacity?",
                                    Faunus::molecules.at(molid).name);
                return false;
            }
        }
        return true;
    };

    auto molecular_groups = reaction.getProducts().second | rv::filter(ReactionData::not_implicit_group) |
                            rv::filter(ReactionData::is_molecular_group);

    return ranges::cpp20::all_of(molecular_groups, can_create);
}

bool ReactionValidator::canProduceAtomicGroups(const ReactionData& reaction) const {
    namespace rv = ranges::cpp20::views;
    auto can_expand = [&](auto key_value) {
        const auto [molid, number_to_insert] = key_value;
        auto mollist = spc.findMolecules(molid, Space::Selection::ALL);
        if (number_to_insert > 0 && range_size(mollist) > 0) {
            if (mollist.begin()->size() + number_to_insert > mollist.begin()->capacity()) {
                faunus_logger->warn("atomic molecule {} is full; increase capacity?", Faunus::molecules.at(molid).name);
                return false;
            }
        }
        return true;
    };

    auto atomic_groups = reaction.getProducts().second | rv::filter(ReactionData::not_implicit_group) |
                         rv::filter(ReactionData::is_atomic_group);

    return ranges::cpp20::all_of(atomic_groups, can_expand);
}

// ----------------------------------------------

/**
 * Randomly assign a new mass center and random orientation
 */
void MolecularGroupDeActivator::setPositionAndOrientation(Group& group) const {
    auto& geometry = spc.geometry;

    // translate to random position within simulation cell
    Point new_mass_center;
    geometry.randompos(new_mass_center, random); // place COM randomly in simulation box
    Point displacement = geometry.vdist(new_mass_center, group.massCenter()->get());
    group.translate(displacement, geometry.getBoundaryFunc());

    // generate random orientation
    const auto rotation_angle = 2.0 * pc::pi * (random() - 0.5); // -pi to pi
    const auto random_unit_vector = randomUnitVector(random);
    Eigen::Quaterniond quaternion(Eigen::AngleAxisd(rotation_angle, random_unit_vector));
    group.rotate(quaternion, geometry.getBoundaryFunc());
}

MolecularGroupDeActivator::ChangeAndBias
MolecularGroupDeActivator::activate(Group& group, GroupDeActivator::OptionalInt num_particles) {
    assert(group.isMolecular());                                      // must be a molecule group
    assert(group.empty());                                            // must be inactive
    group.activate(group.inactive().begin(), group.inactive().end()); // activate all particles
    assert(not group.empty());

    if (num_particles && num_particles.value() != group.capacity()) {
        throw std::runtime_error("only full activation allowed");
    }

    setPositionAndOrientation(group);

    assert(spc.geometry.sqdist(group.mass_center,
                               Geometry::massCenter(group.begin(), group.end(), spc.geometry.getBoundaryFunc(),
                                                    -group.mass_center)) < 1e-9);

    Change::GroupChange group_change;                    // describes the changed - used for energy evaluation
    group_change.group_index = spc.getGroupIndex(group); // index* of moved group
    group_change.all = true;                             // all atoms in group were moved
    group_change.internal = true;
    group_change.relative_atom_indices.resize(group.capacity()); // list of changed atom index
    std::iota(group_change.relative_atom_indices.begin(), group_change.relative_atom_indices.end(), 0);
    return {group_change, -getBondEnergy(group)};
}

MolecularGroupDeActivator::ChangeAndBias
MolecularGroupDeActivator::deactivate(Group& group, GroupDeActivator::OptionalInt num_particles) {
    assert(group.isMolecular());
    assert(!group.empty());
    assert(group.size() == group.capacity());

    if (num_particles && num_particles.value() != group.capacity()) {
        throw std::runtime_error("only full activation allowed");
    }
    group.unwrap(spc.geometry.getDistanceFunc()); // when in storage, remove PBC
    group.deactivate(group.begin(), group.end()); // deactivate whole group
    assert(group.empty());

    Change::GroupChange change_data; // describes the change
    change_data.internal = true;
    change_data.group_index = spc.getGroupIndex(group);
    change_data.all = true;                                     // all atoms in group were moved
    change_data.relative_atom_indices.resize(group.capacity()); // list of changed atom index
    std::iota(change_data.relative_atom_indices.begin(), change_data.relative_atom_indices.end(), 0);

    return {change_data, getBondEnergy(group)};
}

MolecularGroupDeActivator::MolecularGroupDeActivator(Space& spc, Random& random, bool apply_bond_bias)
    : spc(spc)
    , random(random)
    , apply_bond_bias(apply_bond_bias) {}

double MolecularGroupDeActivator::getBondEnergy(const Group& group) const {
    double energy = 0.0;
    if (apply_bond_bias) {
        auto bonds = group.traits().bonds | ranges::views::transform(&Potential::BondData::clone);
        ranges::cpp20::for_each(bonds, [&](auto bond) {
            bond->shiftIndices(spc.getFirstParticleIndex(group));
            bond->setEnergyFunction(spc.particles);
            energy += bond->energyFunc(spc.geometry.getDistanceFunc());
        });
    }
    return energy;
}

// ------------------------------------------

AtomicGroupDeActivator::AtomicGroupDeActivator(Space& spc, Space& old_spc, Random& random)
    : spc(spc)
    , old_spc(old_spc)
    , slump(random) {}

GroupDeActivator::ChangeAndBias AtomicGroupDeActivator::activate(Group& group,
                                                                 GroupDeActivator::OptionalInt number_to_insert) {
    if (!group.isAtomic() || !number_to_insert || number_to_insert.value() + group.size() > group.capacity()) {
        throw std::runtime_error("atomic group expansion bug");
    }
    Change::GroupChange change_data;
    change_data.group_index = spc.getGroupIndex(group);
    change_data.internal = true;
    change_data.dNatomic = true;
    for (int i = 0; i < number_to_insert.value(); i++) {
        group.activate(group.end(), group.end() + 1); // activate one particle
        auto last_atom = group.end() - 1;
        spc.geometry.randompos(last_atom->pos, slump);  // give it a random position
        spc.geometry.getBoundaryFunc()(last_atom->pos); // apply PBC if needed
        change_data.relative_atom_indices.push_back(std::distance(group.begin(), last_atom)); // index relative to group
    }
    change_data.sort();
    return {change_data, 0.0};
}

/**
 * Reduce an atomic group by `number_to_delete` particles. The deleted particles are
 * picked by random; moved the the end of the group; then deactivated. In order for
 * the Hamiltonian to pick up the energy change, particles in the reference Space
 * (`old_group`) are swapped to the same index, albeit not deactivated.
 *
 * @warning Directly modifying the groups in spc and old_spc might interfere with
 *          a future neighbour list implementation.
 */
GroupDeActivator::ChangeAndBias AtomicGroupDeActivator::deactivate(Group& group,
                                                                   GroupDeActivator::OptionalInt number_to_delete) {
    if (!group.isAtomic() || !number_to_delete || number_to_delete.value() > group.size()) {
        throw std::runtime_error("atomic group expansion bug");
    }
    Change::GroupChange change_data; // describes what has changed
    change_data.group_index = spc.getGroupIndex(group);
    change_data.internal = true;
    change_data.dNatomic = true;

    const auto& old_group = old_spc.groups.at(change_data.group_index);

    for (int i = 0; i < number_to_delete.value(); i++) {
        auto particle_to_delete = slump.sample(group.begin(), group.end());
        auto last_particle = group.end() - 1;
        const auto dist = std::distance(particle_to_delete, group.end());

        if (std::distance(particle_to_delete, last_particle) > 1) {
            std::iter_swap(particle_to_delete, last_particle);                     // Shuffle back to end in trial ...
            std::iter_swap(old_group.end() - dist - i, old_group.end() - (1 + i)); // ...and in old group
        }
        const auto deactivated_particle_index = std::distance(group.begin(), last_particle);
        change_data.relative_atom_indices.push_back(deactivated_particle_index);
        group.deactivate(last_particle, group.end()); // deactivate one particle at the time
    }
    change_data.sort();
    return {change_data, 0.0};
}

} // namespace Faunus::Move
