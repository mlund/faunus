#include "bonds.h"
#include "speciation.h"
#include "random.h"
#include "aux/iteratorsupport.h"
#include <algorithm>
#include <ranges>
#include <range/v3/view/sample.hpp>
#include <range/v3/range/conversion.hpp>
#include <doctest/doctest.h>

namespace Faunus::Speciation {
/**
 * @brief Internal exception to carry errors preventing the speciation move.
 */
struct SpeciationMoveException : public std::exception
{
};

// ----------------------------------

void ReactionDirectionRatio::AcceptanceData::update(ReactionData::Direction direction, bool accept)
{
    if (direction == ReactionData::Direction::RIGHT) {
        right += static_cast<double>(accept);
    }
    else {
        left += static_cast<double>(accept);
    }
}

void ReactionDirectionRatio::to_json(json& j) const
{
    auto& _j = j["reactions"] = json::object();
    for (const auto& [reaction, data] : acceptance) {
        _j[reaction->getReactionString()] = {{"attempts", data.left.size() + data.right.size()},
                                             {"acceptance -->", data.right.avg()},
                                             {"acceptance <--", data.left.avg()}};
    }
}

ReactionDirectionRatio::AcceptanceData&
ReactionDirectionRatio::operator[](ReactionDirectionRatio::reaction_iterator iter)
{
    return acceptance[iter];
}

// -----------------------------------------

ReactionValidator::ReactionValidator(const Space& spc)
    : spc(spc)
{
}

/**
 * @param reaction Reaction to check
 * @return If true, it is guarantied that reactants/products are available and no further checks are
 * needed
 */
bool ReactionValidator::isPossible(const ReactionData& reaction) const
{
    return canReduceImplicitGroups(reaction) && canSwapAtoms(reaction) &&
           canReduceMolecularGroups(reaction) && canProduceMolecularGroups(reaction) &&
           canReduceAtomicGroups(reaction) && canProduceAtomicGroups(reaction);
}

bool ReactionValidator::canReduceImplicitGroups(const ReactionData& reaction) const
{
    auto has_enough = [&](const auto& key_value) {
        const auto& [molid, number_to_delete] = key_value;
        return spc.getImplicitReservoir().at(molid) >= number_to_delete;
    };
    auto implicit_reactants =
        reaction.getReactants().second | std::views::filter(ReactionData::is_implicit_group);
    return std::ranges::all_of(implicit_reactants, has_enough);
}

bool ReactionValidator::canSwapAtoms(const ReactionData& reaction) const
{
    namespace rv = std::views;
    if (reaction.containsAtomicSwap()) {
        auto reactive_atoms =
            reaction.getReactants().first | rv::filter(ReactionData::not_implicit_atom);
        for (const auto& [atomid, number_to_swap] : reactive_atoms) {
            auto particles = spc.findAtoms(atomid) | rv::take(number_to_swap);
            if (range_size(particles) != number_to_swap) {
                return false;
            }
        }
    }
    return true;
}

bool ReactionValidator::canReduceMolecularGroups(const ReactionData& reaction) const
{
    namespace rv = std::views;
    auto selection = reaction.only_neutral_molecules ? Space::Selection::ACTIVE_NEUTRAL
                                                     : Space::Selection::ACTIVE;

    auto can_reduce = [&](auto key_value) { // molecular reactant (non-atomic)
        const auto [molid, number_to_delete] = key_value;
        if (number_to_delete > 0) {
            auto groups = spc.findMolecules(molid, selection) | rv::take(number_to_delete);
            if (range_size(groups) != number_to_delete) {
                if (Faunus::molecules.at(molid).activity > 0.0) {
                    faunus_logger->warn("all grand canonical {} molecules have been deleted; "
                                        "increase system volume?",
                                        Faunus::molecules.at(molid).name);
                }
                return false;
            }
        }
        return true;
    };

    auto molecular_groups = reaction.getReactants().second |
                            rv::filter(ReactionData::not_implicit_group) |
                            rv::filter(ReactionData::is_molecular_group);

    return std::ranges::all_of(molecular_groups, can_reduce);
}

bool ReactionValidator::canReduceAtomicGroups(const ReactionData& reaction) const
{
    namespace rv = std::views;
    auto can_reduce = [&](auto key_value) {
        const auto [molid, number_to_delete] = key_value;
        if (number_to_delete > 0) {
            auto groups = spc.findMolecules(molid, Space::Selection::ALL) |
                          ranges::to<Space::ConstGroupRefVector>;
            if (groups.size() != 1) {
                return false;
            }
            const Group& group = groups.front();
            if ((int)group.size() < number_to_delete) {
                faunus_logger->warn("atomic group {} is depleted; increase simulation volume?",
                                    Faunus::molecules.at(molid).name);
                return false;
            }
        }
        return true;
    };

    auto atomic_groups = reaction.getReactants().second |
                         rv::filter(ReactionData::not_implicit_group) |
                         rv::filter(ReactionData::is_atomic_group);
    return std::ranges::all_of(atomic_groups, can_reduce);
}

bool ReactionValidator::canProduceMolecularGroups(const ReactionData& reaction) const
{
    namespace rv = std::views;
    auto selection = reaction.only_neutral_molecules ? Space::Selection::INACTIVE_NEUTRAL
                                                     : Space::Selection::INACTIVE;

    auto can_create = [&](auto key_value) {
        const auto [molid, number_to_insert] = key_value;
        if (number_to_insert > 0) {
            auto inactive = spc.findMolecules(molid, selection) | rv::take(number_to_insert);
            if (range_size(inactive) != number_to_insert) {
                faunus_logger->warn("maximum number of {} molecules reached; increase capacity?",
                                    Faunus::molecules.at(molid).name);
                return false;
            }
        }
        return true;
    };

    auto molecular_groups = reaction.getProducts().second |
                            rv::filter(ReactionData::not_implicit_group) |
                            rv::filter(ReactionData::is_molecular_group);

    return std::ranges::all_of(molecular_groups, can_create);
}

bool ReactionValidator::canProduceAtomicGroups(const ReactionData& reaction) const
{
    namespace rv = std::views;
    auto can_expand = [&](auto key_value) {
        const auto [molid, number_to_insert] = key_value;
        auto groups = spc.findMolecules(molid, Space::Selection::ALL) |
                      ranges::to<Space::ConstGroupRefVector>;
        if (groups.empty()) {
            return false;
        }
        const Group& group = groups.front();
        if (group.size() + number_to_insert > group.capacity()) {
            faunus_logger->warn("atomic molecule {} is full; increase capacity?",
                                group.traits().name);
            return false;
        }
        return true;
    };

    auto atomic_groups = reaction.getProducts().second |
                         rv::filter(ReactionData::not_implicit_group) |
                         rv::filter(ReactionData::is_atomic_group);

    return std::ranges::all_of(atomic_groups, can_expand);
}

// ----------------------------------------------

/**
 * Randomly assign a new mass center and random orientation
 */
void MolecularGroupDeActivator::setPositionAndOrientation(Group& group) const
{
    auto& geometry = spc.geometry;

    // translate to random position within simulation cell
    Point new_mass_center;
    geometry.randompos(new_mass_center, random); // place COM randomly in simulation box
    const Point displacement = geometry.vdist(new_mass_center, group.massCenter()->get());
    group.translate(displacement, geometry.getBoundaryFunc());

    // generate random orientation
    const auto rotation_angle = 2.0 * pc::pi * (random() - 0.5); // -pi to pi
    const auto random_unit_vector = randomUnitVector(random);
    const Eigen::Quaterniond quaternion(Eigen::AngleAxisd(rotation_angle, random_unit_vector));
    group.rotate(quaternion, geometry.getBoundaryFunc());
}

/**
 * Activates a number of particles in an atomic group. Before calling this, make sure that there's
 * sufficient capacity.
 *
 * @param group Group to affect
 * @param num_particles Number of particles to expand with
 */
MolecularGroupDeActivator::ChangeAndBias
MolecularGroupDeActivator::activate(Group& group, GroupDeActivator::OptionalInt num_particles)
{
    assert(group.isMolecular());                                      // must be a molecule group
    assert(group.empty());                                            // must be inactive
    group.activate(group.inactive().begin(), group.inactive().end()); // activate all particles
    assert(not group.empty());

    if (num_particles && num_particles.value() != group.capacity()) {
        throw std::runtime_error("only full activation allowed");
    }

    setPositionAndOrientation(group);

    assert(
        spc.geometry.sqdist(group.mass_center, Geometry::massCenter(group.begin(), group.end(),
                                                                    spc.geometry.getBoundaryFunc(),
                                                                    -group.mass_center)) < 1e-9);

    Change::GroupChange group_change; // describes the changed - used for energy evaluation
    group_change.group_index = spc.getGroupIndex(group); // index* of moved group
    group_change.all = true;                             // all atoms in group were moved
    group_change.internal = true;
    group_change.relative_atom_indices.resize(group.capacity()); // list of changed atom index
    std::iota(group_change.relative_atom_indices.begin(), group_change.relative_atom_indices.end(),
              0);
    return {group_change, -getBondEnergy(group)};
}

MolecularGroupDeActivator::ChangeAndBias
MolecularGroupDeActivator::deactivate(Group& group, GroupDeActivator::OptionalInt num_particles)
{
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
    std::iota(change_data.relative_atom_indices.begin(), change_data.relative_atom_indices.end(),
              0);

    return {change_data, getBondEnergy(group)};
}

MolecularGroupDeActivator::MolecularGroupDeActivator(Space& spc, Random& random,
                                                     bool apply_bond_bias)
    : spc(spc)
    , random(random)
    , apply_bond_bias(apply_bond_bias)
{
}

double MolecularGroupDeActivator::getBondEnergy(const Group& group) const
{
    double energy = 0.0;
    if (apply_bond_bias) {
        auto bonds = group.traits().bonds | std::views::transform(&pairpotential::BondData::clone);
        std::ranges::for_each(bonds, [&](auto bond) {
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
    , random(random)
{
}

GroupDeActivator::ChangeAndBias
AtomicGroupDeActivator::activate(Group& group, GroupDeActivator::OptionalInt number_to_insert)
{
    if (!group.isAtomic() || !number_to_insert ||
        number_to_insert.value() + group.size() > group.capacity()) {
        throw std::runtime_error("atomic group expansion bug");
    }
    Change::GroupChange change_data;
    change_data.group_index = spc.getGroupIndex(group);
    change_data.internal = true;
    change_data.dNatomic = true;
    for (int i = 0; i < number_to_insert.value(); i++) {
        group.activate(group.end(), group.end() + 1); // activate one particle
        auto last_atom = group.end() - 1;
        spc.geometry.randompos(last_atom->pos, random); // give it a random position
        spc.geometry.getBoundaryFunc()(last_atom->pos); // apply PBC if needed
        change_data.relative_atom_indices.push_back(
            std::distance(group.begin(), last_atom)); // index relative to group
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
GroupDeActivator::ChangeAndBias
AtomicGroupDeActivator::deactivate(Group& group, GroupDeActivator::OptionalInt number_to_delete)
{
    if (!group.isAtomic() || !number_to_delete || number_to_delete.value() > group.size()) {
        throw std::runtime_error("atomic group expansion bug");
    }
    Change::GroupChange change_data; // describes what has changed
    change_data.group_index = spc.getGroupIndex(group);
    change_data.internal = true;
    change_data.dNatomic = true;

    const auto& old_group = old_spc.groups.at(change_data.group_index);

    auto delete_particle = [&](const auto i) {
        auto particle_to_delete = random.sample(group.begin(), group.end());
        auto last_particle = group.end() - 1;
        const auto dist = std::distance(particle_to_delete, group.end());

        if (std::distance(particle_to_delete, last_particle) > 1) {
            std::iter_swap(particle_to_delete, last_particle); // Shuffle back to end in trial ...
            std::iter_swap(old_group.end() - dist - i,
                           old_group.end() - (1 + i)); // ...and in old group
        }
        const auto deactivated_particle_index = std::distance(group.begin(), last_particle);
        change_data.relative_atom_indices.push_back(deactivated_particle_index);
        group.deactivate(last_particle, group.end()); // deactivate one particle at the time
    };

    std::ranges::for_each(std::views::iota(0, number_to_delete), delete_particle);
    change_data.sort();
    return {change_data, 0.0};
}

} // namespace Faunus::Speciation

// ----------------------------------

namespace Faunus::move {

void SpeciationMove::_to_json(json& j) const
{
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
void SpeciationMove::atomicSwap(Change& change)
{
    if (!reaction->containsAtomicSwap()) {
        return;
    }
    const auto& atomic_products = reaction->getProducts().first;
    const auto& atomic_reactants = reaction->getReactants().first;

    assert(atomic_products.size() == 1 and atomic_reactants.size() == 1);

    auto atomlist = spc.findAtoms(atomic_reactants.begin()->first); // search all active molecules
    auto& target_particle =
        *random_internal.sample(atomlist.begin(), atomlist.end()); // target particle to swap
    auto& group = *spc.findGroupContaining(target_particle);       // find enclosing group

    auto& group_change = change.groups.emplace_back();
    group_change.relative_atom_indices.push_back(group.getParticleIndex(target_particle));
    group_change.group_index = spc.getGroupIndex(group);
    group_change.internal = true;
    group_change.dNswap = true;

    const auto new_atomid = atomic_products.begin()->first; // atomid of new atom type
    swapParticleProperties(target_particle, new_atomid);
}

/**
 * @brief Swap particle to another type
 * @param particle Existing particle to be swapped out
 * @param new_atomid The new atomtype to swap in
 *
 * This will keep original positions and internal orientation of particles, while swapping
 * the atomid and other properties related to the new atom type
 *
 * @todo This could make use of policies to customize how particles should be updated. Split to
 * helper class.
 */
void SpeciationMove::swapParticleProperties(Particle& particle, const int new_atomid)
{
    Particle new_particle(atoms.at(new_atomid), particle.pos); // new particle with old position
    if (new_particle.hasExtension() ||
        particle.hasExtension()) { // keep also other orientational data
        auto& source_ext = new_particle.getExt();
        auto& target_ext = particle.getExt();
        if (target_ext.isCylindrical() || source_ext.isCylindrical()) {
            source_ext.setDirections(new_particle.traits().sphero_cylinder, target_ext.scdir,
                                     target_ext.patchdir);
        }
        if (target_ext.isDipolar() || source_ext.isQuadrupolar()) {
            throw std::runtime_error("dipolar/quadrupolar properties not yet implemented");
        }
    }
    particle = new_particle; // deep copy
}

void SpeciationMove::activateProducts(Change& change)
{
    activateAtomicGroups(change);
    activateMolecularGroups(change);
}

void SpeciationMove::activateMolecularGroups(Change& change)
{
    namespace rv = std::views;
    auto selection = reaction->only_neutral_molecules ? Space::Selection::INACTIVE_NEUTRAL
                                                      : Space::Selection::INACTIVE;

    auto molecular_products = reaction->getProducts().second |
                              rv::filter(ReactionData::not_implicit_group) |
                              rv::filter(ReactionData::is_molecular_group);

    for (auto [molid, number_to_insert] : molecular_products) {
        auto inactive = spc.findMolecules(molid, selection) |
                        ranges::views::sample(number_to_insert, random_internal.engine) |
                        ranges::to<std::vector<std::reference_wrapper<Group>>>;

        std::ranges::for_each(inactive, [&](auto& group) {
            auto [change_data, bias] = molecular_group_bouncer->activate(group);
            change.groups.emplace_back(change_data);
            bias_energy += bias;
        });
    }
}

void SpeciationMove::activateAtomicGroups(Change& change)
{
    namespace rv = std::views;
    auto atomic_products = reaction->getProducts().second |
                           rv::filter(ReactionData::not_implicit_group) |
                           rv::filter(ReactionData::is_atomic_group);

    for (auto [molid, number_to_insert] : atomic_products) {
        auto groups = spc.findMolecules(molid, Space::Selection::ALL);
        auto [change_data, bias] =
            atomic_group_bouncer->activate(*groups.begin(), number_to_insert);
        change.groups.emplace_back(change_data);
        bias_energy += bias;
    }
}

/**
 * @brief Deactivate all atomic and molecular reactants.
 */
void SpeciationMove::deactivateReactants(Change& change)
{
    deactivateAtomicGroups(change);
    deactivateMolecularGroups(change);
}

void SpeciationMove::deactivateMolecularGroups(Change& change)
{
    namespace rv = std::views;
    auto nonzero_stoichiometric_coeff = [](auto key_value) { return key_value.second > 0; };
    auto selection = (reaction->only_neutral_molecules) ? Space::Selection::ACTIVE_NEUTRAL
                                                        : Space::Selection::ACTIVE;

    auto molecular_reactants =
        reaction->getReactants().second | rv::filter(ReactionData::not_implicit_group) |
        rv::filter(nonzero_stoichiometric_coeff) | rv::filter(ReactionData::is_molecular_group);

    for (const auto& [molid, number_to_delete] : molecular_reactants) {
        auto groups = spc.findMolecules(molid, selection) |
                      ranges::views::sample(number_to_delete, random_internal.engine) |
                      ranges::to<std::vector<std::reference_wrapper<Group>>>;
        std::for_each(groups.begin(), groups.end(), [&](auto& group) {
            auto [change_data, bias] = molecular_group_bouncer->deactivate(group);
            change.groups.emplace_back(change_data);
            bias_energy += bias;
        });
    }
}

void SpeciationMove::deactivateAtomicGroups(Change& change)
{
    namespace rv = std::views;
    auto nonzero_stoichiometric_coeff = [](auto key_value) { return key_value.second > 0; };

    auto atomic_reactants =
        reaction->getReactants().second | rv::filter(ReactionData::not_implicit_group) |
        rv::filter(nonzero_stoichiometric_coeff) | rv::filter(ReactionData::is_atomic_group);

    for (auto [molid, number_to_delete] : atomic_reactants) {
        auto groups = spc.findMolecules(molid, Space::Selection::ALL);
        assert(range_size(groups) == 1);
        auto& group = *groups.begin();
        auto [change_data, bias] = atomic_group_bouncer->deactivate(group, number_to_delete);
        change.groups.push_back(change_data);
        bias_energy += bias;
    }
}

TEST_CASE("[Faunus] Speciation - Ranges::sample")
{
    using ranges::views::sample;
    std::vector<int> vec = {1, 2, 3, 4};
    auto take_nothing = vec | sample(0);
    auto take_less = vec | sample(2);
    auto take_all = vec | sample(4);
    auto take_too_much = vec | sample(10);

    CHECK_EQ(range_size(take_nothing), 0);
    CHECK_EQ(range_size(take_less), 2);
    CHECK_EQ(range_size(take_all), 4);
    CHECK_EQ(range_size(take_too_much), 4);
}

TEST_CASE("[Faunus] ReactionDirectionRatio")
{
    using doctest::Approx;
    using Faunus::Speciation::ReactionDirectionRatio;

    // Set up minimal atoms/molecules/reactions for the test
    Faunus::atoms = R"([{"a": {"r": 1.0}}])"_json.get<decltype(atoms)>();
    Faunus::molecules =
        R"([{"A": {"atomic": false, "activity": 0.5}}, {"B": {"atomic": false, "activity": 0.5}}])"_json
            .get<decltype(molecules)>();
    Faunus::reactions = R"([{"A = B": {"lnK": -1.0}}])"_json.get<decltype(reactions)>();

    ReactionDirectionRatio ratio;
    auto it = Faunus::reactions.begin();

    // Access creates entry
    auto& data = ratio[it];
    CHECK(data.right.empty());
    CHECK(data.left.empty());

    // Update RIGHT direction: 3 accepted, 2 rejected = 5 attempts, ratio 3/5
    data.update(ReactionData::Direction::RIGHT, true);
    data.update(ReactionData::Direction::RIGHT, true);
    data.update(ReactionData::Direction::RIGHT, true);
    data.update(ReactionData::Direction::RIGHT, false);
    data.update(ReactionData::Direction::RIGHT, false);

    CHECK_EQ(data.right.size(), 5);
    CHECK_EQ(data.right.avg(), Approx(3.0 / 5.0));

    // Update LEFT direction: 2 accepted, 1 rejected = 3 attempts, ratio 2/3
    data.update(ReactionData::Direction::LEFT, true);
    data.update(ReactionData::Direction::LEFT, true);
    data.update(ReactionData::Direction::LEFT, false);

    CHECK_EQ(data.left.size(), 3);
    CHECK_EQ(data.left.avg(), Approx(2.0 / 3.0));

    // Verify to_json output
    json j;
    ratio.to_json(j);
    CHECK(j.contains("reactions"));
    CHECK(j["reactions"].contains("A = B"));
    CHECK_EQ(j["reactions"]["A = B"]["attempts"], 8);
    CHECK_EQ(j["reactions"]["A = B"]["acceptance -->"].get<double>(), Approx(3.0 / 5.0));
    CHECK_EQ(j["reactions"]["A = B"]["acceptance <--"].get<double>(), Approx(2.0 / 3.0));
}

TEST_CASE("[Faunus] ReactionValidator")
{
    using namespace Faunus;

    pc::temperature = 298.15_K;

    SUBCASE("atomic GCMC feasibility")
    {
        // Two atom types, two atomic groups
        Faunus::atoms = R"([
            {"c": {"r": 1.0, "mw": 1.0}},
            {"d": {"r": 1.0, "mw": 1.0}}
        ])"_json.get<decltype(atoms)>();

        Faunus::molecules = R"([
            {"C": {"atomic": true, "atoms": ["c"]}},
            {"D": {"atomic": true, "atoms": ["d"]}}
        ])"_json.get<decltype(molecules)>();

        // Reaction: D = C  (remove one d, insert one c)
        Faunus::reactions = R"([{"D = C": {"lnK": 0.0}}])"_json.get<decltype(reactions)>();
        auto& rxn = Faunus::reactions.front();
        rxn.setDirection(ReactionData::Direction::RIGHT);

        Space spc;
        spc.geometry = R"({"type": "cuboid", "length": 50})"_json;

        // Group D: 10 total d-particles, 5 active
        json j_insert = json::array();
        j_insert.push_back({{"D", {{"N", 10}, {"inactive", 5}}}});
        j_insert.push_back({{"C", {{"N", 5}, {"inactive", 5}}}});
        InsertMoleculesInSpace::insertMolecules(j_insert, spc);

        Speciation::ReactionValidator validator(spc);

        // 5 active d's and 5 inactive c capacity => feasible
        CHECK(validator.isPossible(rxn));

        // Deplete D to 0 active
        spc.groups.at(0).resize(0);
        CHECK_FALSE(validator.isPossible(rxn));

        // Refill D, fill C to capacity
        spc.groups.at(0).resize(10);
        spc.groups.at(1).resize(5); // C now full
        CHECK_FALSE(validator.isPossible(rxn));
    }

    SUBCASE("molecular GCMC feasibility")
    {
        Faunus::atoms = R"([
            {"ow": {"r": 1.5, "mw": 16.0}}
        ])"_json.get<decltype(atoms)>();

        Faunus::molecules = R"([{
            "M": {
                "activity": 0.1,
                "structure": [{"ow": [0, 0, 0]}]
            }
        }])"_json.get<decltype(molecules)>();

        // Reaction: = M (pure insertion in RIGHT direction)
        Faunus::reactions = R"([{"= M": {"lnK": 0.0}}])"_json.get<decltype(reactions)>();
        auto& rxn = Faunus::reactions.front();
        rxn.setDirection(ReactionData::Direction::RIGHT);

        Space spc;
        spc.geometry = R"({"type": "cuboid", "length": 50})"_json;

        // 4 molecular groups: 2 active, 2 inactive
        json j_insert = json::array();
        j_insert.push_back({{"M", {{"N", 4}, {"inactive", 2}}}});
        InsertMoleculesInSpace::insertMolecules(j_insert, spc);
        CHECK_EQ(spc.groups.size(), 4);

        Speciation::ReactionValidator validator(spc);

        // 2 inactive groups available => can insert
        CHECK(validator.isPossible(rxn));

        // Activate all groups (no more inactive capacity)
        for (auto& g : spc.groups) {
            if (g.empty()) {
                g.activate(g.inactive().begin(), g.inactive().end());
            }
        }
        CHECK_FALSE(validator.isPossible(rxn));
    }

    SUBCASE("atom swap feasibility")
    {
        Faunus::atoms = R"([
            {"a": {"r": 1.0, "mw": 1.0}},
            {"b": {"r": 1.0, "mw": 1.0}}
        ])"_json.get<decltype(atoms)>();

        Faunus::molecules = R"([
            {"A": {"atomic": true, "atoms": ["a"]}},
            {"B": {"atomic": true, "atoms": ["b"]}}
        ])"_json.get<decltype(molecules)>();

        // Reaction: a + A = b + B (swap a->b plus molecular exchange)
        Faunus::reactions = R"([{"a + A = b + B": {"lnK": 0.0}}])"_json.get<decltype(reactions)>();
        auto& rxn = Faunus::reactions.front();
        rxn.setDirection(ReactionData::Direction::RIGHT);

        Space spc;
        spc.geometry = R"({"type": "cuboid", "length": 50})"_json;

        json j_insert = json::array();
        j_insert.push_back({{"A", {{"N", 5}, {"inactive", 2}}}});
        j_insert.push_back({{"B", {{"N", 5}, {"inactive", 2}}}});
        InsertMoleculesInSpace::insertMolecules(j_insert, spc);

        Speciation::ReactionValidator validator(spc);

        // 3 active 'a' particles, capacity for b => feasible
        CHECK(validator.isPossible(rxn));

        // Empty A group of all active atoms
        spc.groups.at(0).resize(0);
        CHECK_FALSE(validator.isPossible(rxn));
    }

    SUBCASE("implicit group feasibility")
    {
        Faunus::atoms = R"([
            {"x": {"r": 1.0, "mw": 1.0}}
        ])"_json.get<decltype(atoms)>();

        Faunus::molecules = R"([
            {"I": {"implicit": true, "atoms": ["x"]}},
            {"X": {"atomic": true, "atoms": ["x"]}}
        ])"_json.get<decltype(molecules)>();

        // Reaction: I + X = (consume one implicit I and one atomic X)
        Faunus::reactions = R"([{"I + X = ": {"lnK": 0.0}}])"_json.get<decltype(reactions)>();
        auto& rxn = Faunus::reactions.front();
        rxn.setDirection(ReactionData::Direction::RIGHT);

        Space spc;
        spc.geometry = R"({"type": "cuboid", "length": 50})"_json;

        // Insert atomic group X with some particles
        json j_insert = json::array();
        j_insert.push_back({{"X", {{"N", 5}}}});
        InsertMoleculesInSpace::insertMolecules(j_insert, spc);

        // Set implicit reservoir
        auto molid_I = Faunus::findName(molecules, "I")->id();
        spc.getImplicitReservoir()[molid_I] = 10;

        Speciation::ReactionValidator validator(spc);
        CHECK(validator.isPossible(rxn));

        // Deplete reservoir
        spc.getImplicitReservoir()[molid_I] = 0;
        CHECK_FALSE(validator.isPossible(rxn));
    }
}

TEST_CASE("[Faunus] AtomicGroupDeActivator")
{
    using namespace Faunus;

    pc::temperature = 298.15_K;

    Faunus::atoms = R"([
        {"na": {"r": 1.9, "mw": 23.0}}
    ])"_json.get<decltype(atoms)>();

    Faunus::molecules = R"([
        {"salt": {"atomic": true, "atoms": ["na"]}}
    ])"_json.get<decltype(molecules)>();

    auto make_space = [](int total, int inactive) {
        Space spc;
        spc.geometry = R"({"type": "cuboid", "length": 50})"_json;
        json j_insert = json::array();
        j_insert.push_back({{"salt", {{"N", total}, {"inactive", inactive}}}});
        InsertMoleculesInSpace::insertMolecules(j_insert, spc);
        return spc;
    };

    SUBCASE("activate expands group")
    {
        auto spc = make_space(10, 5);     // 10 total, 5 active
        auto old_spc = make_space(10, 5); // identical
        Random random;

        Speciation::AtomicGroupDeActivator bouncer(spc, old_spc, random);

        auto& group = spc.groups.front();
        CHECK_EQ(group.size(), 5);
        CHECK_EQ(group.capacity(), 10);

        auto [change, bias] = bouncer.activate(group, 3);

        CHECK_EQ(group.size(), 8);
        CHECK(change.dNatomic);
        CHECK(change.internal);
        CHECK_EQ(change.relative_atom_indices.size(), 3);
        CHECK_EQ(bias, 0.0);

        // Verify new particles are within geometry bounds
        for (auto index : change.relative_atom_indices) {
            CHECK_FALSE(spc.geometry.collision(group.begin()[index].pos));
        }
    }

    SUBCASE("deactivate contracts group")
    {
        auto spc = make_space(10, 5);     // 10 total, 5 active
        auto old_spc = make_space(10, 5); // identical
        Random random;

        Speciation::AtomicGroupDeActivator bouncer(spc, old_spc, random);

        auto& group = spc.groups.front();
        CHECK_EQ(group.size(), 5);

        auto [change, bias] = bouncer.deactivate(group, 2);

        CHECK_EQ(group.size(), 3);
        CHECK(change.dNatomic);
        CHECK(change.internal);
        CHECK_EQ(change.relative_atom_indices.size(), 2);
        CHECK_EQ(bias, 0.0);
    }

    SUBCASE("error on invalid input")
    {
        auto spc = make_space(10, 5);
        auto old_spc = make_space(10, 5);
        Random random;

        Speciation::AtomicGroupDeActivator bouncer(spc, old_spc, random);
        auto& group = spc.groups.front();

        // Activate more than capacity: should throw
        CHECK_THROWS(bouncer.activate(group, 6));

        // Deactivate more than active size: should throw
        CHECK_THROWS(bouncer.deactivate(group, 6));
    }
}

TEST_CASE("[Faunus] MolecularGroupDeActivator")
{
    using namespace Faunus;

    pc::temperature = 298.15_K;

    Faunus::atoms = R"([
        {"ow": {"r": 1.5, "mw": 16.0, "q": -0.8}},
        {"hw": {"r": 1.0, "mw": 1.0, "q": 0.4}}
    ])"_json.get<decltype(atoms)>();

    Faunus::molecules = R"([{
        "water": {
            "structure": [
                {"ow": [0, 0, 0]},
                {"hw": [1, 0, 0]},
                {"hw": [0, 1, 0]}
            ]
        }
    }])"_json.get<decltype(molecules)>();

    SUBCASE("activate inactive molecular group")
    {
        Space spc;
        spc.geometry = R"({"type": "cuboid", "length": 50})"_json;

        // Insert 1 molecular group, fully inactive
        json j_insert = json::array();
        j_insert.push_back({{"water", {{"N", 1}, {"inactive", 1}}}});
        InsertMoleculesInSpace::insertMolecules(j_insert, spc);

        CHECK_EQ(spc.groups.size(), 1);
        auto& group = spc.groups.front();
        CHECK(group.empty());
        CHECK_EQ(group.capacity(), 3);

        Random random;
        Speciation::MolecularGroupDeActivator bouncer(spc, random, false);

        auto [change, bias] = bouncer.activate(group);

        CHECK_EQ(group.size(), 3); // fully activated
        CHECK(change.all);
        CHECK(change.internal);
        CHECK_EQ(change.relative_atom_indices.size(), 3);
        CHECK_EQ(bias, 0.0); // no bond bias

        // Mass center should be within geometry bounds
        CHECK_FALSE(spc.geometry.collision(group.mass_center));
    }

    SUBCASE("deactivate active molecular group")
    {
        Space spc;
        spc.geometry = R"({"type": "cuboid", "length": 50})"_json;

        // Insert 1 active molecular group
        json j_insert = json::array();
        j_insert.push_back({{"water", {{"N", 1}}}});
        InsertMoleculesInSpace::insertMolecules(j_insert, spc);

        auto& group = spc.groups.front();
        CHECK_EQ(group.size(), 3);

        Random random;
        Speciation::MolecularGroupDeActivator bouncer(spc, random, false);

        auto [change, bias] = bouncer.deactivate(group);

        CHECK(group.empty());
        CHECK(change.all);
        CHECK(change.internal);
        CHECK_EQ(change.relative_atom_indices.size(), 3);
        CHECK_EQ(bias, 0.0); // no bond bias
    }

    SUBCASE("partial activation throws")
    {
        Space spc;
        spc.geometry = R"({"type": "cuboid", "length": 50})"_json;

        json j_insert = json::array();
        j_insert.push_back({{"water", {{"N", 1}, {"inactive", 1}}}});
        InsertMoleculesInSpace::insertMolecules(j_insert, spc);

        Random random;
        Speciation::MolecularGroupDeActivator bouncer(spc, random, false);

        auto& group = spc.groups.front();
        CHECK(group.empty());

        // Requesting partial activation (1 of 3) should throw
        CHECK_THROWS(bouncer.activate(group, 1));
    }

    SUBCASE("partial deactivation throws")
    {
        Space spc;
        spc.geometry = R"({"type": "cuboid", "length": 50})"_json;

        json j_insert = json::array();
        j_insert.push_back({{"water", {{"N", 1}}}});
        InsertMoleculesInSpace::insertMolecules(j_insert, spc);

        Random random;
        Speciation::MolecularGroupDeActivator bouncer(spc, random, false);

        auto& group = spc.groups.front();
        CHECK_EQ(group.size(), 3);

        // Requesting partial deactivation (1 of 3) should throw
        CHECK_THROWS(bouncer.deactivate(group, 1));
    }
}

TEST_CASE("[Faunus] MolecularGroupDeActivator - bond bias")
{
    using namespace Faunus;
    using doctest::Approx;

    pc::temperature = 298.15_K;

    Faunus::atoms = R"([
        {"a1": {"r": 1.5, "mw": 16.0}},
        {"a2": {"r": 1.0, "mw": 1.0}}
    ])"_json.get<decltype(atoms)>();

    // Molecule with a harmonic bond between atoms 0 and 1
    Faunus::molecules = R"([{
        "dimer": {
            "structure": [
                {"a1": [0, 0, 0]},
                {"a2": [2, 0, 0]}
            ],
            "bondlist": [{"harmonic": {"index": [0, 1], "k": 5.0, "req": 1.0}}]
        }
    }])"_json.get<decltype(molecules)>();

    SUBCASE("deactivate returns positive bond energy bias")
    {
        Space spc;
        spc.geometry = R"({"type": "cuboid", "length": 50})"_json;

        json j_insert = json::array();
        j_insert.push_back({{"dimer", {{"N", 1}}}});
        InsertMoleculesInSpace::insertMolecules(j_insert, spc);

        auto& group = spc.groups.front();
        CHECK_EQ(group.size(), 2);

        Random random;
        Speciation::MolecularGroupDeActivator bouncer(spc, random, true);

        auto [change, bias] = bouncer.deactivate(group);

        CHECK(group.empty());
        // Deactivation bias should be positive bond energy (in kT units)
        CHECK_GT(bias, 0.0);
    }

    SUBCASE("activate returns negative bond energy bias")
    {
        Space spc;
        spc.geometry = R"({"type": "cuboid", "length": 50})"_json;

        json j_insert = json::array();
        j_insert.push_back({{"dimer", {{"N", 1}, {"inactive", 1}}}});
        InsertMoleculesInSpace::insertMolecules(j_insert, spc);

        auto& group = spc.groups.front();
        CHECK(group.empty());

        Random random;
        Speciation::MolecularGroupDeActivator bouncer(spc, random, true);

        auto [change, bias] = bouncer.activate(group);

        CHECK_EQ(group.size(), 2);
        // After random placement, bond energy will generally be non-zero
        // Activation bias should be negative of bond energy
        CHECK_NE(bias, 0.0);
        CHECK_LT(bias, 0.0); // negative bond energy
    }
}

TEST_CASE("[Faunus] AtomicGroupDeActivator - old_spc synchronization")
{
    using namespace Faunus;

    pc::temperature = 298.15_K;

    Faunus::atoms = R"([
        {"na": {"r": 1.9, "mw": 23.0}}
    ])"_json.get<decltype(atoms)>();

    Faunus::molecules = R"([
        {"salt": {"atomic": true, "atoms": ["na"]}}
    ])"_json.get<decltype(molecules)>();

    auto make_space = [](int total, int inactive) {
        Space spc;
        spc.geometry = R"({"type": "cuboid", "length": 50})"_json;
        json j_insert = json::array();
        j_insert.push_back({{"salt", {{"N", total}, {"inactive", inactive}}}});
        InsertMoleculesInSpace::insertMolecules(j_insert, spc);
        return spc;
    };

    SUBCASE("old_spc particles rearranged on deactivate")
    {
        auto spc = make_space(10, 5);     // 10 total, 5 active
        auto old_spc = make_space(10, 5); // identical

        // Give each particle in old_spc a unique position for tracking
        auto& old_group = old_spc.groups.front();
        for (int i = 0; i < (int)old_group.capacity(); i++) {
            old_group.begin()[i].pos = Point(i, 0, 0);
        }

        // Copy same positions to spc
        auto& group = spc.groups.front();
        for (int i = 0; i < (int)group.capacity(); i++) {
            group.begin()[i].pos = Point(i, 0, 0);
        }

        Random random;
        Speciation::AtomicGroupDeActivator bouncer(spc, old_spc, random);

        auto [change, bias] = bouncer.deactivate(group, 2);

        // spc group should have shrunk
        CHECK_EQ(group.size(), 3);

        // old_spc group should keep original size but have particles rearranged
        // to match the swap pattern in spc
        CHECK_EQ(old_group.size(), 5); // old_spc not resized

        // The deactivated indices should point to valid positions
        for (auto idx : change.relative_atom_indices) {
            CHECK_LT(idx, 5); // indices within original active range
        }
    }
}

TEST_CASE("[Faunus] swapParticleProperties")
{
    using namespace Faunus;

    pc::temperature = 298.15_K;

    Faunus::atoms = R"([
        {"typeA": {"r": 2.0, "mw": 10.0, "q": 1.0, "sigma": 3.0}},
        {"typeB": {"r": 3.0, "mw": 20.0, "q": -1.0, "sigma": 4.0}}
    ])"_json.get<decltype(atoms)>();

    SUBCASE("basic swap preserves position and updates properties")
    {
        const Point original_pos(1.5, 2.5, 3.5);
        Particle particle(atoms.at(0), original_pos);
        CHECK_EQ(particle.id, 0);
        CHECK_EQ(particle.charge, 1.0);
        CHECK_EQ(particle.pos, original_pos);

        move::SpeciationMove::swapParticleProperties(particle, 1);

        // Position must be preserved
        CHECK_EQ(particle.pos, original_pos);
        // Properties should come from the new atom type
        CHECK_EQ(particle.id, 1);
        CHECK_EQ(particle.charge, -1.0);
    }
}

TEST_CASE("[Faunus] ReactionValidator - reverse direction")
{
    using namespace Faunus;

    pc::temperature = 298.15_K;

    // Two atom types in atomic groups
    Faunus::atoms = R"([
        {"c": {"r": 1.0, "mw": 1.0}},
        {"d": {"r": 1.0, "mw": 1.0}}
    ])"_json.get<decltype(atoms)>();

    Faunus::molecules = R"([
        {"C": {"atomic": true, "atoms": ["c"]}},
        {"D": {"atomic": true, "atoms": ["d"]}}
    ])"_json.get<decltype(molecules)>();

    // Reaction: D = C
    Faunus::reactions = R"([{"D = C": {"lnK": 0.0}}])"_json.get<decltype(reactions)>();
    auto& rxn = Faunus::reactions.front();

    Space spc;
    spc.geometry = R"({"type": "cuboid", "length": 50})"_json;

    // D: 5 active out of 5 (no inactive capacity)
    // C: 0 active out of 5 (5 inactive capacity)
    json j_insert = json::array();
    j_insert.push_back({{"D", {{"N", 5}}}});
    j_insert.push_back({{"C", {{"N", 5}, {"inactive", 5}}}});
    InsertMoleculesInSpace::insertMolecules(j_insert, spc);

    Speciation::ReactionValidator validator(spc);

    // RIGHT: consume D, produce C => feasible (5 D's available, 5 C capacity)
    rxn.setDirection(ReactionData::Direction::RIGHT);
    CHECK(validator.isPossible(rxn));

    // LEFT: consume C, produce D => infeasible (0 active C's, D is full)
    rxn.setDirection(ReactionData::Direction::LEFT);
    CHECK_FALSE(validator.isPossible(rxn));
}

TEST_CASE("[Faunus] ReactionValidator - combined atomic and molecular reaction")
{
    using namespace Faunus;

    pc::temperature = 298.15_K;

    Faunus::atoms = R"([
        {"a": {"r": 1.0, "mw": 1.0}},
        {"m": {"r": 1.5, "mw": 16.0}}
    ])"_json.get<decltype(atoms)>();

    Faunus::molecules = R"([
        {"A": {"atomic": true, "atoms": ["a"]}},
        {"M": {"activity": 0.1, "structure": [{"m": [0, 0, 0]}]}}
    ])"_json.get<decltype(molecules)>();

    // Reaction: A + M = (consume one atomic A and one molecular M)
    Faunus::reactions = R"([{"A + M = ": {"lnK": 0.0}}])"_json.get<decltype(reactions)>();
    auto& rxn = Faunus::reactions.front();
    rxn.setDirection(ReactionData::Direction::RIGHT);

    Space spc;
    spc.geometry = R"({"type": "cuboid", "length": 50})"_json;

    json j_insert = json::array();
    j_insert.push_back({{"A", {{"N", 5}}}});
    j_insert.push_back({{"M", {{"N", 3}, {"inactive", 1}}}});
    InsertMoleculesInSpace::insertMolecules(j_insert, spc);

    Speciation::ReactionValidator validator(spc);

    // 5 active A particles and 2 active M groups => feasible to consume both
    CHECK(validator.isPossible(rxn));

    // Deplete A group
    spc.groups.at(0).resize(0);
    CHECK_FALSE(validator.isPossible(rxn));

    // Refill A, deactivate all M groups
    spc.groups.at(0).resize(5);
    for (size_t i = 1; i < spc.groups.size(); i++) {
        auto& g = spc.groups.at(i);
        if (!g.empty()) {
            g.deactivate(g.begin(), g.end());
        }
    }
    CHECK_FALSE(validator.isPossible(rxn));
}

void SpeciationMove::_move(Change& change)
{
    if (Faunus::reactions.empty()) {
        return;
    }
    bias_energy = 0.0;
    try {
        setRandomReactionAndDirection();
        if (reaction_validator.isPossible(*reaction)) {
            atomicSwap(change);
            deactivateReactants(change);
            activateProducts(change);
            std::sort(change.groups.begin(),
                      change.groups.end()); // change groups *must* be sorted!
            if (change) {
                change.matter_change = true;
                updateGroupMassCenters(change);
            }
        }
    }
    catch (Speciation::SpeciationMoveException&) {
        change.clear();
    }
}

void SpeciationMove::setRandomReactionAndDirection()
{
    reaction = random_internal.sample(reactions.begin(), reactions.end());
    reaction->setRandomDirection(random_internal);
}

/**
 * Speciation move may induce a change in molecular mass centers.
 * Curently updated for the following moves:
 *
 * - Swap moves (as the swapped atom may have a different mass)
 */
void SpeciationMove::updateGroupMassCenters(const Change& change) const
{
    namespace rv = ranges::cpp20::views;

    auto atomic_or_swap = [](const Change::GroupChange& c) { return c.dNatomic || c.dNswap; };
    auto to_group = [&](const Change::GroupChange& c) -> Group& {
        return spc.groups.at(c.group_index);
    };
    auto has_mass_center = [](Group& group) { return group.massCenter().has_value(); };
    auto groups = change.groups | rv::filter(atomic_or_swap) | rv::transform(to_group) |
                  rv::filter(has_mass_center);

    std::ranges::for_each(groups, [&](Group& group) {
        group.updateMassCenter(spc.geometry.getBoundaryFunc(), group.massCenter().value());
    });
}

/**
 * The acceptance/rejection of the move is affected by the equilibrium constant,
 * but unaffected by the change in internal bond energy
 */
double SpeciationMove::bias([[maybe_unused]] Change& change, [[maybe_unused]] double old_energy,
                            [[maybe_unused]] double new_energy)
{
    return reaction->freeEnergy() + bias_energy;
}

void SpeciationMove::_accept([[maybe_unused]] Change& change)
{
    namespace rv = std::views;
    direction_ratio[reaction].update(reaction->getDirection(), true);

    auto implicit_reactants =
        reaction->getReactants().second | rv::filter(ReactionData::is_implicit_group);
    for (const auto& [molid, nu] : implicit_reactants) {
        spc.getImplicitReservoir()[molid] -= nu;
        average_reservoir_size[molid] += spc.getImplicitReservoir().at(molid);
    }

    auto implicit_products =
        reaction->getProducts().second | rv::filter(ReactionData::is_implicit_group);
    for (const auto& [molid, nu] : implicit_products) {
        spc.getImplicitReservoir()[molid] += nu;
        average_reservoir_size[molid] += spc.getImplicitReservoir().at(molid);
    }
}

void SpeciationMove::_reject([[maybe_unused]] Change& change)
{
    namespace rv = std::views;
    direction_ratio[reaction].update(reaction->getDirection(), false);

    auto implicit_reactants =
        reaction->getReactants().second | rv::filter(ReactionData::is_implicit_group);
    for (auto [molid, nu] : implicit_reactants) {
        average_reservoir_size[molid] += spc.getImplicitReservoir().at(molid);
    }

    auto implicit_products =
        reaction->getProducts().second | rv::filter(ReactionData::is_implicit_group);
    for (auto [molid, nu] : implicit_products) {
        average_reservoir_size[molid] += spc.getImplicitReservoir().at(molid);
    }
}

SpeciationMove::SpeciationMove(Space& spc, Space& old_spc, std::string_view name,
                               std::string_view cite)
    : Move(spc, name, cite)
    , random_internal(slump)
    , reaction_validator(spc)
{
    molecular_group_bouncer =
        std::make_unique<Speciation::MolecularGroupDeActivator>(spc, random_internal, true);
    atomic_group_bouncer =
        std::make_unique<Speciation::AtomicGroupDeActivator>(spc, old_spc, random_internal);
}

SpeciationMove::SpeciationMove(Space& spc, Space& old_spc)
    : SpeciationMove(spc, old_spc, "rcmc", "doi:10/fqcpg3")
{
}

void SpeciationMove::_from_json([[maybe_unused]] const json& j) {}

} // namespace Faunus::move
