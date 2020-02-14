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

void SpeciationMove::setOther(Tspace &ospc) { otherspc = &ospc; }

/*
 * This function is only performing checks
 */
bool SpeciationMove::checkInsertProducts(ReactionData &reaction) {
    // Check whether it is possible to insert products (are there any inactive ones?)

    // rewrite to refactored ReactionData
    // auto [atomic_products, molecular_products] = reaction.getProducts();
    // std::cout << atomic_products.size() << " " << molecular_products.size() << std::endl;
    // for (auto [molid, number_to_insert] : molecular_products) {
    // }

    for (auto [molid, number_to_insert] : reaction.moleculesToAdd(forward)) { // Additional checks
        auto molecule_list = spc.findMolecules(molid, Tspace::ALL);
        if (molecules[molid].atomic) {
            if (range_size(molecule_list) != 1) // There can be only one
                throw std::runtime_error("Bad definition: One group per atomic molecule!");
            auto git = molecule_list.begin();
            if ((git->size() + number_to_insert) > git->capacity()) { // Assure that there are atoms enough in the group
                faunus_logger->warn("molecule {} has reached its maximum capacity", Faunus::molecules[molid].name);
                return false;                                 // Slip out the back door
            }
        } else {
            if (neutral) // Only neutral molecules react
                molecule_list = spc.findMolecules(molid, Tspace::INACTIVE_NEUTRAL);
            else
                molecule_list = spc.findMolecules(molid, Tspace::INACTIVE);
            if (range_size(molecule_list) < number_to_insert) {
                faunus_logger->warn("molecule {} has reached its maximum capacity", Faunus::molecules[molid].name);
                return false; // Not possible to perform change, escape through the back door
            }
        }
    }
    return true;
}

bool SpeciationMove::swapReaction(Change &change, ReactionData &reaction) {
    // Perform a swap reaction, e.g., the (de)protonation of an atomic species

    if (reaction.swap) {
        auto [atomic_products, molecular_products] = reaction.getProducts();
        auto [atomic_reactants, molecular_reactants] = reaction.getReactants();
        assert(atomic_products.size() == 1 and atomic_reactants.size() == 1);
    }

    if (reaction.swap == true) {
        auto m1 = reaction.atomsToAdd(not forward); // Swap checks
        auto m2 = reaction.atomsToAdd(forward);
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

/*
 * Reduce an atomic group by `number_to_delete` particles. The deleted particles are
 * picked by random; moved the the end of the group; then deactivated. In order for
 * the Hamiltonian to pick up the energy change, particles in the reference Space
 * (`old_target`) are swapped to the same index, albeit not deactivated.
 */
Change::data SpeciationMove::contractAtomicGroup(Space::Tgroup &target, Space::Tgroup &old_target,
                                                 int number_to_delete) {
    assert(target.atomic);
    Change::data d;

    if ((int)target.size() - number_to_delete >= 0) {
        d.index = &target - &spc.groups.front(); // index of moved group
        d.internal = true;
        d.dNatomic = true;
        for (int i = 0; i < number_to_delete; i++) {                       // deactivate m.second m.first atoms
            auto random_atom = slump.sample(target.begin(), target.end()); // iterator to random atom
            auto last_atom = target.end() - 1;                             // iterator to last atom
            int dist = std::distance(random_atom, target.end());           // distance to random atom from end
            // Shuffle back to end, both in trial and old target
            if (std::distance(random_atom, last_atom) > 1) {
                std::iter_swap(random_atom, last_atom);
                std::iter_swap(old_target.end() - dist - i, old_target.end() - (1 + i));
            }
            d.atoms.push_back(std::distance(target.begin(), last_atom));
            target.deactivate(last_atom, target.end());
        }
        assert(std::is_sorted(d.atoms.begin(), d.atoms.end())); // redundant
    } else {
        faunus_logger->warn("atomic group {} is depleted; increase simulation volume?",
                            Faunus::molecules[target.id].name);
    }
    return d;
}

/*
 * Deactivate an active molecular group. If there are internal bonds, the total
 * bond energy is stored and used to avoid that the bond energy affect acceptance.
 */
Change::data SpeciationMove::deactivateMolecularGroup(Space::Tgroup &target) {
    assert(target.atomic == false);             // group must be molecular
    assert(target.size() == target.capacity()); // group must be active

    target.unwrap(spc.geo.getDistanceFunc()); // when in storage, remove PBC

    // We store the bonded energy of the deactivated molecule
    // The change in bonded energy should not affect the acceptance/rejection of the move
    for (auto &bond : Faunus::molecules.at(target.id).bonds) {
        auto bond_clone = bond->clone();
        bond_clone->shift(std::distance(spc.p.begin(), target.begin()));
        Potential::setBondEnergyFunction(bond_clone, spc.p);
        bond_energy += bond_clone->energy(spc.geo.getDistanceFunc());
    }

    target.deactivate(target.begin(), target.end()); // deactivate whole group

    Change::data d;
    d.internal = true;
    d.index = &target - &spc.groups.front();      // index of moved group
    d.all = true;                                 // all atoms in group were moved
    for (int i = 0; i < target.capacity(); i++) { // crete list w. all atom index
        d.atoms.push_back(i);
    }
    return d;
}

void SpeciationMove::deactivateReactants(Change &change, std::vector<ReactionData>::iterator reaction) {
    // Deactivate reactants
    for (auto [molid, number_to_delete] : reaction->moleculesToAdd(not forward)) { // Delete

        if (molecules[molid].atomic) { // The reagents is an atom
            auto target = spc.findMolecules(molid, Tspace::ALL).begin();
            auto other_target = otherspc->findMolecules(molid, Tspace::ALL).begin();
            auto change_data = contractAtomicGroup(*target, *other_target, number_to_delete);
            if (not change_data.atoms.empty()) {
                change.groups.push_back(change_data);
            }
        } else { // The reactant is a molecule
            // if neutral, only neutral molecules react
            auto selection = (neutral) ? Tspace::ACTIVE_NEUTRAL : Tspace::ACTIVE;
            auto molecule_list = spc.findMolecules(molid, selection);
            if (range_size(molecule_list) >= number_to_delete) {
                for (int i = 0; i < number_to_delete; i++) {
                    auto target = slump.sample(molecule_list.begin(), molecule_list.end());
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

/**
 * Expand a molecular group by `number_to_insert` particles by activating inactive
 * particles at the end of the group. The activated particles are assigned new
 * random positions, guaranteed to fall within the simulation box.
 *
 * If the capacity of the group will be exceeded, a warning is issued and
 * the returned Change::data object will be empty.
 */
Change::data SpeciationMove::expandAtomicGroup(Space::Tgroup &target, int number_to_insert) {
    assert(target.atomic);
    Change::data d;
    if (target.size() + number_to_insert <= target.capacity()) {
        d.index = &target - &spc.groups.front();
        d.internal = true;
        d.dNatomic = true;
        for (int i = 0; i < number_to_insert; i++) {
            target.activate(target.end(), target.end() + 1);
            auto last_atom = target.end() - 1;
            spc.geo.randompos(last_atom->pos, slump);
            spc.geo.getBoundaryFunc()(last_atom->pos);
            d.atoms.push_back(std::distance(target.begin(), last_atom)); // Index of particle rel. to group
        }
    } else {
        faunus_logger->warn("atomic group {} is full; you should probably expand it",
                            Faunus::molecules[target.id].name);
    }
    return d;
}

/*
 * Activate an inactive molecule and assign a new random position and orientation.
 * If the molecule has internal bonds, the bond-energy is calculated to ensure that
 * the bond-energy does not affect the insertion acceptance.
 */
Change::data SpeciationMove::activateMolecularGroup(Space::Tgroup &target) {
    assert(target.size() == 0); // must be inactive

    target.activate(target.inactive().begin(), target.inactive().end());
    Point cm = target.cm;
    spc.geo.randompos(cm, slump);
    target.translate(cm, spc.geo.getBoundaryFunc());
    Point u = ranunit(slump);
    Eigen::Quaterniond Q(Eigen::AngleAxisd(2 * pc::pi * (slump() - 0.5), u));
    target.rotate(Q, spc.geo.getBoundaryFunc());

    assert(spc.geo.sqdist(target.cm, Geometry::massCenter(target.begin(), target.end(), spc.geo.getBoundaryFunc(),
                                                          -target.cm)) < 1e-9);

    // We store the bonded energy of the activated molecule
    // The change in bonded energy should not affect the acceptance/rejection of the move
    for (auto &bond : Faunus::molecules[target.id].bonds) {
        auto bondclone = bond->clone();
        bondclone->shift(std::distance(spc.p.begin(), target.begin()));
        Potential::setBondEnergyFunction(bondclone, spc.p);
        bond_energy -= bondclone->energy(spc.geo.getDistanceFunc());
    }

    Change::data d;
    d.index = &target - &spc.groups.front(); // Integer *index* of moved group
    d.all = true;                            // All atoms in group were moved
    d.internal = true;
    d.atoms.reserve(target.capacity());
    for (size_t i = 0; i < target.capacity(); i++)
        d.atoms.push_back(i);
    return d;
}

void SpeciationMove::activateProducts(Change &change, std::vector<ReactionData>::iterator reaction) {
    for (auto [molid, number_to_insert] : reaction->moleculesToAdd(forward)) {
        const MoleculeData &molecule_data = molecules.at(molid);
        auto mollist = spc.findMolecules(molid, Tspace::ALL);

        if (molecule_data.atomic) { // The product is an atom
            Change::data d = expandAtomicGroup(*mollist.begin(), number_to_insert);
            change.groups.push_back(d);
        } else {           // The product is a molecule
            if (neutral) { // Only neutral molecules react
                mollist = spc.findMolecules(molid, Tspace::INACTIVE_NEUTRAL);
            } else {
                mollist = spc.findMolecules(molid, Tspace::INACTIVE);
            }
            for (int N = 0; N < number_to_insert; N++) {
                auto target = slump.sample(mollist.begin(), mollist.end());
                auto change_data = activateMolecularGroup(*target);
                change.groups.push_back(change_data); // Add to list of moved groups
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
        auto random_reaction_direction = static_cast<ReactionData::Direction>((char)slump.range(0, 1));
        // reaction->setDirection(random_reaction_direction);
        trialprocess = &(*reaction);

        if (reaction->empty(forward)) // Enforce canonic constraint if invoked
            return;              // Out of material, slip out the back door
        if (checkInsertProducts(*reaction) == false)
            return;
        if (swapReaction(change, *reaction) == false)
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
