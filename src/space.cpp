#include "space.h"
#include "io.h"
#include "aux/iteratorsupport.h"
#include "spdlog/spdlog.h"
#include <iostream>
#include "aux/eigensupport.h"

namespace Faunus {

bool Change::data::operator<(const Faunus::Change::data &a) const { return index < a.index; }

void Change::clear() {
    dV = false;
    all = false;
    dN = false;
    moved2moved = true;
    groups.clear();
    assert(empty());
}
bool Change::empty() const {
    if (dV || all || dN || !groups.empty()) {
        return false;
    } else {
        return true;
    }
}
Change::operator bool() const { return not empty(); }

std::vector<int> Change::touchedParticleIndex(const std::vector<Group<Particle>> &group_vector) {
    std::vector<int> atom_indexes;                                   // atom index rel. to first particle in system
    for (const auto &changed : groups) {                             // loop over changed groups
        auto begin_first = group_vector.front().begin();             // first particle, first group
        auto begin_current = group_vector.at(changed.index).begin(); // first particle, current group
        auto offset = std::distance(begin_first, begin_current);     // abs. distance from first particle
        atom_indexes.reserve(atom_indexes.size() + changed.atoms.size());
        for (auto index : changed.atoms) {          // atom index relative to group
            atom_indexes.push_back(index + offset); // atom index relative to first
        }
    }
    return atom_indexes;
}

bool Change::sanityCheck(Space &spc) const {
    // make sure "atoms" belong to the "index"th group
    bool rc = true;
    for (auto &d : groups) {
        auto &g = spc.groups.at(d.index);
        for (size_t atom_index : d.atoms)     // all atoms must be within `g`
            if (atom_index >= g.capacity()) { // atom_index is w. respect to group.begin()
                size_t first = std::distance(spc.p.begin(), g.begin());
                size_t last = std::distance(spc.p.begin(), g.trueend()) - 1;
                faunus_logger->debug("atom {} is outside capacity of group '{}' ({}-{})", atom_index + first,
                                     molecules.at(g.id).name, first, last);
                rc = false;
            }
    }
    return rc;
}

void to_json(json &j, const Change::data &d) {
    j = {{"all", d.all},           {"internal", d.internal}, {"dNswap", d.dNswap},
         {"dNatomic", d.dNatomic}, {"index", d.index},       {"atoms", d.atoms}};
}

void to_json(json &j, const Change &c) {
    j = {{"dV", c.dV}, {"all", c.all}, {"dN", c.dN}, {"moved2moved", c.moved2moved}, {"groups", c.groups}};
}

void Space::clear() {
    p.clear();
    groups.clear();
}

/**
 * The following is considered:
 *
 * - `groups` vector is expanded with a new group at the end
 * - if the space particle vector is relocated, all existing group
 *   iterators are updated to reflect the new memory positions
 * - for molecular groups, the mass center is calculated and set
 *
 * @param molid Molecule id of inserted data
 * @param particles Particles to generate group from
 */
void Space::push_back(int molid, const ParticleVector &particles) {
    if (!particles.empty()) {
        auto original_begin = p.begin();                       // used to detect if `p` is relocated
        p.insert(p.end(), particles.begin(), particles.end()); // insert particle into space
        if (p.begin() != original_begin) {                     // update group iterators if `p` is relocated
            std::for_each(groups.begin(), groups.end(),
                          [&](auto &group) { group.relocate(original_begin, p.begin()); });
        }
        Tgroup group(p.end() - particles.size(), p.end()); // create a group
        const auto &moldata = Faunus::molecules.at(molid);
        group.id = molid;
        group.compressible = moldata.compressible;
        group.atomic = moldata.atomic;
        if (group.isMolecular()) {
            if (particles.size() != moldata.atoms.size()) {
                faunus_logger->error("{} require {} atoms but {} were provided", moldata.name, moldata.atoms.size(),
                                     particles.size());
                throw std::runtime_error("particle size mismatch");
            }
            group.cm = Geometry::massCenter(group.begin(), group.end(), geo.getBoundaryFunc(), -particles.begin()->pos);
        } else {
            if (particles.size() % moldata.atoms.size() != 0) {
                throw std::runtime_error("indivisible by atomic group size: "s + moldata.name);
            }
        }
        groups.push_back(group);
    }
}

/**
 * @param other Space to copy from
 * @param change Change object describing the changes beteeen the two Space objects
 *
 * Copy data from another Space according to Change object. The other space *must* be populated
 * in the exact same way, i.e. must have the same molecules and particles. In DEBUG mode, several
 * assertions are included to ensure this is true. Copied data includes:
 *
 * - geometry
 * - groups
 * - particles
 * - implicit molecules
 */
void Space::sync(const Space &other, const Change &change) {
    if (&other != this && !change.empty()) {
        assert(!groups.empty());
        assert(p.size() == other.p.size());
        assert(groups.size() == other.groups.size());
        assert(implicit_reservoir.size() == other.implicit_reservoir.size());
        if (change.dV or change.all) {
            geo = other.geo; // copy simulation geometry
        }
        if (change.all) {                                                   // deep copy *everything*
            implicit_reservoir = other.implicit_reservoir;                  // copy all implicit molecules
            p = other.p;                                                    // copy all positions
            groups = other.groups;                                          // copy all groups
            assert(p.begin() != other.p.begin());                           // check deep copy problem
            assert(groups.front().begin() != other.groups.front().begin()); // check deep copy problem
        } else {
            for (const auto &changed : change.groups) {             // look over changed groups
                auto &group = groups.at(changed.index);             // old group
                auto &other_group = other.groups.at(changed.index); // new group
                assert(group.id == other_group.id);
                if (group.traits().isImplicit()) { // the molecule is implicit
                    implicit_reservoir[group.id] = other.implicit_reservoir.at(group.id);
                } else if (changed.all) {
                    group = other_group;            // copy everything
                } else {                            // copy only a subset
                    group.shallowcopy(other_group); // copy group data but *not* particles
                    for (auto i : changed.atoms) {  // loop over atom index (rel. to group)
                        group[i] = other_group[i];  // deep copy select particles
                    }
                }
            }
        }
    }
}

Point Space::scaleVolume(double Vnew, Geometry::VolumeMethod method) {
    for (auto &g : groups) // remove periodic boundaries
        if (not g.atomic)
            g.unwrap(geo.getDistanceFunc());

    Point scale = geo.setVolume(Vnew, method);

    for (auto &g : groups) {
        if (not g.empty()) {
            if (g.isAtomic()) { // scale all atoms
                for (auto &i : g) {
                    i.pos = i.pos.cwiseProduct(scale);
                    geo.boundary(i.pos);
                }
            } else { // scale mass center and translate
                Point oldcm = g.cm;
                if (g.compressible) {
                    for (auto &i : g) {
                        i.pos = i.pos.cwiseProduct(scale);
                        geo.boundary(i.pos);
                    }
                    g.cm = Geometry::massCenter(g.begin(), g.end(), geo.getBoundaryFunc(), -oldcm);
                    geo.boundary(g.cm);
                } else {
                    g.cm = g.cm.cwiseProduct(scale);
                    geo.boundary(g.cm);
                    Point delta = g.cm - oldcm;
                    for (auto &i : g) {
                        i.pos += delta;
                        geo.boundary(i.pos);
                    }
#ifndef NDEBUG
                    Point recalc_cm = Geometry::massCenter(g.begin(), g.end(), geo.getBoundaryFunc(), -g.cm);
                    double cm_error = std::fabs(geo.sqdist(g.cm, recalc_cm));
                    if (cm_error > 1e-6) {
                        std::ostringstream o;
                               o  << "error: " << cm_error << std::endl
                                  << "scale: " << scale.transpose() << std::endl
                                  << "delta: " << delta.transpose() << " norm = " << delta.norm() << std::endl
                                  << "|o-n|: " << geo.vdist(oldcm, g.cm).norm() << std::endl
                                  << "oldcm: " << oldcm.transpose() << std::endl
                                  << "newcm: " << g.cm.transpose() << std::endl
                                  << "actual cm: " << recalc_cm.transpose() << std::endl;
                        faunus_logger->error(o.str());
                        assert(false);
                    }
#endif
                }
            }
        }
    }
    double Vold = geo.getVolume();
    // if isochoric, the volume is constant
    if (method == Geometry::ISOCHORIC)
        Vold = std::pow(Vold, 1. / 3.);

    for (auto f : scaleVolumeTriggers)
        f(*this, Vold, Vnew);

    return scale;
}

json Space::info() {
    json j = {{"number of particles", p.size()}, {"number of groups", groups.size()}, {"geometry", geo}};
    auto &_j = j["groups"];
    for (auto &i : groups) {
        auto &name = molecules.at(i.id).name;
        json tmp, d = i;
        d.erase("cm");
        d.erase("id");
        d.erase("atomic");
        auto ndx = i.to_index(p.begin());
        if (not i.empty())
            d["index"] = {ndx.first, ndx.second};
        // d["index"] = std::to_string(ndx.first)+"-"+std::to_string(ndx.second);
        tmp[name] = d;
        _j.push_back(tmp);
    }
    auto &_j2 = j["reactionlist"];
    for (auto &i : Faunus::reactions) {
        json tmp, d = i;
        _j2.push_back(i);
    }
    return j;
}

Space::Tgvec::iterator Space::randomMolecule(int molid, Random &rand, Space::Selection sel) {
    auto m = findMolecules(molid, sel);
    if (not ranges::cpp20::empty(m))
        return groups.begin() + (&*rand.sample(m.begin(), m.end()) - &*groups.begin());
    return groups.end();
}
const std::map<int, int> &Space::getImplicitReservoir() const { return implicit_reservoir; }
std::map<int, int> &Space::getImplicitReservoir() { return implicit_reservoir; }

void to_json(json &j, Space &spc) {
    typedef typename Space::Tpvec Tpvec;
    j["geometry"] = spc.geo;
    j["groups"] = spc.groups;
    j["particles"] = spc.p;
    j["reactionlist"] = reactions;
    j["implicit_reservoir"] = spc.getImplicitReservoir();
}
void from_json(const json &j, Space &spc) {
    typedef typename Space::Tpvec Tpvec;
    using namespace std::string_literals;

    try {
        if (atoms.empty())
            atoms = j.at("atomlist").get<decltype(atoms)>();
        if (molecules.empty())
            molecules = j.at("moleculelist").get<decltype(molecules)>();
        if (reactions.empty())
            if (j.count("reactionlist") > 0)
                reactions = j.at("reactionlist").get<decltype(reactions)>();

        spc.clear();
        spc.geo = j.at("geometry");

        if (j.count("groups") == 0) {
            InsertMoleculesInSpace::insertMolecules(j.at("insertmolecules"), spc);
        } else {
            spc.p = j.at("particles").get<Tpvec>();
            if (!spc.p.empty()) {
                auto begin = spc.p.begin();
                Space::Tgroup g(begin, begin);
                for (auto &i : j.at("groups")) {
                    g.begin() = begin;
                    from_json(i, g);
                    spc.groups.push_back(g);
                    begin = g.trueend();
                }
                if (begin != spc.p.end())
                    throw std::runtime_error("load error");
            }
        }

        if (auto it = j.find("implicit_reservoir"); it != j.end()) {
            assert(it->is_array());
            for (auto vec : *it) {
                assert(vec.is_array() && vec.size() == 2);
                spc.getImplicitReservoir()[vec[0]] = vec[1];
            }
            faunus_logger->trace("{} implicit molecules loaded from json", it->size());
        }

        // check correctness of molecular mass centers
        for (auto &group : spc.groups) {
            if (!group.empty() && group.isMolecular()) {
                if (spc.geo.sqdist(group.cm, Geometry::massCenter(group.begin(), group.end(), spc.geo.getBoundaryFunc(),
                                                                  -group.cm)) > 1e-9) {
                    throw std::runtime_error("mass center mismatch");
                }
            }
        }
    } catch (std::exception &e) {
        throw std::runtime_error("error while constructing Space from JSON: "s + e.what());
    }
}
getActiveParticles::const_iterator::const_iterator(const Space &spc,
                                                   getActiveParticles::const_iterator::Tparticle_iter it)
    : spc(spc), particle_iter(it) {
    groups_iter = spc.groups.begin();
}
getActiveParticles::const_iterator getActiveParticles::const_iterator::operator++() { // advance particles and groups
    if (++particle_iter == groups_iter->end()) {
        do {
            if (++groups_iter == spc.groups.end())
                return *this;
        } while (groups_iter->empty());
        particle_iter = groups_iter->begin();
    }
    return *this;
}
getActiveParticles::const_iterator getActiveParticles::begin() const { return const_iterator(spc, spc.p.begin()); }
getActiveParticles::const_iterator getActiveParticles::end() const {
    return spc.groups.empty() ? begin() : const_iterator(spc, spc.groups.back().end());
}
size_t getActiveParticles::size() const {
    return std::accumulate(spc.groups.begin(), spc.groups.end(), 0,
                           [](size_t sum, const auto &g) { return sum + g.size(); });
}
getActiveParticles::getActiveParticles(const Space &spc) : spc(spc) {}

namespace SpaceFactory {

/**
 * @param space Space to insert into (will be overwritten)
 * @param num_particles Number of salt pairs to insert
 * @param geometry Geometry to use
 *
 * Create a system with two atom types, "Na" and "Cl", forming
 * an atomic molecule, "salt". N salt pairs a randomly inserted
 */
void makeNaCl(Space &space, int num_particles, const Geometry::Chameleon &geometry) {
    pc::temperature = 298.15_K;
    space.geo = geometry;

    Faunus::atoms = R"([
             { "Na": { "sigma": 3.8, "eps": 0.1, "q": 1.0 } },
             { "Cl": { "sigma": 4.0, "eps": 0.05, "q": -1.0 } }
             ])"_json.get<decltype(atoms)>();

    Faunus::molecules = R"([
                { "salt": {"atomic": true, "atoms": ["Na", "Cl"] } }
            ])"_json.get<decltype(molecules)>();

    json j = json::array();
    j.push_back({{"salt", {{"N", num_particles}}}});
    InsertMoleculesInSpace::insertMolecules(j, space);
}

} // namespace SpaceFactory

/**
 * @param moldata Molecule type to insert. Must be `atomic`.
 * @param spc Space to insert into (at the end)
 * @param num_molecules Number of molecules to insert into the *same* group
 * @param num_inactive_molecules Number of molecules to declare inactive
 * @throws If number if inactive is higher than number of active molecules
 *
 * The atoms in the atomic groups are repeated `num_molecules` times, then inserted
 * into Space.
 */
void InsertMoleculesInSpace::insertAtomicGroups(MoleculeData &moldata, Space &spc, int num_molecules,
                                                int num_inactive_molecules) {
    assert(moldata.atomic);
    ParticleVector repeated_particles;
    repeated_particles.reserve(num_molecules * moldata.atoms.size()); // prepare memory
    for (size_t i = 0; i < num_molecules; i++) {                      // repeat insertion into the same atomic group
        auto particles = moldata.getRandomConformation(spc.geo, spc.p);
        repeated_particles.insert(repeated_particles.end(), particles.begin(), particles.end());
    }

    spc.push_back(moldata.id(), repeated_particles); // create new group in Space

    if (num_inactive_molecules > num_molecules) {
        throw std::runtime_error("too many inactive molecules requested");
    } else if (num_inactive_molecules > 0) {
        int num_active_atoms = (num_molecules - num_inactive_molecules) * moldata.atoms.size();
        spc.groups.back().resize(num_active_atoms); // deactivate atoms in molecule
    }
}

void InsertMoleculesInSpace::insertMolecularGroups(MoleculeData &moldata, Space &spc, int num_molecules,
                                                   int num_inactive) {
    assert(moldata.atomic == false);
    for (size_t i = 0; i < num_molecules; i++) { // insert molecules
        spc.push_back(moldata.id(), moldata.getRandomConformation(spc.geo, spc.p));
    }
    if (num_inactive > num_molecules) {
        throw std::runtime_error("too many inactive molecules requested");
    } else { // deactivate groups, starting from the back
        std::for_each(spc.groups.rbegin(), spc.groups.rbegin() + num_inactive, [&](auto &group) {
            group.unwrap(spc.geo.getDistanceFunc()); // make molecules whole (remove PBC) ...
            group.resize(0);                         // ... and then deactivate
        });
    }
}

/**
 * @brief Set positions for num_molecules last groups in Space
 * @param spc Space to insert into
 * @param int num_molecules Number of groups set affect
 * @param particles Position vector for all particles in the N groups
 * @param offset Translate positions by this offset
 * @throws if num_molecules doesn't match, or if positions are outside simulation cell
 *
 * Sets particle positions in space using a given input
 * particle vector. The vector can span several *identical* groups.
 * The number of groups affected must be given in order to update their mass-centers.
 * Only *positions* are affected.
 */
void InsertMoleculesInSpace::setPositionsForTrailingGroups(Space &spc, int num_molecules,
                                                           const Faunus::ParticleVector &particles,
                                                           const Point &offset) {
    assert(spc.groups.size() >= num_molecules);
    if (particles.size() != num_molecules * (spc.groups.rbegin())->traits().atoms.size()) {
        throw std::runtime_error("number of particles doesn't match groups");
    } else {
        // update positions in space, starting from the back
        std::transform(particles.rbegin(), particles.rend(), spc.p.rbegin(), spc.p.rbegin(), [&](auto &src, auto &dst) {
            dst.pos = src.pos + offset;       // shift by offset
            if (spc.geo.collision(dst.pos)) { // check if position is inside simulation volume
                throw std::runtime_error("positions outside box");
            }
            return dst;
        });
        // update mass-centers on modified groups
        std::for_each(spc.groups.rbegin(), spc.groups.rbegin() + num_molecules, [&](auto &g) {
            g.cm = Geometry::massCenter(g.begin(), g.end(), spc.geo.getBoundaryFunc(), -g.begin()->pos);
        });
    }
}

/**
 * @brief Insert implicit groups into Space
 * @param moldata Molecule to insert. Must be implicit.
 * @param spc Space to insert into.
 * @param num_molecules Nunber of implicit molecules to insert
 */
void InsertMoleculesInSpace::insertImplicitGroups(const MoleculeData &moldata, Space &spc, int num_molecules) {
    assert(moldata.isImplicit());
    spc.getImplicitReservoir()[moldata.id()] = num_molecules;
}

/**
 * @brief Insert molecules into Space based on JSON input
 * @param json_array JSON array
 * @param spc Space to insert into
 */
void InsertMoleculesInSpace::insertMolecules(const json &json_array, Space &spc) {
    spc.clear();
    assert(spc.geo.getVolume() > 0);
    if (!json_array.is_array()) {
        throw ConfigurationError("syntax error in insertmolecule");
    }
    for (auto &obj : json_array) { // loop over array of molecules
        if (!obj.is_object() || obj.size() != 1) {
            throw ConfigurationError("syntax error in insertmolecule");
        }
        for (auto &[molname, properties] : obj.items()) {
            if (auto moldata = findName(Faunus::molecules, molname); moldata != Faunus::molecules.end()) {
                int num_molecules = 0; // number of groups to insert
                if (auto it = properties.find("N"); it != properties.end()) {
                    num_molecules = it->get<int>();
                } else {
                    double concentration = properties.at("molarity").get<double>() * 1.0_molar;
                    num_molecules = std::round(concentration * spc.geo.getVolume());
                    if (concentration > pc::epsilon_dbl) {
                        double rel_error = (concentration - num_molecules / spc.geo.getVolume()) / concentration;
                        if (rel_error > 0.01) {
                            faunus_logger->warn("{}: initial molarity differs by {}% from target value", molname,
                                                rel_error * 100);
                        }
                    }
                }
                if (not moldata->isImplicit() and num_molecules < 1) {
                    throw ConfigurationError(molname + ": at least one molecule required. Concentration too low?");
                }
                int num_inactive = getNumberOfInactiveMolecules(properties, num_molecules);

                double molar_concentration = (num_molecules - num_inactive) / spc.geo.getVolume() / 1.0_molar;
                if (moldata->isImplicit()) {
                    faunus_logger->info("adding {} implicit {} molecules --> {} mol/l", num_molecules, molname,
                                        molar_concentration);
                } else {
                    faunus_logger->info("adding {} {} molecules --> {} mol/l ({} inactive)", num_molecules, molname,
                                        molar_concentration, num_inactive);
                }

                if (moldata->atomic) {
                    insertAtomicGroups(*moldata, spc, num_molecules, num_inactive);
                } else if (moldata->isImplicit()) {
                    insertImplicitGroups(*moldata, spc, num_molecules);
                } else {
                    insertMolecularGroups(*moldata, spc, num_molecules, num_inactive);
                    if (auto filename = properties.value("positions", ""s); !filename.empty()) {
                        Space::Tpvec particles; // positions loaded from file
                        if (loadStructure(filename, particles, false)) {
                            faunus_logger->info("{}: loaded position file {}", molname, filename);
                            Point offset = properties.value("translate", Point(0, 0, 0));
                            setPositionsForTrailingGroups(spc, num_molecules, particles, offset);
                        } else {
                            throw ConfigurationError("error loading positions from '" + filename + "'");
                        }
                    }
                }
            } else {
                throw ConfigurationError("cannot insert undefined molecule '" + molname + "'");
            }
        }
    }
}
/**
 * @param j Input json object
 * @param number_of_molecules Total number of molecules
 * @return Number of molecules to be inactive (always smaller than `number_of_molecules`)
 * @throws If Inactive molecules is higher than `number_of_molecules`
 *
 * Looks for key "inactive" and if:
 * - boolean true: all molecules are inactive
 * - number: number of molecules to declare inactive
 * - no `inactive` key found, return zero
 */
int InsertMoleculesInSpace::getNumberOfInactiveMolecules(const json &j, int number_of_molecules) {
    int number_of_inactive_molecules = 0; // number of inactive molecules
    if (auto it = j.find("inactive"); it != j.end()) {
        if (it->is_boolean()) {
            if (*it) {
                number_of_inactive_molecules = number_of_molecules; // all molecules are inactive
            }
        } else if (it->is_number_integer()) {
            number_of_inactive_molecules = *it; // a subset are inactive
            if (number_of_inactive_molecules > number_of_molecules) {
                throw ConfigurationError("too many inactive particles requested");
            }
        }
    }
    return number_of_inactive_molecules;
}

} // namespace Faunus