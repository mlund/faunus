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
    if (dV == false)
        if (all == false)
            if (groups.empty())
                if (dN == false)
                    return true;
    return false;
}
Change::operator bool() const { return not empty(); }

std::vector<int> Change::touchedParticleIndex(const std::vector<Group<Particle>> &group_vector) {
    std::vector<int> atom_indexes; // atom index rel. to first particle in system
    auto begin = group_vector.front().begin();
    for (auto &d : groups)     // loop over changes groups
        for (auto i : d.atoms) // atom index rel. to group
            atom_indexes.push_back(i + std::distance(begin, group_vector.at(d.index).begin()));
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

void Space::push_back(int molid, const Space::Tpvec &in) {
    if (!in.empty()) {
        auto oldbegin = p.begin();
        p.insert(p.end(), in.begin(), in.end());
        if (p.begin() != oldbegin) { // update group iterators if `p` is relocated
            for (auto &g : groups) {
                g.relocate(oldbegin, p.begin());
            }
        }
        Tgroup g(p.end() - in.size(), p.end());
        g.id = molid;
        g.compressible = molecules.at(molid).compressible;
        g.atomic = molecules.at(molid).atomic;

        if (g.atomic == false) {
            g.cm = Geometry::massCenter(in.begin(), in.end(), geo.getBoundaryFunc(), -in.begin()->pos);
            Point cm = Geometry::massCenter(g.begin(), g.end(), geo.getBoundaryFunc(), -g.cm);
            if (geo.sqdist(g.cm, cm) > 1e-6)
                throw std::runtime_error("space: mass center error upon insertion. "s + molecules.at(molid).name +
                                         " molecule too large?\n");
        }

        groups.push_back(g);
        assert(groups.back().begin() == g.begin());
        assert(in.size() == groups.back().capacity());
    }
}

void Space::sync(const Space &other, const Change &change) {
    assert(&other != this);
    assert(p.begin() != other.p.begin());

    if (change.dV or change.all)
        geo = other.geo;

    // deep copy *everything*
    if (change.all) {
        p = other.p; // copy all positions
        assert(p.begin() != other.p.begin() && "deep copy problem");
        groups = other.groups;
        implicit_reservoir = other.implicit_reservoir;

        if (not groups.empty())
            if (groups.front().begin() == other.p.begin())
                for (auto &i : groups)
                    i.relocate(other.p.begin(), p.begin());
    } else {
        for (auto &m : change.groups) {

            auto &g = groups.at(m.index);            // old group
            auto &gother = other.groups.at(m.index); // new group
            assert(gother.id == g.id);

            if (Faunus::molecules[gother.id].isImplicit()) {
                assert(other.implicit_reservoir.count(gother.id));
                implicit_reservoir[g.id] = other.implicit_reservoir.at(gother.id);
            }

            g.shallowcopy(gother); // copy group data but *not* particles

            if (m.all) // copy all particles
                std::copy(gother.begin(), gother.trueend(), g.begin());
            else // copy only a subset
                for (auto i : m.atoms)
                    *(g.begin() + i) = *(gother.begin() + i);
        }
    }
    assert(p.size() == other.p.size());
    assert(p.begin() != other.p.begin());
}

Point Space::scaleVolume(double Vnew, Geometry::VolumeMethod method) {
    for (auto &g : groups) // remove periodic boundaries
        if (not g.atomic)
            g.unwrap(geo.getDistanceFunc());

    Point scale = geo.setVolume(Vnew, method);

    for (auto &g : groups) {
        if (not g.empty()) {
            if (g.atomic) { // scale all atoms
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
        for (auto &i : spc.groups)
            if (not i.empty() and not i.atomic)
                if (spc.geo.sqdist(i.cm, Geometry::massCenter(i.begin(), i.end(), spc.geo.getBoundaryFunc(), -i.cm)) >
                    1e-9)
                    throw std::runtime_error("mass center mismatch");
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
getActiveParticles::const_iterator getActiveParticles::end() const { return spc.groups.empty() ? begin() : const_iterator(spc, spc.groups.back().end()); }
size_t getActiveParticles::size() const {
    return std::accumulate(spc.groups.begin(), spc.groups.end(), 0,
                           [](size_t sum, const auto &g) { return sum + g.size(); });
}
getActiveParticles::getActiveParticles(const Space &spc) : spc(spc){}

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
 * @param inactive Set to true if the molecule should be *inactive*
 */
void InsertMoleculesInSpace::insertAtomicGroups(MoleculeData &moldata, Space &spc, int num_molecules, bool inactive) {
    assert(moldata.atomic == true);
    typename Space::Tpvec p;
    p.reserve(num_molecules * moldata.atoms.size()); // prepare memory
    while (num_molecules-- > 0) {                    // repeat insertion into the same atomic group
        auto particles = moldata.getRandomConformation(spc.geo, spc.p);
        p.insert(p.end(), particles.begin(), particles.end());
    }
    spc.push_back(moldata.id(), p);
    if (inactive) {
        spc.groups.back().resize(0);
    }
}

void InsertMoleculesInSpace::insertMolecularGroups(MoleculeData &moldata, Space &spc, int num_molecules,
                                                   bool inactive) {
    assert(moldata.atomic == false);
    while (num_molecules-- > 0) { // insert molecules
        spc.push_back(moldata.id(), moldata.getRandomConformation(spc.geo, spc.p));
        if (inactive) {
            spc.groups.back().unwrap(spc.geo.getDistanceFunc());
            spc.groups.back().resize(0);
        }
    }
}

/**
 * @brief Set positions for num_molecules last groups in Space
 * @param spc Space to insert into
 * @param int num_molecules Number of groups set affect
 * @param particles Position vector for all particles in the N groups
 * @param offset Translate positions by this offset
 */
void InsertMoleculesInSpace::setPositionsForTrailingGroups(Space &spc, int num_molecules,
                                                           const Faunus::ParticleVector &particles,
                                                           const Point &offset) {
    assert(spc.groups.size() >= num_molecules);
    if (particles.size() == num_molecules * (spc.groups.end() - num_molecules)->traits().atoms.size()) {
        auto target = spc.p.end() - particles.size(); // iterator to first particle to modify
        assert(target >= spc.p.begin());
        for (auto &i : particles) {       // loop over particles
            Point pos = i.pos + offset;   // shift by offset
            if (spc.geo.collision(pos)) { // check if position is inside simulation volume
                throw std::runtime_error("group positions outside box");
            }
            target->pos = pos; // copy position to space
            target++;          // advance to next target particle in space
        }
        // update mass-centers on modified groups
        for (auto it = spc.groups.end() - num_molecules; it != spc.groups.end(); ++it) {
            it->cm = Geometry::massCenter(it->begin(), it->end(), spc.geo.getBoundaryFunc(), -it->begin()->pos);
        }
    } else {
        throw std::runtime_error("mismatch in number of atoms");
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
                bool inactive = properties.value("inactive", false); // active or not?

                {
                    std::string state;
                    if (moldata->isImplicit()) {
                        state = "implicit";
                    } else {
                        state = (inactive) ? "inactive" : "active";
                    }
                    faunus_logger->info("adding {} {} {} molecules --> {} mol/l", num_molecules, state, molname,
                                        num_molecules / spc.geo.getVolume() / 1.0_molar);
                }

                if (moldata->atomic) {
                    insertAtomicGroups(*moldata, spc, num_molecules, inactive);
                } else if (moldata->isImplicit()) {
                    insertImplicitGroups(*moldata, spc, num_molecules);
                } else {
                    insertMolecularGroups(*moldata, spc, num_molecules, inactive);
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

} // namespace Faunus
