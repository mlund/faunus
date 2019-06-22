#pragma once

#include <set>
#include "core.h"
#include "particle.h"
#include "random.h"

namespace Faunus {

namespace Potential {
struct BondData;
}

namespace Geometry {
struct GeometryBase;
}

class MoleculeData;

/**
 * @brief Random position and orientation - typical for rigid bodies
 *
 * Molecule inserters take care of generating molecules
 * for insertion into space and can be used in Grand Canonical moves,
 * Widom analysis, and for generating initial configurations.
 * Inserters will not actually insert anything, but rather
 * return a particle vector with proposed coordinates.
 *
 * All inserters are function objects, expecting
 * a geometry, particle vector, and molecule data.
 */
struct RandomInserter {
    Point dir = {1, 1, 1};     //!< Scalars for random mass center position. Default (1,1,1)
    Point offset = {0, 0, 0};  //!< Added to random position. Default (0,0,0)
    bool rotate = true;        //!< Set to true to randomly rotate molecule when inserted. Default: true
    bool keeppos = false;      //!< Set to true to keep original positions (default: false)
    bool allowoverlap = false; //!< Set to true to skip container overlap check
    int maxtrials = 2e4;       //!< Maximum number of container overlap checks
    int confindex = -1;        //!< Index of last used conformation

    ParticleVector operator()(Geometry::GeometryBase &geo, const ParticleVector &, MoleculeData &mol);
};

/**
 * Possible structure for molecular conformations
 */
struct Conformation {
    std::vector<Point> positions;
    std::vector<double> charges;

    bool empty() const;

    ParticleVector &toParticleVector(ParticleVector &p) const; // copy conformation into particle vector
};

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("[Faunus] Conformation") {
    ParticleVector p(1);
    Conformation c;
    CHECK(c.empty());

    c.positions.push_back({1, 2, 3});
    c.charges.push_back(0.5);
    CHECK(not c.empty());

    c.toParticleVector(p);

    CHECK(p[0].pos == Point(1, 2, 3));
    CHECK(p[0].charge == 0.5);
}
#endif

/**
 * @brief General properties for molecules
 */
class MoleculeData {
  private:
    int _id = -1;

  public:
    typedef std::function<ParticleVector(Geometry::GeometryBase &, const ParticleVector &, MoleculeData &)>
        TinserterFunc;

    TinserterFunc inserterFunctor = nullptr; //!< Function for insertion into space

    int &id();             //!< Type id
    const int &id() const; //!< Type id
    void createMolecularConformations(SingleUseJSON &); //!< Add conformations if appropriate

    std::string name;            //!< Molecule name
    std::string structure;       //!< Structure file (pqr|aam|xyz)
    bool atomic = false;         //!< True if atomic group (salt etc.)
    bool rotate = true;          //!< True if molecule should be rotated upon insertion
    bool keeppos = false;        //!< Keep original positions of `structure`
    bool keepcharges = true;     //!< Set to true to keep charges in PQR file (default: true)
    bool rigid = false;          //!< True if particle should be considered as rigid
    double activity = 0;         //!< Chemical activity (mol/l)
    Point insdir = {1, 1, 1};    //!< Insertion directions
    Point insoffset = {0, 0, 0}; //!< Insertion offset

    std::vector<std::shared_ptr<Potential::BondData>> bonds;
    std::vector<int> atoms;                    //!< Sequence of atoms in molecule (atom id's)
    WeightedDistribution<ParticleVector> conformations; //!< Conformations of molecule

    MoleculeData();

    /** @brief Specify function to be used when inserting into space.
     *
     * By default a random position and orientation is generator and overlap
     * with container is avoided.
     */
    void setInserter(const TinserterFunc &ifunc);

    /**
     * @brief Get random conformation that fits in container
     * @param geo Geometry
     * @param otherparticles Typically `spc.p` is insertion depends on other particle
     *
     * By default the molecule is placed at a random position and orientation with
     * no container overlap using the `RandomInserter` class. This behavior can
     * be changed by specifying another inserter using `setInserter()`.
     */
    ParticleVector getRandomConformation(Geometry::GeometryBase &geo, ParticleVector otherparticles = ParticleVector());

    void loadConformation(const std::string &file, bool keepcharges);
}; // end of class

void to_json(json &j, const MoleculeData &a);

void from_json(const json &j, MoleculeData &a);

void from_json(const json &j, std::vector<MoleculeData> &v);

// global instance of molecule vector
extern std::vector<MoleculeData> molecules;

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("[Faunus] MoleculeData") {
    using doctest::Approx;

    json j = R"(
            { "moleculelist": [
                { "B": {"activity":0.2, "atomic":true, "insdir": [0.5,0,0], "insoffset": [-1.1, 0.5, 10], "atoms":["A"] } },
                { "A": { "atomic":false } }
            ]})"_json;

    molecules = j["moleculelist"].get<decltype(molecules)>(); // fill global instance
    auto &v = molecules;                                      // reference to global molecule vector

    CHECK(v.size() == 2);
    CHECK(v.back().id() == 1);
    CHECK(v.back().name == "A"); // alphabetic order in std::map
    CHECK(v.back().atomic == false);

    MoleculeData m = json(v.front()); // moldata --> json --> moldata

    CHECK(m.name == "B");
    CHECK(m.id() == 0);
    CHECK(m.activity == Approx(0.2_molar));
    CHECK(m.atomic == true);
    CHECK(m.insdir == Point(0.5, 0, 0));
    CHECK(m.insoffset == Point(-1.1, 0.5, 10));
}
#endif

/*
 * @brief General properties of reactions
 */
class ReactionData {
  public:
    typedef std::map<int, int> Tmap;

    std::vector<std::string> _reag, _prod;

    Tmap _reagid_m; // Molecular change, groups. Atomic as Groupwise
    Tmap _reagid_a; // Atomic change, equivalent of swap/titration
    Tmap _prodid_m;
    Tmap _prodid_a;
    // Tmap _Reac, _Prod;

    bool canonic = false; //!< Finite reservoir
    bool swap = false;    //!< True if swap move
    int N_reservoir;      //!< Number of molecules in finite reservoir
    double lnK = 0;       //!< Natural logarithm of molar eq. const.
    double pK = 0;        //!< -log10 of molar eq. const.
    std::string name;     //!< Name of reaction
    std::string formula;  //!< Chemical formula
    double weight;        //!< Statistical weight to be given to reaction in speciation

    bool empty(bool forward) const;

    std::vector<int> participatingMolecules() const; //!< Returns molids of participating molecules

    bool containsMolecule(int molid) const; //!< True of molecule id is part of process

    const Tmap &Molecules2Add(bool forward) const; //!< Map for addition depending on direction

    const Tmap &Atoms2Add(bool forward) const; //!< Map for addition depending on direction

    auto findAtomOrMolecule(const std::string &name) const {
        auto it_a = findName(Faunus::atoms, name);
        auto it_m = findName(Faunus::molecules, name);
        if (it_m == Faunus::molecules.end())
            if (it_a == Faunus::atoms.end())
                throw std::runtime_error("unknown species '" + name + "'");
        return std::make_pair(it_a, it_m);
    } //!< Returns pair of iterators to atomlist and moleculelist. One of them points to end().

}; //!< End of class

inline auto parseProcess(const std::string &process) {
    typedef std::vector<std::string> Tvec;
    Tvec v;
    std::string tmp;
    std::istringstream iss(process);
    while (iss >> tmp)
        v.push_back(tmp);
    v.erase(std::remove(v.begin(), v.end(), "+"), v.end());
    auto it = std::find(v.begin(), v.end(), "=");
    if (it == v.end())
        throw std::runtime_error("products and reactants must be separated by '='");
    return std::make_pair(Tvec(v.begin(), it), Tvec(it + 1, v.end()));
} //!< Parse process string to pair of vectors containing reactant/product species

void from_json(const json &j, ReactionData &a);

void to_json(json &j, const ReactionData &a);

extern std::vector<ReactionData> reactions; // global instance

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("[Faunus] ReactionData") {
    using doctest::Approx;

    json j = R"(
            {
                "atomlist" :
                    [ {"a": { "r":1.1 } } ],
                "moleculelist": [
                    { "A": { "atomic":false, "activity":0.2 } },
                    { "B": { "atomic":true, "atoms":["a"] } }
                ],
                "reactionlist": [
                    {"A = B": {"lnK":-10.051, "canonic":true, "N":100 } }
                ]
            } )"_json;

    Faunus::atoms = j["atomlist"].get<decltype(atoms)>();
    molecules = j["moleculelist"].get<decltype(molecules)>(); // fill global instance

    auto &r = reactions; // reference to global reaction list
    r = j["reactionlist"].get<decltype(reactions)>();

    CHECK(r.size() == 1);
    CHECK(r.front().name == "A = B");
    CHECK(r.front().lnK == Approx(-10.051 - std::log(0.2)));
}
#endif

} // namespace Faunus
