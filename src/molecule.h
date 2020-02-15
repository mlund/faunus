#pragma once

#include <set>
#include "core.h"
#include "auxiliary.h"
#include "particle.h"
#include "random.h"
#include "spdlog/spdlog.h"

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
 * Molecule inserters take care of generating molecules for insertion into
 * space and can be used in Grand Canonical moves, Widom analysis, and for
 * generating initial configurations. Inserters will not actually insert
 * anything, but rather return a particle vector with proposed coordinates.
 *
 * All inserters are function objects, expecting a geometry, particle vector,
 * and molecule data.
 */
struct MoleculeInserter {
    virtual ParticleVector operator()(Geometry::GeometryBase &geo, const ParticleVector &, MoleculeData &mol) = 0;
    virtual void from_json(const json&) {};
    virtual void to_json(json&) const {};
    virtual ~MoleculeInserter() = default;
};

void from_json(const json &j, MoleculeInserter &inserter);
void to_json(json &j, const MoleculeInserter &inserter);

struct RandomInserter : public MoleculeInserter {
    Point dir = {1, 1, 1};     //!< Scalars for random mass center position. Default (1,1,1)
    Point offset = {0, 0, 0};  //!< Added to random position. Default (0,0,0)
    bool rotate = true;           //!< Set to true to randomly rotate molecule when inserted. Default: true
    bool keep_positions = false;  //!< Set to true to keep original positions (default: false)
    bool allow_overlap = false;   //!< Set to true to skip container overlap check
    int max_trials = 20'000;      //!< Maximum number of container overlap checks
    int conformation_ndx = -1;    //!< Index of last used conformation

    ParticleVector operator()(Geometry::GeometryBase &geo, const ParticleVector &, MoleculeData &mol) override;
    void from_json(const json &j) override;
    void to_json(json &j) const override;
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

/**
 * @brief Determines if two particles within a group are excluded from mutual nonbonded interactions.
 *
 * Simple and naïve implementation storing all possible pairs in a matrix.
 * @internal Currently not used. MoleculeData can use this class internally.
 */
class ExclusionsSimple {
    bool any_exclusions = false; //!< true if at least one excluded interaction is present
    int size;                    //!< number of particles within the group; matrix = size × size
    //! 1D representation of excluded interaction matrix; the shared pointer saves copying
    //! Note that std::vector<bool> uses packing with performance implications
    std::shared_ptr<std::vector<unsigned char>> excluded_pairs;

  public:
    //! Creates and populates exclusions.
    static ExclusionsSimple create(int atoms_cnt, const std::vector<std::pair<int, int>> &pairs);
    explicit ExclusionsSimple(int size = 0);
    //! @param i, j indices of atoms within molecule with excluded nonbonded interaction
    void add(int i, int j);
    void add(const std::vector<std::pair<int, int>> &pairs);
    //! @param i, j indices of atoms within molecule with excluded nonbonded interaction
    bool isExcluded(int i, int j) const;
    bool empty() const; //!< true if no excluded interactions at all
    friend void from_json(const json &j, ExclusionsSimple &exclusions);
    friend void to_json(json &j, const ExclusionsSimple &exclusions);
};

inline bool ExclusionsSimple::isExcluded(int i, int j) const {
    if (i > j) {
        std::swap(i, j);
    }
    // use the bracket syntax to skip checks to speed-up reading
    return any_exclusions && (*excluded_pairs)[i * size + j];
}

inline bool ExclusionsSimple::empty() const { return !any_exclusions; }

void from_json(const json &j, ExclusionsSimple &exclusions);
void to_json(json &j, const ExclusionsSimple &exclusions);

/**
 * @brief Determines if two particles within a group are excluded from mutual nonbonded interactions.
 *
 * This is a memory optimized implementation especially for linear molecules. The matrix of excluded pairs
 * has dimensions of number of atoms × maximal distance between excluded neighbours (in terms of atom indices
 * within the molecule).
 *
 * @internal MoleculeData uses this class internally.
 */
class ExclusionsVicinity {
    //! count of atoms in the molecule; unmutable
    int atoms_cnt;
    //! maximal possible distance (difference) between indices of excluded particles
    int  max_bond_distance;
    //! 1D representation of excluded interaction matrix; the shared pointer saves copying
    //! Note that std::vector<bool> uses packing with performance implications
    std::shared_ptr<std::vector<unsigned char>> excluded_pairs;
    int toIndex(int i, int j) const;
    std::pair<int, int> fromIndex(int n) const;

  public:
    /** Generates memory-optimal structure from underlying pair list.
     * @param atoms_cnt number of atoms in the molecule
     * @param pairs list of excluded pairs using intramolecular indices
     */
    static ExclusionsVicinity create(int atoms_cnt, const std::vector<std::pair<int, int>> &pairs);
    /**
     *
     * @param atoms_cnt number of atoms in the molecule
     * @param max_difference maximal allowed distance
     */
    explicit ExclusionsVicinity(int atoms_cnt = 0, int max_difference = 0);
    //! @param i, j indices of atoms within molecule with excluded nonbonded interaction
    void add(int i, int j);
    void add(const std::vector<std::pair<int, int>> &pairs);
    //! @param i, j indices of atoms within molecule with excluded nonbonded interaction
    bool isExcluded(int i, int j) const;
    bool empty() const; //!< true if no excluded interactions at all
    // friend void from_json(const json &j, ExclusionsVicinity &exclusions); // not implemented
    friend void to_json(json &j, const ExclusionsVicinity &exclusions);
};

inline bool ExclusionsVicinity::isExcluded(int i, int j) const {
    if (i > j) {
        std::swap(i, j);
    }
    // use bracket syntax to skip checks to speed-up reading
    return (j - i <= max_bond_distance && (*excluded_pairs)[toIndex(i, j)]);
}

inline bool ExclusionsVicinity::empty() const { return max_bond_distance == 0; }

inline int ExclusionsVicinity::toIndex(int i, int j) const { return i * max_bond_distance + (j - i - 1); }

inline std::pair<int, int> ExclusionsVicinity::fromIndex(int n) const {
    int i = n / max_bond_distance;
    int j = n % max_bond_distance + i + 1;
    return std::make_pair(i, j);
}

// void from_json(const json &j, ExclusionsVicinity &exclusions);  // not implemented
void to_json(json &j, const ExclusionsVicinity &exclusions);

/**
 * @brief General properties for molecules
 */
class MoleculeData {
    json json_cfg; //!< data useful only for to_json
    int _id = -1;

  protected:
    ExclusionsVicinity exclusions; //!< Implementation of isPairExcluded;
                                   //!< various implementation can be provided in the future
  public:
    std::shared_ptr<MoleculeInserter> inserter = nullptr; //!< Functor for insertion into space

    int &id();             //!< Type id
    const int &id() const; //!< Type id
    void createMolecularConformations(const json&); //!< Add conformations if appropriate

    std::string name;            //!< Molecule name
    bool atomic = false;         //!< True if atomic group (salt etc.)
    bool compressible = false;   //!< True if compressible group (scales internally upon volume change)
    bool rigid = false;          //!< True if particle should be considered as rigid
    double activity = 0.0;       //!< Chemical activity (mol/l)

    std::vector<int> atoms;                    //!< Sequence of atoms in molecule (atom id's)
    BasePointerVector<Potential::BondData> bonds;
    WeightedDistribution<ParticleVector, float> conformations; //!< Conformations of molecule

    MoleculeData();
    MoleculeData(const std::string &name, ParticleVector particles,
                 const BasePointerVector<Potential::BondData> &bonds);

    bool isPairExcluded(int i, int j) const;

    /** @brief Specify function to be used when inserting into space.
     *
     * By default a random position and orientation is generator and overlap with container is avoided.
     */
    void setInserter(std::shared_ptr<MoleculeInserter> ins);

    /**
     * @brief Get random conformation that fits in container
     * @param geo Geometry
     * @param otherparticles Typically `spc.p` is insertion depends on other particle
     *
     * By default the molecule is placed at a random position and orientation with
     * no container overlap using the `RandomInserter` class. This behavior can
     * be changed by specifying another inserter using `setInserter()`.
     */
    ParticleVector getRandomConformation(Geometry::GeometryBase &geo,
                                         ParticleVector otherparticles = ParticleVector());

    void loadConformation(const std::string &file, bool keep_positions, bool keep_charges);

    friend class MoleculeBuilder;
    friend void to_json(json &j, const MoleculeData &a);
    friend void from_json(const json &j, MoleculeData &a);
}; // end of class

inline bool MoleculeData::isPairExcluded(int i, int j) const { return exclusions.isExcluded(i, j); }

void to_json(json &j, const MoleculeData &a);

void from_json(const json &j, MoleculeData &a);

void from_json(const json &j, std::vector<MoleculeData> &v);

// global instance of molecule vector
extern std::vector<MoleculeData> molecules;

/**
 * @brief Constructs MoleculeData from JSON.
 *
 * The instance is disposable: The method from_json() may be called only once.
 */
class MoleculeBuilder {
    bool is_used = false; //!< from_json may be called only once
    std::string molecule_name; //!< human readable name of the molecule
    ParticleVector particles;  //!< vector of particles; it unpacks to atoms and conformations[0]
    decltype(MoleculeData::bonds) bonds;
    std::vector<std::pair<int, int>> exclusion_pairs;
  protected:
    bool isFasta(const json &j_properties); //!< the structure is read in FASTA format with harmonic bonds
    void readCompoundValues(const json &j); //!< a director method calling the executive methods in the right order
    void readAtomic(const json &j_properties); //!< reads "atoms" : ["Na", "Cl"]
    void readParticles(const json &j_properties); //!< reads "structure" in any format; uses MoleculeStructureReader
    void readBonds(const json &j_properties); //!< reads "bondlist"
    void readFastaBonds(const json &j_properties); //!< makes up harmonic bonds for a FASTA sequence
    void readExclusions(const json &j_properties); //!< reads "exclusionlist" and "excluded_neighbours"
    std::shared_ptr<MoleculeInserter> createInserter(const json &j_properties);
  public:
    //! initialize MoleculeData from JSON; shall be called only once during the instance lifetime
    void from_json(const json &j, MoleculeData &molecule);
};

/**
 * @brief Fills the particle vector from various sources, e.g., files or JSON array.
 */
class MoleculeStructureReader {
    bool read_charges; //!< shall we also read charges when available, e.g., in PQR
  protected:
    //! reads array with atom types and positions
    void readArray(ParticleVector &particles, const json &j_particles);
    //! reads a FASTA sequence with harmonic bond parameters
    void readFasta(ParticleVector &particles, const json &j_fasta);
  public:
    MoleculeStructureReader(bool read_charges = true) : read_charges(read_charges) {};
    //! reads atom types, positions and optionally charges from a file
    void readFile(ParticleVector &particles, const std::string &filename);
    //! a director determining the executive method based on JSON content
    void readJson(ParticleVector &particles, const json &j);
};

/**
 * @brief Generate all possible atom pairs within a given bond distance. Only 1-2 bonds (e.g., harmonic or FENE)
 * are considered.
 *
 * @internal Used by MoleculeBuilder only to generate exclusions from excluded neighbours count.
 */
class NeighboursGenerator {
    //! a path created from 1-2 bonds (e.g., harmonic or FENE) as an ordered list of atoms involved;
    //! atoms are addressed by intramolecular indices
    typedef std::vector<int> AtomList;
    //! atom → list of directly bonded atoms with a 1-2 bond; addressing by indices
    std::map<int, AtomList> bond_map;
    //! paths indexed by a path length (starting with a zero distance for a single atom)
    //! a path is a sequence (without loops) of atoms connected by a 1-2 bond,
    std::vector<std::set<AtomList>> paths;

    typedef decltype(MoleculeData::bonds) BondVector;
    //! a constructor helper: populates bond_map
    void createBondMap(const BondVector &bonds);
    //! a generatePairs helper: populates paths
    void generatePaths(int bonded_distance);
  public:
    //! atom pairs addressed by intramolecular indices
    typedef std::vector<std::pair<int, int>> AtomPairList;
    NeighboursGenerator(const BondVector &bonds);
    /**
     * Append all atom pairs within a given bond distance to the pair list.
     * @param pairs atom pair list
     * @param bond_distance maximal number of 1-2 bonds between atoms to consider
     */
    void generatePairs(AtomPairList &pairs, int bond_distance);
};

/*
 * @brief General properties of reactions
 *
 * @todo
 * - [x] Rename "products" and "reactants" to "left" and "right"
 * - [ ] make all data private
 * - [x] add `setDirection()` / `getDirection()` functions
 * - [x] add `getProducts()`.
 *       Would point to "left" if direction == RIGHT
 *       Would point to "right" if direction == LEFT
 * - [x] add `getReactants()` like above
 * - [ ] remove awkward `bool forward`
 * - [ ] reversing direction would also change sign of pK
 */
class ReactionData {
  public:
    typedef std::map<int, int> Tmap;
    enum class Direction : char { LEFT = 0, RIGHT = 1 };
    std::vector<std::string> left_names, right_names;

    Tmap left_molecules; // Molecular change, groups. Atomic as Groupwise
    Tmap left_atoms;     // Atomic change, equivalent of swap/titration
    Tmap right_molecules;
    Tmap right_atoms;

  private:
    Direction direction = Direction::RIGHT;

  public:
    void setDirection(Direction);
    Direction getDirection() const;
    void reverseDirection(); //!< reverse direction of reaction
    std::pair<const Tmap &, const Tmap &> getProducts() const;
    std::pair<const Tmap &, const Tmap &> getReactants() const;

    bool canonic = false; //!< Finite reservoir
    bool swap = false;    //!< True if swap move
    int N_reservoir;      //!< Number of molecules in finite reservoir
    double lnK = 0;       //!< Natural logarithm of molar eq. const.
    double pK = 0;        //!< -log10 of molar eq. const.
    bool neutral = false; //!< True if only neutral molecules are involved in the reaction
    std::string reaction; //!< Name of reaction
    double weight;        //!< Statistical weight to be given to reaction in speciation

    bool empty() const;

    std::vector<int> participatingMolecules() const; //!< Returns molids of participating molecules
    bool containsMolecule(int molid) const; //!< True of molecule id is part of process

    auto findAtomOrMolecule(const std::string &name) const {
        auto atom_iter = findName(Faunus::atoms, name);
        auto molecule_iter = findName(Faunus::molecules, name);
        if (molecule_iter == Faunus::molecules.end())
            if (atom_iter == Faunus::atoms.end())
                throw std::runtime_error("unknown species '" + name + "'");
        return std::make_tuple(atom_iter, molecule_iter);
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
    return std::make_tuple(Tvec(v.begin(), it), Tvec(it + 1, v.end()));
} //!< Parse process string to pair of vectors containing reactant/product species

void from_json(const json &j, ReactionData &a);

void to_json(json &j, const ReactionData &a);

extern std::vector<ReactionData> reactions; // global instance

} // namespace Faunus
