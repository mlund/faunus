#pragma once

#include "core.h"
#include "auxiliary.h"
#include "particle.h"
#include "random.h"
#include "spdlog/spdlog.h"
#include "rotate.h"
#include <set>

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
    virtual ParticleVector operator()(Geometry::GeometryBase &geo, MoleculeData &mol,
                                      const ParticleVector &other_particles) = 0;
    virtual void from_json(const json& j);
    virtual void to_json(json& j) const;
    virtual ~MoleculeInserter() = default;
};

void from_json(const json &j, MoleculeInserter &inserter);
void to_json(json &j, const MoleculeInserter &inserter);

/**
 * @brief Inserts molecules into random positions in the container
 */
class RandomInserter : public MoleculeInserter {
  private:
    void translateRotateAtomicGroup(const Geometry::GeometryBase& geo, Faunus::QuaternionRotate& rotator,
                                    ParticleVector& particles) const;
    void translateRotateMolecularGroup(const Geometry::GeometryBase& geo, QuaternionRotate& rotator,
                                       ParticleVector& particles) const;

  public:
    Point dir = {1, 1, 1};       //!< Scalars for random mass center position. Default (1,1,1)
    Point offset = {0, 0, 0};    //!< Added to random position. Default (0,0,0)
    bool rotate = true;          //!< Set to true to randomly rotate molecule when inserted. Default: true
    bool keep_positions = false; //!< Set to true to keep original positions (default: false)
    bool allow_overlap = false;  //!< Set to true to skip container overlap check
    int max_trials = 20'000;     //!< Maximum number of container overlap checks

    ParticleVector operator()(Geometry::GeometryBase &geo, MoleculeData &molecule,
                              const ParticleVector &ignored_other_particles = ParticleVector()) override;
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
    void copyTo(ParticleVector &particles) const; //!< Copy conformation into particle vector
};

/**
 * @brief Determines if two particles within a group are excluded from mutual nonbonded interactions.
 *
 * Simple and naïve implementation storing all possible pairs in a matrix.
 * @internal Currently not used. MoleculeData can use this class internally.
 */
class ExclusionsSimple {
    bool any_exclusions = false;     //!< true if at least one excluded interaction is present
    int size;                        //!< number of particles within the group; matrix = size × size
    using char_bool = unsigned char; //!< obs: std::vector<bool> uses packing with performance implications
    std::shared_ptr<std::vector<char_bool>> excluded_pairs; //!< 1D exclusion matrix; shared ptr saves copying

  public:
    using AtomPair = std::pair<int, int>;
    //! Creates and populates exclusions.
    static ExclusionsSimple create(int atoms_cnt, const std::vector<std::pair<int, int>> &pairs);
    explicit ExclusionsSimple(int size = 0);
    //! @param i, j indices of atoms within molecule with excluded nonbonded interaction
    void add(int i, int j);
    void add(const std::vector<AtomPair>& pairs);
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
    return any_exclusions && static_cast<bool>((*excluded_pairs)[i * size + j]);
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
  public:
    using char_bool = unsigned char; //!< obs: std::vector<bool> uses packing with performance implications
    using AtomPair = std::pair<int, int>;

  private:
    int atoms_cnt = 0;         //!< count of atoms in the molecule; unmutable
    int max_bond_distance = 0; //!< max distance (difference) between indices of excluded particles
    std::shared_ptr<std::vector<char_bool>> excluded_pairs; //!< 1D exclusion matrix; shared ptr saves copying

    int toIndex(int i, int j) const;
    AtomPair fromIndex(int n) const;

  public:
    /** Generates memory-optimal structure from underlying pair list.
     * @param atoms_cnt number of atoms in the molecule
     * @param pairs list of excluded pairs using intramolecular indices
     */
    static ExclusionsVicinity create(int atoms_cnt, const std::vector<AtomPair>& pairs);
    /**
     *
     * @param atoms_cnt number of atoms in the molecule
     * @param max_difference maximal allowed distance
     */
    explicit ExclusionsVicinity(int atoms_cnt = 0, int max_difference = 0);

    //! @param i, j indices of atoms within molecule with excluded nonbonded interaction
    void add(int i, int j);
    void add(const std::vector<std::pair<int, int>> &pairs);
    bool isExcluded(int i,
                    int j) const; //!< @param i, j indices of atoms within molecule with excluded nonbonded interaction
    bool empty() const; //!< true if no excluded interactions at all
    friend void to_json(json &j, const ExclusionsVicinity &exclusions);
};

inline bool ExclusionsVicinity::isExcluded(int i, int j) const {
    if (i > j) {
        std::swap(i, j);
    }
    // use bracket syntax to skip checks to speed-up reading
    return (j - i <= max_bond_distance && static_cast<bool>((*excluded_pairs)[toIndex(i, j)]));
}

inline bool ExclusionsVicinity::empty() const { return max_bond_distance == 0; }

inline int ExclusionsVicinity::toIndex(int i, int j) const { return i * max_bond_distance + (j - i - 1); }

inline ExclusionsVicinity::AtomPair ExclusionsVicinity::fromIndex(const int n) const {
    auto i = n / max_bond_distance;
    auto j = n % max_bond_distance + i + 1;
    return {i, j};
}

// void from_json(const json &j, ExclusionsVicinity &exclusions);  // not implemented
void to_json(json &j, const ExclusionsVicinity &exclusions);

/**
 * @brief General properties for molecules
 */
class MoleculeData {
  public:
    using index_type = int;

  private:
    json json_cfg; //!< data useful only for to_json
    index_type _id = 0;
    bool implicit = false; //!< Is molecule implicit and explicitly absent from simulation cell?

  protected:
    ExclusionsVicinity exclusions; //!< Implementation of isPairExcluded;
                                   //!< various implementation can be provided in the future
  public:
    std::shared_ptr<MoleculeInserter> inserter = nullptr; //!< Functor for insertion into space

    index_type& id();                                 //!< Type id
    const index_type& id() const;                     //!< Type id
    void createMolecularConformations(const json& j); //!< Add conformations if appropriate
    void setConformationWeights(const json& j);       //!< Add weights for conformations

    std::string name;            //!< Molecule name
    bool atomic = false;         //!< True if atomic group (salt etc.)
    bool compressible = false;   //!< True if compressible group (scales internally upon volume change)
    bool rigid = false;          //!< True if particle should be considered as rigid
    double activity = 0.0;       //!< Chemical activity (mol/l)

    std::vector<AtomData::index_type> atoms; //!< Sequence of atoms in molecule (atom id's)
    BasePointerVector<Potential::BondData> bonds;
    WeightedDistribution<ParticleVector> conformations; //!< Conformations of molecule
    size_t numConformations() const;                    //!< Number of conformations

    MoleculeData();
    MoleculeData(const std::string& name, const ParticleVector& particles,
                 const BasePointerVector<Potential::BondData>& bonds);

    bool isImplicit() const; //!< Is molecule implicit and explicitly absent from simulation cell?
    bool isPairExcluded(int i, int j) const;
    bool isMolecular() const;
    bool isAtomic() const;

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
    ParticleVector getRandomConformation(Geometry::GeometryBase& geo,
                                         const ParticleVector& otherparticles = ParticleVector());

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
 * @brief An exception to indicate an unknown molecule name in the input.
 */
struct UnknownMoleculeError: public GenericError {
    explicit UnknownMoleculeError(const std::string &molecule_name);
};

/**
 * @brief Finds a molecule by its name in the global Faunus molecules lexicon.
 *
 * The first matching molecule is returned, or an UnknownMoleculeError is thrown when not found.
 *
 * @param name  a molecule name to look for
 * @return a molecule found
 * @throw UnknownMoleculeError  when no molecule found
 */
MoleculeData& findMoleculeByName(const std::string& name);

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
    static bool isFasta(const json& j);     //!< the structure is read in FASTA format with harmonic bonds
    void readCompoundValues(const json &j); //!< a director method calling the executive methods in the right order
    void readAtomic(const json& j);         //!< reads "atoms" : ["Na", "Cl"]
    void readParticles(const json& j);      //!< reads "structure" in any format; uses MoleculeStructureReader
    void readBonds(const json& j);          //!< reads "bondlist"
    void readFastaBonds(const json& j);     //!< makes up harmonic bonds for a FASTA sequence
    void readExclusions(const json& j);     //!< reads "exclusionlist" and "excluded_neighbours"
    std::shared_ptr<MoleculeInserter> createInserter(const json& j);

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
    static void readArray(ParticleVector& particles, const json& j_particles);
    //! reads a FASTA sequence with harmonic bond parameters
    static void readFasta(ParticleVector& particles, const json& j);

  public:
    MoleculeStructureReader(bool read_charges = true);
    //! reads atom types, positions and optionally charges from a file
    void readFile(ParticleVector& particles, const std::string& filename) const;
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
  public:
    using AtomList = std::vector<int>;
    using AtomPair = std::pair<int, int>;
    using AtomPairList = std::vector<AtomPair>;
    using BondVector = decltype(MoleculeData::bonds);

  private:
    std::map<int, AtomList> bond_map; //!< atom → list of directly bonded atoms with a 1-2 bond; addressing by indices
    //! paths indexed by a path length (starting with a zero distance for a single atom)
    //! a path is a sequence (without loops) of atoms connected by a 1-2 bond,
    std::vector<std::set<AtomList>> paths;

    void createBondMap(const BondVector& bonds); //!< a constructor helper: populates bond_map
    void generatePaths(int bonded_distance);     //!< a generatePairs helper: populates paths

  public:
    //! atom pairs addressed by intramolecular indices
    NeighboursGenerator(const BondVector &bonds);
    /**
     * Append all atom pairs within a given bond distance to the pair list.
     * @param pairs atom pair list
     * @param bond_distance maximal number of 1-2 bonds between atoms to consider
     */
    void generatePairs(AtomPairList &pairs, int bond_distance);
};

/**
 * @brief General properties of reactions
 *
 * Placeholder for chemical reactions used in the RCMC move.
 * A reaction has two sides, left and right, each contained one or more
 * atomic / molecular reactants and products. The reaction is associated with
 * an equilibrium constant, `lnK` or `pK`.
 * If the direction is `RIGHT`, the right-hand side species are products and
 * the left-hand side are reactants. Vice versa if the direction is `LEFT`.
 * The direction can be changed with `setDirection()` which also handles
 * sign changes of `lnK` and `pK`.
 * The functions `getProducts()` and `getReactants()` returns a pair with
 * atomic and molecular reactants/products, always reflecting the current
 * direction.
 *
 * @todo
 * - [x] Enable `canonic` and `reservoir_size`
 * - [ ] Enable reservoir size to be given as _molarity_ and _number_
 * - [ ] Merge products and reactant structures using signed stoichiometric numbers
 * - [x] `reservoir_size` should ideally be associated with an implicit molecule
 * - [ ] `direction` should be a member of ReactionData due to possible data races
 */
class ReactionData {
  public:
    using StoichiometryMap = std::map<int, int>; // key = molid; value = stoichiometic coefficient
    enum class Direction : char { LEFT = 0, RIGHT = 1 };

  private:
    friend void from_json(const json &, ReactionData &);
    friend void to_json(json &, const ReactionData &);

    Direction direction = Direction::RIGHT;           //!< Direction of reaction
    std::vector<std::string> left_names, right_names; //!< Names of reactants and products

    StoichiometryMap left_molecules;  //!< Initial reactants (molecules)
    StoichiometryMap right_molecules; //!< Initial products (molecules)
    StoichiometryMap left_atoms;      //!< Initial reactants (atoms)
    StoichiometryMap right_atoms;     //!< Initial products (atoms)

  public:
    void setDirection(Direction);                               //!< Set directions of the process
    Direction getDirection() const;                             //!< Get direction of the process
    void reverseDirection();                                    //!< Reverse direction of reaction
    std::pair<const StoichiometryMap &, const StoichiometryMap &>
    getProducts() const; //!< Pair with atomic and molecular products
    std::pair<const StoichiometryMap &, const StoichiometryMap &>
    getReactants() const; //!< Pair with atomic and molecular reactants

    std::pair<std::set<int>, std::set<int>> getReactantsAndProducts() const;

    bool swap = false;                   //!< True if swap move
    double lnK = 0.0;                    //!< Effective, natural logarithm of molar eq. const.
    double lnK_unmodified = 0.0;         //!< Natural logarithm of molar eq. const. (unmodified as in input)
    bool only_neutral_molecules = false; //!< Only neutral molecules are involved in the reaction
    std::string reaction_str;            //!< Name of reaction

    /**
     * @brief Find atom name or molecule name
     * @returns pair of iterators to atomlist and moleculelist; one of them points to end().
     *
     * Note that molecules and atoms *cannot* have the same name
     */
    static auto findAtomOrMolecule(const std::string& atom_or_molecule_name) {
        auto atom_iter = findName(Faunus::atoms, atom_or_molecule_name);
        auto molecule_iter = findName(Faunus::molecules, atom_or_molecule_name);
        if (molecule_iter == Faunus::molecules.end() and atom_iter == Faunus::atoms.end()) {
            throw std::runtime_error("unknown species '" + atom_or_molecule_name + "'");
        }
        return std::make_pair(atom_iter, molecule_iter);
    }

}; //!< End of class

/**
 * This parses a string containing a reaction, e.g. "A = B + B" and returns
 * a pair with (1) a vector of reactant names and (2) a vector of product names.
 * Reactants and products are split by a `=` sign. All elements in the string
 * must be separated by a white-space.
 */
inline auto parseReactionString(const std::string &process_string) {
    using Tvec = std::vector<std::string>;
    Tvec names; // vector of atom/molecule names
    std::string atom_or_molecule_name;
    std::istringstream iss(process_string);
    while (iss >> atom_or_molecule_name) { // stream all words into vector
        names.push_back(atom_or_molecule_name);
    }

    names.erase(std::remove(names.begin(), names.end(), "+"), names.end());

    auto it = std::find(names.begin(), names.end(), "=");
    if (it == names.end()) {
        throw std::runtime_error("products and reactants must be separated by ' = '");
    }
    return std::make_pair(Tvec(names.begin(), it), Tvec(it + 1, names.end()));
}

void from_json(const json &, ReactionData &);
void to_json(json &, const ReactionData &);

extern std::vector<ReactionData> reactions; // global instance

} // namespace Faunus
