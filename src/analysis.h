#pragma once

#include "particle.h"
#include "molecule.h"
#include "space.h"
#include "io.h"
#include "scatter.h"
#include "reactioncoordinate.h"
#include "aux/timers.h"
#include "aux/table_2d.h"
#include "aux/equidistant_table.h"
#include <aux/sparsehistogram.h>
#include <Eigen/SparseCore>
#include <limits>
#include <memory>
#include <set>
#include <string_view>
#include <string>
#include <vector>

namespace cereal {
class BinaryOutputArchive;
}

namespace Faunus::pairpotential {
class NewCoulombGalore;
}

namespace Faunus::Energy {
class Hamiltonian;
class EnergyTerm;
class Penalty;
} // namespace Faunus::Energy

namespace Faunus::SASA {
class SASABase;
}

namespace Faunus {
/**
 * Adding a new analysis requires the following steps:
 *
 * 1. Write an analysis class derived from `Analysis::Analysisbase`
 * 2. Add the new class to the factory function `Analysis::createAnalysis()`, see below
 * 3. Recommended: add the required json input format to `docs/schema.yml`
 * 4. Recommended: add documentation to `docs/_docs/analysis.md`
 */
namespace analysis {

/**
 * @brief Base class for all analysis functions
 *
 * Analysis classes are used to perform system analysis
 * over the course or a simulations or for post-analysis
 * of an existing trajectory. The base class adds basic
 * functionality such as timing, number of steps, json IO
 * and enforce a common interface.
 */
class Analysis {
  private:
    virtual void _to_json(json&) const;                   //!< provide json information
    virtual void _from_json(const json&);                 //!< setup from json
    virtual void _sample() = 0;                           //!< perform sample event
    virtual void _to_disk();                              //!< save sampled data to disk
    int number_of_steps = 0;                              //!< counter for total number of steps
    int number_of_skipped_steps = 0;                      //!< steps to skip before sampling (do not modify)
    TimeRelativeOfTotal<std::chrono::microseconds> timer; //!< timer to benchmark `_sample()`

  protected:
    const Space& spc;          //!< Instance of Space to analyse
    int sample_interval = 0;   //!< Steps in between each sample point (do not modify)
    int number_of_samples = 0; //!< counter for number of samples

  public:
    const std::string name;        //!< descriptive name
    std::string cite;              //!< url, doi etc. describing the analysis
    void to_json(json& j) const;   //!< JSON report w. statistics, output etc.
    void from_json(const json& j); //!< configure from json object
    void to_disk();                //!< Save data to disk (if defined)
    void sample();                 //!< Increase step count and sample
    int getNumberOfSteps() const;  //!< Number of steps
    Analysis(const Space& spc, std::string_view name);
    Analysis(const Space& spc, std::string_view name, int sample_interval, int number_of_skipped_steps);
    virtual ~Analysis() = default;
};

void to_json(json& j, const Analysis& base);

/**
 * @brief Function to create an instance of an analysis based on it's `name`
 * @param name Name of analysis to create
 * @param j JSON input for analysis
 * @param spc Space the analysis should operate on
 * @param pot Hamiltonian representing the system
 * @return Shared pointer to analysis
 *
 * After writing a new analysis, it must be added to this function in
 * order to be controlled from the main faunus input.
 */
std::unique_ptr<Analysis> createAnalysis(const std::string& name, const json& j, Space& spc, Energy::Hamiltonian& pot);

/**
 * @brief Aggregator class for storing and selecting multiple analysis instances
 *
 * Holds an arbitrary number of analysis instances and selects them at
 * random, based on user input. This is typically called from the
 * main simulation loop.
 */
class CombinedAnalysis : public BasePointerVector<Analysis> {
  public:
    CombinedAnalysis(const json& json_array, Space& spc, Energy::Hamiltonian& pot);
    void sample();
    void to_disk(); //!< prompt all analysis to save to disk if appropriate
};

/**
 * @brief Base class for perturbation analysis
 *
 * This class provides basic data and functions to support Widom Particle Insertion, Virtual Volume Move
 * etc. If a non-empty `filename` is given, a (compressed) output file used for streaming will be opened.
 *
 * @note Constructor throws if non-empty filename cannot to opened for writing
 */
class PerturbationAnalysis : public Analysis {
  private:
    void _to_disk() override;

  protected:
    Space& mutable_space; //!< This reference to space can be changed
    Energy::EnergyTerm& pot;
    std::string filename;                             //!< output filename (optional)
    std::unique_ptr<std::ostream> stream = nullptr;   //!< output file stream if filename given
    Change change;                                    //!< Change object to describe perturbation
    Average<double> mean_exponentiated_energy_change; //!< < exp(-du/kT) >
    bool collectWidomAverage(double energy_change);   //!< add to exp(-du/kT) incl. safety checks
    PerturbationAnalysis(const std::string& name, Energy::EnergyTerm& pot, Space& spc,
                         const std::string& filename = ""s);
    double meanFreeEnergy() const; //!< Average perturbation free energy, `-ln(<exp(-du/kT)>)`
};

/**
 * @brief Sample and save reaction coordinates to a file
 */
class FileReactionCoordinate : public Analysis {
  private:
    Average<double> mean_reaction_coordinate;
    const std::string filename;
    std::unique_ptr<std::ostream> stream = nullptr;
    const std::unique_ptr<ReactionCoordinate::ReactionCoordinateBase> reaction_coordinate;

    void _to_json(json& j) const override;
    void _sample() override;
    void _to_disk() override;
    FileReactionCoordinate(const Space& spc, const std::string& filename,
                           std::unique_ptr<ReactionCoordinate::ReactionCoordinateBase> reaction_coordinate);

  public:
    FileReactionCoordinate(const json& j, const Space& spc);
};

/**
 * @brief Tracks displacements of particle positions
 *
 * This tracks particle position displacements relative to a (dynamic) reference position. For
 * the first particle only, a trajectory file can be saved to disk.
 * Positions are tracked across periodic boundaries to enable large displacements relative
 * to initial reference positions.
 *
 * @warning Will not work if `molid` is grand canonical
 */
class AtomicDisplacement : public Analysis {
  protected:
    MoleculeData::index_type molid;

  private:
    using average_type = Average<double>;
    std::unique_ptr<std::ostream> single_position_stream; //!< Stream x, y, z, displacement for first position only
    std::vector<Point> reference_positions;               //!< Positions used as references
    std::vector<Point> previous_positions;                //!< Positions from previous analysis event
    std::vector<Eigen::Vector3i> cell_indices;            //!< Tracks in which unit cell the particles are in
    std::vector<average_type> mean_squared_displacement;  //!< Mean squared displacement for each position
    double max_possible_displacement; //!< If any displacement is larger than this, assume unit cell jump

    SparseHistogram<double> displacement_histogram;                 //!< P(r) where r is distance from reference
    std::string displacement_histogram_filename;                    //!< Name of P(r) histogram file
    int reference_reset_interval = std::numeric_limits<int>::max(); //!< Renew reference at given interval

    virtual PointVector getPositions() const;                        //!< Extract current positions from `molid`
    void resetReferencePosition(const Point& position, int index);   //!< Store current positions as reference
    Point getOffset(const Point& diff, Eigen::Vector3i& cell) const; //!< Offset to other cells
    void sampleDisplacementFromReference(const Point& position, int index);
    void _sample() override;
    void _to_json(json& j) const override;
    void _to_disk() override;

  public:
    AtomicDisplacement(const json& j, const Space& spc, std::string_view name = "displacement");
};

/**
 * @brief Tracks displacements of molecular group mass centers
 *
 * Based on `AtomicDisplacement`.
 */
class MassCenterDisplacement : public AtomicDisplacement {
  private:
    PointVector getPositions() const override; //!< Extracts mass centers from all active `molid` groups
  public:
    MassCenterDisplacement(const json& j, const Space& spc, std::string_view name = "displacement_com");
};

/**
 * @brief Excess chemical potential of molecules
 *
 * @todo Migrate `absolute_z_coords` into new `MoleculeInserter` policy
 */
class WidomInsertion : public PerturbationAnalysis {
    std::shared_ptr<MoleculeInserter> inserter; //!< Insertion method
    int number_of_insertions;                   //!< Number of insertions per sample event
    MoleculeData::index_type molid;             //!< Molecule id
    bool absolute_z_coords = false;             //!< Apply abs() on all inserted z coordinates?

    void selectGhostGroup(); //!< Select inactive group to act as group particle
    void updateGroup(Space::GroupType& group, const ParticleVector& particles);
    void _sample() override; //!< Called for each sample event
    void _to_json(json& j) const override;
    void _from_json(const json& j) override;

  public:
    WidomInsertion(const json& j, Space& spc, Energy::Hamiltonian& pot);
};

/**
 * @brief Samples the electric potential correlation as a function of distance, <Φ₁Φ₂>(r)
 */
class PotentialCorrelation : public Analysis {
  private:
    std::unique_ptr<pairpotential::NewCoulombGalore> coulomb;     //!< Class for calculating the potential
    Equidistant2DTable<double, Average<double>> mean_correlation; //!< <Φ₁Φ₂>(r)
    std::string filename;                                         //!< Filename of output <Φ₁Φ₂>(r) table
    std::pair<double, double> range;                              //!< Distance range
    double dr = 0;                                                //!< Distance resolution
    unsigned int calculations_per_sample_event = 1;
    void _to_json(json& j) const override;
    void _sample() override;
    void _to_disk() override;

  public:
    PotentialCorrelation(const json& j, const Space& spc);
};

/**
 * @brief Samples the electric potential at arbitrary positions in the simulation box
 */
class ElectricPotential : public Analysis {
  public:
    enum class Policies { FIXED, RANDOM_WALK, RANDOM_WALK_NO_OVERLAP, INVALID };

  private:
    unsigned int max_overlap_trials = 100; //!< Maximum number of overlap checks before bailing out
    double histogram_resolution = 0.01;    //!< Potential resolution
    unsigned int calculations_per_sample_event = 1;
    std::string file_prefix; //!< Output filename prefix for potential histogram and correlation
    struct Target {
        Point position;                                               //!< Target position
        Average<double> mean_potential;                               //!< mean potential at position
        std::unique_ptr<SparseHistogram<double>> potential_histogram; //!< Histogram of observed potentials
    };
    std::vector<Target> targets;                              //!< List of target points where to sample the potential
    Policies policy;                                          //!< Policy to apply to targets before each sample event
    std::unique_ptr<pairpotential::NewCoulombGalore> coulomb; //!< Class for calculating the potential
    Average<double> mean_potential_correlation;               //!< Correlation between targets, <phi1 x phi2 x ... >
    SparseHistogram<double> potential_correlation_histogram;  //!< Distribution of correlations, P(<phi1 x phi2 x ... >)
    void getTargets(const json& j);                           //!< Get user defined target positions
    void setPolicy(const json& j);                            //!< Set user defined position setting policy
    double calcPotentialOnTarget(const Target& target);       //!< Evaluate net potential of target position
    bool overlapWithParticles(const Point& position) const; //!< Check if position is within the radius of any particle
    std::function<void()> applyPolicy;                      //!< Lambda for position setting policy
    json output_information;                                //!< json output generated during construction
    void _to_json(json& j) const override;
    void _sample() override;
    void _to_disk() override;

  public:
    ElectricPotential(const json& j, const Space& spc);
};

NLOHMANN_JSON_SERIALIZE_ENUM(ElectricPotential::Policies,
                             {
                                 {ElectricPotential::Policies::INVALID, nullptr},
                                 {ElectricPotential::Policies::FIXED, "fixed"},
                                 {ElectricPotential::Policies::RANDOM_WALK_NO_OVERLAP, "no_overlap"},
                                 {ElectricPotential::Policies::RANDOM_WALK, "random_walk"},
                             })

/**
 * @brief Density of atom along axis
 *
 * Calculates the summed density of `atoms` in spherical, cylindrical or planar shells around
 * `origo` which by default is the center of the simulation box
 */
class AtomProfile : public Analysis {
    Equidistant2DTable<double, double> table;         //!< x = distance; y = count or charge
    double dr;                                        //!< resolution for `table`
    std::vector<std::string> names;                   //!< atom names to analyse
    std::set<AtomData::index_type> atom_id_selection; //!< atom ids to analyse
    std::string file;                                 //!< output filename
    Point origin = {0.0, 0.0, 0.0};
    Eigen::Vector3i dir = {1, 1, 1};
    bool count_charge = false;
    int center_of_mass_atom_id = -1; // center at COM of id_com atoms?

    double distanceToOrigin(const Point& position) const;
    void _from_json(const json& j) override;
    void _to_json(json& j) const override;
    void _to_disk() override;
    void _sample() override;

  public:
    AtomProfile(const json& j, const Space& spc);
};

/**
 * @brief Measures the density of atoms along z axis
 */
class SlicedDensity : public Analysis {
    double dz = 0.0;                         //!< Resolution of `histogram` along z
    Table2D<double, unsigned int> histogram; // N(z)
    std::vector<std::string> atom_names;
    std::vector<AtomData::index_type> atom_ids; //!< atom id's to analyse
    std::string file;
    int center_of_mass_atom_id = -1; // center at COM of id_com atoms?

    void _from_json(const json& j) override;
    void _to_json(json& j) const override;
    void _sample() override;
    void _to_disk() override;

  public:
    SlicedDensity(const json& j, const Space& spc);
};

/**
 * Abstract base class for analysing atomic and molecular densities
 */
class DensityBase : public Analysis {
  protected:
    using Table = Equidistant2DTable<unsigned int, double>;
    using id_type = size_t;
    std::map<id_type, Average<double>> mean_density;
    std::map<id_type, std::string_view> names; // id <-> name database
    void _to_disk() override;
    void _sample() override;
    void _to_json(json& j) const override;
    void writeTable(std::string_view name, Table& table);

  private:
    virtual std::map<id_type, int> count() const = 0;
    std::map<MoleculeData::index_type, Table> probability_density;
    Average<double> mean_cubic_root_of_volume;
    Average<double> mean_volume;
    Average<double> mean_inverse_volume;
    double updateVolumeStatistics();

  public:
    template <RequireNamedElements Range>
    DensityBase(const Space& spc, const Range& atoms_or_molecules, std::string_view name)
        : Analysis(spc, name) {
        for (const auto& data : atoms_or_molecules) {
            names[data.id()] = data.name;
            probability_density[data.id()].setResolution(1, 0);
        }
    }
};

/**
 * @brief Analysis of molecular group densities
 */
class MoleculeDensity final : public DensityBase {
  private:
    std::map<id_type, int> count() const override;

  public:
    MoleculeDensity(const json& j, const Space& spc);
};

/**
 * @brief Analysis of single atom densities
 */
class AtomDensity final : public DensityBase {
  private:
    std::map<id_type, Table> atomswap_probability_density;
    void _sample() override;
    void _to_disk() override;
    std::map<id_type, int> count() const override;

  public:
    AtomDensity(const json& j, const Space& spc);
};

/**
 * @brief Generates lambda to calculate properties between groups
 *
 * This looks for the `property` key in the json object. Supported names are
 * `energy`, `com_distance`
 */
std::function<double(const Group&, const Group&)> createGroupGroupProperty(const json& j, const Space& spc,
                                                                           Energy::Hamiltonian& hamiltonian);

std::function<bool(double)> createValueFilter(const std::string& name, const json& j, bool throw_on_error);

/**
 * Base class for streaming pair properties to a sparse matrix. To further
 * reduce the file size, value can be filtered using `value_filter` which by
 * default is disabled (i.e. always returns true). The json object is searched
 * for:
 *
 * - `file`   = output stream file name
 * - `filter` = ExprTk expression for filter
 */
class PairMatrixAnalysis : public Analysis {
  private:
    std::unique_ptr<std::ostream> matrix_stream; //!< Pair matrix output stream
    std::string filename;                        //!< cluster fraction filename
    virtual void setPairMatrix() = 0;            //!< Fills in `pair_matrix` with sampled values
    void _to_json(json& j) const override;
    void _from_json(const json& j) override;
    void _sample() override;
    void _to_disk() override;

  protected:
    const Space& spc;
    std::function<bool(double)> value_filter; //!< Used to filter values generated by `property`
    Eigen::SparseMatrix<double> pair_matrix;  //!< Matrix of pair properties

  public:
    PairMatrixAnalysis(const json& j, const Space& spc);
};

/**
 * @brief Samples a matrix of group properties
 *
 * This calculates a user-defined property between a set of
 * selected groups and outputs a matrix as a function of steps.
 */
class GroupMatrixAnalysis : public PairMatrixAnalysis {
  private:
    std::function<double(const Group&, const Group&)> property; //!< The group-group property to calculate
    std::vector<int> group_indices;                             //!< Selected groups (active and inactive)
    void setPairMatrix() override;                              //!< Fills in `pair_matrix` with sampled values
  public:
    GroupMatrixAnalysis(const json& j, const Space& spc, Energy::Hamiltonian& pot);
};

/**
 * @todo `atom_mean_charges` could be calculated from `atom_histograms`
 */
class ChargeFluctuations : public Analysis {
  private:
    typename decltype(Faunus::molecules)::const_iterator mol_iter; //!< selected molecule type
    using AtomHistogram = std::map<AtomData::index_type, int>;     //!< key = atom id; value = counts
    std::vector<AtomHistogram> atom_histograms;                    //!< one element for each atom in molecule
    std::vector<AverageStdev<double>> atom_mean_charges;           //!< average charges of atomic indexes
    std::string filename;                                          //!< name of PQR file with average charges
    bool verbose = true;                                           //!< set to true for more output

    ParticleVector averageChargeParticles(const Space::GroupType& group);
    std::vector<double> getMeanCharges() const;
    std::vector<double> getChargeStandardDeviation() const;
    std::vector<std::string> getPredominantParticleNames() const;

    void _sample() override;
    void _to_json(json& j) const override;
    void _to_disk() override;

  public:
    ChargeFluctuations(const json& j, const Space& spc);
}; // Fluctuations of atomic charges

/**
 * @brief Molecular multipole moments and their fluctuations
 */
class Multipole : public Analysis {
    struct Data {
        Average<double> charge;
        Average<double> charge_squared;
        Average<double> dipole_moment;
        Average<double> dipole_moment_squared;
    }; // Average sample moment for a molecule
    std::map<MoleculeData::index_type, Data>
        average_moments; //!< Molecular moments and their fluctuations. Key = molid.
    void _sample() override;
    void _to_json(json& j) const override;

  public:
    Multipole(const json& j, const Space& spc);
};

/**
 * @brief Save system energy to disk
 */
class SystemEnergy : public Analysis {
  private:
    Energy::Hamiltonian& hamiltonian;
    std::string file_name;
    std::string separator;
    std::unique_ptr<std::ostream> output_stream;
    std::vector<double> calculateEnergies() const;
    Average<double> mean_energy;
    Average<double> mean_squared_energy;
    Table2D<double, double> energy_histogram; // Density histograms
    double initial_energy = 0.0;

    void createOutputStream();
    void normalize();
    void _sample() override;
    void _to_json(json& j) const override;
    void _from_json(const json& j) override;
    void _to_disk() override;

  public:
    SystemEnergy(const json& j, const Space& spc, Energy::Hamiltonian& hamiltonian);
};

/**
 * @brief Checks if system is sane. If not, abort program.
 */
class SanityCheck : public Analysis {
  private:
    const double mass_center_tolerance = 1.0e-6;
    void _sample() override;
    void checkGroupsCoverParticles();                         //!< Groups must exactly contain all particles in `p`
    void checkMassCenter(const Space::GroupType& group);      //!< check if molecular mass centers are correct
    void checkWithinContainer(const Space::GroupType& group); //!< check if particles are inside container

  public:
    SanityCheck(const json& j, const Space& spc);
};

/**
 * @brief Save simulation state and particle coordinates to disk
 *
 * Using a variety of formats (pqr, gro, xyz, aam, state) this
 * stores the simulation state to disk. The most complete format
 * is `.state` which stores information about groups, particle
 * positions, random number state etc.
 *
 * - if sample interval = -1, analysis is run only once at the simulation end.
 * - if sample interval >= 0, analysis is performed as every nstep.
 * - if `use_numbered_files` = true (default) files are labelled with the step count
 */
class SaveState : public Analysis {
  private:
    std::function<void(const std::string&)> writeFunc = nullptr;
    bool save_random_number_generator_state = false;
    bool use_numbered_files = true;
    bool convert_hexagonal_prism_to_cuboid = false;
    std::string filename;

    void _to_json(json& j) const override;
    void _sample() override;
    void saveAsCuboid(const std::string& filename, const Space& spc, StructureFileWriter& writer) const;
    void saveJsonStateFile(const std::string& filename, const Space& spc) const;
    void saveBinaryJsonStateFile(const std::string& filename, const Space& spc) const;

  public:
    SaveState(json j, const Space& spc);
    ~SaveState() override;
    void setWriteFunction(const Space& spc);
};

/**
 * @brief Base class for distribution functions etc.
 */
class PairFunction : public Analysis {
  protected:
    using index_type = size_t;
    int dimensions = 3;                           //!< dimensions to use when normalizing
    index_type id1 = 0;                           //!< molecule or atom id
    index_type id2 = 0;                           //!< molecule or atom id
    double dr = 0;                                //!< distance resolution
    std::string name1;                            //!< atom or molecule
    std::string name2;                            //!< atom or molecule name
    std::string file;                             //!< output filename
    Equidistant2DTable<double, double> histogram; //!< distance histogram
    Average<double> mean_volume;                  //!< average volume (angstrom^3)
    Eigen::Vector3i slicedir = {0, 0, 0};
    double thickness = 0;

  private:
    void _from_json(const json& j) override;
    void _to_json(json& j) const override;
    void _to_disk() override;
    double volumeElement(double r) const;

  public:
    PairFunction(const Space& spc, const json& j, const std::string_view name);
};

/**
 * @brief Atomic radial distribution function, g(r)
 */
class AtomRDF : public PairFunction {
  private:
    void _sample() override;
    void sampleDistance(const Particle& particle1, const Particle& particle2);
    void sampleDifferent(); //!< particle types are different (id1!=id2)
    void sampleIdentical(); //!< particle types are identical (id1==id2)

  public:
    AtomRDF(const json& j, const Space& spc);
};

/**
 * @brief Same as `AtomRDF` but for molecules. Identical input
 */
class MoleculeRDF : public PairFunction {
  private:
    void _sample() override;
    void sampleDistance(const Group& group_i, const Group& group_j);
    void sampleDifferent(); //!< group types are different (id1!=id2)
    void sampleIdentical(); //!< group types are identical (id1==id2)

  public:
    MoleculeRDF(const json& j, const Space& spc);
};

/**
 * @todo Is this class justified? Messy file handling
 */
class PairAngleFunction : public PairFunction {
  private:
    std::string correlation_filename;

  protected:
    Equidistant2DTable<double, Average<double>> average_correlation_vs_distance;

  private:
    void _from_json(const json& j) override;
    void _to_disk() override;

  public:
    PairAngleFunction(const Space& spc, const json& j, const std::string& name);
};

/** @brief Dipole-dipole correlation function, <\boldsymbol{\mu}(0)\cdot\boldsymbol{\mu}(r)> */
class AtomDipDipCorr : public PairAngleFunction {
    void _sample() override;

  public:
    AtomDipDipCorr(const json& j, const Space& spc);
};

/**
 * @brief Write XTC trajectory file
 *
 * This will save both active and inactive particles as the XTC format must have the
 * same number of particles on all frames. See `QRtraj` for a solution how to disable
 * inactive groups in e.g. VMD.
 */
class XTCtraj : public Analysis {
  private:
    using index_type = MoleculeData::index_type;
    std::vector<index_type> group_ids;      //!< group ids to save to disk (empty = all)
    std::vector<std::size_t> group_indices; //!< indices of groups to save (active AND inactive)
    std::unique_ptr<XTCWriter> writer;
    void _to_json(json& j) const override;
    void _sample() override;
    XTCtraj(const Space& spc, const std::string& filename, const std::vector<std::string>& molecule_names);

  public:
    XTCtraj(const json& j, const Space& spc);
};

/**
 * @brief Excess pressure using virtual volume move
 */
class VirtualVolumeMove : public PerturbationAnalysis {
    Geometry::VolumeMethod volume_scaling_method = Geometry::VolumeMethod::ISOTROPIC;
    double volume_displacement = 0.0;
    void _sample() override;
    void _from_json(const json& j) override;
    void _to_json(json& j) const override;
    void sanityCheck(double old_energy);
    void writeToFileStream(const Point& scale, double energy_change) const;

  public:
    VirtualVolumeMove(const json& j, Space& spc, Energy::EnergyTerm& pot);
};

/**
 * @brief Create histogram of molecule conformation id
 */
class MolecularConformationID : public Analysis {
    MoleculeData::index_type molid;        //!< molecule id to sample
    std::map<int, unsigned int> histogram; //!< key = conformation id; value = count
    void _sample() override;
    void _to_json(json& j) const override;

  public:
    MolecularConformationID(const json& j, const Space& spc);
};

/**
 * @brief Virtual translation move to calculate force
 *
 * Displace a single molecule of `molid` with `dL` in the
 * direction `dir` and measure the free energy of the process
 * using dA=-kT*ln<exp(-dU)> and the resulting force, -dA/dL
 *
 * @todo Does this work with Ewald summation? k-vectors must be refreshed.
 */
class VirtualTranslate : public PerturbationAnalysis {
    MoleculeData::index_type molid; //!< molid to operate on
    Point perturbation_direction = {0.0, 0.0, 1.0};
    double perturbation_distance = 0.0;

    void _sample() override;
    void _from_json(const json& j) override;
    void _to_json(json& j) const override;
    double momentarilyPerturb(Space::GroupType& group);
    void writeToFileStream(double energy_change) const;

  public:
    VirtualTranslate(const json& j, Space& spc, Energy::EnergyTerm& pot);
};

/**
 * @brief Multipolar decomposition between groups as a function of separation
 * @date Malmo 2014
 * @todo Add option to use charge center instead of mass center
 */
class MultipoleDistribution : public Analysis {
    struct Data {
        using average_type = Average<double>;
        average_type exact;
        average_type ion_ion;
        average_type ion_dipole;
        average_type ion_quadrupole;
        average_type dipole_dipole;
        average_type dipole_dipole_correlation;
    };
    std::map<MoleculeData::index_type, Data> mean_energy; //!< Energy distributions
    std::vector<std::string> names;                       //!< Molecule names (len=2)
    std::vector<MoleculeData::index_type> ids;            //!< Molecule ids (len=2)
    std::string filename;                                 //!< output file name
    double dr = 0.0;                                      //!< distance resolution

    double g2g(const Space::GroupType& group1,
               const Space::GroupType& group2); //<! exact ion-ion energy between particles
    void save() const;                          //!< save to disk
    void _sample() override;
    void _to_json(json& j) const override;
    void _to_disk() override;

  public:
    MultipoleDistribution(const json& j, const Space& spc);

}; // end of multipole distribution

/**
 * @brief Sample scattering intensity
 */
class ScatteringFunction : public Analysis {
  private:
    enum class Schemes { DEBYE, EXPLICIT_PBC, EXPLICIT_IPBC }; // three different schemes
    Schemes scheme = Schemes::DEBYE;
    bool mass_center_scattering;                        //!< scatter from mass center, only?
    bool save_after_sample = false;                     //!< if true, save average S(q) after each sample point
    std::string filename;                               //!< output file name
    std::vector<Point> scatter_positions;               //!< vector of scattering points
    std::vector<MoleculeData::index_type> molecule_ids; //!< Molecule ids
    std::vector<std::string> molecule_names;            //!< Molecule names corresponding to `molecule_ids`
    using Tformfactor = Scatter::FormFactorUnity<double>;

    std::unique_ptr<Scatter::DebyeFormula<Tformfactor>> debye;
    std::unique_ptr<Scatter::StructureFactorPBC<>> explicit_average_pbc;
    std::unique_ptr<Scatter::StructureFactorIPBC<>> explicit_average_ipbc;
    void _sample() override;
    void _to_disk() override;
    void _to_json(json& j) const override;

  public:
    ScatteringFunction(const json& j, const Space& spc);
};

/*
 * @brief Sample and save gyration eigenvalues of all particles having the same id
 */
class AtomInertia : public Analysis {
  private:
    std::string filename;
    AtomData::index_type atom_id;
    std::ofstream output_stream;

    Point compute();
    void _to_json(json& j) const override;
    void _sample() override;
    void _to_disk() override;

  public:
    AtomInertia(const json& j, const Space& spc);
};

/**
 * @brief Sample and save the eigenvalues of the inertia tensor for a range of indexes within a molecule
 */
class InertiaTensor : public Analysis {
  private:
    std::string filename;                    //!< file to stream to
    std::unique_ptr<std::ostream> stream;    //!< file output stream
    MoleculeData::index_type group_index;    //!< Group to analyse
    std::vector<size_t> particle_range;      //!< range of indexes within the group
    std::pair<Point, Point> compute() const; //!< Compute eigen values and principal axis
    void _to_json(json& j) const override;
    void _sample() override;
    void _to_disk() override;

  public:
    InertiaTensor(const json& j, const Space& spc);
};

/*
 * @brief Sample and save charge, dipole, and quadrupole moments for a range of indexes within a molecule
 */
class MultipoleMoments : public Analysis {
  private:
    std::string filename;
    std::vector<AtomData::index_type> particle_range; //!< range of indexes within the group
    MoleculeData::index_type group_index;
    bool use_molecular_mass_center = true; //!< Moments w.r.t. the COM of the whole molecule (instead of the subgroup)
    std::ofstream output_stream;

    struct Data {
        double charge = 0.0;                // total charge
        Point dipole_moment{0.0, 0.0, 0.0}; // dipole vector
        Point eivals, eivec, center;        // quadrupole eigenvalues and major axis
    } __attribute__((aligned(128)));

    Data calculateMultipoleMoment() const;
    void _to_json(json& j) const override;
    void _sample() override;
    void _to_disk() override;

  public:
    MultipoleMoments(const json& j, const Space& spc);
};

/**
 * @brief Analysis of polymer shape - radius of gyration, shape etc.
 * @date November, 2011
 *
 * The analysis includes:
 * - gyration tensor and radius of gyration
 * - histogram of Rg
 * - end-to-end distance
 * - shape anisotropy
 */
class PolymerShape : public Analysis {
    struct AverageData {
        using average_type = Average<double>;
        average_type gyration_radius_squared;
        average_type gyration_radius;
        average_type end_to_end_squared;
        average_type shape_factor_squared;
        average_type aspherity;
        average_type acylindricity;
        average_type relative_shape_anisotropy;
    }; //!< Placeholder class for average properties

    AverageData data;               //!< Stores all averages
    MoleculeData::index_type molid; //!< Molecule id to analyse
    std::unique_ptr<std::ostream> tensor_output_stream = nullptr;
    Equidistant2DTable<double, unsigned int> gyration_radius_histogram;

    void _to_json(json& j) const override;
    void _sample() override;
    void _to_disk() override;

  public:
    PolymerShape(const json& j, const Space& spc);
};

/**
 * @brief Trajectory with charge and radius, only, for all (active, inactive) particles
 *
 * For use with VMD to visualize charge fluctuations and grand canonical ensembles. Inactive
 * particles have zero charge and radius. If the `filename` ends with `.gz` a GZip compressed
 * file is created.
 */
class QRtraj : public Analysis {
  protected:
    std::function<void()> write_to_file;            //!< Write a single frame to stream
    std::unique_ptr<std::ostream> stream = nullptr; //!< Output stream

  private:
    std::string filename;    //!< Output filename
    void _sample() override; //!< Samples one frame and outputs to stream
    void _to_json(json& j) const override;
    void _to_disk() override;

  public:
    QRtraj(const json& j, const Space& spc, const std::string& name = "qrtraj");
};

/**
 * Generate PSCs text trajectory containing
 *
 * 1. number of particles (incl. inactive)
 * 2. step number and box dimensions
 * 3. midpoint position; direction; patch direction
 *
 * This can be used to generate VMD visualisation using a conversion
 * script (see scripts/ folder)
 */
class PatchySpheroCylinderTrajectory : public QRtraj {
  public:
    PatchySpheroCylinderTrajectory(const json& j, const Space& spc);
};

/**
 * @brief Trajectory with full Space information
 *
 * The following are saved in (compressed) binary form:
 *
 * - all particle properties (id, position, charge, dipole etc.)
 * - all group properties (id, size, capacity etc.)
 *
 * If zlib compression is enabled the file size
 * is reduced by roughly a factor of two.
 *
 * @todo Geometry information; update z-compression detection
 */
class SpaceTrajectory : public Analysis {
  private:
    const Space::GroupVector& groups; // reference to all groups
    std::string filename;
    std::unique_ptr<std::ostream> stream;
    std::unique_ptr<cereal::BinaryOutputArchive> archive;
    void _sample() override;
    void _to_json(json& j) const override;
    void _to_disk() override;
    bool useCompression() const; //!< decide from filename if zlib should be used

  public:
    SpaceTrajectory(const json& j, const Space& spc);
};

class SamplingPolicyBase;
/**
 * @brief Analysis class for sasa calculations
 *
 * Holds a policy which defines type of sampling to be used and
 * is chosen by the user input
 */
class SASAAnalysis : public Analysis {
    using index_type = Faunus::AtomData::index_type;
    using count_type = size_t;
    using table_type = Equidistant2DTable<double, count_type>;

  public:
    enum class Policies { ATOMIC, MOLECULAR, ATOMS_IN_MOLECULE, INVALID };

  private:
    struct Averages {
        using average_type = Average<double>;
        average_type area;
        average_type area_squared;
    };                     //!< Placeholder class for average properties
    Averages average_data; //!< Stores all averages for the selected molecule

    std::unique_ptr<std::ostream> output_stream; //!< output stream

    double probe_radius; //!< radius of the probe sphere
    int slices_per_atom; //!< number of slices of each sphere in SASA calculation

    std::string filename;                         //!< output file name
    table_type sasa_histogram;                    //!< histogram of sasa values
    std::unique_ptr<Faunus::SASA::SASABase> sasa; //!< sasa object for calculating solute areas
    std::unique_ptr<SamplingPolicyBase> policy;   //!< policy specyfing how sampling will be performed

    virtual void _to_json(json& j) const override;
    virtual void _from_json(const json& input) override;
    virtual void _to_disk() override;
    virtual void _sample() override;

    void setPolicy(const Policies);
    void takeSample(const double area);

    friend class SamplingPolicyBase;

  public:
    SASAAnalysis(const json& j, const Space& spc);
    SASAAnalysis(double probe_radius, int slices_per_atom, double resolution, Policies policy, const Space& spc);
};

NLOHMANN_JSON_SERIALIZE_ENUM(SASAAnalysis::Policies, {{SASAAnalysis::Policies::ATOMIC, "atomic"},
                                                      {SASAAnalysis::Policies::MOLECULAR, "molecular"},
                                                      {SASAAnalysis::Policies::ATOMS_IN_MOLECULE, "atoms_in_molecule"},
                                                      {SASAAnalysis::Policies::INVALID, nullptr}})

/** @brief abstract base class for different SASA sampling policies*/
class SamplingPolicyBase {

  protected:
    template <typename TBegin, typename TEnd>
    void sampleIndividualSASA(TBegin first, TEnd last, SASAAnalysis& analysis);
    template <typename TBegin, typename TEnd> void sampleTotalSASA(TBegin first, TEnd last, SASAAnalysis& analysis);

  public:
    virtual void sample(const Space& spc, SASAAnalysis& analysis) = 0;
    virtual void to_json(json& input) const = 0;
    virtual void from_json(const json& input) = 0;

    SamplingPolicyBase() = default;
    virtual ~SamplingPolicyBase() = default;
};

/** @brief SASA sampling policy which samples individually atoms selected by atom type name*/
class AtomicPolicy : public SamplingPolicyBase {

    AtomData::index_type atom_id; //!< id of atom type to be sampled
    std::string atom_name;        //!< name of the atom type to be sampled

    void sample(const Space& spc, SASAAnalysis& analysis) override;
    void to_json(json& input) const override;
    void from_json(const json& input) override;

  public:
    AtomicPolicy() = default;
};

/** @brief SASA sampling policy which samples individually molecules selected by molecules name*/
class MolecularPolicy : public SamplingPolicyBase {

    MoleculeData::index_type molecule_id; //!< id of molecule to be sampled
    std::string molecule_name;            //!< name of the molecule to be sampled

    void sample(const Space& spc, SASAAnalysis& analysis) override;
    void to_json(json& input) const override;
    void from_json(const json& input) override;

  public:
    MolecularPolicy() = default;
};

/** @brief SASA sampling policy which samples as a whole atom in a selected molecule
 *          if multiple atoms are selected (either by atom names or by indices in a selected molecule)
 *          it samples sum of their SASAs in a given molecule
 *          if single atom name is selected, it samples just the selected atom SASA in a selected molecule
 *          */
class AtomsInMoleculePolicy : public SamplingPolicyBase {

    MoleculeData::index_type molecule_id; //!< id of molecule to be sampled
    std::set<size_t> selected_indices;    //!< selected indices of atoms in the chosen molecule
    std::string molecule_name;            //!<  selected molecule name to be sampled
    std::set<std::string> atom_names;     //!< selected names of atoms in the chosen molecule

    void sample(const Space& spc, SASAAnalysis& analysis) override;
    void to_json(json& input) const override;
    void from_json(const json& input) override;

  public:
    AtomsInMoleculePolicy() = default;
};

/** @brief Example analysis */
template <class T, class Enable = void> struct _analyse {
    void sample(T&) { std::cout << "not a dipole!" << std::endl; } //!< Sample
};                                                                 // primary template

/** @brief Example analysis */
template <class T> struct _analyse<T, typename std::enable_if<std::is_base_of<Dipole, T>::value>::type> {
    void sample(T&) { std::cout << "dipole!" << std::endl; } //!< Sample
};                                                           // specialized template

class SavePenaltyEnergy : public Analysis {
  private:
    const std::string filename; //!< Base file name
    int filenumber = 0;         //!< Counter for each saved file
    std::shared_ptr<Energy::Penalty> penalty_energy;
    void _sample() override;

  public:
    SavePenaltyEnergy(const json& j, const Space& spc, const Energy::Hamiltonian& pot);
};

} // namespace analysis

} // namespace Faunus
