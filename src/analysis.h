#pragma once

#include "space.h"
#include "io.h"
#include "scatter.h"
#include "reactioncoordinate.h"
#include "aux/timers.h"
#include "aux/table_2d.h"
#include "aux/equidistant_table.h"
#include "aux/pairwise_iterator.h"
#include <range/v3/to_container.hpp>
#include <set>

namespace cereal {
class BinaryOutputArchive;
}

namespace Faunus {

namespace Energy {
class Hamiltonian;
class Energybase;
} // namespace Energy

namespace Analysis {

/**
 * @brief Base class for all analysis functions
 *
 * Analysis classes are used to perform system analysis
 * over the course or a simulations or for post-analysis
 * of an existing trajectory. The base class adds basic
 * functionality such as timing, number of steps, json IO
 * and enforce a common interface.
 */
class Analysisbase {
  private:
    virtual void _to_json(json &) const;                  //!< provide json information
    virtual void _from_json(const json &);                //!< setup from json
    virtual void _sample() = 0;                           //!< perform sample event
    virtual void _to_disk();                              //!< save sampled data to disk
    int number_of_steps = 0;                              //!< counter for total number of steps
    int number_of_skipped_steps = 0;                      //!< steps to skip before sampling (do not modify)
    TimeRelativeOfTotal<std::chrono::microseconds> timer; //!< time to benchmark `_sample()`

  protected:
    int sample_interval = 0;   //!< Steps in between each sample point (do not modify)
    int number_of_samples = 0; //!< counter for number of samples

  public:
    std::string name;             //!< descriptive name
    std::string cite;             //!< url, doi etc. describing the analysis
    void to_json(json &) const;   //!< JSON report w. statistics, output etc.
    void from_json(const json &); //!< configure from json object
    void to_disk();               //!< Save data to disk (if defined)
    void sample();                //!< Increase step count and sample
    int getNumberOfSteps() const; //!< Number of steps
    virtual ~Analysisbase() = default;
};

void to_json(json &, const Analysisbase &);

/*
 * @brief Sample and save reaction coordinates to a file
 */
class FileReactionCoordinate : public Analysisbase {
  private:
    Average<double> avg;
    std::string type, filename;
    std::unique_ptr<std::ostream> stream = nullptr;
    std::shared_ptr<ReactionCoordinate::ReactionCoordinateBase> rc = nullptr;

    void _to_json(json &j) const override;
    void _sample() override;
    void _to_disk() override;

  public:
    FileReactionCoordinate(const json &, Space &);
};

/**
 * @brief Excess chemical potential of molecules
 *
 * @todo While `inserter` is currently limited to random
 * insertion, the code is designed for arbitrary insertion
 * schemes inheriting from `MoleculeInserter`.
 */
class WidomInsertion : public Analysisbase {
    Space &space;
    Energy::Hamiltonian &hamiltonian;           //!< Potential energy method
    std::shared_ptr<MoleculeInserter> inserter; //!< Insertion method
    int number_of_insertions;                   //!< Number of insertions per sample event
    int molid;                                  //!< Molecule id
    bool absolute_z_coords = false;             //!< Apply abs() on all inserted z coordinates?
    Average<double> exponential_average;        //!< Widom average, <exp(-dU/kT)>
    Change change;

    void selectGhostGroup(); //!< Select inactive group to act as group particle
    void _sample() override; //!< Called for each sample event
    void _to_json(json &) const override;
    void _from_json(const json &) override;

  public:
    WidomInsertion(const json &, Space &, Energy::Hamiltonian &);
};

/**
 * @brief Density of atom along axis
 *
 * Calculates the summed density of `atoms` in spherical, cylindrical or planar shells around
 * `origo` which by default is the center of the simulation box
 */
class AtomProfile : public Analysisbase {
    Space &spc;
    Equidistant2DTable<double, double> tbl;
    std::vector<std::string> names; // atom names to analyse
    std::set<int> ids;              // atom ids to analyse
    std::string file;               // output filename
    Point ref = {0, 0, 0};
    Eigen::Vector3i dir = {1, 1, 1};
    double dr; // radial resolution
    bool count_charge = false;
    std::string atom_com;
    int id_com = -1; // center at COM of id_com atoms?

    void _from_json(const json &j) override;
    void _to_json(json &j) const override;
    void _to_disk() override;
    void _sample() override;

  public:
    AtomProfile(const json &j, Space &spc);
};

/**
 * @brief Measures the density of atoms along z axis
 */
class SlicedDensity : public Analysisbase {
    Space &spc;
    Table2D<double, unsigned int> N; // N(z)
    std::vector<std::string> names;
    std::vector<int> ids;
    std::string file;
    double dz;
    std::string atom_com;
    int id_com = -1; // center at COM of id_com atoms?

    void _from_json(const json &j) override;
    void _to_json(json &j) const override;
    void _sample() override;
    void _to_disk() override;

  public:
    SlicedDensity(const json &j, Space &spc);
};

/**
 * @brief Analysis of particle densities
 */
class Density : public Analysisbase {
    Space &spc;
    typedef Equidistant2DTable<unsigned int, double> Ttable; // why double?

    std::map<int, Ttable> swpdhist; // Probability density of swapping atoms
    std::map<int, Ttable> atmdhist; // Probability density of atomic molecules
    std::map<int, Ttable> moldhist; // Probability density of polyatomic molecules
    std::map<int, Average<double>> rho_mol, rho_atom;
    std::map<int, int> Nmol, Natom;
    Average<double> Lavg, Vavg, invVavg;

    // int capacity_limit = 10; // issue warning if capacity get lower than this

    void _sample() override;
    void _to_json(json &) const override;
    void _to_disk() override;

  public:
    Density(const json &, Space &);
};

class ChargeFluctuations : public Analysisbase {
  private:
    Space &spc;
    typename decltype(Faunus::molecules)::const_iterator mol_iter; // selected molecule type

    std::vector<std::map<int, int>> idcnt; // populations of types of atomic indexes
    std::vector<Average<double>> charge;   // average charges of atomic indexes
    std::string file;                      // name of PQR file with average charges
    bool verbose;                          // set to true for more output

    void _sample() override;

    void _to_json(json &) const override;

    /**
     * @brief Saves average PQR file to disk if `pqrfile` input is given
     */
    void _to_disk() override;

  public:
    ChargeFluctuations(const json &j, Space &spc);
}; // Fluctuations of atomic charges

/**
 * @brief Molecular multipole moments and their fluctuations
 */
class Multipole : public Analysisbase {
    const Space &spc;
    struct Data {
        Average<double> charge;
        Average<double> charge_squared;
        Average<double> dipole_moment;
        Average<double> dipole_moment_squared;
    };                                   // Average sample moment for a molecule
    std::map<int, Data> average_moments; //!< Molecular moments and their fluctuations. Key = molid.
    void _sample() override;
    void _to_json(json &) const override;

  public:
    Multipole(const json &, const Space &);
};

/**
 * @brief Save system energy to disk
 */
class SystemEnergy : public Analysisbase {
  private:
    std::string file_name, separator = " ";
    std::unique_ptr<std::ostream> output_stream = nullptr;
    std::function<std::vector<double>()> energyFunc;
    Average<double> mean_energy, mean_squared_energy;
    std::vector<std::string> names_of_energy_terms;
    Table2D<double, double> energy_histogram; // Density histograms
    double initial_energy = 0.0;

    void normalize();
    void _sample() override;
    void _to_json(json &) const override;
    void _from_json(const json &) override;
    void _to_disk() override;

  public:
    SystemEnergy(const json &, Energy::Hamiltonian &);
};

/**
 * @brief Checks if system is sane. If not, abort program.
 */
class SanityCheck : public Analysisbase {
  private:
    Space &spc;
    void _sample() override;

  public:
    SanityCheck(const json &, Space &);
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
class SaveState : public Analysisbase {
  private:
    std::function<void(const std::string &)> writeFunc = nullptr;
    bool save_random_number_generator_state = false;
    bool use_numbered_files = true;
    bool convert_hexagonal_prism_to_cuboid = false;
    std::string filename;
    void _to_json(json &) const override;
    void _sample() override;

  public:
    SaveState(json, Space &);
    ~SaveState();
};

/**
 * @brief Base class for distribution functions etc.
 */
class PairFunctionBase : public Analysisbase {
  protected:
    int dim = 3;            // dimentions to use when normalizing
    int id1 = -1, id2 = -1; // particle id (mol or atom)
    double dr = 0;          // distance resolution
    Eigen::Vector3i slicedir = {0, 0, 0};
    double thickness = 0;
    Equidistant2DTable<double, double> hist;
    std::string name1, name2; // atom/molecule names
    std::string file;         // output filename
    double Rhypersphere = -1; // Radius of 2D hypersphere
    Average<double> V;        // average volume (angstrom^3)

  private:
    void _from_json(const json &) override;
    void _to_json(json &) const override;
    void _to_disk() override;

  public:
    PairFunctionBase(const json &);
};

class PairAngleFunctionBase : public PairFunctionBase {
  protected:
    Equidistant2DTable<double, Average<double>> hist2;

  private:
    void _from_json(const json &) override;
    void _to_disk() override;

  public:
    PairAngleFunctionBase(const json &);
};

/** @brief Atomic radial distribution function, g(r) */
class AtomRDF : public PairFunctionBase {
    Space &spc;

    void _sample() override;

  public:
    AtomRDF(const json &, Space &);
};

/** @brief Same as `AtomRDF` but for molecules. Identical input. */
class MoleculeRDF : public PairFunctionBase {
    Space &spc;
    void _sample() override;
  public:
    MoleculeRDF(const json &, Space &);
};

/** @brief Dipole-dipole correlation function, <\boldsymbol{\mu}(0)\cdot\boldsymbol{\mu}(r)> */
class AtomDipDipCorr : public PairAngleFunctionBase {
    Space &spc;
    void _sample() override;
  public:
    AtomDipDipCorr(const json &, Space &);
};

/** @brief Write XTC trajectory file */
class XTCtraj : public Analysisbase {
    std::vector<int> molids;        // molecule ids to save to disk
    std::vector<std::string> names; // molecule names of above
    std::function<bool(Particle &)> filter; // function to filter molecule ids
    Space &spc;
    std::shared_ptr<XTCWriter> writer;

    void _to_json(json &) const override;
    void _from_json(const json &) override;
    void _sample() override;
  public:
    XTCtraj(const json &j, Space &s);
};

/**
 * @brief Excess pressure using virtual volume move
 */
class VirtualVolume : public Analysisbase {
    Space &spc;
    Geometry::VolumeMethod volume_scaling_method = Geometry::ISOTROPIC;
    std::string filename;                                  // output filename (optional)
    std::unique_ptr<std::ostream> output_stream = nullptr; // output file stream
    double dV;                                             // volume perturbation
    Change change;
    Energy::Energybase &pot;
    Average<double> mean_exponentiated_energy_change; // < exp(-du/kT) >

    void _sample() override;
    void _from_json(const json &) override;
    void _to_json(json &) const override;
    void _to_disk() override;

  public:
    VirtualVolume(const json &, Space &, Energy::Energybase &);
};

/**
 * @brief Pressure analysis using the virial theorem
 *
 * This calculates the excess pressure tensor defined as
 * @f[
 * \mathcal{P} = \frac{1}{3V}\left <
 * \sum_{i}^{N-1} \sum_{j=i+1}^N \mathbf{r}_{ij} \otimes \mathbf{f}_{ij}
 * \right >_{NVT}
 * @f]
 *
 * @todo Under construction. Was in Faunus v1 but later abandoned
 */
template <typename forcefunctor> class VirialPressure : public Analysisbase {
    using Tgroup = Space::Tgroup;
    Space &spc;
    Tensor pressure_tensor;
    std::vector<int> user_excluded_molids;          //!< User defined molids excluded from internal pressure
    std::vector<int> groups_with_internal_pressure; //!< Index to `spc.groups[]`

    void _from_json(const json &) override;

    void _to_json(json &j) const override {
        j["no_internal"] =
            user_excluded_molids |
            ranges::cpp20::views::transform([](auto molid) { return Faunus::molecules.at(molid).name; }) |
            ranges::to<std::vector<std::string>>();
    }
    void _to_disk() override;

    Tensor distance_x_force(const Particle &particle1, const Particle &particle2) const {
        const Point distance = spc.geo.vdist(particle1.pos, particle2.pos);
        const Point force = forcefunctor(particle1, particle2);
        return distance * force.transpose();
    }; // todo: get this from Hamiltonian...

    Tensor group_to_group(const Tgroup &group1, const Tgroup &group2) const {
        Tensor pressure_tensor;
        pressure_tensor.setZero();
        for (const auto &particle_i : group1) {
            for (const auto &particle_j : group2) {
                pressure_tensor += distance_x_force(particle_i, particle_j);
            }
        }
        return pressure_tensor;
    }

    Tensor group_internal(const Tgroup &group) const {
        Tensor pressure_tensor;
        pressure_tensor.setZero();
        for (auto particle1 = group.begin(); particle1 != group.end(); ++particle1) {
            for (auto particle2 = particle1; ++particle2 != group.end();) {
                pressure_tensor += distance_x_force(*particle1, *particle2);
            }
        }
        return pressure_tensor;
    }

    void _sample() override {
        // contributions from internal pressure
        for (const auto index : groups_with_internal_pressure) {
            pressure_tensor += group_internal(spc.groups.at(index));
        }
        // contributions from group-group interactions
        for (const auto [group_i, group_j] : PairwiseIterator::internal_pairs(spc.groups)) {
            pressure_tensor += group_to_group(group_i, group_j);
        }
    }

  public:
    VirialPressure(const json &j, Space &spc, Energy::Energybase &pot) {
        pressure_tensor.setZero();
        const auto excluded_molecules = j.value("no_internal", std::vector<std::string>());
        user_excluded_molids = Faunus::names2ids(Faunus::molecules, excluded_molecules);
        std::sort(user_excluded_molids.begin(), user_excluded_molids.end());

        // Find groups with internal pressure (not user excluded ⋀ not rigid ⋁ marked compressible)
        auto excluded = [&](const auto &group) {
            return std::binary_search(user_excluded_molids.begin(), user_excluded_molids.end(), group.id);
        };
        int index = 0;
        for (const Space::Tgroup &group : spc.groups) {
            if (excluded(group)) {
                faunus_logger->debug("{}: excluding {} ({}) from internal pressure", name, group.traits().name, index);
            } else {
                if (not group.traits().rigid or group.traits().compressible) {
                    groups_with_internal_pressure.push_back(index);
                }
            }
            index++;
        }
    }
};

/**
 * @brief Create histogram of molecule conformation id
 */
class MolecularConformationID : public Analysisbase {
    Space &spc;
    int molid;                             //!< molecule id to sample
    std::map<int, unsigned int> histogram; //!< key is conformation id; value is count
    void _sample() override;
    void _to_json(json &j) const override;

  public:
    MolecularConformationID(const json &j, Space &spc);
};

/**
 * @brief Virtual translation move to calculate force
 *
 * Displace a single molecule of `molid` with `dL` in the
 * direction `dir` and measure the free energy of the process
 * using dA=-kT*ln<exp(-dU)> and the resulting force, -dA/dL
 */
class VirtualTranslate : public Analysisbase {
    Change change;         //!< Change object for energy calc.
    Change::data data;     //!< Change data for molecule
    std::string file;      //!< output filename
    int molid;             //!< molid to operate on
    Point dir = {0, 0, 1}; //!< direction to move
    double dL = 0;         //!< distance perturbation
    Energy::Energybase &pot;
    Space &spc;
    Average<double> average_exp_du; //!< <exp(-du/kT)>
    std::ofstream output_file;      // output filestream

    void _sample() override;
    void _from_json(const json &) override;
    void _to_json(json &) const override;
    void _to_disk() override;

  public:
    VirtualTranslate(const json &, Space &, Energy::Energybase &);
};

/**
 * @brief Multipolar decomposition between groups as a function of separation
 * @date Malmo 2014
 * @todo Add option to use charge center instead of mass center
 */
class MultipoleDistribution : public Analysisbase {
    typedef typename Space::Tgroup Tgroup;

    struct data {
        Average<double> exact, ii, id, iq, dd, mucorr;
    };

    std::vector<std::string> names; //!< Molecule names (len=2)
    std::vector<int> ids;           //!< Molecule ids (len=2)
    std::string filename;           //!< output file name
    // int id1, id2;                   //!< pair of molecular id's to analyse
    double dr;             //!< distance resolution
    std::map<int, data> m; //!< Energy distributions
    Space &spc;

    double g2g(const Tgroup &g1, const Tgroup &g2); //<! exact ion-ion energy between particles
    void save() const;                              //!< save to disk
    void _sample() override;
    void _to_json(json &j) const override;
    void _to_disk() override;

  public:
    MultipoleDistribution(const json &j, Space &spc);

}; // end of multipole distribution

/**
 * @brief Sample scattering intensity
 */
class ScatteringFunction : public Analysisbase {
  private:
    enum Schemes { DEBYE, EXPLICIT_PBC, EXPLICIT_IPBC}; // three different schemes
    Schemes scheme = DEBYE;
    Space &spc;
    bool use_com;                   // scatter from mass center, only?
    bool save_after_sample = false; // if true, save average S(q) after each sample point
    std::string filename;           // output file name
    std::vector<Point> p;           // vector of scattering points
    std::vector<int> ids;           // Molecule ids
    std::vector<std::string> names; // Molecule names
    typedef Scatter::FormFactorUnity<double> Tformfactor;

    std::shared_ptr<Scatter::DebyeFormula<Tformfactor>> debye;
    std::shared_ptr<Scatter::StructureFactorPBC<>> explicit_average_pbc;
    std::shared_ptr<Scatter::StructureFactorIPBC<>> explicit_average_ipbc;
    void _sample() override;
    void _to_disk() override;
    void _to_json(json &j) const override;

  public:
    ScatteringFunction(const json &j, Space &spc);
};

/*
 * @brief Sample and save gyration eigenvalues of all particles having the same id
 */
class AtomInertia : public Analysisbase {
  private:
    Space &spc;
    std::string filename;
    int index; // atom id
    std::ofstream file;

    Point compute();
    void _to_json(json &j) const override;
    void _sample() override;
    void _to_disk() override;

  public:
    AtomInertia(const json &j, Space &spc);
};

/*
 * @brief Sample and save the eigenvalues of the inertia tensor for a range of indexes within a molecule
 */
class InertiaTensor : public Analysisbase {
  private:
    Space &spc;
    std::string filename;
    std::vector<size_t> indexes; // range of indexes within the group
    int index; // group indes
    std::ofstream file;
    struct Data {
        Point eivals, eivec; // eigenvalues and principal axis
    };

    Data compute();
    void _to_json(json &j) const override;
    void _sample() override;
    void _to_disk() override;

  public:
    InertiaTensor(const json &j, Space &spc);
};

/*
 * @brief Sample and save charge, dipole and quadrupole moments for a range of indexes within a molecule
 */
class MultipoleMoments : public Analysisbase {
  private:
    Space &spc;
    std::string filename;
    std::vector<size_t> indexes; // range of indexes within the group
    size_t index; // group index
    bool mol_cm;
    std::ofstream file;
    struct Data {
        int q = 0; // total charge
        Point mu {0,0,0}; // dipole vector
        Point eivals, eivec, center; // quadrupole eigenvalues and major axis
    };

    Data compute();
    void _to_json(json &j) const override;
    void _sample() override;
    void _to_disk() override;

  public:
    MultipoleMoments(const json &j, Space &spc);
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
class PolymerShape : public Analysisbase {
    struct AverageData {
        using average_type = Average<double>;
        average_type gyration_radius_squared;
        average_type gyration_radius;
        average_type end_to_end_squared;
        average_type shape_factor_squared;
        average_type aspherity;
        average_type acylindricity;
        average_type relative_shape_anisotropy;
    };                //!< Placeholder class for average polymer properties
    AverageData data; //!< Stores all averages
    Equidistant2DTable<double, unsigned int> gyration_radius_histogram;
    int molid; //!< Molecule id to analyse
    Space &spc;
    std::unique_ptr<std::ostream> tensor_output_stream = nullptr; //!< Output file for tensor

    void _to_json(json &j) const override;
    void _sample() override;
    void _to_disk() override;

  public:
    PolymerShape(const json &j, Space &spc);
};

/**
 * @brief Trajectory with charge and radius, only, for all (active, inactive) particles
 *
 * For use with VMD to visualize charge fluctuations and grand canonical ensembles. Inactive
 * particles have zero charge and radius. If the `filename` ends with `.gz` a GZip compressed
 * file is created.
 */
class QRtraj : public Analysisbase {
  private:
    std::string filename;                           //!< Output filename
    std::unique_ptr<std::ostream> stream = nullptr; //!< Output stream
    std::function<void()> write_to_file;            //!< Write a single frame to stream
    void _sample() override;                        //!< Samples one frame and outputs to stream
    void _to_json(json &j) const override;
    void _to_disk() override;

  public:
    QRtraj(const json &, Space &spc);
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
 * @todo Geometry information
 */
class SpaceTrajectory : public Analysisbase {
  private:
    Space::Tgvec &groups; // reference to all groups
    std::string filename;
    std::unique_ptr<std::ostream> stream;
    std::unique_ptr<cereal::BinaryOutputArchive> archive;
    void _sample() override;
    void _to_json(json &j) const override;
    void _to_disk() override;
    bool useCompression() const; //!< decide from filename if zlib should be used

  public:
    SpaceTrajectory(const json &, Space::Tgvec &);
};

struct CombinedAnalysis : public BasePointerVector<Analysisbase> {
    CombinedAnalysis(const json &j, Space &spc, Energy::Hamiltonian &pot);
    void sample();
    void to_disk(); // prompt all analysis to safe to disk if appropriate
}; //!< Aggregates analysis

/** @brief Example analysis */
template <class T, class Enable = void> struct _analyse {
    void sample(T &) { std::cout << "not a dipole!" << std::endl; } //!< Sample
};                                                                  // primary template

/** @brief Example analysis */
template <class T> struct _analyse<T, typename std::enable_if<std::is_base_of<Dipole, T>::value>::type> {
    void sample(T &) { std::cout << "dipole!" << std::endl; } //!< Sample
};                                                            // specialized template

} // namespace Analysis

} // namespace Faunus
