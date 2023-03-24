#pragma once

#include "geometry.h"
#include "space.h"

namespace Faunus {

namespace Energy {
class Hamiltonian;
class NonbondedBase;
} // namespace Energy
class XTCWriter;

/**
 * @brief Performs a task or "action" on the system, typically before or after a simulation
 */
struct SystemAction {
    virtual void operator()(Space& space, Energy::Hamiltonian& hamiltonian) = 0;
    virtual ~SystemAction() = default;
};

/**
 * @brief Rotate and translate two molecules to explore all poses
 *
 * For each pose, the energy is calculated and streamed to a file
 * along with two quaternions that describe rotations from the initial
 * orientations. A trajectory with the poses can be optionally saved.
 */
class AngularScan : public SystemAction {
    inline static const std::string name = "angular scan"; //!< Name used for logging

    /// @brief Helper class to analyse (free) energies
    class EnergyAnalysis {
        double partition_sum = 0;    //!< Partition function (per COM separation)
        double thermal_energy = 0;   //!< Thermal energy sum for each COM separation run
        Average<double> free_energy; //!< Free energy
      public:
        void clear();                    //!< Zeros all data
        void add(double energy);         //!< Add sample point
        double getFreeEnergy() const;    //!< w = -ln < exp(-energy/kT) >
        double getThermalEnergy() const; //!< <u> = ∑ u * exp(-energy/kT) / Q
        void info() const;               //!< Print to global logger
    };

    /// @brief Helper class for store information about each of the two rigid bodies
    struct Molecule {
        Space::GroupVector::size_type index;                          //!< Group index in `Space::groups`
        std::vector<Point> ref_positions;                             //!< Original reference positions of particles
        void initialize(const Space::GroupVector& groups, int index); //!< Set molecule index and reset
        ParticleVector getRotatedReference(const Space::GroupVector& groups, const Eigen::Quaterniond& q);
    };

    double zmin = 0;                //!< Minimum mass center separation (along z)
    double zmax = 0;                //!< Maximum mass center separation (along z)
    double dz = 0;                  //!< z resolution
    double max_energy = pc::infty;  //!< Skip poses with energies higher than this
    EnergyAnalysis energy_analysis; //!< Helper to calculate the (free) energy

    std::pair<Molecule, Molecule> molecules; //!< The two molecules to scan
    std::unique_ptr<std::ostream> stream;    //!< Output file with poses
    std::unique_ptr<XTCWriter> trajectory;   //!< Output each pose to Gromacs trajectory file
    Geometry::TwobodyAngles angles;          //!< Helper class handling angular space

    /// Calculate energy and report to stream and trajectory
    void report(const Group& group1, const Group& group2, const Eigen::Quaterniond& q1, const Eigen::Quaterniond& q2,
                Energy::NonbondedBase& nonbonded);

  public:
    AngularScan(const json& input, const Space& spc);
    void operator()(Space& space, Energy::Hamiltonian& hamiltonian) override;
};

/// @brief Create single action from JSON input
std::unique_ptr<SystemAction> createAction(std::string_view name, const json& properties, Space& spc);

/// @brief Create vector of actions from JSON list input
std::vector<std::unique_ptr<SystemAction>> createActionList(const json& input, Space& spc);

} // namespace Faunus
