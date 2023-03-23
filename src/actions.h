#pragma once

#include "molecule.h"
#include "geometry.h"
#include "space.h"
#include "io.h"

namespace Faunus {

namespace Energy {
class Hamiltonian;
}

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
 * along with two quaternions that describe the rotation from the initial
 * orientation. A trajectory with the poses can be optionally saved.
 */
class AngularScan : public SystemAction {
    struct Molecule {
        Space::GroupVector::size_type index;                          //!< Group index in `Space::groups`
        std::vector<Point> ref_positions;                             //!< Original reference positions of particles
        void initialize(const Space::GroupVector& groups, int index); //!< Set molecule index and reset
        ParticleVector getRotatedReference(const Space::GroupVector& groups, const Eigen::Quaterniond& q);
    };

    double zmin = 0;                         //!< Minimum mass center separation (along z)
    double zmax = 0;                         //!< Maximum mass center separation (along z)
    double dz = 0;                           //!< z resolution
    double max_energy = pc::infty;           //!< Skip poses with energies higher than this
    std::pair<Molecule, Molecule> molecules; //!< The two molecules to scan
    std::unique_ptr<std::ostream> stream;    //!< Output file with poses
    std::unique_ptr<XTCWriter> trajectory;   //!< Output each pose to Gromacs trajectory file
    Geometry::TwobodyAnglesState angles;     //!< Helper class handling angular space
                                             // void initialize(const Space& spc, int molecule_index);

  public:
    AngularScan(const json& input, const Space& spc);
    void operator()(Space& space, Energy::Hamiltonian& hamiltonian) override;
};

/// @brief Create single action from JSON input
std::unique_ptr<SystemAction> createAction(std::string_view name, const json& properties, Space& spc);

/// @brief Create vector of actions from JSON list input
std::vector<std::unique_ptr<SystemAction>> createActionList(const json& input, Space& spc);

} // namespace Faunus
