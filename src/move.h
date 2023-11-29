#pragma once

#include "average.h"
#include "mpicontroller.h"
#include "molecule.h"
#include "geometry.h"
#include "space.h"
#include "io.h"
#include "smart_montecarlo.h"
#include "aux/timers.h"
#include "aux/sparsehistogram.h"
#include <range/v3/view/filter.hpp>
#include <range/v3/view/iota.hpp>
#include <range/v3/view/indirect.hpp>
#include <range/v3/algorithm/count_if.hpp>
#include <optional>

namespace Faunus {

namespace Speciation {
class GroupDeActivator;
}

namespace Energy {
class Hamiltonian;
}

namespace move {

class MoveCollection;

/**
 * @brief Base class for all moves (MC, Langevin, ...)
 *
 * This will propagate the system and return a `Change` object that describes which
 * parts of the system that was updated. This object is later used to calculate the energy
 * of the change, but this is not the responsibility of the move classes.
 *
 * @todo Privatize and rename members
 */
class Move {
  private:
    virtual void _move(Change& change) = 0;                    //!< Perform move and modify change object
    virtual void _accept(Change& change);                      //!< Call after move is accepted
    virtual void _reject(Change& change);                      //!< Call after move is rejected
    virtual void _to_json(json& j) const = 0;                  //!< Extra info for report if needed
    virtual void _from_json(const json& j) = 0;                //!< Extra info for report if needed
    TimeRelativeOfTotal<std::chrono::microseconds> timer;      //!< Timer for whole move
    TimeRelativeOfTotal<std::chrono::microseconds> timer_move; //!< Timer for _move() only

    friend MoveCollection;
    unsigned long number_of_accepted_moves = 0;
    unsigned long number_of_rejected_moves = 0;
    unsigned int sweep_interval = 1; //!< Run interval for defused moves (with weight = 0)
    const std::string cite;          //!< Reference, preferable a short-doi, e.g. "doi:10/b9jq"

  protected:
    const std::string name; //!< Name of move
    Space& spc;             //!< Space to operate on
    int repeat = 1;
    unsigned long number_of_attempted_moves = 0; //!< Counter for total number of move attempts

  public:
    static Random slump; //!< Shared for all moves

    void from_json(const json& j);
    void to_json(json& j) const; //!< JSON report w. statistics, output etc.
    void move(Change& change);   //!< Perform move and modify given change object
    void accept(Change& change);
    void reject(Change& change);
    void setRepeat(int repeat);
    virtual double bias(Change& change, double old_energy,
                        double new_energy); //!< Extra energy not captured by the Hamiltonian
    Move(Space& spc, std::string_view name, std::string_view cite);
    inline virtual ~Move() = default;
    bool isStochastic() const; //!< True if move should be called stochastically
    const std::string& getName() const;
};

void from_json(const json&, Move&); //!< Configure any move via json
void to_json(json&, const Move&);

/**
 * @brief Replay simulation from a trajectory
 *
 * Particles' positions are updated in every step based on coordinates read from the trajectory. Currently only
 * XTCReader is supported.
 */
class ReplayMove : public Move {
    std::unique_ptr<XTCReader> reader = nullptr; //!< trajectory reader
    TrajectoryFrame frame;                       //!< recently read frame (w/o coordinates)
    bool end_of_trajectory = false;              //!< flag raised when end of trajectory was reached
    // FIXME resolve always accept / always reject on the Faunus level
    const double force_accept = -1e12; //!< a very negative value of energy difference to force-accept the move

    void _move(Change& change) override;
    void _to_json(json&) const override;
    void _from_json(const json&) override;
    double bias(Change&, double, double) override;

  protected:
    using Move::spc;
    ReplayMove(Space& spc, std::string name, std::string cite);

  public:
    explicit ReplayMove(Space& spc);
};

/**
 * @brief Swap the charge of a single atom
 */
class AtomicSwapCharge : public Move {
    int molid = -1;
    double pKa, pH;
    Average<double> msqd; // mean squared displacement
    double _sqd, _bias;   // squared displament
    std::string molname;  // name of molecule to operate on
    Change::GroupChange cdata;

    void _to_json(json&) const override;
    void _from_json(const json&) override; //!< Configure via json object
    typename ParticleVector::iterator randomAtom();
    void _move(Change& change) override;
    double bias(Change&, double, double) override; //!< adds extra energy change not captured by the Hamiltonian
    void _accept(Change&) override;
    void _reject(Change&) override;

  protected:
    using Move::spc;
    AtomicSwapCharge(Space& spc, std::string name, std::string cite);

  public:
    explicit AtomicSwapCharge(Space& spc);
};

/**
 * @brief Translate and rotate a molecular group
 */
class AtomicTranslateRotate : public Move {
    ParticleVector::const_iterator latest_particle;        //!< Iterator to last moved particle
    const Energy::Hamiltonian& hamiltonian;                //!< Reference to Hamiltonian
    std::map<int, SparseHistogram<>> energy_histogram;     //!< Energy histogram (value) for each particle type (key)
    double energy_resolution = 0.0;                        //!< Resolution of sampled energy histogram
    double latest_displacement_squared;                    //!< temporary squared displacement
    void sampleEnergyHistogram();                          //!< Update energy histogram based on latest move
    void saveHistograms();                                 //!< Write histograms for file
    void checkMassCenter(Space::GroupType& group) const;   //!< Perform test to see if the move violates PBC
    void groupToDisk(const Space::GroupType& group) const; //!< Save structure to disk in case of failure

  protected:
    int molid = -1;                           //!< Molecule id to move
    Point directions = {1, 1, 1};             //!< displacement directions
    Average<double> mean_square_displacement; //!< mean squared displacement
    std::string molecule_name;                //!< name of molecule to operate on
    Change::GroupChange cdata;                //!< Data for change object

    void _to_json(json&) const override;
    void _from_json(const json&) override; //!< Configure via json object
    ParticleVector::iterator randomAtom(); //!< Select random particle to move

    virtual void translateParticle(ParticleVector::iterator, double); //!< translate single particle
    void _move(Change&) override;
    void _accept(Change&) override;
    void _reject(Change&) override;

    AtomicTranslateRotate(Space& spc, const Energy::Hamiltonian& hamiltonian, std::string name, std::string cite);

  public:
    AtomicTranslateRotate(Space& spc, const Energy::Hamiltonian& hamiltonian);
    ~AtomicTranslateRotate() override;
};

/**
 * @brief Translate and rotate an atom on a 2D hypersphere-surface
 * @todo under construction
 */
/*
template<typename Tspace>
   class Atomic2dTranslateRotate : public AtomicTranslateRotate {
       protected:
           typedef AtomicTranslateRotate base;
           using base::spc;

           void translateParticle(Space::Tpvec::iterator p, double dp) override {
               auto &g = spc.groups[base::cdata.index];
               Point oldpos = p->pos;

               Point rtp = xyz2rtp(p->pos); // Get the spherical coordinates of the particle
               double slump_theta = dp * (base::slump() - 0.5); // Get random theta-move
               double slump_phi = dp * (base::slump() - 0.5);   // Get random phi-move

               double scalefactor_theta = spc.geo.getRadius() * sin(rtp.z()); // Scale-factor for theta
               double scalefactor_phi = spc.geo.getRadius();                  // Scale-factor for phi

               Point theta_dir = Point(-sin(rtp.y()), cos(rtp.y()), 0); // Unit-vector in theta-direction
               Point phi_dir = Point(cos(rtp.y()) * cos(rtp.z()), sin(rtp.y()) * cos(rtp.z()),
                                     -sin(rtp.z())); // Unit-vector in phi-direction
               Point xyz = oldpos + scalefactor_theta * theta_dir * slump_theta +
                           scalefactor_phi * phi_dir * slump_phi; // New position
               p->pos = spc.geo.getRadius() * xyz / xyz.norm();   // Convert to cartesian coordinates

               spc.geo.boundary(p->pos);
               base::_sqd = spc.geo.sqdist(oldpos, p->pos); // squared displacement
               if (not g.atomic) {                          // recalc mass-center for non-molecular groups
                   g.cm = Geometry::massCenter(g.begin(), g.end(), spc.geo.getBoundaryFunc(), -g.cm);

#ifndef NDEBUG
                   Point cmbak = g.cm;                             // backup mass center
                   g.translate(-cmbak, spc.geo.getBoundaryFunc()); // translate to {0,0,0}
                   double should_be_zero = spc.geo.sqdist({0, 0, 0}, Geometry::massCenter(g.begin(),
g.end())); if (should_be_zero > 1e-6) throw std::runtime_error("atomic move too large"); else g.translate(cmbak,
spc.geo.getBoundaryFunc()); #endif
               }
           }

       public:
           Atomic2dTranslateRotate(Tspace &spc, std::string name, std::string cite) : base(spc, name, cite) {}
           explicit Atomic2dTranslateRotate(Tspace &spc) : Atomic2dTranslateRotate(spc, "transrot 2d", "") {}
   };*/

/**
 * @brief Translate and rotate a molecular group
 */
class TranslateRotate : public Move {
  protected:
    using OptionalGroup = std::optional<std::reference_wrapper<Space::GroupType>>;
    int molid = -1; //!< Molecule ID of the molecule(s) to move
    void _to_json(json& j) const override;
    TranslateRotate(Space& spc, std::string name, std::string cite);

  private:
    Average<double> mean_squared_displacement;
    Average<double> mean_squared_rotation_angle;

    double latest_displacement_squared = 0.0;
    double latest_rotation_angle_squared = 0.0;
    double translational_displacement = 0.0;   //!< User defined displacement parameter
    double rotational_displacement = 0.0;      //!< User defined rotational displacement parameter
    Point translational_direction = {1, 1, 1}; //!< User defined directions along x, y, z
    Point fixed_rotation_axis = {0, 0, 0};     //!< Axis of rotation. 0,0,0 == random.

    virtual OptionalGroup findRandomMolecule();
    double translateMolecule(Space::GroupType& group);
    double rotateMolecule(Space::GroupType& group);
    void checkMassCenter(const Space::GroupType& group) const; // sanity check of move

    void _from_json(const json& j) override;
    void _move(Change& change) override;
    void _accept(Change&) override;
    void _reject(Change&) override;

  public:
    explicit TranslateRotate(Space& spc);
};

/**
 * @brief Smart Monte Carlo version of molecular translation and rotation
 *
 * The important steps in modifying the original TranslateRotate is to
 * replace:
 *
 * 1. `findRandomMolecule()` which performs the group selection
 * 2. `bias()` as we need to correct for the asymmetric sampling
 */
class SmarterTranslateRotate : public TranslateRotate {
  private:
    SmarterMonteCarlo::MoveSupport<Group> smartmc;
    void _to_json(json& j) const override;
    OptionalGroup findRandomMolecule() override;
    double bias(Change& change, double old_energy, double new_energy) override;

  public:
    SmarterTranslateRotate(Space& spc, const json& j);
};

/**
 * @brief Move that will swap conformation of a molecule
 *
 * This will swap between different molecular conformations
 * as defined in `MoleculeData` with `traj` and `weight`.
 * If defined, the weight
 * distribution is respected, otherwise all conformations
 * have equal intrinsic weight. Upon insertion, the new conformation
 * is randomly oriented and placed on top of the mass-center of
 * an exising molecule. That is, there is no mass center movement.
 *
 * @todo Add feature to align molecule on top of an exiting one
 */
class ConformationSwap : public Move {
  public:
    enum class CopyPolicy { ALL, POSITIONS, CHARGES, PATCHES, INVALID }; //!< What to copy from conformation library
  private:
    CopyPolicy copy_policy;
    RandomInserter inserter;
    int molid = -1; //!< Molecule ID to operate on
    void copyConformation(ParticleVector& source_particle, ParticleVector::iterator destination) const;
    void _to_json(json& j) const override;
    void _from_json(const json& j) override;
    void _move(Change& change) override;
    void setRepeat();                   //!< Set move repeat
    void checkConformationSize() const; //!< Do conformations fit simulation cell?
    void checkMassCenterDrift(const Point& old_mass_center, const ParticleVector& particles); //!< Check for CM drift
    void registerChanges(Change& change, const Space::GroupType& group) const;                //!< Update change object
    ConformationSwap(Space& spc, const std::string& name, const std::string& cite);

  public:
    explicit ConformationSwap(Space& spc);
}; // end of conformation swap move

NLOHMANN_JSON_SERIALIZE_ENUM(ConformationSwap::CopyPolicy, {{ConformationSwap::CopyPolicy::INVALID, nullptr},
                                                            {ConformationSwap::CopyPolicy::ALL, "all"},
                                                            {ConformationSwap::CopyPolicy::PATCHES, "patches"},
                                                            {ConformationSwap::CopyPolicy::POSITIONS, "positions"},
                                                            {ConformationSwap::CopyPolicy::CHARGES, "charges"}})

class VolumeMove : public Move {
  protected:
    Geometry::VolumeMethod volume_scaling_method = Geometry::VolumeMethod::ISOTROPIC;
    double logarithmic_volume_displacement_factor = 0.0;
    double old_volume = 0.0;
    double new_volume = 0.0;
    void _move(Change& change) override;
    void _from_json(const json& j) override;

  private:
    Average<double> mean_volume;
    Average<double> mean_square_volume_change;

    virtual void setNewVolume();
    void _to_json(json& j) const override;
    void _accept(Change& change) override;
    void _reject(Change& change) override;

  public:
    VolumeMove(Space& spc, std::string_view name);
    explicit VolumeMove(Space& spc);
}; // end of VolumeMove

/**
 * @brief Displaces charge on a single atom
 */
class ChargeMove : public Move {
  private:
    Average<double> mean_squared_charge_displacement;
    Change::GroupChange group_change;

    void _move(Change& change) override;
    void _accept(Change&) override;
    void _reject(Change&) override;
    void _from_json(const json& j) override;
    virtual double getChargeDisplacement(const Particle& particle) const;
    ChargeMove(Space& spc, std::string_view name, std::string_view cite);

  protected:
    double max_charge_displacement = 0.0;
    double charge_displacement = 0.0;
    ParticleVector::size_type particle_index;
    void _to_json(json& j) const override;

  public:
    explicit ChargeMove(Space& spc);
};

/**
 * @brief As `ChargeMove` but performs displacements in squared charge
 *
 * The displacement is q' = sqrt(q^2 + dq^2) and the following bias
 * must be added to enforce symmetry: u/kT = ln( |q'/q| )
 */
class QuadraticChargeMove : public ChargeMove {
  private:
    Average<double> mean_bias;
    void _to_json(json& j) const override;
    double getChargeDisplacement(const Particle& particle) const override;
    double bias(Change& change, double old_energy, double new_energy) override;

  public:
    explicit QuadraticChargeMove(Space& spc);
};

/**
 * @brief Transfers charge between two molecules
 */
class ChargeTransfer : public Move {
  private:
    Average<double> msqd; // mean squared displacement
    double dq = 0, deltaq = 0;

    struct moldata {
        double charges = 0;
        double moves = 0;
        int numOfAtoms = 0;
        int id = 0;
        std::string molname;
        std::vector<double> min, max;
        std::vector<double> molrange;
        std::vector<double> ratio;
        std::vector<double> changeQ;
        Change::GroupChange cdata;
    };

    moldata mol1, mol2;

    double sumTemp = 0;
    int i = 0;

    void _to_json(json& j) const override;
    void _from_json(const json& j) override;
    void _move(Change& change) override;
    void _accept(Change&) override;
    void _reject(Change&) override;

  protected:
    using Move::spc;
    ChargeTransfer(Space& spc, std::string name, std::string cite);

  public:
    explicit ChargeTransfer(Space& spc);
};

/**
 * @brief QuadrantJump translates a molecule to another quadrant
 * considering as the origin the center of the box or the center of mass
 * of a range of atomic indexes specified by "index": [start:stop].
 */
class QuadrantJump : public Move {
  private:
    int molid = -1;
    Point dir = {1, 1, 1};
    std::vector<size_t> index;
    double _sqd;          // squared displacement
    Average<double> msqd; // mean squared displacement

    void _to_json(json& j) const override;
    void _from_json(const json& j) override; //!< Configure via json object
    void _move(Change& change) override;
    void _accept(Change&) override { msqd += _sqd; }
    void _reject(Change&) override { msqd += 0; }

  protected:
    using Move::spc;
    QuadrantJump(Space& spc, std::string name, std::string cite);

  public:
    explicit QuadrantJump(Space& spc);
};

#ifdef ENABLE_MPI

/**
 * @brief Helper class for the Gibbs ensemble generalized for multi-component systems
 *
 * This is the base used for volume and matter exchange in order to determine phase
 * co-existence. Additional information:
 *
 * - [_Phase equilibria by simulation in the Gibbs ensemble_](https://dx.doi.org/10/cvzgw9)
 * - Frenkel and Smith, 2nd Ed., Chapter 8
 */
class GibbsEnsembleHelper {
  public:
    using VectorOfMolIds = std::vector<MoleculeData::index_type>;
    const MPI::Controller& mpi;
    int partner_rank = -1;                                          //!< Either rank 0 or 1
    double total_volume = 0;                                        //!< Total volume of both cells
    int total_num_particles = 0;                                    //! Total number of particles in both cells
    VectorOfMolIds molids;                                          //!< Molecule id's to exchange. Must be molecular.
    std::pair<int, int> currentNumParticles(const Space& spc) const;  //!< Current number of particles in cell 1 and 2
    std::pair<double, double> currentVolumes(const Space& spc) const; //!< Current volumes in cell 1 and 2
    double exchange(double value) const;                            //!< MPI exchange a double with partner
    GibbsEnsembleHelper(const Space& spc, const MPI::Controller& mpi, const VectorOfMolIds& molids);
};

/**
 * @brief Volume exchange move for the Gibbs ensemble (doi:10/cvzgw9)
 */
class GibbsVolumeMove : public VolumeMove {
  private:
    MPI::Controller& mpi;
    std::unique_ptr<GibbsEnsembleHelper> gibbs;
    bool direct_volume_displacement = true; //!< True if direct displacement in V; false if lnV displacement
    void setNewVolume() override;
    void _from_json(const json& j) override;
    bool volumeTooExtreme() const; //!< Check if volume is too small or too large

  protected:
    void _move(Change& change) override;

  public:
    GibbsVolumeMove(Space& spc, MPI::Controller& mpi);
    double bias(Change& change, double old_energy, double new_energy) override;
};

/**
 * @brief Matter (particle) exchange move for the Gibbs ensemble (doi:10/cvzgw9)
 *
 * - Each of the two cells runs as a separate MPI process
 * - Multi-component systems are supported; just add a move instance for each species
 */
class GibbsMatterMove : public Move {
  private:
    bool insert;                    //!< Insert or delete particle?
    MoleculeData::index_type molid; //!< Molid to insert or delete
    MPI::Controller& mpi;
    std::unique_ptr<GibbsEnsembleHelper> gibbs;
    std::unique_ptr<Speciation::GroupDeActivator> molecule_bouncer;
    void _from_json(const json& j) override;
    void _to_json(json& j) const override;

  protected:
    void _move(Change& change) override;

  public:
    GibbsMatterMove(Space& spc, MPI::Controller& mpi);
    double bias(Change& change, double old_energy, double new_energy) override;
};

/**
 * @brief Class for parallel tempering (aka replica exchange) using MPI
 *
 * Although not completely correct, the recommended way of performing a temper move
 * is to do `N` Monte Carlo passes with regular moves and then do a tempering move.
 * This is because the MPI nodes must be in sync and if you have a system where
 * the random number generator calls are influenced by the Hamiltonian we could
 * end up in a deadlock.
 *
 * @date Lund 2012, 2018
 */
class ParallelTempering : public Move {
  private:
    const MPI::Controller& mpi;
    MPI::ExchangeParticles exchange_particles; //!< Helper class to exchange particles
    std::unique_ptr<MPI::Partner> partner;     //!< Policy for finding MPI partners
    Geometry::VolumeMethod volume_scaling_method = Geometry::VolumeMethod::ISOTROPIC; //!< How to scale volumes
    std::map<MPI::Partner::PartnerPair, Average<double>> acceptance_map;              //!< Exchange statistics
    Random slump; // static instance of Random (shared for all in ParallelTempering)
    void _to_json(json& j) const override;
    void _from_json(const json& j) override;
    void _move(Change& change) override;
    void _accept(Change& change) override;
    void _reject(Change& change) override;
    double bias(Change& change, double uold, double unew) override; //!< Energy change in partner replica
    double exchangeEnergy(double energy_change);                    //!< Exchange energy with partner
    void exchangeState(Change& change);                             //!< Exchange positions, charges, volume etc.
    void exchangeGroupSizes(Space::GroupVector& groups, int partner_rank);

  public:
    explicit ParallelTempering(Space& spc, const MPI::Controller& mpi);
};

#endif

/**
 * Factory for creating instances of fully constructed moves based on names (string)
 *
 * @returns unique pointer to move
 * @throw if invalid name or input parameters
 */
std::unique_ptr<Move> createMove(const std::string& name, const json& properties, Space& spc,
                                 Energy::Hamiltonian& hamiltonian);

/**
 * @brief Class storing a list of MC moves with their probability weights and
 * randomly selecting one.
 */
class MoveCollection {
  private:
    unsigned int number_of_moves_per_sweep;                //!< Sum of all weights
    BasePointerVector<Move> moves;                         //!< list of moves
    std::vector<double> repeats;                           //!< list of repeats (weights) for `moves`
    std::discrete_distribution<unsigned int> distribution; //!< Probability distribution for `moves`
    using move_iterator = decltype(moves.vec)::iterator;   //!< Iterator to move pointer
    move_iterator sample();                                //!< Pick move from a weighted, random distribution

  public:
    MoveCollection(const json& list_of_moves, Space& spc, Energy::Hamiltonian& hamiltonian, Space &old_spc);
    void addMove(std::shared_ptr<Move>&& move);                     //!< Register new move with correct weight
    const BasePointerVector<Move>& getMoves() const;                //!< Get list of moves
    friend void to_json(json& j, const MoveCollection& propagator); //!< Generate json output

    /**
     * Generates a range of repeated, randomized move pointers guaranteed to be valid.
     * All moves with `repeat=0` are excluded. Effectively each registered move is called
     * with a probability proportional to it's `repeat` value. The random picking of moves is repeated
     * $\sum repeat_i$ times so that running over the range constitutes a complete MC "sweep".
     *
     * @returns Range of valid move pointers
     */
    auto repeatedStochasticMoves() {
        auto is_valid_and_stochastic = [&](auto move) { return move < moves.end() && (*move)->isStochastic(); };
        return ranges::cpp20::views::iota(0U, number_of_moves_per_sweep) |
               ranges::cpp20::views::transform([&]([[maybe_unused]] auto count) { return sample(); }) |
               ranges::cpp20::views::filter(is_valid_and_stochastic) | ranges::views::indirect; // dereference iterator
    }

    /**
     * @brief Range of moves excluded from the `sample()` algorithm above due to ZERO weight
     *
     * Moves added with zero weight are excluded from the `sample()` function but can be
     * accessed through this function. This is used to run these moves at a fixed frequency
     * Monte Carlo sweep frequency. Used by e.g. parallel tempering that in the current
     * implementation must be run at fixed intervals due to MPI concurrency.
     *
     * @param sweep_number Current sweep count to decide if move should be included based on `sweep_interval`
     * @returns Range of valid move pointers to be run at given sweep number, i.e. non-stochastically
     */
    auto constantIntervalMoves(const unsigned int sweep_number = 1) {
        auto is_static_and_time_to_sample = [&, sweep_number](const auto& move) {
            return (!move->isStochastic()) && (sweep_number % move->sweep_interval == 0);
        };
        return moves | ranges::cpp20::views::filter(is_static_and_time_to_sample);
    }
};

void to_json(json& j, const MoveCollection& propagator);

} // namespace move
} // namespace Faunus
