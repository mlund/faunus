#pragma once

#include "average.h"
#include "mpicontroller.h"
#include "molecule.h"
#include "geometry.h"
#include "space.h"
#include "io.h"
#include "aux/timers.h"
#include <range/v3/view/filter.hpp>
#include <optional>

namespace Faunus {

namespace Energy {
class Hamiltonian;
}

namespace Move {

class MoveBase {
  private:
    virtual void _move(Change& change) = 0;                    //!< Perform move and modify change object
    virtual void _accept(Change& change);                      //!< Call after move is accepted
    virtual void _reject(Change& change);                      //!< Call after move is rejected
    virtual void _to_json(json& j) const = 0;                  //!< Extra info for report if needed
    virtual void _from_json(const json& j) = 0;                //!< Extra info for report if needed
    TimeRelativeOfTotal<std::chrono::microseconds> timer;      //!< Timer for whole move
    TimeRelativeOfTotal<std::chrono::microseconds> timer_move; //!< Timer for _move() only
  protected:
    Space& spc;                                  //!< Space to operate on
    unsigned long number_of_attempted_moves = 0; //!< Counter for total number of move attempts
    unsigned long number_of_accepted_moves = 0;
    unsigned long number_of_rejected_moves = 0;

  public:
    static Random slump; //!< Shared for all moves
    std::string name;    //!< Name of move
    std::string cite;    //!< Reference, preferable a short-doi, e.g. "doi:10/b9jq"
    int repeat = 1;      //!< How many times the move should be repeated per sweep
    size_t steps_between_samples = 1; //!< Run interval for defused moves (with weight = 0)

    void from_json(const json& j);
    void to_json(json& j) const; //!< JSON report w. statistics, output etc.
    void move(Change& change);   //!< Perform move and modify given change object
    void accept(Change& change);
    void reject(Change& change);
    virtual double bias(Change& change, double old_energy,
                        double new_energy); //!< adds extra energy change not captured by the Hamiltonian
    MoveBase(Space& spc, const std::string& name, const std::string& cite);
    inline virtual ~MoveBase() = default;
};

void from_json(const json &, MoveBase &); //!< Configure any move via json
void to_json(json &, const MoveBase &);

/**
 * @brief Replay simulation from a trajectory
 *
 * Particles' positions are updated in every step based on coordinates read from the trajectory. Currently only
 * XTCReader is supported.
 */
class ReplayMove : public MoveBase {
    std::unique_ptr<XTCReader> reader = nullptr; //!< trajectory reader
    TrajectoryFrame frame;                       //!< recently read frame (w/o coordinates)
    bool end_of_trajectory = false;              //!< flag raised when end of trajectory was reached
    // FIXME resolve always accept / always reject on the Faunus level
    const double force_accept = -1e12;           //!< a very negative value of energy difference to force-accept the move

    void _move(Change &change) override;
    void _to_json(json &) const override;
    void _from_json(const json &) override;
    double bias(Change &, double, double) override;

  protected:
    using MoveBase::spc;
    ReplayMove(Space &spc, std::string name, std::string cite);

  public:
    explicit ReplayMove(Space &spc);
};

/**
 * @brief Swap the charge of a single atom
 */
class AtomicSwapCharge : public MoveBase {
    int molid = -1;
    double ln10 = log(10);
    double pKa, pH;
    Average<double> msqd; // mean squared displacement
    double _sqd, _bias;   // squared displament
    std::string molname;  // name of molecule to operate on
    Change::data cdata;

    void _to_json(json &) const override;
    void _from_json(const json &) override; //!< Configure via json object
    typename ParticleVector::iterator randomAtom();
    void _move(Change &change) override;
    double bias(Change &, double, double) override; //!< adds extra energy change not captured by the Hamiltonian
    void _accept(Change &) override;
    void _reject(Change &) override;

  protected:
    using MoveBase::spc;
    AtomicSwapCharge(Space &spc, std::string name, std::string cite);

  public:
    explicit AtomicSwapCharge(Space &spc);
};

/**
 * @brief Histogram for an arbitrary set of values using a sparse memory layout (map)
 *
 * Builds a histogram by binning given values to a specified resolution. Values are stored
 * in a memory efficient map-structure with log(N) lookup complexity.
 */
template <typename T = double> class SparseHistogram {
    T resolution;
    std::map<int, unsigned int> data;

  public:
    explicit SparseHistogram(T resolution) : resolution(resolution) {}
    void add(const T value) {
        if (std::isfinite(value)) {
            data[static_cast<int>(std::round(value / resolution))]++;
        } else {
            faunus_logger->warn("histogram: skipping inf/nan number");
        }
    }
    friend auto& operator<<(std::ostream& stream, const SparseHistogram& histogram) {
        std::for_each(histogram.data.begin(), histogram.data.end(), [&](const auto& sample) {
            stream << fmt::format("{:.6E} {}\n", T(sample.first) * histogram.resolution, sample.second);
        });
        return stream;
    }
};

/**
 * @brief Translate and rotate a molecular group
 */
class AtomicTranslateRotate : public MoveBase {
    Space::Tpvec::const_iterator latest_particle;      //!< Iterator to last moved particle
    const Energy::Hamiltonian& hamiltonian;            //!< Reference to Hamiltonian
    std::map<int, SparseHistogram<>> energy_histogram; //!< Energy histogram (value) for each particle type (key)
    double energy_resolution = 0.0;                    //!< Resolution of sampled energy histogram
    double latest_displacement_squared;                //!< temporary squared displacement
    void sampleEnergyHistogram();                      //!< Update energy histogram based on latest move
    void saveHistograms();                             //!< Write histograms for file
    void checkMassCenter(Space::Tgroup& group) const;  //!< Perform test to see if the move violates PBC
    void groupToDisk(const Space::Tgroup& group) const; //!< Save structure to disk in case of failure

  protected:
    int molid = -1;                           //!< Molecule id to move
    Point directions = {1, 1, 1};             //!< displacement directions
    Average<double> mean_square_displacement; //!< mean squared displacement
    std::string molecule_name;                //!< name of molecule to operate on
    Change::data cdata;                       //!< Data for change object

    void _to_json(json &) const override;
    void _from_json(const json &) override; //!< Configure via json object
    ParticleVector::iterator randomAtom();  //!< Select random particle to move

    virtual void translateParticle(ParticleVector::iterator, double); //!< translate single particle
    void _move(Change &) override;
    void _accept(Change &) override;
    void _reject(Change &) override;

    AtomicTranslateRotate(Space& spc, const Energy::Hamiltonian& hamiltonian, std::string name, std::string cite);

  public:
    AtomicTranslateRotate(Space& spc, const Energy::Hamiltonian& hamiltonian);
    ~AtomicTranslateRotate();
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
class TranslateRotate : public MoveBase {
  protected:
    int molid = -1; //!< Molecule ID of the molecule(s) to move
    Average<double> mean_squared_displacement;
    Average<double> mean_squared_rotation_angle;

    double latest_displacement_squared = 0.0;
    double latest_rotation_angle_squared = 0.0;
    double translational_displacement = 0.0;   //!< User defined displacement parameter
    double rotational_displacement = 0.0;      //!< User defined rotational displacement parameter
    Point translational_direction = {1, 1, 1}; //!< User defined directions along x, y, z
    Point fixed_rotation_axis = {0, 0, 0};     //!< Axis of rotation. 0,0,0 == random.

    std::optional<std::reference_wrapper<Space::Tgroup>> findRandomMolecule() const;
    double translateMolecule(Space::Tgroup &group);
    double rotateMolecule(Space::Tgroup &group);
    void checkMassCenter(const Space::Tgroup& group) const; // sanity check of move

    void _to_json(json &j) const override;
    void _from_json(const json &j) override;
    void _move(Change &change) override;
    void _accept(Change &) override;
    void _reject(Change &) override;

    TranslateRotate(Space &spc, std::string name, std::string cite);

  public:
    explicit TranslateRotate(Space &spc);
};

/**
 * @brief Move that preferentially displaces molecules within a specified region around a specified atom type
 * Idea based on the chapter 'Smarter Monte Carlo' in 'Computer Simulation of Liquids' by Allen & Tildesley (p. 317)
 * The current region implemented is an ellipsoid for which you specify radii along length and width of ellipsoid
 *
 */

class SmartTranslateRotate : public MoveBase {
  protected:
    typedef typename Space::Tpvec Tpvec;
    using MoveBase::spc;

    int molid = -1, refid1 = -1, refid2 = -1; // molecule to displace, reference atoms 1 and 2 defining geometry
    unsigned long cnt;
    double dptrans = 0, dprot = 0;
    double p = 1; // initializing probability that a molecule outside geometry is kept as selected molecule
    double r_x = 0,
           r_y = 0; // defining lengths of perpendicular radii defining the ellipsoid (or sphere if a and b are equal)
    double _sqd;    // squared displacement
    Average<double> msqd, countNin_avg, countNin_avgBlocks, countNout_avg,
        countNout_avgBlocks; // mean squared displacement and particle counters

    double cosTheta, theta;            // geometrical variables
    double x, y;                       // x and y coordinate relative to center of geometry of chosen molecule
    double coord, coordNew, coordTemp; // normalized coordinates to decide if molecule is inside or outside geometry
    double randNbr;
    double _bias = 0, rsd = 0.01, Nin, countNin, countNout, Ntot = 0,
           cntInner = 0; // bias to add when crossing boundary between in and out, counters keeping track of molecules
                         // inside, outside geomtry etc...

    Point dir = {1, 1, 1};
    Point cylAxis = {0, 0, 0}; // axis/vector connecting the two reference atoms
    Point origo = {0, 0, 0};
    Point molV = {0, 0, 0}; // coordinate vector of chosen molecule

    bool findBias = true;

    void _to_json(json &j) const override;
    void _from_json(const json &j) override; //!< Configure via json object
    void _move(Change &change) override;
    double bias(Change &, double, double) override;
    void _accept(Change &) override { msqd += _sqd; }
    void _reject(Change &) override { msqd += 0; }

    SmartTranslateRotate(Space &spc, std::string name, std::string cite);

  public:
    explicit SmartTranslateRotate(Space &spc);
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
class ConformationSwap : public MoveBase {
  public:
    enum CopyPolicy { ALL, POSITIONS, CHARGES, INVALID }; //!< What to copy from conformation library
  private:
    CopyPolicy copy_policy;
    RandomInserter inserter;
    int molid = -1; //!< Molecule ID to operate on
    void copyConformation(ParticleVector& source_particle, ParticleVector::iterator destination) const;
    void _to_json(json &j) const override;
    void _from_json(const json &j) override;
    void _move(Change &change) override;
    void setRepeat();                   //!< Set move repeat
    void checkConformationSize() const; //!< Do conformations fit simulation cell?
    void checkMassCenterDrift(const Point& old_mass_center, const ParticleVector& particles); //!< Check for CM drift
    void registerChanges(Change& change, const Space::Tgroup& group) const;                   //!< Update change object
    ConformationSwap(Space& spc, const std::string& name, const std::string& cite);

  public:
    explicit ConformationSwap(Space &spc);
}; // end of conformation swap move

NLOHMANN_JSON_SERIALIZE_ENUM(ConformationSwap::CopyPolicy, {{ConformationSwap::INVALID, nullptr},
                                                            {ConformationSwap::ALL, "all"},
                                                            {ConformationSwap::POSITIONS, "positions"},
                                                            {ConformationSwap::CHARGES, "charges"}})

class VolumeMove : public MoveBase {
  private:
    Geometry::VolumeMethod volume_scaling_method = Geometry::VolumeMethod::ISOTROPIC;
    Average<double> mean_volume;
    Average<double> mean_square_volume_change;
    double old_volume = 0.0;
    double new_volume = 0.0;
    double logarithmic_volume_displacement_factor = 0.0;

    void _to_json(json &j) const override;
    void _from_json(const json &j) override;
    void _move(Change &change) override;
    void _accept(Change& change) override;
    void _reject(Change& change) override;

  public:
    explicit VolumeMove(Space &spc);
}; // end of VolumeMove

/**
 * @brief Displaces charge on a single atom
 */
class ChargeMove : public MoveBase {
  private:
    typedef typename Space::Tpvec Tpvec;
    Average<double> msqd; // mean squared displacement
    double dq = 0, deltaq = 0;
    int atomIndex;
    Change::data cdata;

    void _to_json(json &j) const override;
    void _from_json(const json &j) override;
    void _move(Change &change) override;
    void _accept(Change &) override;
    void _reject(Change &) override;

  protected:
    using MoveBase::spc;
    ChargeMove(Space &spc, std::string name, std::string cite);

  public:
    ChargeMove(Space &spc);
};

/**
 * @brief Transfers charge between two molecules
 */
class ChargeTransfer : public MoveBase {
  private:
    typedef typename Space::Tpvec Tpvec;
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
        Change::data cdata;
    };

    moldata mol1, mol2;

    double sumTemp = 0;
    int i = 0;

    void _to_json(json &j) const override;
    void _from_json(const json &j) override;
    void _move(Change &change) override;
    void _accept(Change &) override;
    void _reject(Change &) override;

  protected:
    using MoveBase::spc;
    ChargeTransfer(Space &spc, std::string name, std::string cite);

  public:
    explicit ChargeTransfer(Space &spc);
};

/**
 * @brief QuadrantJump translates a molecule to another quadrant
 * considering as the origin the center of the box or the center of mass
 * of a range of atomic indexes specified by "index": [start:stop].
 */
class QuadrantJump : public MoveBase {
  private:
    typedef typename Space::Tpvec Tpvec;
    int molid = -1;
    Point dir = {1, 1, 1};
    std::vector<size_t> index;
    double _sqd;          // squared displacement
    Average<double> msqd; // mean squared displacement

    void _to_json(json &j) const override;
    void _from_json(const json &j) override; //!< Configure via json object
    void _move(Change &change) override;
    void _accept(Change &) override { msqd += _sqd; }
    void _reject(Change &) override { msqd += 0; }

  protected:
    using MoveBase::spc;
    QuadrantJump(Space &spc, std::string name, std::string cite);

  public:
    explicit QuadrantJump(Space &spc);
};

#ifdef ENABLE_MPI

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
class ParallelTempering : public MoveBase {
  private:
    Geometry::VolumeMethod volume_scaling_method = Geometry::VolumeMethod::ISOTROPIC; //!< How to scale volumes
    double very_small_volume = 1e-9;
    MPI::MPIController &mpi;
    std::unique_ptr<ParticleVector> partner_particles;
    Random random;
    int partner = -1;              //!< Exchange replica (partner)
    enum extradata { VOLUME = 0 }; //!< Structure of extra data to send
    std::map<std::string, Average<double>> acceptance_map;

    MPI::FloatTransmitter float_transmitter;                       //!< Class for transmitting floats over MPI
    MPI::ParticleTransmitter<ParticleVector> particle_transmitter; //!< Class for transmitting particles over MPI

    void findPartner(); //!< Find replica to exchange with
    bool goodPartner(); //!< Is partner valid?
    void _to_json(json &j) const override;
    void _move(Change &change) override;
    double exchangeEnergy(double energy_change);              //!< Exchange energy with partner
    void exchangeState(Change &change);                       //!< Exchange positions, charges, volume etc.
    double bias(Change &, double uold, double unew) override; //!< Energy change in partner replica
    std::string id() const;                                   //!< Unique string to identify set of partners
    void _accept(Change &) override;
    void _reject(Change &) override;
    void _from_json(const json &j) override;

  public:
    ParallelTempering(Space &spc, MPI::MPIController &mpi);
    ~ParallelTempering();
};

#endif

/**
 * @brief Class storing a list of MC moves with their probability weights and
 * randomly selecting one.
 */
class Propagator {
  private:
    int _repeat;
    std::discrete_distribution<> distribution;
    BasePointerVector<MoveBase> _moves; //!< list of moves
    std::vector<double> _weights;       //!< list of weights for each move
    void addWeight(double weight = 1);

  public:
    Propagator() = default;
    Propagator(const json &j, Space &spc, Energy::Hamiltonian &pot, MPI::MPIController &mpi);
    auto repeat() const -> decltype(_repeat) { return _repeat; }
    auto moves() const -> const decltype(_moves) & { return _moves; };
    auto sample() {
        if (!_moves.empty()) {
            assert(_weights.size() == _moves.size());
#ifdef ENABLE_MPI
            //!< Avoid parallel processes to get out of sync
            //!< Needed for replica exchange or parallel tempering
            int offset = distribution(MPI::mpi.random.engine);
#else
            int offset = distribution(Move::MoveBase::slump.engine);
#endif
            return _moves.begin() + offset;
        }
        return _moves.end();
    } //!< Pick move from a weighted, random distribution
    auto end() { return _moves.end(); }

    friend void to_json(json &j, const Propagator &propagator);

    /**
     * @brief Range of moves excluded from the `sample()` algorithm above due to ZERO weight
     *
     * Moves added with zero weight are excluded from the `sample()` function but can be
     * accessed through this function. This is used to run these _hidden_ moves exactly
     * once per Monte Carlo sweep.
     */
    auto defusedMoves() {
        return _moves | ranges::cpp20::views::filter([&](auto move) { return move->repeat == 0; });
    }
};

void to_json(json &j, const Propagator &propagator);

} // namespace Move


} // namespace Faunus
