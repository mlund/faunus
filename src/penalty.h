#pragma once

#include "mpicontroller.h"
#include "externalpotential.h"
#include "reactioncoordinate.h"
#include "aux/table_1d.h"

namespace Faunus::Energy {

/**
 * `sum_of_energy_increments` is the total change of updating the energy function. If
 * not handled this will appear as an energy drift (which it is!). To
 * avoid this, this term is added to the energy but since it's the
 * same in both the trial and old state energies it will not affect
 * MC move acceptance.
 */
class Penalty : public EnergyTerm {
  private:
    const Space& spc;
    std::string penalty_function_filename;
    std::string histogram_filename;
    bool overwrite_penalty = true; // overwrite input penalty function?
    void loadPenaltyFunction(const std::string& filename);
    void toDisk(); //!< Save penalty function and histogram to disk
    void initializePenaltyFunction(const json& j);
    void to_json(json& j) const override;
    virtual void updatePenalty(const std::vector<double>& coordinate);

  protected:
    typedef typename std::shared_ptr<ReactionCoordinate::ReactionCoordinateBase> Tcoord;
    bool verbose = false;                      //!< kæft op?
    bool avoid_energy_drift = true;            //!< avoid energy drift when upgrading penalty function
    size_t number_of_reaction_coordinates = 0; //!< number of reaction coordinate
    size_t update_counter = 0;                 //!< number of calls to `sync()`
    size_t number_of_steps_between_updates;    //!< update frequency [steps]
    size_t penalty_function_exchange_counter = 0;
    size_t samplings;
    double sum_of_energy_increments = 0;                //!< total energy change of updating penalty function
    double energy_increment_scaling_factor = 1.0;       //!< scaling factor for f0
    double energy_increment = 0.0;                      //!< penalty increment
    std::vector<Tcoord> reaction_coordinates_functions; //!< vector of reaction coordinate functions (length = 1 or 2)
    std::vector<double> latest_coordinate;              //!< latest reaction coordinate (length = 1 or 2)

    Table<unsigned int> histogram;      //!< count how often a reaction coordinate is visited
    Table<double> penalty_energy;       //!< penalty energy as a function of coordinates
    void logBarrierInformation() const; //!< Add barrier information to output log

  public:
    Penalty(const json& j, const Space& spc);
    virtual ~Penalty(); //!< destruct and save to disk (!)
    double energy(const Change& change) override;
    void sync(EnergyTerm* other, const Change& change) override;
    void streamPenaltyFunction(std::ostream &stream) const;
    void streamHistogram(std::ostream &stream) const;
};

#ifdef ENABLE_MPI
/**
 * @brief Penalty function with MPI exchange
 */
class PenaltyMPI : public Penalty {
  private:
    const MPI::Controller& mpi;
    Eigen::VectorXi weights;                                     //!< array w. minimum histogram counts
    Eigen::VectorXd buffer;                                      //!< receive buffer for penalty functions
    void updatePenalty(const std::vector<double>& coordinate) override; //!< Average penalty function across all nodes
    void averagePenaltyFunctions();                              //!< Average penalty functions over all MPI nodes
  public:
    PenaltyMPI(const json& j, Space& spc, const MPI::Controller& mpi);
};
#endif

} // namespace Faunus::Energy
