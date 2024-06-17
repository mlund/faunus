#include "penalty.h"
#include "space.h"
#include <spdlog/spdlog.h>

namespace Faunus::Energy {

Penalty::Penalty(const json& j, const Space& spc) : spc(spc) {
    name = "penalty";
    overwrite_penalty = j.value("overwrite", true);
    energy_increment = j.at("f0").get<double>();
    verbose = !j.value("quiet", true);
    number_of_steps_between_updates = j.at("update").get<size_t>();
    samplings = j.value("samplings", 1);
    avoid_energy_drift = j.value("nodrift", true);
    penalty_function_filename = j.at("file").get<std::string>();
    histogram_filename = j.value("histogram", "penalty-histogram.dat");

    energy_increment_scaling_factor = j.at("scale").get<double>();
    if (energy_increment_scaling_factor < 0 or energy_increment_scaling_factor > 1) {
        throw std::runtime_error("`scale` must be in the interval [0:1]");
    }

    initializePenaltyFunction(j.at("coords"));
    loadPenaltyFunction(MPI::prefix + penalty_function_filename);
}
void Penalty::initializePenaltyFunction(const json& j) {
    if (!j.is_array()) {
        throw ConfigurationError("array of reaction coordinates required");
    }
    std::vector<double> resolutions, minimum_values, maximum_values;
    for (const auto& i : j) {
        auto& reaction_coordinate =
            reaction_coordinates_functions.emplace_back(ReactionCoordinate::createReactionCoordinate(i, spc));
        if (reaction_coordinate->minimum_value >= reaction_coordinate->maximum_value ||
            reaction_coordinate->resolution <= 0) {
            throw ConfigurationError("min<max and resolution>0 required for penalty reaction coordinate");
        }
        resolutions.push_back(reaction_coordinate->resolution);
        minimum_values.push_back(reaction_coordinate->minimum_value);
        maximum_values.push_back(reaction_coordinate->maximum_value);
    }
    number_of_reaction_coordinates = reaction_coordinates_functions.size();
    if (number_of_reaction_coordinates < 1 or number_of_reaction_coordinates > 2) {
        throw std::runtime_error("exactly one or two coordinates required");
    }

    latest_coordinate.resize(number_of_reaction_coordinates, 0);
    histogram.reInitializer(resolutions, minimum_values, maximum_values);
    penalty_energy.reInitializer(resolutions, minimum_values, maximum_values);
}

void Penalty::loadPenaltyFunction(const std::string& filename) {
    if (std::ifstream stream(filename); stream) {
        faunus_logger->info("Loading penalty function {}", filename);
        std::string ignore;
        stream >> ignore >> energy_increment >> samplings >> penalty_function_exchange_counter; // header line
        for (int row = 0; row < penalty_energy.rows(); row++) {
            for (int col = 0; col < penalty_energy.cols(); col++) {
                if (stream.eof()) {
                    throw std::runtime_error("penalty file dimension mismatch");
                }
                stream >> penalty_energy(row, col);
            }
        }
    }
}

/**
 * @todo Foul to let constructor be responsible for i/o...
 */
Penalty::~Penalty() { toDisk(); }

/**
 * @brief Stream penalty function, offset with the minimum observed energy
 */
void Penalty::streamPenaltyFunction(std::ostream& stream) const {
    stream.precision(16);
    stream << fmt::format("# {} {} {}\n", energy_increment, samplings, penalty_function_exchange_counter)
           << penalty_energy.array() - penalty_energy.minCoeff() << "\n";
}

void Penalty::streamHistogram(std::ostream& stream) const {
    stream << histogram;
}

void Penalty::toDisk() {
    if (overwrite_penalty) {
        if (std::ofstream stream(MPI::prefix + penalty_function_filename); stream) {
            streamPenaltyFunction(stream);
        }
    }
    if (std::ofstream stream(MPI::prefix + histogram_filename); stream) {
        stream << histogram << "\n";
    }
}

void Penalty::to_json(json& j) const {
    j["file"] = penalty_function_filename;
    j["scale"] = energy_increment_scaling_factor;
    j["update"] = number_of_steps_between_updates;
    j["nodrift"] = avoid_energy_drift;
    j["histogram"] = histogram_filename;
    j["f0_final"] = energy_increment;
    j["overwrite"] = overwrite_penalty;
    if (penalty_function_exchange_counter > 0) {
        j["mpi exchanges"] = penalty_function_exchange_counter;
    }
    auto& coordinates_j = j["coords"] = json::array();
    for (const auto& reaction_coordinate : reaction_coordinates_functions) {
        coordinates_j.emplace_back(*reaction_coordinate); // `ReactionCoordinateBase` --> `json`
    }
}
double Penalty::energy(const Change& change) {
    double energy = 0.0;
    if (change) {
        for (size_t i = 0; i < number_of_reaction_coordinates; i++) {
            latest_coordinate[i] = reaction_coordinates_functions[i]->operator()();
            if (not reaction_coordinates_functions[i]->inRange(latest_coordinate[i])) {
                return pc::infty; // coordinate outside allowed range -> infinite energy
            }
        }
        penalty_energy.to_index(latest_coordinate);
        energy = penalty_energy[latest_coordinate];
    }
    // reaching here, `latest_coordinate` always reflects the current reaction coordinate
    if (avoid_energy_drift) {
        return energy - sum_of_energy_increments;
    }
    return energy;
}

/**
 * @todo: If this is called before `energy()`, the latest_coordinate
 * is never calculated and causes undefined behavior
 */
void Penalty::updatePenalty(const std::vector<double>& coordinate) {
    update_counter++;
    if (update_counter % number_of_steps_between_updates == 0 and energy_increment > 0.0) {
        if (histogram.minCoeff() >= (int)samplings) {
            const auto minimum_penalty_energy = penalty_energy.minCoeff();    // define min. penalty energy...
            penalty_energy = penalty_energy.array() - minimum_penalty_energy; // ...to zero
            sum_of_energy_increments -= minimum_penalty_energy;
            if (verbose) {
                logBarrierInformation();
            }
            energy_increment *= energy_increment_scaling_factor; // reduce penalty energy
            samplings = std::ceil(double(samplings) / energy_increment_scaling_factor);
            histogram.setZero(); // reset histogram
        }
    }
    latest_coordinate = coordinate;
    histogram[coordinate]++;
    penalty_energy[coordinate] += energy_increment;
    sum_of_energy_increments += energy_increment;
}
void Penalty::logBarrierInformation() const {
    faunus_logger->info("energy barriers: penalty function = {} kT, histogram = {} kT", penalty_energy.maxCoeff(),
                        std::log(double(histogram.maxCoeff()) / histogram.minCoeff()));
}

/**
 *  Called when a move is accepted or rejected, as well as when initializing the system
 *  @todo: this doubles the MPI communication
 */
void Penalty::sync(EnergyTerm* other, [[maybe_unused]] const Change& change) {
    auto* other_penalty = dynamic_cast<decltype(this)>(other);
    if (other_penalty == nullptr) {
        throw std::runtime_error("error in Penalty::sync - please report");
    }
    updatePenalty(other_penalty->latest_coordinate);
    other_penalty->updatePenalty(other_penalty->latest_coordinate); // keep update_counter and samplings in sync
    assert(samplings == other_penalty->samplings);
    assert(latest_coordinate == other_penalty->latest_coordinate);
    assert(update_counter == other_penalty->update_counter);
    assert(sum_of_energy_increments == other_penalty->sum_of_energy_increments);
}

#ifdef ENABLE_MPI

PenaltyMPI::PenaltyMPI(const json& j, Space& spc, const MPI::Controller& mpi) : Penalty(j, spc), mpi(mpi) {
    weights.resize(mpi.world.size());
}

/**
 * @note These are generic MPI calls previously used:
 *
 * MPI_Allgather(&least_sampled_in_histogram, 1, MPI_INT, weights.data(), 1, MPI_INT, mpi.comm);
 */
void PenaltyMPI::updatePenalty(const std::vector<double>& coordinate) {
    const auto old_penalty_energy = penalty_energy[coordinate];
    update_counter++;
    if (update_counter % number_of_steps_between_updates == 0 and energy_increment > 0.0) {
        int least_sampled_in_histogram = histogram.minCoeff(); // if > 0 --> all RC's visited
        mpi.world.allgather(least_sampled_in_histogram, weights.data());

        // if at least one walker has sampled full RC space at least `samplings` times
        if (weights.maxCoeff() > samplings) { // change to minCoeff()?
            averagePenaltyFunctions();
            if (verbose && least_sampled_in_histogram > 0) {
                logBarrierInformation();
            }
            energy_increment *= energy_increment_scaling_factor; // reduce penalty energy
            samplings = std::ceil(samplings / energy_increment_scaling_factor);
            histogram.setZero();
        }
    }
    latest_coordinate = coordinate;
    histogram[coordinate]++;
    penalty_energy[coordinate] += energy_increment;
    sum_of_energy_increments += penalty_energy[coordinate] - old_penalty_energy;
}

/**
 * 1. Master gathers penalty functions from all ranks
 * 2. Master calculates an average, shifter so that minimum energy is at zero
 * 3. Master broadcasts average to all nodes
 * 4. When done, all penalty functions shall be identical
 */
void PenaltyMPI::averagePenaltyFunctions() {
    penalty_function_exchange_counter += 1;
    const auto layout = mpl::contiguous_layout<double>(penalty_energy.size());
    buffer.resize(penalty_energy.size() * mpi.world.size());
    mpi.world.gather(mpi.master_rank, penalty_energy.data(), layout, buffer.data(), layout);
    if (mpi.isMaster()) {
        penalty_energy.setZero();
        for (int i = 0; i < mpi.world.size(); i++) {
            const auto offset = i * penalty_energy.size();
            penalty_energy +=
                Eigen::Map<Eigen::MatrixXd>(buffer.data() + offset, penalty_energy.rows(), penalty_energy.cols());
        }
        penalty_energy = (penalty_energy.array() - penalty_energy.minCoeff()) / static_cast<double>(mpi.world.size());
    }
    mpi.world.bcast(mpi.master_rank, penalty_energy.data(), layout);
}

#endif
} // namespace Faunus::Energy
