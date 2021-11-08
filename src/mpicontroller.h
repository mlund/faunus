#pragma once

#include <string>

#ifdef ENABLE_MPI
#include "random.h"
#include "core.h"
#include "particle.h"

#include <vector>
#include <fstream>
#include <range/v3/view/iota.hpp>
#include <range/v3/view/filter.hpp>
#include <range/v3/algorithm/for_each.hpp>
#include <mpl/mpl.hpp>

// Expose classes for MPI serialization
MPL_REFLECTION(Faunus::Point, x(), y(), z())
MPL_REFLECTION(Faunus::Particle, id, charge, pos)
#endif

namespace Faunus::MPI {

/**
 * @brief Filename prefix for MPI related input and output (empty if no MPI)
 *
 * Used to generate rank-based I/O for MPI processes. If the number
 * of ranks is above 2, `prefix = "mpi{rank}."` or otherwise empty.
 * It's a good habbit to append `MPI::prefix` to analysis output in e.g. the `Analysis`
 * namespace.
 */
extern std::string prefix;

#ifdef ENABLE_MPI

/**
 * @brief Main container for MPI related functionality
 *
 * A global instance is available.
 */
class Controller {
  private:
    std::ofstream stream; //!< Redirect stdout to here for rank-based file output
  public:
    Controller();
    const mpl::communicator& world; //!< MPI World communicator
    Random random;                  //!< Random number generator for MPI calls
    const int master_rank = 0;      //!< MPI rank of master process
    bool isMaster() const;          //!< Determines if current rank is the master
    std::ostream& cout();           //!< Rank specific (file) stream
    void to_json(json& j) const;
};

extern Controller mpi;

/**
 * @brief Split N items into nproc parts
 *
 * This returns a pair with the first and last
 * item for the current rank.
 */
template <class T = int> std::pair<T, T> splitEven(const mpl::communicator& communicator, T N) {
    static_assert(std::is_integral_v<T>());
    auto M = static_cast<T>(communicator.size());
    auto i = static_cast<T>(communicator.rank());
    auto beg = (N * i) / M;
    auto end = (N * i + N) / M - 1;
    return {beg, end};
}

/**
 * @brief Reduced sum
 *
 * Each rank sends "local" to master who sums them up.
 * Master sends back (broadcasts) sum to all ranks.
 */
double reduceDouble(const mpl::communicator& communicator, double local);

/**
 * @brief Class for serializing particle vector
 *
 * This is used to serialize select information from a particle vector into a
 * continuous block of memory (std::vector<double>) for use with MPI
 * communication.
 */
class ParticleBuffer {
  private:
    int packet_size = 0; //!< Number of doubles per particle
  public:
    /**
     * @brief Particle information to be copied
     *
     * XYZ -> positions.
     * XYZQ -> positions, charge.
     * XYZQI -> positions, charge, atom id.
     */
    enum class Format { XYZ, XYZQ, XYZQI };
    Format format = Format::XYZQI; //!< Data format to send/receive - default is XYZQ
    void setFormat(const std::string& format_string);
    std::vector<double> buffer;
    void copyParticlesToBuffer(const ParticleVector& particles); //!< Copy source particle vector to send buffer
    void copyBufferToParticles(ParticleVector& particles);       //!< Copy receive buffer to target particle vector
    int packetSize() const;                                      //!< Number of doubles per particle
};

/**
 * @brief Sum tables computed by parallel processes
 *
 * @details Slave processes send histograms to the master. The master computes the
 * average and sends it back to the slaves. Ttable can be Table, Table2D or Table3D in auxiliary.h.
 *
 */
template <class Ttable> void avgTables(const mpl::communicator& communicator, Ttable& table, int& size) {
    std::vector<double> send_buffer; // data to be sent
    std::vector<double> recv_buffer; // buffer for recieving data
    if (!mpi.isMaster()) {
        send_buffer = table.hist2buf(size);
        recv_buffer.resize(send_buffer.size());
        const auto tag = mpl::tag_t(0);
        const auto layout = mpl::contiguous_layout<double>(send_buffer.size());
        communicator.sendrecv(send_buffer.data(), layout, mpi.master_rank, tag, recv_buffer.data(), layout,
                              mpi.master_rank, tag);
        table.buf2hist(recv_buffer);
    } else {
        send_buffer = table.hist2buf(size);
        recv_buffer.resize(size);
        auto slaves = ranges::cpp20::views::iota(0, communicator.size()) |
                      ranges::cpp20::views::filter([&](auto rank) { return rank != mpi.master_rank; });

        ranges::cpp20::for_each(slaves, [&](auto rank) {
            communicator.recv(recv_buffer, rank);
            send_buffer.insert(send_buffer.end(), recv_buffer.begin(), recv_buffer.end());
        });

        table.buf2hist(send_buffer);
        send_buffer = table.hist2buf(size);
        ranges::cpp20::for_each(slaves, [&](auto rank) { communicator.send(send_buffer, rank); });
    }
}
#endif

} // namespace Faunus::MPI
