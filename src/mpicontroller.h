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
    int masterRank() const;         //!< MPI rank of master process
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
    T M = communicator.size();
    T i = communicator.rank();
    T beg = (N * i) / M;
    T end = (N * i + N) / M - 1;
    return std::pair<T, T>(beg, end);
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
  public:
    enum dataformat { XYZ = 3, XYZQ = 4, XYZQI = 5 };
    void setFormat(dataformat);
    void setFormat(const std::string&);
    dataformat getFormat() const;
    std::vector<double> buffer;
    void copyParticlesToBuffer(const ParticleVector& particles); //!< Copy source particle vector to send buffer
    void copyBufferToParticles(ParticleVector& particles);       //!< Copy receive buffer to target particle vector
  private:
    dataformat format = XYZQI; //!< Data format to send/receive - default is XYZQ
};

/*
 * @brief Sum tables computed by parallel processes
 *
 * @details Slave processes send histograms to the master. The master computes the
 * average and sends it back to the slaves. Ttable can be Table, Table2D or Table3D in auxiliary.h.
 *
 */
template <class Ttable>
void avgTables(const mpl::communicator& communicator, Ttable& table, int& size) {
    if (!mpi.isMaster()) {
        std::vector<double> sendBuf = table.hist2buf(size); // data to be sent
        std::vector<double> recvBuf(sendBuf.size());        // buffer for recieving data

        const auto tag = mpl::tag_t(0);
        const auto layout = mpl::contiguous_layout<double>(sendBuf.size());
        communicator.sendrecv(sendBuf.data(), layout, mpi.masterRank(), tag, recvBuf.data(), layout, mpi.masterRank(),
                              tag);
        table.buf2hist(recvBuf);
    } else {
        std::vector<double> sendBuf = table.hist2buf(size);
        std::vector<double> recvBuf(size);
        auto slaves = ranges::cpp20::views::iota(0, communicator.size()) |
                      ranges::cpp20::views::filter([&](auto rank) { return rank != mpi.masterRank(); });

        ranges::cpp20::for_each(slaves, [&](auto rank) {
            communicator.recv(recvBuf, rank);
            sendBuf.insert(sendBuf.end(), recvBuf.begin(), recvBuf.end());
        });

        table.buf2hist(sendBuf);
        sendBuf = table.hist2buf(size);

        ranges::cpp20::for_each(slaves, [&](auto rank) { communicator.send(sendBuf, rank); });
    }
}
#endif

} // namespace Faunus::MPI
