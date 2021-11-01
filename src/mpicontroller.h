#pragma once
#include "random.h"
#include "core.h"
#include "particle.h"

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <cstdio>

#ifdef ENABLE_MPI
#include "mpl/mpl.hpp"
// Expose classes for MPI serialization
MPL_REFLECTION(Faunus::Point, x(), y(), z())
MPL_REFLECTION(Faunus::Particle, id, charge, pos)
#endif

namespace Faunus::MPI {

extern std::string prefix; //!< Filneme prefix for MPI related input and output (empty if no MPI)

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
    const mpl::communicator& world_comm; //!< MPI World communicator
    Random random;                      //!< Random number generator for MPI calls
    int masterRank() const;             //!< MPI rank of master process
    bool isMaster() const;              //!< Determines if current rank is the master
    std::ostream& cout();               //!< Rank specific (file) stream
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
 * @brief Class for transmitting floating point arrays over MPI
 */
class FloatTransmitter {
  public:
    std::vector<double> swapf(const mpl::communicator& communicator, std::vector<double>&,
                              int); //!< Swap data with another process
    void sendf(const mpl::communicator& communicator, std::vector<double>&, int dst);     //!< Send vector of floats
    void recvf(const mpl::communicator& communicator, int src, std::vector<double>& dst); //!< Receive vector of floats
};

/**
 * @brief Class for sending/receiving particle vectors over MPI.
 *
 * This will take a particle vector and send selected information though MPI. It is
 * possible to send only coordinates using the dataformat `XYZ` or, if charges should be
 * send too, `XYZQ`.
 *
 * Besides particle data it is possible to send extra floats by adding
 * these to the `sendExtra` vector; received extras will be stored in `recvExtra`. Before
 * transmitting extra data, make sure that `recvExtra` and `sendExtra` have the
 * same size.
 *
 *     int dst_rank = 1;
 *     floatp extra1 = 2.34, extra2 = -1.23
 *     Tpvec myparticles(200); // we have 200 particles
 *
 *     Faunus::MPI::MPIController mpi;
 *     Faunus::MPI::ParticleTransmitter pt;
 *
 *     pt.sendExtra.addGroup(extra1);
 *     pt.sendExtra.addGroup(extra2);
 *
 *     pt.send(mpi, myparticles, dst_rank);
 *     pt.waitsend();
 *
 * @date Lund 2012
 * @todo This class needs refactoring for clarity and better safety
 *
 */
template <typename Tpvec> class ParticleTransmitter : public FloatTransmitter {
  public:
    enum dataformat { XYZ = 3, XYZQ = 4, XYZQI = 5 };
    std::vector<double> sendExtra; //!< Put extra data to send here.
    std::vector<double> recvExtra; //!< Received extra data will be stored here
    ParticleTransmitter();
    void send(const mpl::communicator& communicator, const Tpvec&, int); //!< Send particle vector to another node
    void recv(const mpl::communicator& communicator, int, Tpvec&);       //!< Receive particle vector from another node
    void waitrecv();
    void setFormat(dataformat);
    void setFormat(const std::string&);
    dataformat getFormat() const;

  private:
    dataformat format; //!< Data format to send/receive - default is XYZQ
    std::vector<double> sendBuf;
    std::vector<double> recvBuf;
    Tpvec* dstPtr;               //!< pointer to receiving particle vector
    void pvec2buf(const Tpvec&); //!< Copy source particle vector to send buffer
    void buf2pvec(Tpvec&);       //!< Copy receive buffer to target particle vector
};

template <typename Tpvec> ParticleTransmitter<Tpvec>::ParticleTransmitter() { setFormat(XYZQI); }

template <typename Tpvec> void ParticleTransmitter<Tpvec>::setFormat(dataformat d) { format = d; }

template <typename Tpvec> void ParticleTransmitter<Tpvec>::setFormat(const std::string& s) {
    setFormat(XYZQI);
    if (s == "XYZQ") {
        setFormat(XYZQ);
    }
    if (s == "XYZ") {
        setFormat(XYZ);
    }
}

template <typename Tpvec>
typename ParticleTransmitter<Tpvec>::dataformat ParticleTransmitter<Tpvec>::getFormat() const {
    return format;
}

/*!
 * \param mpi MPI controller to use
 * \param src Source particle vector
 * \param dst Destination node
 */
template <typename Tpvec>
void ParticleTransmitter<Tpvec>::send(const mpl::communicator& communicator, const Tpvec& src, int dst) {
    assert(dst >= 0 && dst < communicator.size() && "Invalid MPI destination");
    pvec2buf(src);
    FloatTransmitter::sendf(communicator, sendBuf, dst);
}

template <typename Tpvec> void ParticleTransmitter<Tpvec>::pvec2buf(const Tpvec& src) {
    sendBuf.clear();
    sendBuf.reserve(src.size() * 5 + sendExtra.size());
    for (auto& p : src) {
        sendBuf.push_back(p.pos.x());
        sendBuf.push_back(p.pos.y());
        sendBuf.push_back(p.pos.z());
        if (format == XYZQ) {
            sendBuf.push_back(p.charge);
        }
        if (format == XYZQI) {
            sendBuf.push_back(p.charge);
            sendBuf.push_back(static_cast<double>(p.id));
        }
    }
    for (auto i : sendExtra) {
        sendBuf.push_back(i);
    }
}

/*!
 * \param mpi MPI controller to use
 * \param src Source node
 * \param dst Destination particle vector
 */
template <typename Tpvec>
void ParticleTransmitter<Tpvec>::recv(const mpl::communicator& communicator, int src, Tpvec& dst) {
    assert(src >= 0 && src < communicator.size() && "Invalid MPI source");
    dstPtr = &dst; // save a pointer to the destination particle vector
    if (format == XYZ) {
        recvBuf.resize(3 * dst.size());
    }
    if (format == XYZQ) {
        recvBuf.resize(4 * dst.size());
    }
    if (format == XYZQI) {
        recvBuf.resize(5 * dst.size());
    }

    // resize to fit extra data (if any)
    recvExtra.resize(sendExtra.size());
    recvBuf.resize(recvBuf.size() + recvExtra.size());

    FloatTransmitter::recvf(communicator, src, recvBuf);
}

template <typename Tpvec> void ParticleTransmitter<Tpvec>::waitrecv() {
    buf2pvec(*dstPtr);
}

template <typename Tpvec> void ParticleTransmitter<Tpvec>::buf2pvec(Tpvec& dst) {
    int i = 0;
    for (auto& particle : dst) {
        particle.pos.x() = recvBuf[i++];
        particle.pos.y() = recvBuf[i++];
        particle.pos.z() = recvBuf[i++];
        if (format == XYZQ) {
            particle.charge = recvBuf[i++];
        }
        if (format == XYZQI) {
            particle.charge = recvBuf[i++];
            particle.id = (int)recvBuf[i++];
        }
    }
    for (auto& x : recvExtra) {
        x = recvBuf[i++];
    }
    assert((size_t)i == recvBuf.size());
    if ((size_t)i != recvBuf.size()) {
        std::cerr << "Particle transmitter says: !!!!!!!!!!!" << std::endl;
    }
}

/*
 * @brief Sum tables computed by parallel processes
 *
 * @details Slave processes send histograms to the master. The master computes the
 * average and sends it back to the slaves. Ttable can be Table, Table2D or Table3D in auxiliary.h.
 *
 */
template <class Ttable>
void avgTables(const mpl::communicator& communicator, FloatTransmitter& ft, Ttable& table, int& size) {
    if (!mpi.isMaster()) {
        std::vector<double> sendBuf = table.hist2buf(size);
        std::vector<double> recvBuf = ft.swapf(communicator, sendBuf, mpi.masterRank());
        table.buf2hist(recvBuf);
    }
    if (mpi.isMaster()) {
        std::vector<double> sendBuf = table.hist2buf(size);
        std::vector<double> recvBuf(size);
        for (int i = 0; i < communicator.size(); ++i) {
            if (i != mpi.masterRank()) {
                ft.recvf(communicator, i, recvBuf);
                sendBuf.insert(sendBuf.end(), recvBuf.begin(), recvBuf.end());
            }
        }
        table.buf2hist(sendBuf);
        sendBuf = table.hist2buf(size);
        for (int i = 0; i < communicator.size(); ++i) {
            if (i != mpi.masterRank()) {
                ft.sendf(communicator, sendBuf, i);
                //ft.waitsend();
            }
        }
    }
}
#endif

} // namespace Faunus::MPI
