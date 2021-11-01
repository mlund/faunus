#include "mpicontroller.h"
#include "core.h"
#include <vector>

namespace Faunus::MPI {

std::string prefix;

#ifdef ENABLE_MPI
double reduceDouble(const mpl::communicator& communicator, double local) {
    double sum = 0.0;
    communicator.allreduce(mpl::plus<double>(), local, sum);
    return sum;
}
int Controller::masterRank() const { return 0; }

bool Controller::isMaster() const { return world_comm.rank() == masterRank(); }

Controller::Controller() : world_comm(mpl::environment::comm_world()) {
    if (world_comm.size() > 1) {
        prefix = fmt::format("mpi{}.", world_comm.rank());
        stream.open((prefix + "stdout"));
    } else {
        prefix.clear();
    }
}

void Controller::to_json(json& j) const {
    j = {{"rank", world_comm.rank()}, {"nproc", world_comm.size()}, {"prefix", prefix}, {"master", masterRank()}};
}

std::ostream& Controller::cout() {
    return stream.is_open() ? stream : std::cout;
}

void FloatTransmitter::sendf(const mpl::communicator& communicator, std::vector<double>& src, int dst) {
    communicator.send(src, dst);
    // Generic: MPI_Issend(&src[0], src.size(), MPI_DOUBLE, dst, tag, mpi.comm, &sendReq);
}

void FloatTransmitter::recvf(const mpl::communicator& communicator, int src, std::vector<double>& dst) {
    communicator.recv(dst, src);
    // Generic: MPI_Irecv(&dst[0], dst.size(), MPI_DOUBLE, src, tag, mpi.comm, &recvReq);
}

/**
 * This will send a vector of floats and at the same time wait for the destination process
 * to send back another vector of the same size.
 *
 * @param communicator MPI communicator to use
 * @param src Vector to send
 * @param dst Node to send/receive to/from
 * @return Received data
 */
std::vector<double> FloatTransmitter::swapf(const mpl::communicator& communicator, std::vector<double>& src, int dst) {
    std::vector<double> buffer(src.size());
    auto tag = mpl::tag_t(0);
    auto layout = mpl::contiguous_layout<double>(src.size());
    communicator.sendrecv(src.data(), layout, dst, tag, buffer.data(), layout, dst, tag);
    return buffer;
}

Controller mpi; //!< Global instance of MPI controller

#endif

} // namespace Faunus::MPI
