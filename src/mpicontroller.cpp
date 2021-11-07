#include "mpicontroller.h"
#include "core.h"
#include <vector>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/join.hpp>
#include <tuple>

namespace Faunus::MPI {

std::string prefix;

#ifdef ENABLE_MPI
double reduceDouble(const mpl::communicator& communicator, double local) {
    double sum = 0.0;
    communicator.allreduce(mpl::plus<double>(), local, sum);
    return sum;
}
int Controller::masterRank() const { return 0; }

bool Controller::isMaster() const { return world.rank() == masterRank(); }

Controller::Controller() : world(mpl::environment::comm_world()) {
    if (world.size() > 1) {
        prefix = fmt::format("mpi{}.", world.rank());
        stream.open((prefix + "stdout"));
    } else {
        prefix.clear();
    }
}

void Controller::to_json(json& j) const {
    j = {{"rank", world.rank()}, {"nproc", world.size()}, {"prefix", prefix}, {"master", masterRank()}};
}

std::ostream& Controller::cout() {
    return stream.is_open() ? stream : std::cout;
}

void ParticleBuffer::setFormat(dataformat d) { format = d; }

void ParticleBuffer::setFormat(const std::string& format) {
    if (format == "XYZQ") {
        setFormat(XYZQ);
    } else if (format == "XYZ") {
        setFormat(XYZ);
    } else {
        setFormat(XYZQI);
    }
}

typename ParticleBuffer::dataformat ParticleBuffer::getFormat() const { return format; }

void ParticleBuffer::copyParticlesToBuffer(const ParticleVector& particles) {
    auto expand = [](const Particle& particle) {
        return std::vector<const double*>({&particle.charge, &particle.charge});
    };
    auto expand2 = [](const Particle& particle) -> std::vector<std::reference_wrapper<const double>> {
        return {particle.charge};
    };
    auto expand3 = [](const Particle& particle) -> std::vector<double> {
        return {particle.charge};
    };

    using namespace ranges::cpp20::views;
    auto flatten = particles | transform(expand2) | ranges::view::join;
    switch (format) {
    case XYZ:
        buffer.resize(3 * particles.size());
        break;
    case XYZQ:
        buffer.resize(4 * particles.size());
        break;
    case XYZQI:
        buffer.resize(5 * particles.size());
        break;
    }
    size_t i = 0;
    for (const auto& particle : particles) {
        buffer.at(i++) = particle.pos.x();
        buffer.at(i++) = particle.pos.y();
        buffer.at(i++) = particle.pos.z();
        if (format == XYZQ) {
            buffer.at(i++) = particle.charge;
        } else if (format == XYZQI) {
            buffer.at(i++) = particle.charge;
            buffer.at(i++) = static_cast<double>(particle.id);
        }
    }
    if (i != buffer.size()) {
        throw std::runtime_error("buffer mismatch");
    }

}

void ParticleBuffer::copyBufferToParticles(ParticleVector& particles) {
    size_t i = 0;
    for (auto& particle : particles) {
        particle.pos.x() = buffer.at(i++);
        particle.pos.y() = buffer.at(i++);
        particle.pos.z() = buffer.at(i++);
        if (format == XYZQ) {
            particle.charge = buffer.at(i++);
        } else if (format == XYZQI) {
            particle.charge = buffer.at(i++);
            particle.id = static_cast<AtomData::index_type>(buffer.at(i++));
        }
    }
    if (i != buffer.size()) {
        throw std::runtime_error("buffer <-> particle mismatch");
    }
}

Controller mpi; //!< Global instance of MPI controller

#endif

} // namespace Faunus::MPI
