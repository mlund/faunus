#include "mpicontroller.h"
#include <range/v3/algorithm/copy.hpp>
#include <iostream>

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
    std::function<void(const Particle&, decltype(buffer)::iterator&)> copy_to_buffer;
    switch (format) {
    case XYZ: // copy x, y, z
        buffer.resize(3 * particles.size());
        copy_to_buffer = [&](auto& particle, auto& destination) {
            ranges::cpp20::copy(particle.pos, destination);
            std::advance(destination, 3);
        };
        break;
    case XYZQ: // copy x, y, z, charge
        buffer.resize(4 * particles.size());
        copy_to_buffer = [&](auto& particle, auto& destination) {
            ranges::cpp20::copy(particle.pos, destination);
            std::advance(destination, 3);
            *(destination++) = particle.charge;
        };
        break;
    case XYZQI: // copy x, y, z, charge, id
        buffer.resize(5 * particles.size());
        copy_to_buffer = [&](auto& particle, auto& destination) {
            ranges::cpp20::copy(particle.pos, destination);
            std::advance(destination, 3);
            *(destination++) = particle.charge;
            *(destination++) = static_cast<double>(particle.id);
        };
        break;
    }
    auto destination = buffer.begin(); // set *after* buffer resize
    ranges::cpp20::for_each(particles, std::bind(copy_to_buffer, std::placeholders::_1, std::ref(destination)));
    if (destination != buffer.end()) {
        throw std::runtime_error("buffer mismatch");
    }
}

void ParticleBuffer::copyBufferToParticles(ParticleVector& particles) {
    std::function<void(Particle&)> copy_to_particle;
    auto source = buffer.begin();
    switch (format) {
    case XYZ:
        copy_to_particle = [&](auto& particle) {
            std::copy(source, source + 3, particle.pos.begin());
            std::advance(source, 3);
        };
        break;
    case XYZQ:
        copy_to_particle = [&](auto& particle) {
            std::copy(source, source + 3, particle.pos.begin());
            std::advance(source, 3);
            particle.charge = *source++;
        };
        break;
    case XYZQI:
        copy_to_particle = [&](auto& particle) {
            std::copy(source, source + 3, particle.pos.begin());
            std::advance(source, 3);
            particle.charge = *source++;
            particle.id = static_cast<AtomData::index_type>(*source++);
        };
        break;
    }
    ranges::cpp20::for_each(particles, copy_to_particle);
    if (source != buffer.end()) {
        throw std::runtime_error("buffer mismatch");
    }
}

Controller mpi; //!< Global instance of MPI controller

#endif

} // namespace Faunus::MPI
