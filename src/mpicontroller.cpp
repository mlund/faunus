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

bool Controller::isMaster() const { return world.rank() == master_rank; }

Controller::Controller() : world(mpl::environment::comm_world()) {
    if (world.size() > 1) {
        prefix = fmt::format("mpi{}.", world.rank());
        stream.open((prefix + "stdout"));
    } else {
        prefix.clear();
    }
}

void Controller::to_json(json& j) const {
    j = {{"rank", world.rank()}, {"nproc", world.size()}, {"prefix", prefix}, {"master", master_rank}};
}

std::ostream& Controller::cout() { return stream.is_open() ? stream : std::cout; }

void ParticleBuffer::setFormat(const std::string& format_string) {
    if (format_string == "XYZQI") {
        format = Format::XYZQI;
    } else if (format_string == "XYZQ") {
        format = Format::XYZQ;
    } else if (format_string == "XYZ") {
        format = Format::XYZ;
    } else {
        throw std::runtime_error("unknown format: " + format_string);
    }
}

void ParticleBuffer::copyParticlesToBuffer(const ParticleVector& particles) {
    std::function<void(const Particle&, decltype(buffer)::iterator&)> copy_to_buffer;
    switch (format) {
    case Format::XYZ: // copy x, y, z
        packet_size = 3;
        copy_to_buffer = [&](auto& particle, auto& destination) {
            ranges::cpp20::copy(particle.pos, destination);
            std::advance(destination, 3);
        };
        break;
    case Format::XYZQ: // copy x, y, z, charge
        packet_size = 4;
        copy_to_buffer = [&](auto& particle, auto& destination) {
            ranges::cpp20::copy(particle.pos, destination);
            std::advance(destination, 3);
            *(destination++) = particle.charge;
        };
        break;
    case Format::XYZQI: // copy x, y, z, charge, id
        packet_size = 5;
        copy_to_buffer = [&](auto& particle, auto& destination) {
            ranges::cpp20::copy(particle.pos, destination);
            std::advance(destination, 3);
            *(destination++) = particle.charge;
            *(destination++) = static_cast<double>(particle.id);
        };
        break;
    }
    buffer.resize(packet_size * particles.size());
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
    case Format::XYZ:
        copy_to_particle = [&](auto& particle) {
            std::copy(source, source + 3, particle.pos.begin());
            std::advance(source, 3);
        };
        break;
    case Format::XYZQ:
        copy_to_particle = [&](auto& particle) {
            std::copy(source, source + 3, particle.pos.begin());
            std::advance(source, 3);
            particle.charge = *source++;
        };
        break;
    case Format::XYZQI:
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
int ParticleBuffer::packetSize() const { return packet_size; }

Controller mpi; //!< Global instance of MPI controller

#endif

} // namespace Faunus::MPI
