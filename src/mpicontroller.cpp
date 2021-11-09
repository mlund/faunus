#include "mpicontroller.h"
#include <range/v3/algorithm/copy.hpp>
#include <algorithm>
#include <iostream>

namespace Faunus::MPI {

std::string prefix;

#ifdef ENABLE_MPI

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

void ParticleBuffer::setFormat(ParticleBuffer::Format data_format) {
    format = data_format;
    switch (format) {
    case Format::XYZ:
        packet_size = 3;
        from_particle = [&](auto& particle, auto& destination) {
            ranges::cpp20::copy(particle.pos, destination);
            std::advance(destination, 3);
        };
        to_particle = [&](auto& source, auto& particle) {
            std::copy(source, source + 3, particle.pos.begin());
            std::advance(source, 3);
        };
        break;
    case Format::XYZQ:
        packet_size = 4;
        from_particle = [&](auto& particle, auto& destination) {
            ranges::cpp20::copy(particle.pos, destination);
            std::advance(destination, 3);
            *(destination++) = particle.charge;
        };
        to_particle = [&](auto& source, auto& particle) {
            std::copy(source, source + 3, particle.pos.begin());
            std::advance(source, 3);
            particle.charge = *source++;
        };
        break;
    case Format::XYZQI:
        packet_size = 5;
        from_particle = [&](auto& particle, auto& destination) {
            ranges::cpp20::copy(particle.pos, destination);
            std::advance(destination, 3);
            *(destination++) = particle.charge;
            *(destination++) = static_cast<double>(particle.id);
        };
        to_particle = [&](auto& source, auto& particle) {
            std::copy(source, source + 3, particle.pos.begin());
            std::advance(source, 3);
            particle.charge = *source++;
            particle.id = static_cast<AtomData::index_type>(*source++);
        };
        break;
    default:
        throw std::runtime_error("unknown format");
    }
}

void ParticleBuffer::copyToBuffer(const ParticleVector& particles) {
    buffer.resize(packet_size * particles.size());
    auto destination = buffer.begin(); // set *after* buffer resize
    ranges::cpp20::for_each(particles, std::bind(from_particle, std::placeholders::_1, std::ref(destination)));
    if (destination != buffer.end()) {
        throw std::runtime_error("buffer mismatch");
    }
}

void ParticleBuffer::copyFromBuffer(ParticleVector& particles) {
    if (buffer.size() != particles.size() * packet_size) {
        throw std::out_of_range("particles out of range");
    }
    auto source = buffer.begin();
    ranges::cpp20::for_each(particles, std::bind(to_particle, std::ref(source), std::placeholders::_1));
    assert(source == buffer.end());
}

ParticleBuffer::Format ParticleBuffer::getFormat() const { return format; }
ParticleBuffer::ParticleBuffer() { setFormat(Format::XYZQI); }
ParticleBuffer::buffer_iterator ParticleBuffer::begin() { return buffer.begin(); }
ParticleBuffer::buffer_iterator ParticleBuffer::end() { return buffer.end(); }

Controller mpi; //!< Global instance of MPI controller

#endif

} // namespace Faunus::MPI
