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

Partner::Partner(PartnerPolicy policy) : policy(policy) {}

bool Partner::isValid(const mpl::communicator& mpi, int partner) {
    return (partner >= 0 && partner < mpi.size() && partner != mpi.rank());
}

Partner::PartnerPair Partner::getPair(const mpl::communicator& mpi) const {
    if (rank.has_value()) {
        // note `std::minmax(a,b)` takes _references_; the initializer list (used here) takes a _copy_
        return std::minmax({mpi.rank(), rank.value()});
    }
    throw std::runtime_error("bad partner");
}

OddEvenPartner::OddEvenPartner() : Partner(PartnerPolicy::ODDEVEN) {}

/**
 * If true is returned, a valid partner was found
 */
bool OddEvenPartner::generate(const mpl::communicator& mpi, Random& random) {
    int rank_increment = static_cast<bool>(random.range(0, 1)) ? 1 : -1;
    if (mpi.rank() % 2 == 0) { // even replica
        rank = mpi.rank() + rank_increment;
    } else { // odd replica
        rank = mpi.rank() - rank_increment;
    }
    if (!isValid(mpi, rank.value())) {
        rank = std::nullopt;
    }
    return rank.has_value();
}

std::unique_ptr<Partner> createMPIPartnerPolicy(PartnerPolicy policy) {
    switch (policy) {
    case PartnerPolicy::ODDEVEN:
        return std::make_unique<OddEvenPartner>();
    default:
        throw std::runtime_error("unknown policy");
    }
}

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

/**
 * @brief Perform the exchange
 * @param mpi MPI Controller
 * @param partner_rank MPI partner to exchange with
 * @param particles The particles to send to partner_rank
 * @return Reference to recieved particles from partner MPI process
 * @todo This involves a lot of copying...
 */
const ParticleVector& ExchangeParticles::operator()(const Controller& mpi, int partner_rank,
                                                    const ParticleVector& particles) {
    if (!partner_particles) {
        partner_particles = std::make_unique<ParticleVector>();
    }
    partner_particles->resize(particles.size());

    particle_buffer.copyToBuffer(particles); // particle data -> vector of doubles
    mpi.world.sendrecv_replace(particle_buffer.begin(), particle_buffer.end(), partner_rank, mpl::tag_t(0),
                               partner_rank, mpl::tag_t(0));
    particle_buffer.copyFromBuffer(*partner_particles);
    return *partner_particles;
}

/**
 * @param comm MPI communicator
 * @param partner_rank Destination and source
 * @param particles Particle vector to send/recieve
 */
void ExchangeParticles::replace(const mpl::communicator& comm, int partner_rank, ParticleVector& particles) {
    particle_buffer.copyToBuffer(particles); // particle data -> vector of doubles
    comm.sendrecv_replace(particle_buffer.begin(), particle_buffer.end(), partner_rank, mpl::tag_t(0),
                               partner_rank, mpl::tag_t(0));
    particle_buffer.copyFromBuffer(particles);
}

ParticleBuffer::Format ExchangeParticles::getFormat() const { return particle_buffer.getFormat(); }

void ExchangeParticles::setFormat(ParticleBuffer::Format format) { particle_buffer.setFormat(format); }

bool exchangeVolume(const Controller& mpi, int partner_rank, Geometry::GeometryBase& geometry,
                    Geometry::VolumeMethod& volume_scaling_method) {
    const auto old_volume = geometry.getVolume();
    auto new_volume = 0.0;
    mpi.world.sendrecv(old_volume, partner_rank, mpl::tag_t(0), new_volume, partner_rank, mpl::tag_t(0));
    if (new_volume <= pc::epsilon_dbl) {
        mpi.world.abort(1);
    }
    if (std::fabs(new_volume - old_volume) > pc::epsilon_dbl) {
        geometry.setVolume(new_volume, volume_scaling_method);
        return true;
    }
    return false;
}

/**
 * @param comm MPI communicator
 * @param random Random number object to test. Will be propagated.
 * @return True if all random number engines are in sync, i.e. returns the same values
 */
bool checkRandomEngineState(const mpl::communicator& comm, Random& random) {
    double random_number = random();
    std::vector<double> buffer(comm.size(), 0.0);
    comm.gather(0, random_number, buffer.data());
    return std::adjacent_find(buffer.begin(), buffer.end(), std::not_equal_to<>()) == buffer.end();
}

Controller mpi; //!< Global instance of MPI controller

#endif

} // namespace Faunus::MPI
