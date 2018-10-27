#include "particle.h"

namespace Faunus {
    void ParticlePropertyBase::rotate(const Eigen::Quaterniond &q, const Eigen::Matrix3d &) {}

    void Radius::to_json(json &j) const { j["r"] = radius; }

    void Radius::from_json(const json &j) { radius = j.value("r", 0.0); }

    void Charge::to_json(json &j) const { j["q"] = charge; }

    void Charge::from_json(const json &j) { charge = j.value("q", 0.0); }

    void Dipole::rotate(const Eigen::Quaterniond &q, const Eigen::Matrix3d &) {
        mu = q * mu;
    }

    void Dipole::to_json(json &j) const {
        j["mu"] = mu;
        j["mulen"] = mulen;
    }

    void Dipole::from_json(const json &j) {
        mu = j.value("mu", Point(1,0,0) );
        mulen = j.value("mulen", mulen);
    }

    void Quadrupole::rotate(const Eigen::Quaterniond &q, const Eigen::Matrix3d &m) { Q.rotate(m); }

    void Quadrupole::to_json(json &j) const { j["Q"] = Q; }

    void Quadrupole::from_json(const json &j) { Q = j.value("Q", Q); }

    void Cigar::rotate(const Eigen::Quaterniond &q, const Eigen::Matrix3d &) {
        scdir = q * scdir;
    }

    void Cigar::to_json(json &j) const {
        j["scdir"] = scdir;
        j["sclen"] = sclen;
    }

} // end of Faunus namespace
