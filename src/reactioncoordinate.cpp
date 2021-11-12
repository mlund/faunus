#include <doctest/doctest.h>
#include "space.h"
#include "reactioncoordinate.h"
#include "average.h"
#include "multipole.h"
#include <Eigen/Dense>
#include "aux/eigensupport.h"
#include <range/v3/view/filter.hpp>
#include <range/v3/algorithm/count_if.hpp>
#include <range/v3/view/join.hpp>
#include "spdlog/spdlog.h"

namespace Faunus::ReactionCoordinate {

ReactionCoordinateBase::ReactionCoordinateBase(const json& j) {
    auto range = j.value("range", std::vector<double>({0.0, 0.0}));
    if (range.size() != 2 || range[0] > range[1]) {
        throw std::runtime_error(name + ": 'range' requires [min, max>=min]");
    }
    minimum_value = range[0];
    maximum_value = range[1];
    resolution = j.value("resolution", 0.5); // @todo require resolution user input w. at().
}

void ReactionCoordinateBase::_to_json([[maybe_unused]] json& j) const {}

double ReactionCoordinateBase::operator()() {
    assert(function != nullptr);
    return function();
}

bool ReactionCoordinateBase::inRange(double coord) const { return (coord >= minimum_value && coord <= maximum_value); }

void to_json(json& j, const ReactionCoordinateBase& reaction_coordinate) {
    assert(!reaction_coordinate.name.empty());
    auto& _j = j[reaction_coordinate.name];
    _j = {{"range", {reaction_coordinate.minimum_value, reaction_coordinate.maximum_value}},
          {"resolution", reaction_coordinate.resolution}};
    reaction_coordinate._to_json(_j);
}

} // namespace

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("[Faunus] ReactionCoordinateBase") {
    using doctest::Approx;
    Faunus::ReactionCoordinate::ReactionCoordinateBase c(R"({"range":[-1.5, 2.1], "resolution":0.2})"_json);
    CHECK(c.minimum_value == Approx(-1.5));
    CHECK(c.maximum_value == Approx(2.1));
    CHECK(c.resolution == Approx(0.2));
    CHECK(c.inRange(-1.5) == true);
    CHECK(c.inRange(-1.51) == false);
    CHECK(c.inRange(2.11) == false);
    CHECK(c.inRange(2.1) == true);
}
#endif

namespace Faunus::ReactionCoordinate {

/**
 * Factory function for all known penalty functions. The json instance
 * must be an object of exact size one, for example:
 *
 *     atom: {resolution: 0.1, ... }
 */
std::unique_ptr<ReactionCoordinateBase> createReactionCoordinate(const json& j, const Space& spc) {
    try {
        const auto& [key, j_params] = jsonSingleItem(j);
        try {
            if (key == "atom") {
                return std::make_unique<AtomProperty>(j_params, spc);
            }
            if (key == "molecule") {
                return std::make_unique<MoleculeProperty>(j_params, spc);
            }
            if (key == "system") {
                return std::make_unique<SystemProperty>(j_params, spc);
            }
            throw ConfigurationError("unknown reaction coordinate");
        } catch (std::exception& e) {
            usageTip.pick(fmt::format("coords=[{}]", key));
            throw ConfigurationError("'{}': {}", key, e.what());
        }
    } catch (std::exception& e) { throw ConfigurationError("reaction coordinate: {}", e.what()).attachJson(j); }
}

void SystemProperty::_to_json(json &j) const { j["property"] = property; }

SystemProperty::SystemProperty(const json &j, const Space &spc) : ReactionCoordinateBase(j) {
    name = "system";
    property = j.at("property").get<std::string>();
    if (property == "V") {
        function = [&geometry = spc.geometry]() { return geometry.getVolume(); };
    } else if (property == "Lx") {
        function = [&geometry = spc.geometry]() { return geometry.getLength().x(); };
    } else if (property == "Ly") {
        function = [&geometry = spc.geometry]() { return geometry.getLength().y(); };
    } else if (property == "Lz" or property == "height") {
        function = [&geometry = spc.geometry]() { return geometry.getLength().z(); };
    } else if (property == "radius") {
        if (spc.geometry.type == Geometry::Variant::CUBOID or spc.geometry.type == Geometry::Variant::SLIT) {
            faunus_logger->warn("`radius` coordinate unavailable for geometry");
        } else {
            function = [&geometry = spc.geometry]() { return 0.5 * geometry.getLength().x(); };
        }
    } else if (property == "Q") { // system net charge
        function = [&groups = spc.groups]() {
            auto charges = groups | ranges::cpp20::views::join | ranges::cpp20::views::transform(&Particle::charge);
            return std::accumulate(charges.begin(), charges.end(), 0.0);
        };
    } else if (property == "mu") { // system dipole moment
        function = [&groups = spc.groups]() {
            auto particles = groups | ranges::cpp20::views::join;
            return Faunus::dipoleMoment(particles.begin(), particles.end()).norm();
        };
    } else if (property == "mu_x") { // system dipole moment
        function = [&groups = spc.groups]() {
            auto particles = groups | ranges::cpp20::views::join;
            return Faunus::dipoleMoment(particles.begin(), particles.end()).x();
        };
    } else if (property == "mu_y") { // system dipole moment
        function = [&groups = spc.groups]() {
            auto particles = groups | ranges::cpp20::views::join;
            return Faunus::dipoleMoment(particles.begin(), particles.end()).y();
        };
    } else if (property == "mu_z") { // system dipole moment
        function = [&groups = spc.groups]() {
            auto particles = groups | ranges::cpp20::views::join;
            return Faunus::dipoleMoment(particles.begin(), particles.end()).z();
        };
    } else if (property == "N") { // number of particles
        function = [&groups = spc.groups]() {
            auto sizes = groups | ranges::cpp20::views::transform(&Space::GroupType::size);
            return static_cast<double>(std::accumulate(sizes.begin(), sizes.end(), 0U));
        };
    }
    if (function == nullptr) {
        usageTip.pick("coords=[system]");
        throw ConfigurationError("{}: unknown property '{}'", name, property);
    }
}

void AtomProperty::_to_json(json &j) const {
    j["property"] = property;
    j["index"] = index;
    if (dir.squaredNorm() > 1e-9) {
        j["dir"] = dir;
    }
}

AtomProperty::AtomProperty(const json &j, const Space &spc) : ReactionCoordinateBase(j) {
    name = "atom";
    index = j.at("index");
    if (index >= spc.particles.size()) {
        throw ConfigurationError("invalid index");
    }
    property = j.at("property").get<std::string>();
    if (property == "x") {
        function = [&particle = spc.particles.at(index)]() { return particle.pos.x(); };
    } else if (property == "y") {
        function = [&particle = spc.particles.at(index)]() { return particle.pos.y(); };
    } else if (property == "z") {
        function = [&particle = spc.particles.at(index)]() { return particle.pos.z(); };
    } else if (property == "R") {
        function = [&particle = spc.particles.at(index)]() { return particle.pos.norm(); };
    } else if (property == "q") {
        function = [&particle = spc.particles.at(index)]() { return particle.charge; };
    } else if (property == "N") {
        function = [&groups = spc.groups, id = index]() {
            auto particles = groups | ranges::cpp20::views::join;
            return static_cast<double>(
                ranges::cpp20::count_if(particles, [&](const auto& particle) { return particle.id == id; }));
        };
    }

    if (function == nullptr) {
        usageTip.pick("coords=[atom]");
        throw ConfigurationError("{}: unknown property '{}'", name, property);
    }
}

void MoleculeProperty::_to_json(json& j) const {
    j["property"] = property;
    j["index"] = index;
    if (direction.squaredNorm() > 1e-9) {
        j["dir"] = direction;
    }
    if (indexes.size() >= 2) {
        j["indexes"] = indexes;
    }
}

MoleculeProperty::MoleculeProperty(const json &j, const Space &spc) : ReactionCoordinateBase(j) {
    name = "molecule";
    index = j.value("index", 0);
    if (index >= spc.groups.size()) {
        throw ConfigurationError("invalid index");
    }
    auto b = spc.geometry.getBoundaryFunc();
    const auto& group = spc.groups.at(index);

    property = j.at("property").get<std::string>();

    if (property == "active") { // if molecule is active (1) or not (0)
        function = [&group]() { return static_cast<double>(!group.empty()); };
    } else if (property == "confid") {
        function = [&group]() { return static_cast<double>(group.conformation_id); };
    } else if (property == "com_x") {
        function = [&group]() { return group.mass_center.x(); };
    } else if (property == "com_y") {
        function = [&group]() { return group.mass_center.y(); };
    } else if (property == "com_z") {
        function = [&group]() { return group.mass_center.z(); };
    } else if (property == "N") {
        function = [&group]() { return static_cast<double>(group.size()); };
    } else if (property == "Q") {
        function = [&group]() { return monopoleMoment(group.begin(), group.end()); };
    } else if (property == "mu_x") {
        function = [&group, b]() { return dipoleMoment(group.begin(), group.end(), b).x(); };
    } else if (property == "mu_y") {
        function = [&group, b]() { return dipoleMoment(group.begin(), group.end(), b).y(); };
    } else if (property == "mu_z") {
        function = [&group, b]() { return dipoleMoment(group.begin(), group.end(), b).z(); };
    } else if (property == "mu") {
        function = [&group, b]() { return dipoleMoment(group.begin(), group.end(), b).norm(); };
    } else if (property == "end2end") {
        function = [&spc, &group]() {
            return std::sqrt(spc.geometry.sqdist(group.begin()->pos, std::prev(group.end())->pos));
        };
    } else if (property == "Rg") {
        selectGyrationRadius(spc);
    } else if (property == "muangle") {
        selectDipoleAngle(j, spc, b);
    } else if (property == "atomatom") {
        selectAtomAtomDistance(j, spc);
    } else if (property == "cmcm_z") {
        selectMassCenterDistanceZ(j, spc);
    } else if (property == "cmcm") {
        selectMassCenterDistance(j, spc);
    } else if (property == "L/R") {
        selectLengthOverRadiusRatio(j, spc);
    } else if (property == "mindist") {
        selectMinimumGroupDistance(j, spc);
    } else if (property == "Rinner") {
        selectRinner(j, spc);
    } else if (property == "angle") {
        selectAngleWithVector(j, spc);
    }
    if (function == nullptr) {
        usageTip.pick("coords=[molecule]");
        throw ConfigurationError("{}: unknown or impossible property property '{}'", name, property);
    }
}
void MoleculeProperty::selectLengthOverRadiusRatio(const json& j, const Space& spc) {
    direction = j.at("dir");
    indexes = j.value("indexes", decltype(indexes)());
    assert(indexes.size() == 2 && "An array of 2 indexes should be specified.");
    function = [&spc, &dir = direction, i = indexes[0], j = indexes[1]]() {
        Average<double> mean_radius_j;
        Average<double> Rin;
        Average<double> Rout;
        auto particles_i = spc.findAtoms(i);
        auto mass_center_i =
            Geometry::massCenter(particles_i.begin(), particles_i.end(), spc.geometry.getBoundaryFunc());
        for (const auto& particle : spc.findAtoms(j)) {
            mean_radius_j += spc.geometry.vdist(particle.pos, mass_center_i).cwiseProduct(dir.cast<double>()).norm();
        }
        const auto Rjavg = mean_radius_j.avg();
        for (const auto& particle_i : particles_i) {
            const auto radial_distance =
                spc.geometry.vdist(particle_i.pos, mass_center_i).cwiseProduct(dir.cast<double>()).norm();
            if (radial_distance < Rjavg) {
                Rin += radial_distance;
            } else if (radial_distance > Rjavg) {
                Rout += radial_distance;
            }
        }
        return 2.0 * spc.geometry.getLength().z() / (Rin.avg() + Rout.avg());
    };
}
void MoleculeProperty::selectMassCenterDistanceZ(const json& j, const Space& spc) {
    indexes = j.value("indexes", decltype(indexes)());
    assert(indexes.size() > 1 && "An array of 2 or 4 indexes should be specified.");
    if (indexes.size() == 4) {
        function = [&spc, i = indexes[0], j = indexes[1] + 1, k = indexes[2], l = indexes[3] + 1]() {
            auto cm1 = Geometry::massCenter(spc.particles.begin() + i, spc.particles.begin() + j,
                                            spc.geometry.getBoundaryFunc());
            auto cm2 = Geometry::massCenter(spc.particles.begin() + k, spc.particles.begin() + l,
                                            spc.geometry.getBoundaryFunc());
            return spc.geometry.vdist(cm1, cm2).z();
        };
    } else if (indexes.size() == 2) {
        function = [&geometry = spc.geometry, &group1 = spc.groups.at(indexes[0]),
                    &group2 = spc.groups.at(indexes[1])]() {
            return geometry.vdist(group1.mass_center, group2.mass_center).z();
        };
    }
}
void MoleculeProperty::selectAtomAtomDistance(const json& j, const Space& spc) {
    direction = j.at("dir");
    indexes = j.at("indexes").get<decltype(indexes)>();
    if (indexes.size() != 2) {
        throw std::runtime_error("exactly two indices expected");
    }
    function = [&spc, &dir = direction, i = indexes[0], j = indexes[1]]() {
        return spc.geometry.vdist(spc.particles.at(i).pos, spc.particles.at(j).pos)
            .cwiseProduct(dir.cast<double>())
            .norm();
    };
}
void MoleculeProperty::selectGyrationRadius(const Space& spc) {
    function = [&spc, &group = spc.groups.at(index)]() {
        assert(group.size() > 1);
        auto S = Geometry::gyration(group.begin(), group.end(), group.mass_center, spc.geometry.getBoundaryFunc());
        return sqrt(S.trace()); // S.trace() == S.eigenvalues().sum() but faster
    };
}
void MoleculeProperty::selectDipoleAngle(const json& j, const Space& spc, Geometry::BoundaryFunction& b) {
    direction = j.at("dir").get<Point>().normalized();
    if (spc.groups.at(index).isMolecular()) {
        function = [&group = spc.groups.at(index), b, &dir = direction]() {
            const auto dot_product = dipoleMoment(group.begin(), group.end(), b).dot(dir);
            return acos(dot_product) * 180.0 / pc::pi;
        };
    }
}

void MoleculeProperty::selectMassCenterDistance(const json& j, const Space& spc) {
    direction = j.at("dir");
    indexes = j.value("indexes", decltype(indexes)());
    assert(indexes.size() > 1 && "An array of 2 or 4 indexes should be specified.");
    if (indexes.size() == 4) {
        function = [&spc, dir = direction, i = indexes[0], j = indexes[1] + 1, k = indexes[2], l = indexes[3] + 1]() {
            auto cm1 = Geometry::massCenter(spc.particles.begin() + i, spc.particles.begin() + j,
                                            spc.geometry.getBoundaryFunc());
            auto cm2 = Geometry::massCenter(spc.particles.begin() + k, spc.particles.begin() + l,
                                            spc.geometry.getBoundaryFunc());
            return spc.geometry.vdist(cm1, cm2).cwiseProduct(dir.cast<double>()).norm();
        };
    } else if (indexes.size() == 2) {
        function = [&spc, dir = direction, i = indexes[0], j = indexes[1]]() {
            auto& cm1 = spc.groups.at(i).mass_center;
            auto& cm2 = spc.groups.at(j).mass_center;
            return spc.geometry.vdist(cm1, cm2).cwiseProduct(dir.cast<double>()).norm();
        };
    }
}
void MoleculeProperty::selectMinimumGroupDistance(const json& j, const Space& spc) {
    indexes = j.value("indexes", decltype(indexes)());
    assert(indexes.size() == 2 && "An array of 2 indexes should be specified.");
    function = [&spc, i = indexes[0], j = indexes[1]]() {
        auto minimum_distance_squared = spc.geometry.getLength().norm();
        for (const auto& particle_i : spc.findAtoms(i)) {
            for (const auto& particle_j : spc.findAtoms(j)) {
                const auto distance_squared = spc.geometry.sqdist(particle_i.pos, particle_j.pos);
                minimum_distance_squared = std::min(minimum_distance_squared, distance_squared);
            }
        }
        return sqrt(minimum_distance_squared);
    };
}
void MoleculeProperty::selectRinner(const json& j, const Space& spc) {
    direction = j.at("dir");
    indexes = j.value("indexes", decltype(indexes)());
    if (indexes.size() != 3) {
        throw ConfigurationError("An array of at least 3 indexes should be specified.");
    }
    function = [&spc, &dir = direction, i = indexes[0], j = indexes[1], k = indexes[2], l = indexes[3]]() {
        Average<double> mean_radius_j;
        Average<double> mean_radius_i;
        auto slicei = spc.findAtoms(i);
        const auto mass_center = Geometry::massCenter(slicei.begin(), slicei.end(), spc.geometry.getBoundaryFunc());
        for (const auto& particle : spc.findAtoms(j)) {
            mean_radius_j += spc.geometry.vdist(particle.pos, mass_center).cwiseProduct(dir.cast<double>()).norm();
        }
        const auto mean_radius = mean_radius_j.avg();
        for (const auto& particle : spc.activeParticles()) {
            if ((particle.id == k) or (particle.id == l)) {
                const auto radius =
                    spc.geometry.vdist(particle.pos, mass_center).cwiseProduct(dir.cast<double>()).norm();
                if (radius < mean_radius) {
                    mean_radius_i += radius;
                }
            }
        }
        return mean_radius_i.avg();
    };
}
void MoleculeProperty::selectAngleWithVector(const json& j, const Space& spc) {
    direction = j.at("dir").get<Point>().normalized();
    if (spc.groups.at(index).isMolecular()) {
        function = [&spc, &dir = direction, &group = spc.groups.at(index)]() {
            auto& cm = group.mass_center;
            auto S = Geometry::gyration(group.begin(), group.end(), cm, spc.geometry.getBoundaryFunc());
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> esf(S);
            Point eivals = esf.eigenvalues();
            ptrdiff_t i_eival;
            eivals.minCoeff(&i_eival);
            Point vec = esf.eigenvectors().col(i_eival).real();
            const auto cosine = vec.dot(dir);
            const auto angle = acos(fabs(cosine)) * 180.0 / pc::pi;
            return angle;
        };
    }
}
} // namespace Faunus::ReactionCoordinate
