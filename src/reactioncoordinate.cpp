#include <doctest/doctest.h>
#include "space.h"
#include "reactioncoordinate.h"
#include "average.h"
#include "multipole.h"
#include <Eigen/Dense>
#include "aux/eigensupport.h"
#include <range/v3/view/filter.hpp>
#include <range/v3/algorithm/count_if.hpp>
#include <range/v3/algorithm/for_each.hpp>
#include <range/v3/view/join.hpp>
#include "spdlog/spdlog.h"

namespace Faunus::ReactionCoordinate {

ReactionCoordinateBase::ReactionCoordinateBase(const json& j)
{
    auto range = j.value("range", std::vector<double>({0.0, 0.0}));
    if (range.size() != 2 || range[0] > range[1]) {
        throw std::runtime_error(name + ": 'range' requires [min, max>=min]");
    }
    minimum_value = range[0];
    maximum_value = range[1];
    resolution = j.value("resolution", 0.5); // @todo require resolution user input w. at().
}

void ReactionCoordinateBase::_to_json([[maybe_unused]] json& j) const {}

double ReactionCoordinateBase::operator()()
{
    assert(function != nullptr);
    return function();
}

bool ReactionCoordinateBase::inRange(double coord) const
{
    return (coord >= minimum_value && coord <= maximum_value);
}

void to_json(json& j, const ReactionCoordinateBase& reaction_coordinate)
{
    assert(!reaction_coordinate.name.empty());
    auto& _j = j[reaction_coordinate.name];
    _j = {{"range", {reaction_coordinate.minimum_value, reaction_coordinate.maximum_value}},
          {"resolution", reaction_coordinate.resolution}};
    reaction_coordinate._to_json(_j);
}

} // namespace Faunus::ReactionCoordinate

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("[Faunus] ReactionCoordinateBase")
{
    using doctest::Approx;
    Faunus::ReactionCoordinate::ReactionCoordinateBase c(
        R"({"range":[-1.5, 2.1], "resolution":0.2})"_json);
    CHECK_EQ(c.minimum_value, Approx(-1.5));
    CHECK_EQ(c.maximum_value, Approx(2.1));
    CHECK_EQ(c.resolution, Approx(0.2));
    CHECK_EQ(c.inRange(-1.5), true);
    CHECK_EQ(c.inRange(-1.51), false);
    CHECK_EQ(c.inRange(2.11), false);
    CHECK_EQ(c.inRange(2.1), true);
}
#endif

namespace Faunus::ReactionCoordinate {

/**
 * Factory function for all known penalty functions. The json instance
 * must be an object of exact size one, for example:
 *
 *     atom: {resolution: 0.1, ... }
 */
std::unique_ptr<ReactionCoordinateBase> createReactionCoordinate(const json& j, const Space& spc)
{
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
        }
        catch (std::exception& e) {
            usageTip.pick(fmt::format("coords=[{}]", key));
            throw ConfigurationError("'{}': {}", key, e.what());
        }
    }
    catch (std::exception& e) {
        throw ConfigurationError("reaction coordinate: {}", e.what()).attachJson(j);
    }
}

void SystemProperty::_to_json(json& j) const
{
    j["property"] = property;
}

SystemProperty::SystemProperty(const json& j, const Space& spc)
    : ReactionCoordinateBase(j)
{
    namespace rv = ranges::cpp20::views;
    name = "system";
    property = j.at("property").get<std::string>();
    if (property == "V") {
        function = [&geometry = spc.geometry]() { return geometry.getVolume(); };
    }
    else if (property == "Lx") {
        function = [&geometry = spc.geometry]() { return geometry.getLength().x(); };
    }
    else if (property == "Ly") {
        function = [&geometry = spc.geometry]() { return geometry.getLength().y(); };
    }
    else if (property == "Lz" or property == "height") {
        function = [&geometry = spc.geometry]() { return geometry.getLength().z(); };
    }
    else if (property == "radius") {
        if (spc.geometry.type == Geometry::Variant::CUBOID or
            spc.geometry.type == Geometry::Variant::SLIT) {
            faunus_logger->warn("`radius` coordinate unavailable for geometry");
        }
        else {
            function = [&geometry = spc.geometry]() { return 0.5 * geometry.getLength().x(); };
        }
    }
    else if (property == "Q") { // system net charge
        function = [&spc] {
            auto charges = spc.groups | rv::join | rv::transform(&Particle::charge);
            return std::accumulate(charges.begin(), charges.end(), 0.0);
        };
    }
    else if (property == "mu") { // system dipole moment
        function = [&spc]() {
            auto particles = spc.groups | rv::join;
            return Faunus::dipoleMoment(particles.begin(), particles.end()).norm();
        };
    }
    else if (property == "mu_x") { // system dipole moment
        function = [&spc]() {
            auto particles = spc.groups | rv::join;
            return Faunus::dipoleMoment(particles.begin(), particles.end()).x();
        };
    }
    else if (property == "mu_y") { // system dipole moment
        function = [&spc]() {
            auto particles = spc.groups | rv::join;
            return Faunus::dipoleMoment(particles.begin(), particles.end()).y();
        };
    }
    else if (property == "mu_z") { // system dipole moment
        function = [&spc]() {
            auto particles = spc.groups | rv::join;
            return Faunus::dipoleMoment(particles.begin(), particles.end()).z();
        };
    }
    else if (property == "N") { // number of particles
        function = [&spc]() {
            auto sizes = spc.groups | rv::transform(&Space::GroupType::size);
            return static_cast<double>(std::accumulate(sizes.begin(), sizes.end(), size_t(0)));
        };
    }
    if (function == nullptr) {
        usageTip.pick("coords=[system]");
        throw ConfigurationError("{}: unknown property '{}'", name, property);
    }
}

void AtomProperty::_to_json(json& j) const
{
    j["property"] = property;
    j["index"] = index;
    if (dir.squaredNorm() > 1e-9) {
        j["dir"] = dir;
    }
}

/**
 * @warning For the lambda capture, always capture "Space&" and not unerlying objects like
 * `particles` or `groups`. This is because the memory location of the latter may be modified after
 * the lambda cration, thus leading to undefined dereferencing.
 */
AtomProperty::AtomProperty(const json& j, const Space& spc)
    : ReactionCoordinateBase(j)
{
    name = "atom";
    index = j.at("index");
    if (index >= spc.particles.size()) {
        throw ConfigurationError("invalid index");
    }
    property = j.at("property").get<std::string>();
    if (property == "x") {
        function = [&spc, i = index]() { return spc.particles.at(i).pos.x(); };
    }
    else if (property == "y") {
        function = [&spc, i = index]() { return spc.particles.at(i).pos.y(); };
    }
    else if (property == "z") {
        function = [&spc, i = index]() { return spc.particles.at(i).pos.z(); };
    }
    else if (property == "R") {
        function = [&spc, i = index]() { return spc.particles.at(i).pos.norm(); };
    }
    else if (property == "q") {
        function = [&spc, i = index]() { return spc.particles.at(i).charge; };
    }
    else if (property == "N") {
        function = [&spc, id = index]() {
            return static_cast<double>(
                ranges::cpp20::count_if(spc.activeParticles(), [&](const Particle& particle) {
                    return particle.id == id;
                }));
        };
    }

    if (function == nullptr) {
        usageTip.pick("coords=[atom]");
        throw ConfigurationError("{}: unknown property '{}'", name, property);
    }
}

void MoleculeProperty::_to_json(json& j) const
{
    j["property"] = property;
    j["index"] = index;
    if (direction.squaredNorm() > 1e-9) {
        j["dir"] = direction;
    }
    if (indexes.size() >= 2) {
        j["indexes"] = indexes;
    }
}

/**
 * @warning For the lambda capture, always capture "Space&" and not unerlying objects like
 * `particle` or `group`. This is because the memory location of the latter may be modified after
 * the lambda cration, thus leading to undefined dereferencing.
 */
MoleculeProperty::MoleculeProperty(const json& j, const Space& spc)
    : ReactionCoordinateBase(j)
{
    name = "molecule";
    index = j.value("index", 0);
    if (index >= spc.groups.size()) {
        throw ConfigurationError("invalid index");
    }
    auto b = spc.geometry.getBoundaryFunc();

    property = j.at("property").get<std::string>();

    if (property == "active") { // if molecule is active (1) or not (0)
        function = [&spc, i = index]() { return static_cast<double>(!spc.groups.at(i).empty()); };
    }
    else if (property == "confid") {
        function = [&spc, i = index]() {
            return static_cast<double>(spc.groups.at(i).conformation_id);
        };
    }
    else if (property == "com_x") {
        function = [&spc, i = index]() { return spc.groups.at(i).mass_center.x(); };
    }
    else if (property == "com_y") {
        function = [&spc, i = index]() { return spc.groups.at(i).mass_center.y(); };
    }
    else if (property == "com_z") {
        function = [&spc, i = index]() { return spc.groups.at(i).mass_center.z(); };
    }
    else if (property == "N") {
        function = [&spc, i = index]() { return static_cast<double>(spc.groups.at(i).size()); };
    }
    else if (property == "Q") {
        function = [&spc, i = index]() {
            return monopoleMoment(spc.groups.at(i).begin(), spc.groups.at(i).end());
        };
    }
    else if (property == "mu_x") {
        function = [&spc, i = index, b]() {
            return dipoleMoment(spc.groups.at(i).begin(), spc.groups.at(i).end(), b).x();
        };
    }
    else if (property == "mu_y") {
        function = [&spc, i = index, b]() {
            return dipoleMoment(spc.groups.at(i).begin(), spc.groups.at(i).end(), b).y();
        };
    }
    else if (property == "mu_z") {
        function = [&spc, i = index, b]() {
            return dipoleMoment(spc.groups.at(i).begin(), spc.groups.at(i).end(), b).z();
        };
    }
    else if (property == "mu") {
        function = [&spc, i = index, b]() {
            return dipoleMoment(spc.groups.at(i).begin(), spc.groups.at(i).end(), b).norm();
        };
    }
    else if (property == "end2end") {
        function = [&spc, i = index]() {
            return std::sqrt(spc.geometry.sqdist(spc.groups.at(i).begin()->pos,
                                                 std::prev(spc.groups.at(i).end())->pos));
        };
    }
    else if (property == "Rg") {
        selectGyrationRadius(spc);
    }
    else if (property == "muangle") {
        selectDipoleAngle(j, spc, b);
    }
    else if (property == "atomatom") {
        selectAtomAtomDistance(j, spc);
    }
    else if (property == "cmcm_z") {
        selectMassCenterDistanceZ(j, spc);
    }
    else if (property == "cmcm") {
        selectMassCenterDistance(j, spc);
    }
    else if (property == "L/R") {
        selectLengthOverRadiusRatio(j, spc);
    }
    else if (property == "mindist") {
        selectMinimumGroupDistance(j, spc);
    }
    else if (property == "Rinner") {
        selectRinner(j, spc);
    }
    else if (property == "angle") {
        selectAngleWithVector(j, spc);
    }
    if (function == nullptr) {
        usageTip.pick("coords=[molecule]");
        throw ConfigurationError("{}: unknown or impossible property property '{}'", name,
                                 property);
    }
}

void MoleculeProperty::selectLengthOverRadiusRatio(const json& j, const Space& spc)
{
    direction = j.at("dir");
    indexes = j.value("indexes", decltype(indexes)());
    if (indexes.size() != 2) {
        throw ConfigurationError("An array of 2 indexes should be specified");
    }
    function = [&spc, &dir = direction, i = indexes[0], j = indexes[1]]() {
        Average<double> mean_radius_j;
        Average<double> Rin;
        Average<double> Rout;
        auto particles_i = spc.findAtoms(i);
        Point mass_center_i = Geometry::massCenter(particles_i.begin(), particles_i.end(),
                                                   spc.geometry.getBoundaryFunc());
        for (const Particle& particle : spc.findAtoms(j)) {
            mean_radius_j += spc.geometry.vdist(particle.pos, mass_center_i)
                                 .cwiseProduct(dir.cast<double>())
                                 .norm();
        }
        const auto Rjavg = mean_radius_j.avg();
        for (const auto& particle_i : particles_i) {
            const auto radial_distance = spc.geometry.vdist(particle_i.pos, mass_center_i)
                                             .cwiseProduct(dir.cast<double>())
                                             .norm();
            if (radial_distance < Rjavg) {
                Rin += radial_distance;
            }
            else if (radial_distance > Rjavg) {
                Rout += radial_distance;
            }
        }
        return 2.0 * spc.geometry.getLength().z() / (Rin.avg() + Rout.avg());
    };
}

void MoleculeProperty::selectMassCenterDistanceZ(const json& j, const Space& spc)
{
    indexes = j.value("indexes", decltype(indexes)());
    if (indexes.size() == 4) {
        function = [&spc, i = indexes[0], j = indexes[1] + 1, k = indexes[2],
                    l = indexes[3] + 1]() {
            Point cm1 = Geometry::massCenter(spc.particles.begin() + i, spc.particles.begin() + j,
                                             spc.geometry.getBoundaryFunc());
            Point cm2 = Geometry::massCenter(spc.particles.begin() + k, spc.particles.begin() + l,
                                             spc.geometry.getBoundaryFunc());
            return spc.geometry.vdist(cm1, cm2).z();
        };
    }
    else if (indexes.size() == 2) {
        function = [&spc, i = indexes[0], j = indexes[1]]() {
            return spc.geometry.vdist(spc.groups.at(i).mass_center, spc.groups.at(j).mass_center)
                .z();
        };
    }
    throw ConfigurationError("An array of 2 or 4 indexes should be specified.");
}

void MoleculeProperty::selectAtomAtomDistance(const json& j, const Space& spc)
{
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

void MoleculeProperty::selectGyrationRadius(const Space& spc)
{
    function = [&spc, i = index]() {
        assert(spc.groups.at(i).size() > 1);
        Tensor S = Geometry::gyration(spc.groups.at(i).begin(), spc.groups.at(i).end(),
                                      spc.groups.at(i).mass_center, spc.geometry.getBoundaryFunc());
        return sqrt(S.trace()); // S.trace() == S.eigenvalues().sum() but faster
    };
}

void MoleculeProperty::selectDipoleAngle(const json& j, const Space& spc,
                                         Geometry::BoundaryFunction& b)
{
    direction = j.at("dir").get<Point>().normalized();
    if (spc.groups.at(index).isMolecular()) {
        function = [&spc, i = index, b, &dir = direction]() {
            const auto dot_product =
                dipoleMoment(spc.groups.at(i).begin(), spc.groups.at(i).end(), b).dot(dir);
            return acos(dot_product) * 180.0 / pc::pi;
        };
    }
}

void MoleculeProperty::selectMassCenterDistance(const json& j, const Space& spc)
{
    direction = j.at("dir");
    indexes = j.value("indexes", decltype(indexes)());
    assert(indexes.size() > 1 && "An array of 2 or 4 indexes should be specified.");
    if (indexes.size() == 4) {
        function = [&spc, dir = direction.cast<double>(), i = indexes[0], j = indexes[1] + 1,
                    k = indexes[2], l = indexes[3] + 1]() {
            Point cm1 = Geometry::massCenter(spc.particles.begin() + i, spc.particles.begin() + j,
                                             spc.geometry.getBoundaryFunc());
            Point cm2 = Geometry::massCenter(spc.particles.begin() + k, spc.particles.begin() + l,
                                             spc.geometry.getBoundaryFunc());
            return spc.geometry.vdist(cm1, cm2).cwiseProduct(dir).norm();
        };
    }
    else if (indexes.size() == 2) {
        function = [&spc, dir = direction.cast<double>(), i = indexes[0], j = indexes[1]]() {
            const auto& cm1 = spc.groups.at(i).mass_center;
            const auto& cm2 = spc.groups.at(j).mass_center;
            return spc.geometry.vdist(cm1, cm2).cwiseProduct(dir).norm();
        };
    }
}

void MoleculeProperty::selectMinimumGroupDistance(const json& j, const Space& spc)
{
    indexes = j.value("indexes", decltype(indexes)());
    if (indexes.size() != 2) {
        throw ConfigurationError("indexes must have two elements");
    }
    function = [&spc, i = indexes.at(0), j = indexes.at(1)]() {
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

void MoleculeProperty::selectRinner(const json& j, const Space& spc)
{
    direction = j.at("dir");
    indexes = j.value("indexes", decltype(indexes)());
    if (indexes.size() != 4) {
        throw ConfigurationError("indexes must have four elements");
    }
    function = [&spc, &dir = direction, i = indexes.at(0), j = indexes.at(1), k = indexes.at(2),
                l = indexes.at(3)]() {
        namespace rv = ranges::cpp20::views;
        auto slicei = spc.findAtoms(i);
        const auto mass_center =
            Geometry::massCenter(slicei.begin(), slicei.end(), spc.geometry.getBoundaryFunc());

        // filter and transform lambdas
        auto k_or_l = [k, l](auto& particle) { return (particle.id == k) or (particle.id == l); };
        auto radius = [&, dir = dir.cast<double>()](auto& particle) {
            return spc.geometry.vdist(particle.pos, mass_center).cwiseProduct(dir).norm();
        };
        auto mean = [](auto& radii) {
            Average<double> mean;
            ranges::cpp20::for_each(radii, [&](auto radius) { mean += radius; });
            return mean.avg();
        };

        auto radii_j = spc.findAtoms(j) | rv::transform(radius);
        auto mean_radii_j = mean(radii_j);

        auto radii = spc.activeParticles() | rv::filter(k_or_l) | rv::transform(radius) |
                     rv::filter([mean_radii_j](auto radius) { return radius < mean_radii_j; });
        return mean(radii);
    };
}

void MoleculeProperty::selectAngleWithVector(const json& j, const Space& spc)
{
    direction = j.at("dir").get<Point>().normalized();
    if (spc.groups.at(index).isMolecular()) {
        function = [&spc, &dir = direction, i = index]() {
            auto& group = spc.groups.at(i);
            auto& cm = group.mass_center;
            Tensor S =
                Geometry::gyration(group.begin(), group.end(), cm, spc.geometry.getBoundaryFunc());
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> esf(S);
            const Point& eivals = esf.eigenvalues();
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
