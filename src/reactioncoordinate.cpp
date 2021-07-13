#include <doctest/doctest.h>
#include "space.h"
#include "reactioncoordinate.h"
#include "average.h"
#include "multipole.h"
#include <Eigen/Dense>
#include "aux/eigensupport.h"
#include "spdlog/spdlog.h"

namespace Faunus::ReactionCoordinate {

ReactionCoordinateBase::ReactionCoordinateBase(const json &j) {
    resolution = j.value("resolution", 0.5);
    auto range = j.value("range", std::vector<double>({0, 0}));
    if (range.size() == 2) {
        if (range[0] <= range[1]) {
            minimum_value = range[0];
            maximum_value = range[1];
            return;
        }
    }
    throw std::runtime_error(name + ": 'range' require two numbers: [min, max>=min]");
}

void ReactionCoordinateBase::_to_json(json &) const {}

double ReactionCoordinateBase::normalize(double) const { return 1.; }

double ReactionCoordinateBase::operator()() {
    assert(function != nullptr);
    return function();
}

bool ReactionCoordinateBase::inRange(double coord) const { return (coord >= minimum_value && coord <= maximum_value); }

void to_json(json &j, const ReactionCoordinateBase &rc) {
    assert(!rc.name.empty());
    auto &_j = j[rc.name];
    _j = {{"range", {rc.minimum_value, rc.maximum_value}}, {"resolution", rc.resolution}};
    rc._to_json(_j);
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
std::shared_ptr<ReactionCoordinateBase> createReactionCoordinate(const json& j, Space& spc) {
    std::shared_ptr<ReactionCoordinateBase> reaction_coordinate;
    try {
        const auto& [key, j_params] = jsonSingleItem(j);
        try {
            if (key == "atom") {
                reaction_coordinate = std::make_shared<AtomProperty>(j_params, spc);
            } else if (key == "molecule") {
                reaction_coordinate = std::make_shared<MoleculeProperty>(j_params, spc);
            } else if (key == "system") {
                reaction_coordinate = std::make_shared<SystemProperty>(j_params, spc);
            } else {
                throw ConfigurationError("unknown reaction coordinate");
            }
        } catch (std::exception& e) {
            usageTip.pick(fmt::format("coords=[{}]", key));
            throw ConfigurationError("'{}': {}", key, e.what());
        }
    } catch (std::exception& e) {
        throw ConfigurationError("reaction coordinate: {}", e.what()).attachJson(j);
    }
    return reaction_coordinate;
}

void SystemProperty::_to_json(json &j) const { j["property"] = property; }

SystemProperty::SystemProperty(const json &j, Space &spc) : ReactionCoordinateBase(j) {
    name = "system";
    property = j.at("property").get<std::string>();
    if (property == "V")
        function = [&g = spc.geo]() { return g.getVolume(); };
    else if (property == "Lx")
        function = [&g = spc.geo]() { return g.getLength().x(); };
    else if (property == "Ly")
        function = [&g = spc.geo]() { return g.getLength().y(); };
    else if (property == "Lz" or property == "height")
        function = [&g = spc.geo]() { return g.getLength().z(); };
    else if (property == "radius") {
        if (spc.geo.type == Geometry::CUBOID or spc.geo.type == Geometry::SLIT)
            faunus_logger->warn("`radius` coordinate unavailable for geometry");
        else
            function = [&g = spc.geo]() { return 0.5 * g.getLength().x(); };
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
            size_t number_of_active_particles = 0;
            for (const auto& group : groups) { // loops over groups
                number_of_active_particles += group.size();
            }
            return number_of_active_particles;
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
    if (dir.squaredNorm() > 1e-9)
        j["dir"] = dir;
}

AtomProperty::AtomProperty(const json &j, Space &spc) : ReactionCoordinateBase(j) {
    name = "atom";
    index = j.at("index");
    if (index >= spc.p.size())
        throw ConfigurationError("invalid index");
    property = j.at("property").get<std::string>();
    if (property == "x")
        function = [&p = spc.p, i = index]() { return p[i].pos.x(); };
    else if (property == "y")
        function = [&p = spc.p, i = index]() { return p[i].pos.y(); };
    else if (property == "z")
        function = [&p = spc.p, i = index]() { return p[i].pos.z(); };
    else if (property == "R")
        function = [&p = spc.p, i = index]() { return p[i].pos.norm(); };
    else if (property == "q")
        function = [&p = spc.p, i = index]() { return p[i].charge; };
    else if (property == "N") // number of atom of id=index
        function = [&groups = spc.groups, i = index]() {
            int N_sum = 0;
            for (auto &g : groups) // loops over groups
                for (auto &p : g)  // loops over particles
                    if (p.id == i)
                        N_sum++;
            return N_sum;
        };

    if (function == nullptr) {
        usageTip.pick("coords=[atom]");
        throw ConfigurationError("{}: unknown property '{}'", name, property);
    }
}

void MoleculeProperty::_to_json(json &j) const {
    j["property"] = property;
    j["index"] = index;
    if (dir.squaredNorm() > 1e-9)
        j["dir"] = dir;
    if (indexes.size() >= 2)
        j["indexes"] = indexes;
}

MoleculeProperty::MoleculeProperty(const json &j, Space &spc) : ReactionCoordinateBase(j) {
    name = "molecule";
    index = j.value("index", 0);
    if (index >= spc.groups.size())
        throw ConfigurationError("invalid index");
    auto b = spc.geo.getBoundaryFunc();
    property = j.at("property").get<std::string>();

    if (property == "active") // if molecule is active (1) or not (0)
        function = [&g = spc.groups, i = index]() { return (double)(!g[i].empty()); };
    if (property == "confid")
        function = [&g = spc.groups, i = index]() { return g[i].confid; };
    else if (property == "com_x")
        function = [&g = spc.groups, i = index]() { return g[i].cm.x(); };
    else if (property == "com_y")
        function = [&g = spc.groups, i = index]() { return g[i].cm.y(); };
    else if (property == "com_z")
        function = [&g = spc.groups, i = index]() { return g[i].cm.z(); };
    else if (property == "N")
        function = [&g = spc.groups, i = index]() { return g[i].size(); };
    else if (property == "Q")
        function = [&g = spc.groups, i = index]() { return monopoleMoment(g[i].begin(), g[i].end()); };

    else if (property == "mu_x")
        function = [&g = spc.groups, i = index, b]() { return dipoleMoment(g[i].begin(), g[i].end(), b).x(); };

    else if (property == "mu_y")
        function = [&g = spc.groups, i = index, b]() { return dipoleMoment(g[i].begin(), g[i].end(), b).y(); };

    else if (property == "mu_z")
        function = [&g = spc.groups, i = index, b]() { return dipoleMoment(g[i].begin(), g[i].end(), b).z(); };

    else if (property == "mu")
        function = [&g = spc.groups, i = index, b]() { return dipoleMoment(g[i].begin(), g[i].end(), b).norm(); };

    else if (property == "end2end")
        function = [&spc, i = index]() {
            assert(spc.groups[i].size() > 1);
            return std::sqrt(spc.geo.sqdist(spc.groups[i].begin()->pos, (spc.groups[i].end() - 1)->pos));
        };

    else if (property == "Rg")
        function = [&spc, i = index]() {
            assert(spc.groups[i].size() > 1);
            auto S = Geometry::gyration(spc.groups[i].begin(), spc.groups[i].end(), spc.groups[i].cm, spc.geo.getBoundaryFunc());
            return std::sqrt(S.trace()); // S.trace() == S.eigenvalues().sum() but faster
        };

    else if (property == "muangle") {
        dir = j.at("dir").get<Point>().normalized();
        if (not spc.groups.at(index).atomic) {
            function = [&g = spc.groups, i = index, b, &dir = dir]() {
                Point mu = dipoleMoment(g[i].begin(), g[i].end(), b);
                return std::acos(mu.dot(dir)) * 180 / pc::pi;
            };
        }
    }

    else if (property == "atomatom") {
        dir = j.at("dir");
        indexes = j.at("indexes").get<decltype(indexes)>();
        if (indexes.size() != 2)
            throw std::runtime_error("exactly two indices expected");
        function = [&spc, &dir = dir, i = indexes[0], j = indexes[1]]() {
            auto &pos1 = spc.p.at(i).pos;
            auto &pos2 = spc.p.at(j).pos;
            return spc.geo.vdist(pos1, pos2).cwiseProduct(dir.cast<double>()).norm();
        };
    }

    else if (property == "cmcm_z") {
        indexes = j.value("indexes", decltype(indexes)());
        assert(indexes.size() > 1 && "An array of 2 or 4 indexes should be specified.");
        if (indexes.size() == 4) {
            function = [&spc, dir = dir, i = indexes[0], j = indexes[1] + 1, k = indexes[2], l = indexes[3] + 1]() {
                auto cm1 = Geometry::massCenter(spc.p.begin() + i, spc.p.begin() + j, spc.geo.getBoundaryFunc());
                auto cm2 = Geometry::massCenter(spc.p.begin() + k, spc.p.begin() + l, spc.geo.getBoundaryFunc());
                return spc.geo.vdist(cm1, cm2).z();
            };
        }
        if (indexes.size() == 2) {
            function = [&spc, dir = dir, i = indexes[0], j = indexes[1]]() {
                auto cm1 = spc.groups[i].cm;
                auto cm2 = spc.groups[j].cm;
                return spc.geo.vdist(cm1, cm2).z();
            };
        }
    }

    else if (property == "cmcm") {
        dir = j.at("dir");
        indexes = j.value("indexes", decltype(indexes)());
        assert(indexes.size() > 1 && "An array of 2 or 4 indexes should be specified.");
        if (indexes.size() == 4) {
            function = [&spc, dir = dir, i = indexes[0], j = indexes[1] + 1, k = indexes[2], l = indexes[3] + 1]() {
                auto cm1 = Geometry::massCenter(spc.p.begin() + i, spc.p.begin() + j, spc.geo.getBoundaryFunc());
                auto cm2 = Geometry::massCenter(spc.p.begin() + k, spc.p.begin() + l, spc.geo.getBoundaryFunc());
                return spc.geo.vdist(cm1, cm2).cwiseProduct(dir.cast<double>()).norm();
            };
        }
        if (indexes.size() == 2) {
            function = [&spc, dir = dir, i = indexes[0], j = indexes[1]]() {
                auto cm1 = spc.groups[i].cm;
                auto cm2 = spc.groups[j].cm;
                return spc.geo.vdist(cm1, cm2).cwiseProduct(dir.cast<double>()).norm();
            };
        }
    }

    else if (property == "L/R") {
        dir = j.at("dir");
        indexes = j.value("indexes", decltype(indexes)());
        assert(indexes.size() == 2 && "An array of 2 indexes should be specified.");
        function = [&spc, &dir = dir, i = indexes[0], j = indexes[1]]() {
            Average<double> Rj, Rin, Rout;
            auto slicei = spc.findAtoms(i);
            auto cm = Geometry::massCenter(slicei.begin(), slicei.end(), spc.geo.getBoundaryFunc());
            auto slicej = spc.findAtoms(j);
            for (auto p : slicej)
                Rj += spc.geo.vdist(p.pos, cm).cwiseProduct(dir.cast<double>()).norm();
            double Rjavg = Rj.avg();
            for (auto p : slicei) {
                double d = spc.geo.vdist(p.pos, cm).cwiseProduct(dir.cast<double>()).norm();
                if (d < Rjavg)
                    Rin += d;
                else if (d > Rjavg)
                    Rout += d;
            }
            return 2 * spc.geo.getLength().z() / (Rin.avg() + Rout.avg());
        };
    }

    else if (property == "mindist") {
        indexes = j.value("indexes", decltype(indexes)());
        assert(indexes.size() == 2 && "An array of 2 indexes should be specified.");
        function = [&spc, i = indexes[0], j = indexes[1]]() {
            auto slicei = spc.findAtoms(i);
            auto slicej = spc.findAtoms(j);
            double dmin = spc.geo.getLength().norm();
            for (auto pi : slicei) {
                for (auto pj : slicej) {
                    double d = spc.geo.sqdist(pi.pos, pj.pos);
                    if (d < dmin)
                        dmin = d;
                }
            }
            return std::sqrt(dmin);
        };
    }

    else if (property == "Rinner") {
        dir = j.at("dir");
        indexes = j.value("indexes", decltype(indexes)());
        if (indexes.size() != 3)
            throw ConfigurationError("An array of at least 3 indexes should be specified.");
        function = [&spc, &dir = dir, i = indexes[0], j = indexes[1], k = indexes[2], l = indexes[3]]() {
            Average<double> Rj, Ri;
            auto slicei = spc.findAtoms(i);
            const auto mass_center = Geometry::massCenter(slicei.begin(), slicei.end(), spc.geo.getBoundaryFunc());
            auto slicej = spc.findAtoms(j);
            for (auto p : slicej)
                Rj += spc.geo.vdist(p.pos, mass_center).cwiseProduct(dir.cast<double>()).norm();
            const auto mean_radius = Rj.avg();
            for (const auto &particle : spc.activeParticles()) {
                if ((particle.id==k) or (particle.id==l)) {
                    const auto distance  = spc.geo.vdist(particle.pos, mass_center).cwiseProduct(dir.cast<double>()).norm();
                    if (distance < mean_radius)
                        Ri += distance;
                }
            }
            return Ri.avg();
        };
    }

    else if (property == "angle") {
        dir = j.at("dir").get<Point>().normalized();
        if (not spc.groups.at(index).atomic) {
            function = [&spc, &dir = dir, i = index]() {
                auto &cm = spc.groups[i].cm;
                auto S = Geometry::gyration(spc.groups[i].begin(), spc.groups[i].end(), cm, spc.geo.getBoundaryFunc());
                Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> esf(S);
                Point eivals = esf.eigenvalues();
                std::ptrdiff_t i_eival;
                eivals.minCoeff(&i_eival);
                Point vec = esf.eigenvectors().col(i_eival).real();
                double cosine = vec.dot(dir);
                double angle = std::acos(std::fabs(cosine)) * 180. / pc::pi;
                return angle;
            };
        }
    }
    if (function == nullptr) {
        usageTip.pick("coords=[molecule]");
        throw ConfigurationError("{}: unknown or impossible property property '{}'", name, property);
    }
}
} // namespace Faunus::ReactionCoordinate
