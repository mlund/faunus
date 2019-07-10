#include <Eigen/Dense>
#include "space.h"
#include "reactioncoordinate.h"
#include "average.h"
#include "multipole.h"

namespace Faunus {
namespace ReactionCoordinate {

ReactionCoordinateBase::ReactionCoordinateBase(const json &j) {
    binwidth = j.value("resolution", 0.5);
    auto range = j.value("range", std::vector<double>({0, 0}));
    if (range.size() == 2) {
        if (range[0] <= range[1]) {
            min = range[0];
            max = range[1];
            return;
        }
    }
    throw std::runtime_error(name + ": 'range' require two numbers: [min, max>=min]");
}

void ReactionCoordinateBase::_to_json(json &) const {}

double ReactionCoordinateBase::normalize(double) const { return 1.; }

double ReactionCoordinateBase::operator()() {
    assert(f != nullptr);
    return f();
}

bool ReactionCoordinateBase::inRange(double coord) const { return (coord >= min && coord <= max); }

void to_json(json &j, const ReactionCoordinateBase &rc) {
    assert(!rc.name.empty());
    auto &_j = j[rc.name];
    _j = {{"range", {rc.min, rc.max}}, {"resolution", rc.binwidth}};
    rc._to_json(_j);
}

/**
 * Factory function for all known penalty functions. The json instance
 * must be an object of exact size one, for example:
 *
 *     atom: {resolution: 0.1, ... }
 */
std::shared_ptr<ReactionCoordinateBase> createReactionCoordinate(const json &j, Space &spc) {
    assert(j.is_object());
    assert(j.size() == 1);

    std::shared_ptr<ReactionCoordinateBase> rc;

    for (auto it : j.items()) {
        try {
            if (it.key() == "atom")
                rc = std::make_shared<AtomProperty>(it.value(), spc);
            else if (it.key() == "molecule")
                rc = std::make_shared<MoleculeProperty>(it.value(), spc);
            else if (it.key() == "system")
                rc = std::make_shared<SystemProperty>(it.value(), spc);
            else
                throw std::runtime_error("unknown type");
        } catch (std::exception &e) {
            throw std::runtime_error("error creating reaction coordinate '" + it.key() + "': " + e.what() +
                                     usageTip["coords=[" + it.key() + "]"]);
        }
    }
    return rc;
}

void SystemProperty::_to_json(json &j) const { j["property"] = property; }

SystemProperty::SystemProperty(const json &j, Space &spc) : ReactionCoordinateBase(j) {
    name = "system";
    property = j.at("property").get<std::string>();
    if (property == "V")
        f = [&g = spc.geo]() { return g.getVolume(); };
    else if (property == "Lx")
        f = [&g = spc.geo]() { return g.getLength().x(); };
    else if (property == "Ly")
        f = [&g = spc.geo]() { return g.getLength().y(); };
    else if (property == "Lz" or property == "height")
        f = [&g = spc.geo]() { return g.getLength().z(); };
    else if (property == "radius") {
        if (spc.geo.type == Geometry::CUBOID or spc.geo.type == Geometry::SLIT)
            std::cerr << "`radius` coordinate unavailable for geometry" << endl;
        else
            f = [&g = spc.geo]() { return 0.5 * g.getLength().x(); };
    } else if (property == "Q") // system net charge
        f = [&groups = spc.groups]() {
            double charge_sum = 0;
            for (auto &g : groups) // loops over groups
                for (auto &p : g)  // loop over particles
                    charge_sum += p.charge;
            return charge_sum;
        };
    else if (property == "N") // number of particles
        f = [&groups = spc.groups]() {
            int N_sum = 0;
            for (auto &g : groups) // loops over groups
                N_sum += g.size();
            return N_sum;
        };
    if (f == nullptr)
        throw std::runtime_error(name + ": unknown property '" + property + "'" + usageTip["coords=[system]"]);
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
        throw std::runtime_error("invalid index");
    property = j.at("property").get<std::string>();
    if (property == "x")
        f = [&p = spc.p, i = index]() { return p[i].pos.x(); };
    else if (property == "y")
        f = [&p = spc.p, i = index]() { return p[i].pos.y(); };
    else if (property == "z")
        f = [&p = spc.p, i = index]() { return p[i].pos.z(); };
    else if (property == "R")
        f = [&p = spc.p, i = index]() { return p[i].pos.norm(); };
    else if (property == "q")
        f = [&p = spc.p, i = index]() { return p[i].charge; };
    else if (property == "N") // number of atom of id=index
        f = [&groups = spc.groups, i = index]() {
            int N_sum = 0;
            for (auto &g : groups) // loops over groups
                for (auto &p : g)  // loops over particles
                    if (p.id == i)
                        N_sum++;
            return N_sum;
        };
    if (f == nullptr)
        throw std::runtime_error(name + ": unknown property '" + property + "'" + usageTip["coords=[atom]"]);
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
        throw std::runtime_error("invalid index");
    auto b = spc.geo.getBoundaryFunc();
    property = j.at("property").get<std::string>();

    if (property == "confid")
        f = [&g = spc.groups, i = index]() { return g[i].confid; };
    else if (property == "com_x")
        f = [&g = spc.groups, i = index]() { return g[i].cm.x(); };
    else if (property == "com_y")
        f = [&g = spc.groups, i = index]() { return g[i].cm.y(); };
    else if (property == "com_z")
        f = [&g = spc.groups, i = index]() { return g[i].cm.z(); };
    else if (property == "N")
        f = [&g = spc.groups, i = index]() { return g[i].size(); };
    else if (property == "Q")
        f = [&g = spc.groups, i = index]() { return monopoleMoment(g[i].begin(), g[i].end()); };

    else if (property == "mu_x")
        f = [&g = spc.groups, i = index, b]() { return dipoleMoment(g[i].begin(), g[i].end(), b).x(); };

    else if (property == "mu_y")
        f = [&g = spc.groups, i = index, b]() { return dipoleMoment(g[i].begin(), g[i].end(), b).y(); };

    else if (property == "mu_z")
        f = [&g = spc.groups, i = index, b]() { return dipoleMoment(g[i].begin(), g[i].end(), b).z(); };

    else if (property == "mu")
        f = [&g = spc.groups, i = index, b]() { return dipoleMoment(g[i].begin(), g[i].end(), b).norm(); };

    else if (property == "end2end")
        f = [&spc, i = index]() {
            assert(spc.groups[i].size() > 1);
            return std::sqrt(spc.geo.sqdist(spc.groups[i].begin()->pos, (spc.groups[i].end() - 1)->pos));
        };

    else if (property == "muangle") {
        dir = j.at("dir").get<Point>().normalized();
        if (not spc.groups.at(index).atomic) {
            f = [&g = spc.groups, i = index, b, &dir = dir]() {
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
        f = [&spc, &dir = dir, i = indexes[0], j = indexes[1]]() {
            auto &pos1 = spc.p.at(i).pos;
            auto &pos2 = spc.p.at(j).pos;
            return spc.geo.vdist(pos1, pos2).cwiseProduct(dir.cast<double>()).norm();
        };
    }

    else if (property == "cmcm_z") {
        indexes = j.value("indexes", decltype(indexes)());
        if (indexes.size() != 4)
            throw std::runtime_error("exactly four indices expected");
        f = [&spc, dir = dir, i = indexes[0], j = indexes[1] + 1, k = indexes[2], l = indexes[3] + 1]() {
            auto cm1 = Geometry::massCenter(spc.p.begin() + i, spc.p.begin() + j, spc.geo.getBoundaryFunc());
            auto cm2 = Geometry::massCenter(spc.p.begin() + k, spc.p.begin() + l, spc.geo.getBoundaryFunc());
            return spc.geo.vdist(cm1, cm2).z();
        };
    }

    else if (property == "cmcm") {
        dir = j.at("dir");
        indexes = j.value("indexes", decltype(indexes)());
        assert(indexes.size() == 4 && "An array of 4 indexes should be specified.");
        f = [&spc, dir = dir, i = indexes[0], j = indexes[1] + 1, k = indexes[2], l = indexes[3] + 1]() {
            auto cm1 = Geometry::massCenter(spc.p.begin() + i, spc.p.begin() + j, spc.geo.getBoundaryFunc());
            auto cm2 = Geometry::massCenter(spc.p.begin() + k, spc.p.begin() + l, spc.geo.getBoundaryFunc());
            return spc.geo.vdist(cm1, cm2).cwiseProduct(dir.cast<double>()).norm();
        };
    }

    else if (property == "L/R") {
        dir = j.at("dir");
        indexes = j.value("indexes", decltype(indexes)());
        assert(indexes.size() == 2 && "An array of 2 indexes should be specified.");
        f = [&spc, &dir = dir, i = indexes[0], j = indexes[1]]() {
            Average<double> Rj, Rin, Rout;
            Space::Tgroup g(spc.p.begin(), spc.p.end());
            auto slicei = g.find_id(i);
            auto cm = Geometry::massCenter(slicei.begin(), slicei.end(), spc.geo.getBoundaryFunc());
            auto slicej = g.find_id(j);
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

    else if (property == "Rhex") {
        dir = j.at("dir");
        indexes = j.value("indexes", decltype(indexes)());
        assert(indexes.size() == 2 && "An array of 2 indexes should be specified.");
        f = [&spc, &dir = dir, i = indexes[0], j = indexes[1]]() {
            Average<double> Rj, Ri;
            Space::Tgroup g(spc.p.begin(), spc.p.end());
            auto slicei = g.find_id(i);
            auto cm = Geometry::massCenter(slicei.begin(), slicei.end(), spc.geo.getBoundaryFunc());
            auto slicej = g.find_id(j);
            for (auto p : slicej)
                Rj += spc.geo.vdist(p.pos, cm).cwiseProduct(dir.cast<double>()).norm();
            double Rjavg = Rj.avg();
            for (auto p : slicei) {
                double d = spc.geo.vdist(p.pos, cm).cwiseProduct(dir.cast<double>()).norm();
                if (d < Rjavg)
                    Ri += d;
            }
            return Ri.avg();
        };
    }

    else if (property == "angle") {
        dir = j.at("dir").get<Point>().normalized();
        if (not spc.groups.at(index).atomic) {
            f = [&spc, &dir = dir, i = index]() {
                auto &cm = spc.groups[i].cm;
                auto S = Geometry::gyration(spc.groups[i].begin(), spc.groups[i].end(), spc.geo.getBoundaryFunc(), cm);
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

    if (f == nullptr)
        throw std::runtime_error(name + ": unknown or impossible property '" + property + "'" +
                                 usageTip["coords=[molecule]"]);
}
} // namespace ReactionCoordinate
} // namespace Faunus
