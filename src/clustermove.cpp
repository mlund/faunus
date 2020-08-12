#include "clustermove.h"
#include "aux/eigensupport.h"
#include <algorithm>
#include <range/v3/view/cartesian_product.hpp>

namespace Faunus {
namespace Move {

double Cluster::clusterProbability(const Cluster::Tgroup &group1, const Cluster::Tgroup &group2) const {
    if (!group1.empty() and !group2.empty()) {
        if (spc.geo.sqdist(group1.cm, group2.cm) <= thresholds_squared(group1.id, group2.id)) {
            return 1.0;
        }
    }
    return 0.0;
}
void Cluster::_to_json(json &j) const {
    using namespace u8;
    j = {{"dir", translation_direction},
         {"dp", translation_displacement_factor},
         {"dprot", rotation_displacement_factor},
         {"single_layer", single_layer},
         {"dirrot", rotation_axis},
         {rootof + bracket("r" + squared), std::sqrt(translational_mean_square_displacement.avg())},
         {rootof + bracket(theta + squared) + "/" + degrees,
          std::sqrt(rotational_mean_square_displacement.avg()) / 1.0_deg},
         {bracket("N"), average_cluster_size.avg()},
         {"bias rejection rate", double(bias_rejected) / cnt},
         {"clusterdistribution", cluster_size_distribution}};
    _roundjson(j, 3);

    // print threshold matrix
    auto &_j = j["threshold"];
    for (auto i : molids)
        for (auto j : molids)
            if (i >= j) {
                std::string key = Faunus::molecules[i].name + " " + Faunus::molecules[j].name;
                _j[key] = std::sqrt(thresholds_squared(i, j));
                _roundjson(_j[key], 3);
            }

    // print satellite molecules
    if (not satellites.empty()) {
        auto &_j = j["satellites"];
        _j = json::array();
        for (auto id : satellites)
            _j.push_back(Faunus::molecules[id].name);
    }
}

void Cluster::registerSatellites(const std::vector<std::string> &satellite_names) {
    const auto ids = Faunus::names2ids(Faunus::molecules, satellite_names); // names --> molids
    satellites = std::set<int>(ids.begin(), ids.end());

    auto invalid_satellite = [&](auto molid) { return std::find(molids.begin(), molids.end(), molid) == molids.end(); };
    if (std::any_of(satellites.begin(), satellites.end(), invalid_satellite)) {
        throw std::runtime_error("satellite molecules must be defined in `molecules`");
    }

    if (satellites.size() >= molids.size()) {
        throw std::runtime_error("all molecules cannot be satellites");
    }
}

/**
 * @param j json object corresponding to "threshold" key in user input
 * @throws If invalid input
 *
 * This reads the cluster threshold(s) either as a single number, or
 * as a matrix of molecule-molecule thresholds.
 */
void Cluster::parseThresholds(const json &j) {
    if (j.is_number()) { // threshold is given as a single number
        for (auto [id1, id2] : ranges::views::cartesian_product(molids, molids)) {
            thresholds_squared.set(id1, id2, std::pow(j.get<double>(), 2));
        }
    } else if (j.is_object()) { // threshold is given as pairs of clustering molecules
        const int threshold_combinations = molids.size() * (molids.size() + 1) / 2; // N*(N+1)/2
        if (j.size() != threshold_combinations) {
            throw ConfigurationError(
                "exactly {} molecule pairs must be given in threshold matrix to cover all combinations",
                threshold_combinations);
        }
        for (const auto &[key, value] : j.items()) {
            if (auto name_pair = Faunus::words2vec<std::string>(key); name_pair.size() == 2) {
                auto it1 = Faunus::findName(Faunus::molecules, name_pair[0]);
                auto it2 = Faunus::findName(Faunus::molecules, name_pair[1]);
                if (it1 == Faunus::molecules.end() or it2 == Faunus::molecules.end()) {
                    throw ConfigurationError("unknown threshold molecule(s): [{} {}]", name_pair[0], name_pair[1]);
                }
                thresholds_squared.set(it1->id(), it2->id(), std::pow(value.get<double>(), 2));
            } else {
                throw ConfigurationError("threshold requires exactly two space-separated molecules");
            }
        }
    } else {
        throw ConfigurationError("cluster threshold must be a number or object");
    }
}

void Cluster::_from_json(const json &j) {
    translation_displacement_factor = j.at("dp");
    translation_direction = j.value("dir", Point(1, 1, 1));
    rotation_axis = j.value("dirrot", Point(0, 0, 0)); // predefined axis of rotation
    rotation_axis.normalize();                         // make sure dirrot is a unit-vector
    rotation_displacement_factor = j.at("dprot");
    single_layer = j.value("single_layer", false);
    if (j.count("spread")) {
        faunus_logger->warn("{}: 'spread' is deprecated, use 'single_layer' instead", name);
    }
    molecule_names = j.at("molecules").get<decltype(molecule_names)>(); // molecule names
    molids = Faunus::names2ids(Faunus::molecules, molecule_names);      // names --> molids
    std::sort(molids.begin(), molids.end());
    updateMoleculeIndex();

    repeat = std::max<size_t>(1, molecule_index.size());

    // read satellite ids (molecules NOT to be considered as cluster centers)
    registerSatellites(j.value("satellites", std::vector<std::string>()));
    parseThresholds(j.at("threshold"));
}

/**
 * @return Index of random seed molecule
 * @throws if no seed can be found (no molecules)
 *
 * Satellites are molecules that cannot be used to seed a cluster but can be
 * part of it, i.e. rotated and translated. Here we therefore filter out any
 * satellite molecules, and pick a random molecule index of what remains.
 */
std::size_t Cluster::findSeed() {
    auto avoid_satellites = [&](auto index) { return (satellites.count(spc.groups[index].id) == 0); };
    auto not_satellites = molecule_index | ranges::cpp20::views::filter(avoid_satellites);
    if (ranges::cpp20::empty(not_satellites)) {
        throw std::runtime_error("cannot find cluster seed: out of molecules");
    }
    return *slump.sample(not_satellites.begin(), not_satellites.end());
}

/**
 * Find cluster
 *
 * @param spc Space to operate on
 * @param seed_index Index of seed_index group to evaluate the cluster around
 * @returns Destination vector for group indices of the found cluster
 */
std::vector<size_t> Cluster::findCluster(Space &spc, size_t seed_index) {
    assert(seed_index < spc.p.size());
    std::set<size_t> pool(molecule_index.begin(), molecule_index.end());
    assert(pool.count(seed_index) == 1);

    std::vector<size_t> cluster;
    cluster.reserve(molecule_index.size()); // ensures safe resizing without invalidating iterators
    cluster.push_back(seed_index);          // 'seed_index' is the index of the seed molecule
    pool.erase(seed_index);                 // ...which is already in the cluster and not part of pool

    // cluster search algorithm
    for (auto it1 = cluster.begin(); it1 != cluster.end(); it1++) {
        for (auto it2 = pool.begin(); it2 != pool.end();) {
            double P = clusterProbability(spc.groups.at(*it1), spc.groups.at(*it2)); // probability to cluster
            if (Movebase::slump() <= P) {
                cluster.push_back(*it2); // add to cluster
                it2 = pool.erase(it2);   // erase and advance (c++11)
            } else {
                ++it2;
            }
        }
        if (single_layer) { // stop after one iteration around 'seed_index'
            break;
        }
    }
    std::sort(cluster.begin(), cluster.end()); // required for correct energy evaluation

    // check if cluster is too large
    perform_rotation = true;
    double max = spc.geo.getLength().minCoeff() / 2;
    for (auto i : cluster) {
        for (auto j : cluster) {
            if (j > i) {
                if (spc.geo.sqdist(spc.groups.at(i).cm, spc.groups.at(j).cm) >= max * max) {
                    perform_rotation = false; // skip rotation if cluster larger than half the box length
                }
            }
        }
    }
    assert(std::adjacent_find(cluster.begin(), cluster.end()) == cluster.end()); // check for duplicates
    return cluster;
}

Point Cluster::clusterMassCenter(const std::vector<size_t> &cluster_index) const {
    auto boundary = spc.geo.getBoundaryFunc();
    double mass_sum = 0.0;
    Point mass_center(0, 0, 0);
    Point origin = spc.groups[*cluster_index.begin()].cm;
    for (auto i : cluster_index) { // loop over clustered molecules (index)
        Point r = spc.groups[i].cm - origin;
        boundary(r);
        double mass = spc.groups[i].mass();
        mass_center += mass * r;
        mass_sum += mass;
    }
    mass_center = mass_center / mass_sum + origin;
    boundary(mass_center);
    return mass_center;
}

Eigen::Quaterniond Cluster::setRandomRotation() {
    Eigen::Quaterniond quaternion;
    if (perform_rotation) {
        Point axis = (rotation_axis.count() > 0) ? rotation_axis : ranunit(slump);
        rotation_angle = rotation_displacement_factor * (slump() - 0.5);
        quaternion = Eigen::AngleAxisd(rotation_angle, axis); // quaternion
    } else {
        rotation_angle = 0.0;
    }
    return quaternion;
}

void Cluster::_move(Change &change) {
    using namespace ranges::cpp20;
    auto boundary = spc.geo.getBoundaryFunc();
    _bias = 0;
    perform_rotation = true;
    updateMoleculeIndex();
    if (not molecule_index.empty()) {
        auto seed_index = findSeed();                // find "nuclei" or cluster center
        auto cluster = findCluster(spc, seed_index); // find cluster
        Point COM = clusterMassCenter(cluster);      // cluster mass center
        auto quaternion = setRandomRotation();       // set rotation

        auto rotate_group = [&](auto &group) {
            Geometry::rotate(group.begin(), group.end(), quaternion, boundary, -COM);
            group.cm = group.cm - COM;
            boundary(group.cm);
            group.cm = quaternion * group.cm + COM;
            boundary(group.cm);
        };

        // translate cluster
        translation_displacement = ranunit(slump, translation_direction) * translation_displacement_factor * slump();
        auto cluster_groups =
            cluster | views::transform([&](size_t index) -> Space::Tgroup & { return spc.groups.at(index); });
        for (auto &group : cluster_groups) {
            if (perform_rotation) {
                rotate_group(group);
            }
            group.translate(translation_displacement, boundary);
            Change::data d;
            d.all = true;
            d.index = std::distance(&spc.groups.front(), &group);
            change.groups.push_back(d);
        }
        change.moved2moved = false; // do not calc. internal cluster energy

        /*
         * Reject if cluster composition changes during move
         * Note: this only works for the binary 0/1 probability function currently implemented in `findCluster()`.
         */
        if (auto aftercluster = findCluster(spc, seed_index); aftercluster == cluster) {
            _bias = 0;
        } else {
            _bias = pc::infty; // bias is infinite --> reject
            bias_rejected++;   // count how many time we reject due to bias
        }

        average_cluster_size += cluster.size();      // average cluster size
        cluster_size_distribution[cluster.size()]++; // update cluster size distribution

        if constexpr (false) { // debug code: check if cluster mass center movement matches displacement
            if (_bias == 0) {
                auto newCOM = clusterMassCenter(cluster);             // org. cluster center
                Point d = spc.geo.vdist(COM, newCOM);                 // distance between new and old COM
                double _zero = (d + translation_displacement).norm(); // |d+dp| should ideally be zero...
                if (std::fabs(_zero) > 1e-9) {
                    _bias = pc::infty; // by setting bias=oo the move is rejected
                    faunus_logger->warn("Skipping too large cluster: COM difference = {}", _zero);
                }
            }
        }
    }
}
double Cluster::bias(Change &, double, double) { return _bias; }
void Cluster::_reject(Change &) {
    translational_mean_square_displacement += 0.0;
    rotational_mean_square_displacement += 0.0;
}
void Cluster::_accept(Change &) {
    translational_mean_square_displacement += translation_displacement.squaredNorm();
    rotational_mean_square_displacement += rotation_angle * rotation_angle;
}
Cluster::Cluster(Space &spc) : spc(spc) {
    cite = "doi:10/cj9gnn";
    name = "cluster";
    repeat = -1;
}

/**
 * Search for molecules participating in the cluster move.
 * This should be called for every move event as a
 * grand canonical move may have changed the number
 * of particles.
 *
 * Atomic groups and inactive groups are ignored.
 */
void Cluster::updateMoleculeIndex() {
    assert(not molids.empty());
    molecule_index.clear();
    for (auto &group : spc.groups) {
        if (!group.atomic && group.size() == group.capacity()) {
            if (std::binary_search(molids.begin(), molids.end(), group.id)) {
                molecule_index.push_back(&group - &spc.groups.front());
            }
        }
    }
}

} // namespace Move
} // namespace Faunus
