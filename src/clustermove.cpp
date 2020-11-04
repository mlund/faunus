#include "clustermove.h"
#include "aux/eigensupport.h"
#include <algorithm>
#include <range/v3/view/cartesian_product.hpp>
#include <range/v3/range/conversion.hpp>
#include <set>

namespace Faunus {
namespace Move {

ClusterShapeAnalysis::ClusterShapeAnalysis(bool shape_anisotropy_use_com, const std::string &filename,
                                           bool dump_pqr_files)
    : save_pqr_files(dump_pqr_files), shape_anisotropy_use_com(shape_anisotropy_use_com) {
    if (!filename.empty()) {
        if (stream = IO::openCompressedOutputStream(MPI::prefix + filename); !*stream) { // may be gzip compressed
            throw std::runtime_error("could not create "s + filename);
        } else {
            *stream << "# size seed shape_anisotropy\n";
        }
    }
}
ClusterShapeAnalysis::ClusterShapeAnalysis(const json &j)
    : ClusterShapeAnalysis(j.value("com", true), j.value("file", std::string()), j.value("save_pqr", false)) {}

decltype(ClusterShapeAnalysis::pqr_distribution)::iterator ClusterShapeAnalysis::findPQRstream(size_t cluster_size) {
    if (auto it = pqr_distribution.find(cluster_size); it == pqr_distribution.end()) {
        const std::string file = MPI::prefix + fmt::format("clustersize{}.pqr", cluster_size);
        return pqr_distribution.emplace(cluster_size, std::ofstream(file)).first;
    } else {
        return it;
    }
}

void to_json(json &j, const ClusterShapeAnalysis &shape) {
    j["size distribution"] = shape.size_distribution;
    auto &json_shape = j["shape distribution"] = json::array();
    for (const auto &[size, shape_anisotropy] : shape.shape_distribution) {
        json_shape.push_back({size, shape_anisotropy.avg()});
    }
}

FindCluster::FindCluster(Space &spc, const json &j) : spc(spc) {
    single_layer = j.value("single_layer", false);
    if (j.count("spread")) {
        faunus_logger->error("cluster: 'spread' is deprecated, use 'single_layer' instead");
    }
    molecule_names = j.at("molecules").get<decltype(molecule_names)>(); // molecule names
    molids = Faunus::names2ids(Faunus::molecules, molecule_names);      // names --> molids
    std::sort(molids.begin(), molids.end());
    updateMoleculeIndex();

    // read satellite ids (molecules NOT to be considered as cluster centers)
    registerSatellites(j.value("satellites", std::vector<std::string>()));
    parseThresholds(j.at("threshold"));
}

double FindCluster::clusterProbability(const FindCluster::Tgroup &group1, const FindCluster::Tgroup &group2) const {
    if (!group1.empty() and !group2.empty()) {
        if (spc.geo.sqdist(group1.cm, group2.cm) <= thresholds_squared(group1.id, group2.id)) {
            return 1.0;
        }
    }
    return 0.0;
}

/**
 * Search for molecules participating in the cluster move.
 * This should be called for every move event as a
 * grand canonical move may have changed the number
 * of particles.
 *
 * Atomic groups and inactive groups are ignored.
 */
void FindCluster::updateMoleculeIndex() {
    auto to_group_index = [&](const Tgroup &group) { return &group - &spc.groups.front(); };
    auto matching_molid = [&](const Tgroup &group) {
        if (group.isAtomic() || group.end() != group.trueend()) {
            return false; // skip atomic and incomplete groups
        } else {
            return std::binary_search(molids.begin(), molids.end(), group.id);
        }
    };
    molecule_index = spc.groups | ranges::cpp20::views::filter(matching_molid) |
                     ranges::cpp20::views::transform(to_group_index) | ranges::to<decltype(molecule_index)>();
}

void FindCluster::registerSatellites(const std::vector<std::string> &satellite_names) {
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
void FindCluster::parseThresholds(const json &j) {
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

/**
 * Satellites are molecules that cannot be used to seed a cluster but can be
 * part of it, i.e. rotated and translated. Here we therefore filter out any
 * satellite molecules, and pick a random molecule index of what remains.
 */
std::optional<size_t> FindCluster::findSeed(Random &random) {
    updateMoleculeIndex();
    auto avoid_satellites = [&](auto index) { return (satellites.count(spc.groups[index].id) == 0); };
    auto not_satellites = molecule_index | ranges::cpp20::views::filter(avoid_satellites);
    if (ranges::cpp20::empty(not_satellites)) {
        return std::nullopt;
    } else {
        return *random.sample(not_satellites.begin(), not_satellites.end());
    }
}

/**
 * Find cluster
 *
 * @param seed_index Index of seed_index group to evaluate the cluster around
 * @returns Pair w. vector for group indices in cluster and a bool if the cluster
 *          can be safely rotated in a PBC environment
 */
std::pair<std::vector<size_t>, bool> FindCluster::findCluster(size_t seed_index) {
    assert(seed_index < spc.p.size());
    std::set<size_t> pool(molecule_index.begin(), molecule_index.end());
    assert(pool.count(seed_index) == 1);

    std::vector<size_t> cluster;            // molecular index
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

    // check if cluster is too large to be rotated
    // (larger than half the box length)
    bool safe_to_rotate = true;
    const double max = spc.geo.getLength().minCoeff() / 2;
    for (auto i : cluster) {
        for (auto j : cluster) {
            if (j > i) {
                if (spc.geo.sqdist(spc.groups.at(i).cm, spc.groups.at(j).cm) >= max * max) {
                    safe_to_rotate = false;
                    break;
                }
            }
        }
    }
    assert(std::adjacent_find(cluster.begin(), cluster.end()) == cluster.end()); // check for duplicates
    return {cluster, safe_to_rotate};
}

void to_json(json &j, const FindCluster &cluster) {
    j["single_layer"] = cluster.single_layer;
    auto &_j = j["threshold"];
    for (const auto k : cluster.molids) {
        for (const auto l : cluster.molids) {
            if (k >= l) {
                const std::string key = Faunus::molecules[k].name + " " + Faunus::molecules[l].name;
                _j[key] = std::sqrt(cluster.thresholds_squared(k, l));
                _roundjson(_j[key], 3);
            }
        }
    }

    // print satellite molecules
    if (not cluster.satellites.empty()) {
        auto &satellites_json = j["satellites"] = json::array();
        for (auto id : cluster.satellites) {
            satellites_json.push_back(Faunus::molecules[id].name);
        }
    }
}

void Cluster::_to_json(json &j) const {
    j = {{"dir", translation_direction},
         {"dp", translation_displacement_factor},
         {"dprot", rotation_displacement_factor},
         {"dirrot", rotation_axis},
         {"√⟨r²⟩", std::sqrt(translational_mean_square_displacement.avg())},
         {"√⟨θ²⟩/°", std::sqrt(rotational_mean_square_displacement.avg()) / 1.0_deg},
         {"bias rejection rate", static_cast<double>(bias_rejected) / cnt},
         {"⟨N⟩", average_cluster_size.avg()},
         {"cluster analysis", *shape_analysis},
         {"cluster analysis interval", shape_analysis_interval}};
    Move::to_json(j, *find_cluster);
    _roundjson(j, 3);
}

void Cluster::_from_json(const json &j) {
    translation_displacement_factor = j.at("dp").get<double>();
    translation_direction = j.value("dir", Point(1, 1, 1));
    rotation_axis = j.value("dirrot", Point(0, 0, 0)); // predefined axis of rotation
    rotation_axis.normalize();                         // make sure dirrot is a unit-vector
    rotation_displacement_factor = j.at("dprot");
    find_cluster = std::make_shared<FindCluster>(spc, j);
    repeat = std::max<size_t>(1, find_cluster->molecule_index.size());

    auto analysis_json = j.value("analysis", json::object());
    shape_analysis = std::make_shared<ClusterShapeAnalysis>(analysis_json);
    shape_analysis_interval = analysis_json.value("interval", 10);
}

Point Cluster::clusterMassCenter(const std::vector<size_t> &cluster_index) const {
    auto boundary = spc.geo.getBoundaryFunc();
    double mass_sum = 0.0;
    Point mass_center(0, 0, 0);
    Point origin = spc.groups[*cluster_index.begin()].cm;
    for (auto i : cluster_index) { // loop over clustered molecules (index)
        Point r = spc.groups[i].cm - origin;
        boundary(r);
        const double mass = spc.groups[i].mass();
        mass_center += mass * r;
        mass_sum += mass;
    }
    mass_center = mass_center / mass_sum + origin;
    boundary(mass_center);
    return mass_center;
}

Eigen::Quaterniond Cluster::setRandomRotation() {
    Point axis = (rotation_axis.count() > 0) ? rotation_axis : ranunit(slump);
    rotation_angle = rotation_displacement_factor * (slump() - 0.5);
    return static_cast<Eigen::Quaterniond>(Eigen::AngleAxisd(rotation_angle, axis));
}

void Cluster::_move(Change &change) {
    _bias = 0.0;
    rotation_angle = 0.0;
    translation_displacement.setZero();
    if (auto seed_index = find_cluster->findSeed(slump)) {
        const auto [cluster, safe_to_rotate] = find_cluster->findCluster(seed_index.value());
        auto cluster_groups = cluster | ranges::cpp20::views::transform(
                                            [&](size_t index) -> Space::Tgroup & { return spc.groups.at(index); });

        const auto boundary = spc.geo.getBoundaryFunc();
        const Point COM = clusterMassCenter(cluster); // cluster mass center

        if (cnt % shape_analysis_interval == 0) {
            shape_analysis->sample(cluster_groups, COM, spc);
        }

        if (safe_to_rotate) {
            const auto quaternion = setRandomRotation(); // set rotation
            auto rotate_group = [&](Tgroup &group) {
                Geometry::rotate(group.begin(), group.end(), quaternion, boundary, -COM);
                group.cm = group.cm - COM;
                boundary(group.cm);
                group.cm = quaternion * group.cm + COM;
                boundary(group.cm);
            };
            std::for_each(cluster_groups.begin(), cluster_groups.end(), rotate_group);
        }

        translation_displacement = ranunit(slump, translation_direction) * translation_displacement_factor * slump();
        auto translate_group = [&](Tgroup &group) { group.translate(translation_displacement, boundary); };
        std::for_each(cluster_groups.begin(), cluster_groups.end(), translate_group);

        biasRejectOrAccept(seed_index.value(), cluster);

        average_cluster_size += cluster.size(); // average cluster size
        change = createChangeObject(cluster);

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

/**
 * Reject if cluster composition changes during move.
 * Current implementation works only for the binary 0/1 probability function
 * currently implemented in `findCluster()`.
 */
void Cluster::biasRejectOrAccept(const size_t seed_index, const std::vector<size_t> &cluster_index) {
    auto [aftercluster, safe_to_rotate] = find_cluster->findCluster(seed_index);
    if (aftercluster == cluster_index) {
        _bias = 0.0;
    } else {
        _bias = pc::infty; // bias is infinite --> reject
        bias_rejected++;   // count how many times we reject due to bias
    }
}

/**
 * Update change object with changed molecules
 */
Change Cluster::createChangeObject(const std::vector<size_t> &cluster_index) const {
    Change change;
    change.groups.reserve(cluster_index.size());
    std::for_each(cluster_index.begin(), cluster_index.end(), [&](auto index) { // register changes
        change.groups.emplace_back();
        auto &change_data = change.groups.back();
        change_data.all = true;
        change_data.internal = false;
        change_data.index = index;
    });
    change.moved2moved = false; // do not calc. internal cluster energy
    return change;
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

} // namespace Move
} // namespace Faunus
