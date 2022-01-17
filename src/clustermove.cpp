#include "clustermove.h"
#include "aux/eigensupport.h"
#include <algorithm>
#include <range/v3/view/cartesian_product.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/algorithm/for_each.hpp>
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

FindCluster::FindCluster(const Space& spc, const json& j)
    : spc(spc) {
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

double FindCluster::clusterProbability(const Group& group1, const Group& group2) const {
    if (!group1.empty() and !group2.empty()) {
        if (spc.geometry.sqdist(group1.mass_center, group2.mass_center) <= thresholds_squared(group1.id, group2.id)) {
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
    namespace rv = ranges::cpp20::views;
    auto group_to_index = [&](auto& group) { return &group - &spc.groups.front(); };
    auto matching_molid = [&](auto& group) {
        if (group.isAtomic() || group.end() != group.trueend()) {
            return false; // skip atomic and incomplete groups
        } else {
            return std::binary_search(molids.begin(), molids.end(), group.id);
        }
    };
    molecule_index = spc.groups | rv::filter(matching_molid) | rv::transform(group_to_index) |
                     ranges::to<decltype(molecule_index)>();
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
                const auto molecule1 = findMoleculeByName(name_pair[0]);
                const auto molecule2 = findMoleculeByName(name_pair[1]);
                thresholds_squared.set(molecule1.id(), molecule2.id(), std::pow(value.get<double>(), 2));
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
    auto avoid_satellites = [&](auto index) { return (satellites.count(spc.groups.at(index).id) == 0); };
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
    namespace rv = ranges::cpp20::views;
    assert(seed_index < spc.particles.size());
    std::set<size_t> pool(molecule_index.begin(), molecule_index.end()); // decreasing pool of candidate groups
    assert(pool.count(seed_index) == 1);

    std::vector<size_t> cluster;            // molecular index
    cluster.reserve(molecule_index.size()); // ensures safe resizing without invalidating iterators
    cluster.push_back(seed_index);          // 'seed_index' is the index of the seed molecule
    pool.erase(seed_index);                 // ...which is already in the cluster and not part of pool

    // cluster search algorithm
    for (auto it1 = cluster.begin(); it1 != cluster.end(); it1++) {
        for (auto it2 = pool.begin(); it2 != pool.end();) {
            const auto p = clusterProbability(spc.groups.at(*it1), spc.groups.at(*it2)); // probability to cluster
            if (MoveBase::slump() <= p) {                                                // is group part of cluster?
                cluster.push_back(*it2);                                                 // yes, expand cluster...
                it2 = pool.erase(it2);                                                   // ...and remove from pool
            } else {
                ++it2; // not part of cluster, move along to next group in pool
            }
        }
        if (single_layer) { // stop after one iteration around 'seed_index'
            break;
        }
    }
    std::sort(cluster.begin(), cluster.end()); // required for correct energy evaluation

    // check if cluster is too large to be rotated (larger than half the box length)
    auto mass_centers = cluster | rv::transform([&](auto i) -> const Point& { return spc.groups.at(i).mass_center; });
    const auto max = 0.5 * spc.geometry.getLength().minCoeff();
    bool safe_to_rotate = true;
    for (auto i = mass_centers.begin(); i != mass_centers.end(); ++i) {
        for (auto j = i; ++j != mass_centers.end();) {
            if (spc.geometry.sqdist(*i, *j) < max * max) {
                continue;
            }
            safe_to_rotate = false;
            break;
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
                const auto key = fmt::format("{} {}", Faunus::molecules.at(k).name, Faunus::molecules.at(l).name);
                _j[key] = std::sqrt(cluster.thresholds_squared(k, l));
                roundJSON(_j[key], 3);
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
    j = {{"dir", translate->direction},
         {"dp", translate->displacement_factor},
         {"dprot", rotate->displacement_factor},
         {"dirrot", rotate->direction},
         {"√⟨r²⟩", std::sqrt(translate->mean_square_displacement.avg())},
         {"√⟨θ²⟩/°", std::sqrt(rotate->mean_square_displacement.avg()) / 1.0_deg},
         {"bias rejection rate", static_cast<double>(bias_rejected) / number_of_attempted_moves},
         {"⟨N⟩", average_cluster_size.avg()},
         {"cluster analysis", *shape_analysis},
         {"cluster analysis interval", shape_analysis_interval}};
    Move::to_json(j, *find_cluster);
    roundJSON(j, 3);
}

void Cluster::_from_json(const json &j) {
    translate = std::make_unique<GroupTranslator>(j.at("dp").get<double>(), j.value("dir", Point(1.0, 1.0, 1.0)));
    rotate = std::make_unique<GroupRotator>(j.at("dprot"), j.value("dirrot", Point(0, 0, 0)));
    find_cluster = std::make_unique<FindCluster>(spc, j);
    repeat = std::max<size_t>(1, find_cluster->molecule_index.size());

    auto analysis_json = j.value("analysis", json::object());
    shape_analysis = std::make_unique<ClusterShapeAnalysis>(analysis_json);
    shape_analysis_interval = analysis_json.value("interval", 10);
}

/**
 * @param indices Indices of groups in cluster
 * @return Mass center of the cluster
 */
Point Cluster::clusterMassCenter(const std::vector<size_t>& indices) const {
    assert(!indices.empty());
    namespace rv = ranges::cpp20::views;
    auto groups = indices | rv::transform(index_to_group);
    auto positions = groups | rv::transform(&Group::mass_center);
    auto masses = groups | rv::transform(&Group::mass);
    return Geometry::weightedCenter(positions, masses, spc.geometry.getBoundaryFunc(), -groups.begin()->mass_center);
}

void Cluster::_move(Change& change) {
    clearForMove();
    auto seed_index = find_cluster->findSeed(slump);
    if (!seed_index) {
        return;
    }
    const auto [cluster, safe_to_rotate] = find_cluster->findCluster(seed_index.value());
    const Point cluster_mass_center = clusterMassCenter(cluster); // cluster mass center
    auto groups = cluster | ranges::cpp20::views::transform(index_to_group);

    average_cluster_size += cluster.size(); // average cluster size
    if (number_of_attempted_moves % shape_analysis_interval == 0) {
        shape_analysis->sample(groups, cluster_mass_center, spc);
    }

    using ranges::cpp20::for_each;
    auto boundary = spc.geometry.getBoundaryFunc();
    if (safe_to_rotate) {
        for_each(groups, rotate->getLambda(boundary, cluster_mass_center, slump));
    }
    for_each(groups, translate->getLambda(boundary, slump));

    calculateBias(seed_index.value(), cluster);
    setChange(change, cluster);
}

void Cluster::clearForMove() {
    translate->current_displacement.setZero();
    rotate->current_displacement = 0.0;
    bias_energy = 0.0;
}

/**
 * Reject if cluster composition changes during move.
 * Current implementation works only for the binary 0/1 probability function
 * currently implemented in `findCluster()`.
 */
void Cluster::calculateBias(const size_t seed_index, const std::vector<size_t>& cluster_index) {
    [[maybe_unused]] auto [aftercluster, safe_to_rotate] = find_cluster->findCluster(seed_index);
    if (aftercluster == cluster_index) {
        bias_energy = 0.0;
    } else {
        bias_energy = pc::infty; // bias is infinite --> reject
        bias_rejected++;   // count how many times we reject due to bias
    }
}

/**
 * Update change object with changed molecules
 */
void Cluster::setChange(Change& change, const std::vector<size_t>& group_indices) const {
    change.groups.reserve(group_indices.size());
    std::for_each(group_indices.begin(), group_indices.end(), [&](auto index) { // register changes
        auto& change_data = change.groups.emplace_back();
        change_data.all = true;
        change_data.internal = false;
        change_data.group_index = index;
    });
    change.moved_to_moved_interactions = false; // do not calc. internal cluster energy
}

double Cluster::bias([[maybe_unused]] Change& change, [[maybe_unused]] double old_energy,
                     [[maybe_unused]] double new_energy) {
    return bias_energy;
}

void Cluster::_reject([[maybe_unused]] Change& change) {
    translate->mean_square_displacement += 0.0;
    rotate->mean_square_displacement += 0.0;
}

void Cluster::_accept([[maybe_unused]] Change& change) {
    translate->mean_square_displacement += translate->current_displacement.squaredNorm();
    rotate->mean_square_displacement += std::pow(rotate->current_displacement, 2);
}

Cluster::Cluster(Space& spc, std::string_view name, std::string_view cite)
    : MoveBase(spc, name, cite) {
    repeat = -1;
}

Cluster::Cluster(Space& spc)
    : Cluster(spc, "cluster", "doi:10/cj9gnn") {
    index_to_group = [&](auto index) -> Group& { return this->spc.groups.at(index); };
}

// -------------------------------------------

GroupMover::GroupMover(double displacement_factor, const Point& direction)
    : displacement_factor(displacement_factor)
    , direction(direction) {}

// -------------------------------------------

GroupTranslator::GroupTranslator(double displacement_factor, const Point& direction)
    : GroupMover(displacement_factor, direction) {}

/** Lambda to translate group in cluster. Will set a random displacement */
std::function<void(Group&)> GroupTranslator::getLambda(Geometry::BoundaryFunction boundary, Random& slump) {
    current_displacement = randomUnitVector(slump, direction) * displacement_factor * slump();
    if (std::fabs(displacement_factor) <= pc::epsilon_dbl) { // if zero displacement...
        return [](auto&) {};                                 // ...then return NOP lambda
    }
    return [&, boundary](Group& group) { group.translate(current_displacement, boundary); };
}

// -------------------------------------------

GroupRotator::GroupRotator(double displacement_factor, const Point& rotation_axis)
    : GroupMover(displacement_factor,
                 (rotation_axis.squaredNorm() > pc::epsilon_dbl) ? rotation_axis.normalized() : Point::Zero()) {}

Eigen::Quaterniond GroupRotator::setRandomRotation(Random& random) {
    Point axis = (direction.squaredNorm() > pc::epsilon_dbl) ? direction : randomUnitVector(random);
    current_displacement = displacement_factor * (random() - 0.5);
    return static_cast<Eigen::Quaterniond>(Eigen::AngleAxisd(current_displacement, axis));
}

/** Lambda to rotate a group. Will set a random rotation axis and angle. */
std::function<void(Group&)> GroupRotator::getLambda(Geometry::BoundaryFunction boundary, const Point& rotation_origin,
                                                    Random& random) {
    if (std::fabs(displacement_factor) <= pc::epsilon_dbl) { // if zero displacement...
        return [](auto&) {};                                 // ...then return NOP lambda
    }
    auto quaternion = setRandomRotation(random);
    return [=](Group& group) {
        Geometry::rotate(group.begin(), group.end(), quaternion, boundary, -rotation_origin);
        group.mass_center = group.mass_center - rotation_origin;
        boundary(group.mass_center);
        group.mass_center = quaternion * group.mass_center + rotation_origin;
        boundary(group.mass_center);
    };
}

// -------------------------------------------

} // namespace Move
} // namespace Faunus
