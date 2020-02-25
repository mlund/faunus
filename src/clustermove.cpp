#include "clustermove.h"
#include "aux/eigensupport.h"
#include <algorithm>

namespace Faunus {
namespace Move {

double Cluster::clusterProbability(const Cluster::Tgroup &g1, const Cluster::Tgroup &g2) const {
    if (spc.geo.sqdist(g1.cm, g2.cm) <= group_thresholds(g1.id, g2.id))
        return 1.0;
    return 0.0;
}
void Cluster::_to_json(json &j) const {
    using namespace u8;
    j = {{"dir", dir},
         {"dp", dptrans},
         {"dprot", dprot},
         {"spread", spread},
         {"dirrot", dirrot},
         {rootof + bracket("r" + squared), std::sqrt(msqd.avg())},
         {rootof + bracket(theta + squared) + "/" + degrees, std::sqrt(msqd_angle.avg()) / 1.0_deg},
         {bracket("N"), N.avg()},
         {"bias rejection rate", double(bias_rejected) / cnt},
         {"clusterdistribution", clusterSizeDistribution}};
    _roundjson(j, 3);

    // print threshold matrix
    auto &_j = j["threshold"];
    for (auto i : molids)
        for (auto j : molids)
            if (i >= j) {
                auto str = Faunus::molecules[i].name + " " + Faunus::molecules[j].name;
                _j[str] = std::sqrt(group_thresholds(i, j));
                _roundjson(_j[str], 3);
            }

    // print satellite molecules
    if (not satellites.empty()) {
        auto &_j = j["satellites"];
        _j = json::array();
        for (auto id : satellites)
            _j.push_back(Faunus::molecules[id].name);
    }
}
void Cluster::_from_json(const json &j) {
    dptrans = j.at("dp");
    dir = j.value("dir", Point(1, 1, 1));
    dirrot = j.value("dirrot", Point(0, 0, 0)); // predefined axis of rotation
    dirrot.normalize();                         // make sure dirrot is a unit-vector
    dprot = j.at("dprot");
    spread = j.value("spread", true);
    molecule_names = j.at("molecules").get<decltype(molecule_names)>(); // molecule names
    molids = names2ids(molecules, molecule_names);                      // names --> molids

    updateMoleculeIndex();

    repeat = std::max<size_t>(1, molecule_index.size());

    // read satellite ids (molecules NOT to be considered as cluster centers)
    auto satnames = j.value("satellites", std::vector<std::string>()); // molecule names
    auto vec = names2ids(Faunus::molecules, satnames);                 // names --> molids
    satellites = std::set<int>(vec.begin(), vec.end());

    // make sure the satellites are in `molids`
    for (auto id : satellites)
        if (std::find(molids.begin(), molids.end(), id) == molids.end())
            throw std::runtime_error("satellite molecules must be defined in `molecules`");

    // read cluster thresholds
    if (j.count("threshold") == 1) {
        auto &_j = j.at("threshold");
        // threshold is given as a single number
        if (_j.is_number()) {
            for (auto i : molids)
                for (auto j : molids)
                    if (i >= j)
                        group_thresholds.set(i, j, std::pow(_j.get<double>(), 2));
        }
        // threshold is given as pairs of clustering molecules
        else if (_j.is_object()) {
            for (auto it = _j.begin(); it != _j.end(); ++it) {
                auto v = words2vec<std::string>(it.key());
                if (v.size() == 2) {
                    auto it1 = findName(Faunus::molecules, v[0]);
                    auto it2 = findName(Faunus::molecules, v[1]);
                    if (it1 == Faunus::molecules.end() or it2 == Faunus::molecules.end())
                        throw std::runtime_error("unknown molecule(s): ["s + v[0] + " " + v[1] + "]");
                    group_thresholds.set(it1->id(), it2->id(), std::pow(it.value().get<double>(), 2));
                } else
                    throw std::runtime_error("threshold requires exactly two space-separated molecules");
            }
        } else
            throw std::runtime_error("threshold must be a number or object");
    }
}
void Cluster::findCluster(Space &spc, size_t first, std::set<size_t> &cluster) {
    assert(first < spc.p.size());
    std::set<size_t> pool(molecule_index.begin(), molecule_index.end());
    assert(pool.count(first) > 0);

    cluster.clear();
    cluster.insert(first);
    pool.erase(first);

    size_t n;
    do { // find cluster (not very clever...)
    start:
        n = cluster.size();
        for (size_t i : cluster)
            if (not spc.groups.at(i).empty()) // check if group is inactive
                for (size_t j : pool)
                    if (i != j)
                        if (not spc.groups.at(j).empty()) { // check if group is inactive
                            // probability to cluster
                            double P = clusterProbability(spc.groups.at(i), spc.groups.at(j));
                            if (Movebase::slump() <= P) {
                                cluster.insert(j);
                                pool.erase(j);
				if(spread)
				  goto start; // wow, first goto ever!
                            }
                        }
    } while (cluster.size() != n);

    // check if cluster is too large
    double max = spc.geo.getLength().minCoeff() / 2;
    for (auto i : cluster)
        for (auto j : cluster)
            if (j > i)
                if (spc.geo.sqdist(spc.groups.at(i).cm, spc.groups.at(j).cm) >= max * max)
                    rotate = false; // skip rotation if cluster larger than half the box length
}
void Cluster::_move(Change &change) {
    _bias = 0;
    rotate = true;
    updateMoleculeIndex();
    if (not molecule_index.empty()) {
        std::set<size_t> cluster; // all group index in cluster

        // find "nuclei" or cluster center and exclude any molecule id listed as "satellite".
        size_t index_of_nuclei;
        do {
            index_of_nuclei = *slump.sample(molecule_index.begin(), molecule_index.end()); // random molecule (nuclei)
        } while (satellites.count(spc.groups[index_of_nuclei].id) != 0);

        findCluster(spc, index_of_nuclei, cluster); // find cluster around first

        N += cluster.size();                       // average cluster size
        clusterSizeDistribution[cluster.size()]++; // update cluster size distribution
        Change::data d;
        d.all = true;
        dp = ranunit(slump, dir) * dptrans * slump();

        auto boundary = spc.geo.getBoundaryFunc();

        // lambda function to calculate cluster COM
        auto clusterCOM = [&]() {
            double mass_sum = 0;
            Point cm(0, 0, 0);
            Point O = spc.groups[*cluster.begin()].cm;
            for (auto i : cluster) { // loop over clustered molecules (index)
                auto &g = spc.groups[i];
                Point t = g.cm - O;
                boundary(t);
                double m = g.mass();
                cm += m * t;
                mass_sum += m;
            }
            cm = cm / mass_sum + O;
            boundary(cm);
            return cm;
        };

        Point COM = clusterCOM(); // org. cluster center
        Eigen::Quaterniond Q;
        if (rotate) {
            Point u;
            if (dirrot.count() > 0)
                u = dirrot;
            else
                u = ranunit(slump);
            angle = dprot * (slump() - 0.5);
            Q = Eigen::AngleAxisd(angle, u); // quaternion
        } else
            angle = 0;

        for (auto i : cluster) { // loop over molecules in cluster
            auto &g = spc.groups[i];
            if (rotate) {
                Geometry::rotate(g.begin(), g.end(), Q, boundary, -COM);
                g.cm = g.cm - COM;
                boundary(g.cm);
                g.cm = Q * g.cm + COM;
                boundary(g.cm);
            }
            g.translate(dp, boundary);
            d.index = i;
            change.groups.push_back(d);
        }

        change.moved2moved = false; // do not calc. internal cluster energy

        // Reject if cluster composition changes during move
        // Note: this only works for the binary 0/1 probability function
        // currently implemented in `findCluster()`.

        std::set<size_t> aftercluster;         // all group index in cluster _after_move
        findCluster(spc, index_of_nuclei, aftercluster); // find cluster around first
        if (aftercluster == cluster)
            _bias = 0;
        else {
            _bias = pc::infty; // bias is infinite --> reject
            bias_rejected++;   // count how many time we reject due to bias
        }
#ifndef NDEBUG
        // check if cluster mass center movement matches displacement
        if (_bias == 0) {
            Point newCOM = clusterCOM();          // org. cluster center
            Point d = spc.geo.vdist(COM, newCOM); // distance between new and old COM
            double _zero = (d + dp).norm();       // |d+dp| should ideally be zero...
            if (std::fabs(_zero) > 1e-9) {
                _bias = pc::infty; // by setting bias=oo the move is rejected
                faunus_logger->warn("Skipping too large cluster: COM difference = {}", _zero);
            }
        }
#endif
    }
}
double Cluster::bias(Change &, double, double) { return _bias; }
void Cluster::_reject(Change &) {
    msqd += 0;
    msqd_angle += 0;
}
void Cluster::_accept(Change &) {
    msqd += dp.squaredNorm();
    msqd_angle += angle * angle;
}
Cluster::Cluster(Space &spc) : spc(spc) {
    cite = "doi:10/cj9gnn";
    name = "cluster";
    repeat = -1;
}

/**
 * search for molecules participating in the cluster
 * move. This is called for every move event as a
 * grand canonical move may have changed the number
 * of particles.
 */
void Cluster::updateMoleculeIndex() {
    assert(not molids.empty());
    molecule_index.clear();
    for (auto &g : spc.groups) {          // loop over all groups
        if (not g.atomic)                 // only molecular groups
            if (g.size() == g.capacity()) // only active particles
                if (std::find(molids.begin(), molids.end(), g.id) != molids.end())
                    molecule_index.push_back(&g - &spc.groups.front());
    }
}
} // namespace Move
} // namespace Faunus
