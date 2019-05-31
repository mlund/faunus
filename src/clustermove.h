#pragma once

#include "move.h"

namespace Faunus {
namespace Move {
/**
 * @brief Molecular cluster move
 *
 * @todo fix so that it works w. GC molecules (index are calculating before simulation)
 */
class Cluster : public Movebase {
  private:
    typedef typename Space::Tpvec Tpvec;
    typedef typename Space::Tgroup Tgroup;
    Tspace &spc;
    Average<double> msqd, msqd_angle, N;
    double dptrans = 0, dprot = 0, angle = 0, _bias = 0;
    size_t bias_rejected = 0;
    bool rotate; // true if cluster should be rotated
    bool spread; // true if cluster should spread outside of the first layer of surrounding molecules
    Point dir = {1, 1, 1}, dp;
    std::vector<std::string> names; // names of molecules to be considered
    std::vector<int> ids;           // molecule id's of molecules to be considered
    std::set<int> satellites;       // subset of molecules to cluster, but NOT act as nuclei (cluster centers)
    std::vector<size_t> index;      // index of all possible molecules to be considered
    std::map<size_t, size_t> clusterSizeDistribution; // distribution of cluster sizes
    PairMatrix<double, true> thresholdsq;

    virtual double clusterProbability(const Tgroup &g1, const Tgroup &g2) const;

    void _to_json(json &j) const override;
    void _from_json(const json &j) override;

    /**
     * @param spc Space
     * @param first Index of initial molecule (randomly selected)
     * @param index w. all molecules clustered around first (first included)
     */
    void findCluster(Space &spc, size_t first, std::set<size_t> &cluster);
    void _move(Change &change) override;
    double bias(Change &, double, double) override; //!< adds extra energy change not captured by the Hamiltonian
    void _reject(Change &) override;
    void _accept(Change &) override;

  public:
    Cluster(Space &spc);
};

} // namespace Move
} // namespace Faunus
