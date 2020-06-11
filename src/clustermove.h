#pragma once

#include "move.h"
#include <set>

namespace Faunus {
namespace Move {
/**
 * @brief Molecular cluster move
 */
class Cluster : public Movebase {
  private:
    typedef typename Space::Tpvec Tpvec;
    typedef typename Space::Tgroup Tgroup;
    Tspace &spc;
    Average<double> msqd, msqd_angle, N;
    double dptrans = 0, dprot = 0, angle = 0, _bias = 0;
    size_t bias_rejected = 0;
    bool rotate = true;        // true if cluster should be rotated
    bool single_layer = false; // stop cluster search after first layer of neighbors
    Point dir = {1, 1, 1}, dp;
    Point dirrot = {0, 0, 0}; // predefined axis of rotation
    std::vector<std::string> molecule_names; // names of molecules to be considered
    std::vector<int> molids;                 // molecule id's of molecules to be considered (must be sorted!)
    std::set<int> satellites;           // subset of molecule id's to cluster, but NOT act as nuclei (cluster centers)
    std::vector<size_t> molecule_index; // index of all possible molecules to be considered
    std::map<size_t, size_t> clusterSizeDistribution; // distribution of cluster sizes
    PairMatrix<double, true> group_thresholds;

    void updateMoleculeIndex(); //!< update `molecule_index`

    virtual double clusterProbability(const Tgroup &g1, const Tgroup &g2) const;

    void _to_json(json &j) const override;
    void _from_json(const json &j) override;

    /**
     * @param spc Space
     * @param seed_index Index of initial molecule (randomly selected)
     * @param index w. all molecules clustered around seed_index (seed_index included)
     */
    void findCluster(Space &spc, size_t seed_index, std::vector<size_t> &cluster);
    void _move(Change &change) override;
    double bias(Change &, double, double) override; //!< adds extra energy change not captured by the Hamiltonian
    void _reject(Change &) override;
    void _accept(Change &) override;

  public:
    Cluster(Space &spc);
};

} // namespace Move
} // namespace Faunus
