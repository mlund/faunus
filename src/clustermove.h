#pragma once

#include "move.h"
#include "externalpotential.h"
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
    bool rotate; // true if cluster should be rotated
    bool spread; // true if cluster should spread outside of the first layer of surrounding molecules
    Point dir = {1, 1, 1}, dp;
    Point dirrot = {0, 0, 0}; // predefined axis of rotation
    std::vector<std::string> molecule_names; // names of molecules to be considered
    std::vector<int> molids;                 // molecule id's of molecules to be considered
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

/**
 * @brief Non-rejective cluster translation.
 *
 * This type of move will attempt to translate collective sets of macromolecules that
 * with a symmetric transition matrix (no flow through the clusters).
 * See detailed description [here](http://dx.doi.org/10/fthw8k).
 *
 * Setting the boolean `skipEnergyUpdate` to true (default is false) updates of the
 * total energy are skipped to speed up the move.
 * While this has no influence on the Markov chain it will cause an apparent energy
 * drift. It is recommended that this is enabled only for long production runs after
 * having properly checked that no drifts occur with `skipEnergyUpdate=false`.
 *
 * Keyword     | Description
 * :-----------| :---------------------------------------------
 * `dp`        | Displacement parameter (default: 0)
 * `skipenergy`| Skip energy update, see above (default: false)
 * `prob`      | Runfraction (default: 1.0)
 *
 * @note Requirements for usage:
 * - Compatible only with purely molecular systems
 * - Works only with periodic containers
 * - External potentials are ignored
 *
 * @author Bjoern Persson
 * @date Lund 2009-2010
 * @warning Legacy code under construction
 * @todo Energy calc. before and after can be optimized by only looping over `moved` with
 * `remaining`.
 */
class ClusterTranslateNR : public Movebase {
  private:
    std::vector<int> moved;
    std::vector<int> remaining;
    Space &spc;
    Point dir = {1, 1, 1};             //!< Directions to translate
    double translational_displacement; //!< Displacement parameter [aa]
    Average<double> move_fraction;     //!< Fraction of groups moved

    void _move(Change &) override;
    void _accept(Change &) override;
    void _reject(Change &) override;
    void _to_json(json &) const override;

  public:
    ClusterTranslateNR(Space &, json &);
};

} // namespace Move
} // namespace Faunus
