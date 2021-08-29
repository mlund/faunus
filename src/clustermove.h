#pragma once

#include "move.h"
#include "aux/pairmatrix.h"
#include <map>

namespace Faunus {
namespace Move {

/**
 * @brief Helper class for Cluster move for analysing shape and size of found clusters
 */
class ClusterShapeAnalysis {
  private:
    unsigned long int number_of_samples = 0;
    bool save_pqr_files = false;          //!< Save cluster to PQR files? One file for each size.
    bool shape_anisotropy_use_com = true; //!< true if kappa^2 should be based on mass centers instead of particles
    std::unique_ptr<std::ostream> stream; //!< log current cluster into (compressed) stream
    std::map<size_t, size_t> size_distribution;                                  //!< distribution of cluster sizes
    std::map<size_t, AverageObj<Geometry::ShapeDescriptors>> shape_distribution; //!< average shape
    std::map<size_t, std::ofstream> pqr_distribution;                            //!< Each cluster size is saved to disk
    decltype(pqr_distribution)::iterator findPQRstream(size_t cluster_size);     //!< Create or find PQR files stream

    /** brief Calculates the gyration tensor for a collection of groups */
    template <typename Range>
    Tensor gyrationFromMassCenterPositions(
        const Range &groups, const Point &mass_center_of_groups,
        const Geometry::BoundaryFunction boundary = [](auto &) {}) {
        using namespace ranges::cpp20::views;
        auto positions = groups | transform([](const auto &group) { return group.cm; });
        auto masses = groups | transform([](const auto &group) { return group.mass(); });
        return Geometry::gyration(positions.begin(), positions.end(), masses.begin(), mass_center_of_groups, boundary);
    }

    /** brief Calculates the gyration tensor for a collection of groups */
    template <typename Range>
    Tensor gyrationFromParticlePositions(
        const Range &groups, const Point &mass_center_of_groups,
        const Geometry::BoundaryFunction boundary = [](auto &) {}) {
        using namespace ranges::cpp20::views;
        auto positions = groups | join | transform([](const auto &particle) { return particle.pos; });
        auto masses = groups | join | transform([](const auto &particle) { return particle.traits().mw; });
        return Geometry::gyration(positions.begin(), positions.end(), masses.begin(), mass_center_of_groups, boundary);
    }

    friend void to_json(json &, const ClusterShapeAnalysis &);

  public:
    template <typename Range> void sample(const Range &groups, const Point &mass_center_of_groups, const Space &spc) {
        const auto cluster_size = std::distance(groups.begin(), groups.end());
        Tensor gyration_tensor;
        if (shape_anisotropy_use_com) {
            gyration_tensor = gyrationFromMassCenterPositions(groups, mass_center_of_groups, spc.geo.getBoundaryFunc());
        } else {
            gyration_tensor = gyrationFromParticlePositions(groups, mass_center_of_groups, spc.geo.getBoundaryFunc());
        }
        const auto shape = Geometry::ShapeDescriptors(gyration_tensor);
        size_distribution[cluster_size]++;
        shape_distribution[cluster_size] += shape;
        if (stream) {
            const auto seed_index = &(*groups.begin()) - &spc.groups.at(0); // index of first molecule in cluster
            *stream << fmt::format("{} {} {:.3f}\n", cluster_size, seed_index, shape.relative_shape_anisotropy);
        }
        if (save_pqr_files) {
            auto it = findPQRstream(cluster_size);
            it->second << fmt::format("REMARK   0 Sample number = {}\n", number_of_samples)
                       << fmt::format("REMARK   0 Relative shape anisotropy = {:.3f}\n",
                                      shape.relative_shape_anisotropy);
            PQRWriter pqr;
            pqr.style = PQRWriter::Style::PQR;
            pqr.save(it->second, groups, spc.geo.getLength());
        }
        ++number_of_samples;
    }
    ClusterShapeAnalysis(bool shape_anisotropy_use_com, const std::string &filename, bool dump_pqr_files);
    ClusterShapeAnalysis(const json &j);
};

void to_json(json &j, const ClusterShapeAnalysis &shape);

/**
 * @brief Helper class for Cluster move for finding clusters
 * @todo This is a break-out from `Cluster` and further separation of functionality
 * could be made
 */
class FindCluster {
  private:
    typedef typename Space::Tgroup Tgroup;
    Space &spc;
    bool single_layer = false;               //!< stop cluster search after first layer of neighbors
    std::vector<std::string> molecule_names; //!< names of molecules to be considered
    std::vector<MoleculeData::Tid> molids;   //!< molecule id's of molecules to be considered (must be sorted!)
    std::set<int> satellites; //!< subset of molecule id's to cluster, but NOT act as nuclei (cluster centers)
    PairMatrix<double, true> thresholds_squared; //!< Cluster thresholds for pairs of groups

    void parseThresholds(const json &j); //!< Read thresholds from json input
    double clusterProbability(const Tgroup &group1, const Tgroup &group2) const;
    void registerSatellites(const std::vector<std::string> &); //!< Register satellites
    void updateMoleculeIndex();                                //!< update `molecule_index`
    friend void to_json(json &j, const FindCluster &cluster);

  public:
    std::vector<size_t> molecule_index;             //!< index of all possible molecules to be considered
    std::optional<size_t> findSeed(Random &random); //!< Find first group; exclude satellites
    std::pair<std::vector<size_t>, bool> findCluster(size_t seed_index); //!< Find cluster
    FindCluster(Space &spc, const json &j);
};

void to_json(json &j, const FindCluster &cluster); //!< Serialize to json

/**
 * @brief Molecular cluster move
 */
class Cluster : public MoveBase {
  private:
    typedef typename Space::Tgroup Tgroup;
    std::shared_ptr<FindCluster> find_cluster;
    std::shared_ptr<ClusterShapeAnalysis> shape_analysis;
    Average<double> average_cluster_size;
    Average<double> translational_mean_square_displacement;
    Average<double> rotational_mean_square_displacement;
    double _bias = 0;                           //!< Current bias energy (currently zero or infinity)
    double translation_displacement_factor = 0; //!< User defined scaling of translation
    double rotation_displacement_factor = 0;    //!< User defined scaling of rotation
    double rotation_angle = 0;                  //!< Current rotation angle
    unsigned int bias_rejected = 0;             //!< Number of times rejection occurred due to bias rejection
    unsigned int shape_analysis_interval = 0;   //!< Number of sample events between shape analysis
    Point translation_direction = {1, 1, 1};    //!< Directions which to translate
    Point translation_displacement = {0, 0, 0}; //!< Current translational displacement
    Point rotation_axis = {0, 0, 0};            //! predefined axis of rotation (if zero, generate randomly)

    Point clusterMassCenter(const std::vector<size_t> &) const; //!< Calculates cluster mass center
    Eigen::Quaterniond setRandomRotation();                     //!< Sets random rotation angle and axis
    void _move(Change &change) override;                        //!< Performs the move
    double bias(Change &, double, double) override; //!< adds extra energy change not captured by the Hamiltonian
    void _reject(Change &) override;
    void _accept(Change &) override;
    void _to_json(json &j) const override;
    void _from_json(const json &molid) override;
    Change createChangeObject(const std::vector<size_t> &cluster_index) const;
    void biasRejectOrAccept(const size_t seed_index, const std::vector<size_t> &cluster_index);

  protected:
    using MoveBase::spc;
    Cluster(Space &spc, std::string name, std::string cite);

  public:
    explicit Cluster(Space &spc);
};

} // namespace Move
} // namespace Faunus
