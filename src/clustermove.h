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
        auto positions = groups | transform(&Group::mass_center);
        auto masses = groups | transform(&Group::mass);
        return Geometry::gyration(positions.begin(), positions.end(), masses.begin(), mass_center_of_groups, boundary);
    }

    /** brief Calculates the gyration tensor for a collection of groups */
    template <typename Range>
    Tensor gyrationFromParticlePositions(
        const Range &groups, const Point &mass_center_of_groups,
        const Geometry::BoundaryFunction boundary = [](auto &) {}) {
        using namespace ranges::cpp20::views;
        auto positions = groups | join | transform(&Particle::pos);
        auto masses = groups | join | transform(&Particle::traits) | transform(&AtomData::mw);
        return Geometry::gyration(positions.begin(), positions.end(), masses.begin(), mass_center_of_groups, boundary);
    }

    friend void to_json(json &, const ClusterShapeAnalysis &);

  public:
    template <typename Range> void sample(const Range &groups, const Point &mass_center_of_groups, const Space &spc) {
        const auto cluster_size = std::distance(groups.begin(), groups.end());
        Tensor gyration_tensor;
        if (shape_anisotropy_use_com) {
            gyration_tensor =
                gyrationFromMassCenterPositions(groups, mass_center_of_groups, spc.geometry.getBoundaryFunc());
        } else {
            gyration_tensor =
                gyrationFromParticlePositions(groups, mass_center_of_groups, spc.geometry.getBoundaryFunc());
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
            pqr.save(it->second, groups, spc.geometry.getLength());
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
    const Space& spc;
    bool single_layer = false;               //!< stop cluster search after first layer of neighbors
    std::vector<std::string> molecule_names; //!< names of molecules to be considered
    std::vector<int> molids;                 //!< molecule id's of molecules to be considered (must be sorted!)
    std::set<int> satellites; //!< subset of molecule id's to cluster, but NOT act as nuclei (cluster centers)
    PairMatrix<double, true> thresholds_squared; //!< Cluster thresholds for pairs of groups

    void parseThresholds(const json &j); //!< Read thresholds from json input
    double clusterProbability(const Group& group1, const Group& group2) const;
    void registerSatellites(const std::vector<std::string> &); //!< Register satellites
    void updateMoleculeIndex();                                //!< update `molecule_index`
    friend void to_json(json &j, const FindCluster &cluster);

  public:
    std::vector<size_t> molecule_index;             //!< index of all possible molecules to be considered
    std::optional<size_t> findSeed(Random &random); //!< Find first group; exclude satellites
    std::pair<std::vector<size_t>, bool> findCluster(size_t seed_index); //!< Find cluster
    FindCluster(const Space& spc, const json& j);
};

void to_json(json &j, const FindCluster &cluster); //!< Serialize to json

/** @brief Helper base class for translating and rotating groups */
class GroupMover {
  public:
    Average<double> mean_square_displacement; //!< Must be updated manually
    const double displacement_factor;         //!< Displacement to be scaled by a random number
    const Point direction;                    //!< Directions or axis of translation or rotation
    GroupMover(double displacement_factor, const Point& direction);
};

/** @brief Helper class to translate a group */
class GroupTranslator : public GroupMover {
  public:
    Point current_displacement = Point::Zero(); //!< Currently set displacement
    GroupTranslator(double displacement_factor, const Point& direction = Point::Ones());
    std::function<void(Group&)> getLambda(Geometry::BoundaryFunction boundary, Random& slump);
};

/** @brief Helper class to translate a group
 *
 * If `direction` is zero, then generate random rotation axis
 */
class GroupRotator : public GroupMover {
  private:
    Eigen::Quaterniond setRandomRotation(Random& random);

  public:
    double current_displacement = 0.0; //!< Current displacement angle
    GroupRotator(double displacement_factor, const Point& rotation_axis = Point::Zero());
    std::function<void(Group&)> getLambda(Geometry::BoundaryFunction boundary, const Point& rotation_origin,
                                          Random& random); //!< Lambda to rotate group in cluster
};

/**
 * @brief Molecular cluster move
 */
class Cluster : public MoveBase {
  private:
    std::unique_ptr<FindCluster> find_cluster;
    std::unique_ptr<ClusterShapeAnalysis> shape_analysis;
    std::unique_ptr<GroupTranslator> translate; //!< Helper to translate group in cluster
    std::unique_ptr<GroupRotator> rotate;       //!< Helper to rotate group in cluster

    std::function<Group&(size_t)> index_to_group;
    Average<double> average_cluster_size;
    double bias_energy = 0;                     //!< Current bias energy (currently zero or infinity)
    unsigned int bias_rejected = 0;             //!< Number of times rejection occurred due to bias rejection
    unsigned int shape_analysis_interval = 0;   //!< Number of sample events between shape analysis

    Point clusterMassCenter(const std::vector<size_t>& indices) const; //!< Calculates cluster mass center
    void _move(Change& change) override;                               //!< Performs the move
    double bias(Change& change, double old_energy,
                double new_energy) override; //!< adds extra energy change not captured by the Hamiltonian
    void _reject(Change& change) override;
    void _accept(Change& change) override;
    void _to_json(json& j) const override;
    void _from_json(const json& j) override;
    void setChange(Change& change, const std::vector<size_t>& group_indices) const;
    void calculateBias(size_t seed_index, const std::vector<size_t>& cluster_index);
    void clearForMove(); //!< Clear displacements and bias

  protected:
    Cluster(Space& spc, std::string_view name, std::string_view cite);

  public:
    explicit Cluster(Space &spc);
};

} // namespace Move
} // namespace Faunus
