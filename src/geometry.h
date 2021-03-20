#pragma once

#include "units.h"
#include "core.h"
#include "particle.h"
#include "tensor.h"
#include <Eigen/Geometry>
#include <iostream>
#include <cereal/types/base_class.hpp>

/** @brief Faunus main namespace */
namespace Faunus {

class Random;

/**
 * @brief Simulation geometries and related operations.
 *
 * Other parts of Faunus use directly only Chameleon geometry which serves as an interface. Based on the provided
 * configuration, Chameleon initializes an appropriate concrete implementation, which it encapsulates.
 *
 * To add a new geometry implementation, a class derived from GeometryImplementation is created. Geometry::Variant
 * enum type is extended and an initialization within Chameleon::makeGeometry() is provided. In order to make
 * geometry constructable from a json configuration, the map Chameleon::names is extended. When performance is
 * an issue, inlineable implementation of vdist and boundary can be added into respective methods of Chameleon.
 *
 * All geometry implementation shall be covered by unit tests.
 *
 */
namespace Geometry {

//! Function to apply PBC to a position, i.e. wrap around the borders if applicable for the given container geometry
typedef std::function<void(Point &)> BoundaryFunction;

//! Function to calculate the (minimum) distance between two points depending on contained geometry
typedef std::function<Point(const Point &, const Point &)> DistanceFunction;

//! Geometry variant used for Chameleon
enum Variant { CUBOID = 0, SPHERE, CYLINDER, SLIT, HEXAGONAL, OCTAHEDRON, HYPERSPHERE2D };

//! Various methods of volume scaling, @see GeometryBase::setVolume.
enum VolumeMethod { ISOTROPIC, ISOCHORIC, XY, Z, INVALID };

NLOHMANN_JSON_SERIALIZE_ENUM(VolumeMethod, {{VolumeMethod::INVALID, nullptr},
                                            {VolumeMethod::ISOTROPIC, "isotropic"},
                                            {VolumeMethod::ISOCHORIC, "isochoric"},
                                            {VolumeMethod::XY, "xy"},
                                            {VolumeMethod::Z, "z"}})

enum Coordinates { ORTHOGONAL, ORTHOHEXAGONAL, TRUNC_OCTAHEDRAL, NON3D };
enum Boundary { FIXED, PERIODIC };

/**
 * @brief A structure containing a type of boundary condition in each direction.
 *
 * A stub. It can be extended to fully json-configurable boundary conditions.
 */
struct BoundaryCondition {
    using BoundaryXYZ = Eigen::Matrix<Boundary, 3, 1>;
    // typedef std::pair<std::string, BoundaryXYZ> BoundaryName;
    // static const std::map<std::string, BoundaryXYZ> names; //!< boundary names

    Coordinates coordinates;
    BoundaryXYZ direction;

    template <class Archive> void serialize(Archive &archive) {
        archive(coordinates, direction);
    } //!< Cereal serialisation

    BoundaryCondition(Coordinates coordinates = ORTHOGONAL, BoundaryXYZ boundary = {FIXED, FIXED, FIXED})
        : coordinates(coordinates), direction(boundary){};
};

/**
 * @brief An interface for all geometries.
 */
struct GeometryBase {
    virtual Point setVolume(double, VolumeMethod = ISOTROPIC) = 0; //!< Set volume
    virtual double getVolume(int = 3) const = 0;                   //!< Get volume
    virtual void boundary(Point &) const = 0;                      //!< Apply boundary conditions
    virtual bool collision(const Point &) const = 0;               //!< Overlap with boundaries
    virtual void randompos(Point &, Random &) const = 0;           //!< Generate random position
    virtual Point vdist(const Point &, const Point &) const = 0;   //!< (Minimum) distance between two points
    virtual Point getLength() const = 0;                           //!< Side lengths
    virtual ~GeometryBase();
    virtual void to_json(json &j) const = 0;
    virtual void from_json(const json &j) = 0;

    inline BoundaryFunction getBoundaryFunc() const {
        return [this](Point &i) { boundary(i); };
    } //!< Lambda for applying boundary conditions on a point

    inline DistanceFunction getDistanceFunc() const {
        return [this](const Point &i, const Point &j) { return vdist(i, j); };
    } //!< Lambda for calculating the (minimum) distance vector between two positions

  protected:
    template <typename T = double> inline int anint(T x) const {
        return int(x > 0.0 ? x + 0.5 : x - 0.5);
    } //!< Round to int

}; //!< Base class for all geometries

/**
 * @brief A base class for various geometries implementations.
 */
class GeometryImplementation : public GeometryBase {
  public:
    BoundaryCondition boundary_conditions;

    virtual ~GeometryImplementation();

    //! A unique pointer to a copy of self. To be used in copy constructors.
    virtual std::unique_ptr<GeometryImplementation> clone() const = 0;

    //! Cereal serialisation
    template <class Archive> void serialize(Archive &archive) { archive(boundary_conditions); }
};

/**
 * @brief The cuboid geometry with periodic boundary conditions possibly applied in all three directions.
 */
class Cuboid : public GeometryImplementation {
  protected:
    Point box, box_half, box_inv;

  public:
    Point getLength() const override;
    double getVolume(int dim = 3) const final; // finalized to help the compiler with inlining
    void setLength(const Point &len);          // todo shall be protected
    Point setVolume(double volume, VolumeMethod method = ISOTROPIC) override;
    Point vdist(const Point &a, const Point &b) const override;
    void boundary(Point &a) const override;
    bool collision(const Point &a) const override;
    void randompos(Point &m, Random &rand) const override;
    void from_json(const json &j) override;
    void to_json(json &j) const override;
    Cuboid(const Point &side_length);
    Cuboid();

    std::unique_ptr<GeometryImplementation> clone() const override {
        return std::make_unique<Cuboid>(*this);
    }; //!< A unique pointer to a copy of self.

    //! Cereal serialisation
    template <class Archive> void serialize(Archive &archive) {
        archive(cereal::base_class<GeometryImplementation>(this), box);
    }
};

/**
 * @brief A legacy class for the cuboid geometry with periodic boundary conditions only in xy directions.
 *
 * @deprecated Shall be replaced by Cuboid with a proper periodic boundary set on initialization.
 */
class Slit : public Cuboid {
    using Tbase = Cuboid;

  public:
    Slit(const Point &p);
    Slit(double x, double y, double z);
    Slit(double x = 0.0);

    std::unique_ptr<GeometryImplementation> clone() const override {
        return std::make_unique<Slit>(*this);
    }; //!< A unique pointer to a copy of itself.
};

/**
 * @brief The spherical geometry where no periodic boundary condition could be applied.
 */
class Sphere : public GeometryImplementation {
  protected:
    double radius;

  public:
    Point getLength() const override;
    double getVolume(int dim = 3) const override;
    Point setVolume(double volume, VolumeMethod method = ISOTROPIC) override;
    Point vdist(const Point &a, const Point &b) const override;
    double sqdist(const Point &a, const Point &b) const { return (a - b).squaredNorm(); };
    void boundary(Point &a) const override;
    bool collision(const Point &point) const override;
    void randompos(Point &m, Random &rand) const override;
    void from_json(const json &j) override;
    void to_json(json &j) const override;
    Sphere(double radius = 0.0);

    std::unique_ptr<GeometryImplementation> clone() const override {
        return std::make_unique<Sphere>(*this);
    }; //!< A unique pointer to a copy of self.

    //! Cereal serialisation
    template <class Archive> void serialize(Archive &archive) {
        archive(cereal::base_class<GeometryImplementation>(this), radius);
    }
};

class Hypersphere2d : public Sphere {
  public:
    Point vdist(const Point &a, const Point &b) const override;
    bool collision(const Point &a) const override;
    void randompos(Point &m, Random &rand) const override;
    Hypersphere2d(double radius = 0.0);

    std::unique_ptr<GeometryImplementation> clone() const override {
        return std::make_unique<Hypersphere2d>(*this);
    }; //!< A unique pointer to a copy of self.
};

/**
 * @brief The cylindrical geometry with periodic boundary conditions in z-axis (the height of the cylinder).
 */
class Cylinder : public GeometryImplementation {
  protected:
    double radius, height;

  public:
    Point getLength() const override;
    double getVolume(int dim = 3) const override;
    Point setVolume(double volume, VolumeMethod method = ISOTROPIC) override;
    Point vdist(const Point &a, const Point &b) const override;
    void boundary(Point &a) const override;
    bool collision(const Point &a) const override;
    void randompos(Point &m, Random &rand) const override;
    void from_json(const json &j) override;
    void to_json(json &j) const override;
    Cylinder(double radius = 0.0, double height = 0.0);

    std::unique_ptr<GeometryImplementation> clone() const override {
        return std::make_unique<Cylinder>(*this);
    }; //!< A unique pointer to a copy of self.

    //! Cereal serialisation
    template <class Archive> void serialize(Archive &archive) {
        archive(cereal::base_class<GeometryImplementation>(this), radius, height);
    }
};

/**
 * @brief The hexagonal prism geometry with periodic boundary conditions.
 *
 * The prism is oriented in the coordination system as follows: z height, xy base ⬢ with a shorter length
 * (the diameter of an inscribed circle d = 2r) in x direction, and a longer length (the diameter of a
 * circumscribed circle D = 2R) in y direction.
 */
class HexagonalPrism : public GeometryImplementation {
    //! Change matrices from rhombic to cartesian coordinates and back.
    //! The Y (and Z) axis is identical in both coordination systems, while the X-axis is tilted by -30deg
    //! (i.e., clockwise), forming a 120deg angle with the Y axis in the rhombic coordinates.
    static const Eigen::Matrix3d rhombic2cartesian;
    static const Eigen::Matrix3d cartesian2rhombic;

    Point box; //!< x = inscribed circle diameter, y = circumscribed circle diameter, z = height
    void set_box(double side, double height);

  public:
    Point getLength() const override;
    double getVolume(int dim = 3) const override;
    Point setVolume(double volume, VolumeMethod method = ISOTROPIC) override;
    Point vdist(const Point &a, const Point &b) const override;
    void boundary(Point &a) const override;
    bool collision(const Point &a) const override;
    void randompos(Point &m, Random &rand) const override;
    void from_json(const json &j) override;
    void to_json(json &j) const override;
    HexagonalPrism(double side = 0.0, double height = 0.0);

    std::unique_ptr<GeometryImplementation> clone() const override {
        return std::make_unique<HexagonalPrism>(*this);
    }; //!< A unique pointer to a copy of self.

    //! Cereal serialisation
    template <class Archive> void serialize(Archive &archive) {
        archive(cereal::base_class<GeometryImplementation>(this), box);
    }

    double innerRadius() const; //!< Inner hexagonal radius
    double outerRadius() const; //!< Outer radius / side-length
    double height() const;      //!< Prism height
};

/**
 * @brief The truncated octahedron geoemtry with periodic boundary conditions in all directions.
 */
class TruncatedOctahedron : public GeometryImplementation {
    double side;

  public:
    Point getLength() const override;
    double getVolume(int dim = 3) const override;
    Point setVolume(double volume, VolumeMethod method = ISOTROPIC) override;
    Point vdist(const Point &a, const Point &b) const override;
    void boundary(Point &a) const override;
    bool collision(const Point &a) const override;
    void randompos(Point &pos, Random &rand) const override;
    void from_json(const json &j) override;
    void to_json(json &j) const override;
    TruncatedOctahedron(double side = 0.0);

    std::unique_ptr<GeometryImplementation> clone() const override {
        return std::make_unique<TruncatedOctahedron>(*this);
    }; //!< A unique pointer to a copy of self.

    //! Cereal serialisation
    template <class Archive> void serialize(Archive &archive) {
        archive(cereal::base_class<GeometryImplementation>(this), side);
    }
};

/**
 * @brief Geometry class for spheres, cylinders, cuboids, hexagonal prism, truncated octahedron, slits. It is
 * a wrapper of a concrete geometry implementation.
 *
 * The class re-implements the time-critical functions vdist and boundary for the orthogonal periodic boundary
 * conditions. Hence the call can be inlined by the compiler. That would not be possible otherwise due to the
 * polymorphism of the concrete implementations. Other functions calls are delegated directly to the concrete
 * implementation.
 *
 * Note that the class implements a copy constructor and overloads the assignment operator.
 *
 * @note
 * - [Efficient Coding of the Minimum Image Convention](http://doi.org/kvs)
 * - [Fast Coding of the Minimum Image Convention](http://doi.org/ck2nrd)
 *
 * @todo Implement unit tests
 */
class Chameleon : public GeometryBase {
  private:
    Point len_or_zero = {0, 0, 0}; //!< Box length (if PBC) or zero (if no PBC) in given direction
    Point len, len_half, len_inv;  //!< Cached box dimensions, their half-values, and reciprocal values.
    std::unique_ptr<GeometryImplementation> geometry = nullptr; //!< A concrete geometry implementation.
    Variant _type;                                              //!< Type of concrete geometry.
    std::string _name;                                          //!< Name of concrete geometry, e.g., for json.
    void makeGeometry(const Variant type = CUBOID); //!< Creates and assigns a concrete geometry implementation.
    void _setLength(const Point &l);

  public:
    const Variant &type = _type;     //!< Type of concrete geometry, read-only.
    const std::string &name = _name; //!< Name of concrete geometry, e.g., for json, read-only.
    double getVolume(int dim = 3) const override;
    Point setVolume(double, VolumeMethod = ISOTROPIC) override;
    Point getLength() const override; //!< A minimal containing cubic box.
    // setLength() needed only for Move::ReplayMove (stems from IO::XTCReader).
    void setLength(const Point &);                            //!< Sets the box dimensions.
    void boundary(Point &) const override;                    //!< Apply boundary conditions
    Point vdist(const Point &, const Point &) const override; //!< (Minimum) distance between two points
    double sqdist(const Point &, const Point &) const;        //!< (Minimum) squared distance between two points
    void randompos(Point &, Random &) const override;
    bool collision(const Point &) const override;
    void from_json(const json &) override;
    void to_json(json &j) const override;

    const BoundaryCondition &boundaryConditions() const; //!< Get info on boundary conditions

    static const std::map<std::string, Variant> names; //!< Geometry names.
    typedef std::pair<std::string, Variant> VariantName;

    static VariantName variantName(const std::string &name);

    static VariantName variantName(const json &j);

    Chameleon(const Variant type = CUBOID);
    Chameleon(const GeometryImplementation &geo, const Variant type);

    //! Copy everything, but clone the geometry.
    Chameleon(const Chameleon &geo);

    //! During the assignment copy everything, but clone the geometry.
    Chameleon &operator=(const Chameleon &geo);

    std::shared_ptr<GeometryImplementation> asSimpleGeometry();
};

inline void Chameleon::randompos(Point &m, Random &rand) const {
    assert(geometry);
    geometry->randompos(m, rand);
}

inline bool Chameleon::collision(const Point &a) const {
    assert(geometry);
    return geometry->collision(a);
}

inline void Chameleon::boundary(Point &a) const {
    const auto &boundary_conditions = geometry->boundary_conditions;
    if (boundary_conditions.coordinates == ORTHOGONAL) {
        if (boundary_conditions.direction.x() == PERIODIC) {
            if (std::fabs(a.x()) > len_half.x())
                a.x() -= len.x() * anint(a.x() * len_inv.x());
        }
        if (boundary_conditions.direction.y() == PERIODIC) {
            if (std::fabs(a.y()) > len_half.y())
                a.y() -= len.y() * anint(a.y() * len_inv.y());
        }
        if (boundary_conditions.direction.z() == PERIODIC) {
            if (std::fabs(a.z()) > len_half.z())
                a.z() -= len.z() * anint(a.z() * len_inv.z());
        }
    } else {
        geometry->boundary(a);
    }
}

inline Point Chameleon::vdist(const Point &a, const Point &b) const {
    Point distance;
    const auto &boundary_conditions = geometry->boundary_conditions;
    if (boundary_conditions.coordinates == ORTHOGONAL) {
        distance = a - b;
        if (boundary_conditions.direction.x() == PERIODIC) {
            if (distance.x() > len_half.x())
                distance.x() -= len.x();
            else if (distance.x() < -len_half.x())
                distance.x() += len.x();
        }
        if (boundary_conditions.direction.y() == PERIODIC) {
            if (distance.y() > len_half.y())
                distance.y() -= len.y();
            else if (distance.y() < -len_half.y())
                distance.y() += len.y();
        }
        if (boundary_conditions.direction.z() == PERIODIC) {
            if (distance.z() > len_half.z())
                distance.z() -= len.z();
            else if (distance.z() < -len_half.z())
                distance.z() += len.z();
        }
    } else {
        distance = geometry->vdist(a, b);
    }
    return distance;
}

inline double Chameleon::sqdist(const Point &a, const Point &b) const {
    if (geometry->boundary_conditions.coordinates == ORTHOGONAL) {
        if constexpr (true) {
            Point d((a - b).cwiseAbs());
            return (d - (d.array() > len_half.array()).cast<double>().matrix().cwiseProduct(len_or_zero)).squaredNorm();
        } else { // more readable alternative(?), nearly same speed
            Point d(a - b);
            for (int i = 0; i < 3; ++i) {
                d[i] = std::fabs(d[i]);
                d[i] = d[i] - len_or_zero[i] * static_cast<double>(d[i] > len_half[i]); // casting faster than branching
            }
            return d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
        }
    } else {
        return geometry->vdist(a, b).squaredNorm();
    }
}

void to_json(json &, const Chameleon &);
void from_json(const json &, Chameleon &);

enum class weight { MASS, CHARGE, GEOMETRIC };

/**
 * @brief Calculate Center of particle range using arbitrary weight functions applied to each particle
 * @param begin Begin particle iterator
 * @param end End parti cle iterator
 * @param apply_boundary Boundary function to apply PBC (default: no PBC)
 * @param weight_function Functor return weight for a given particle
 * @param shift Shift by this vector before calculating center, then add again. For PBC removal; default: 0,0,0
 * @return Center position
 * @throws if the sum of weights is zero, thereby hampering normalization
 */
template <typename iterator, typename weightFunc>
Point anyCenter(iterator begin, iterator end, BoundaryFunction apply_boundary, const weightFunc &weight_function,
                const Point &shift = {0, 0, 0}) {
    double weight_sum = 0.0;
    Point center(0.0, 0.0, 0.0);
    std::for_each(begin, end, [&](const auto &particle) {
        const auto weight = weight_function(particle);
        Point pos = particle.pos + shift; // translate
        apply_boundary(pos);
        center += weight * pos;
        weight_sum += weight;
    });
    if (weight_sum > pc::epsilon_dbl) {
        center = center / weight_sum - shift; // translate back
        apply_boundary(center);
        return center;
    } else {
        throw std::runtime_error("cannot calculate center with zero weights");
    }
} //!< Mass, charge, or geometric center of a collection of particles

/**
 * @brief Calculates mass center of range of particles
 * @param begin Begin particle iterator
 * @param end End particle iterator
 * @param apply_boundary Boundary function to apply PBC (default: no PBC)
 * @param shift Shift by this vector before calculating center, then add again. For PBC removal; default: 0,0,0
 * @return Mass center position
 * @throws if the sum of masses is zero, thereby hampering normalization
 */
template <typename iterator>
Point massCenter(
    iterator begin, iterator end, BoundaryFunction apply_boundary = [](Point &) {},
    const Point &shift = {0.0, 0.0, 0.0}) {
    auto particle_mass = [](const auto &particle) { return particle.traits().mw; };
    return anyCenter(begin, end, apply_boundary, particle_mass, shift);
}

/**
 * @param begin Begin iterator
 * @param end End iterator
 * @param displacement Displacement vector
 * @param apply_boundary Boundary function to apply PBC (default: no PBC)
 */
template <typename iterator>
void translate(
    iterator begin, iterator end, const Point &displacement, BoundaryFunction apply_boundary = [](Point &) {}) {
    std::for_each(begin, end, [&](auto &particle) {
        particle.pos += displacement;
        apply_boundary(particle.pos);
    });
}

/**
 * @brief Move the mass center of the particle range to origin (0,0,0) and wrap to PBC
 * @param begin Begin iterator
 * @param end End iterator
 * @param apply_boundary Boundary function to apply PBC (default: none)
 */
template <typename iterator>
void translateToOrigin(
    iterator begin, iterator end, BoundaryFunction apply_boundary = [](Point &) {}) {
    Point cm = massCenter(begin, end, apply_boundary);
    translate(begin, end, -cm, apply_boundary);
}

/**
 * @brief Rotate range of particles using a Quaternion
 * @param begin Begin iterator
 * @param end End iterator
 * @param quaternion Quaternion used for rotation
 * @param apply_boundary Boundary function to apply PBC (default: none)
 * @param shift This value is added before rotation to aid PBC remove (default: 0,0,0)
 * @todo Currently both quaternion and rotation matrix are passed, but one of them should be enough
 *
 * This will rotate both positions and internal coordinates in the particle (dipole moments, tensors etc.)
 */
template <typename iterator>
void rotate(
    iterator begin, iterator end, const Eigen::Quaterniond &quaternion,
    BoundaryFunction apply_boundary = [](Point &) {}, const Point &shift = Point(0.0, 0.0, 0.0)) {
    const auto rotation_matrix = quaternion.toRotationMatrix(); // rotation matrix
    std::for_each(begin, end, [&](auto &particle) {
        particle.rotate(quaternion, rotation_matrix); // rotate internal coordinates
        particle.pos += shift;
        apply_boundary(particle.pos);
        particle.pos = (quaternion * particle.pos) - shift; // rotate positions
        apply_boundary(particle.pos);
    });
} //!< Rotate particle pos and internal coordinates

/*
 * @brief Calculate mass center of cluster of particles in unbounded environment
 *
 * [More info](http://dx.doi.org/10.1080/2151237X.2008.10129266)
 */
template <class Tspace, class GroupIndex>
Point trigoCom(const Tspace &spc, const GroupIndex &groups, const std::vector<int> &dir = {0, 1, 2}) {
    assert(!dir.empty() && dir.size() <= 3);
    Point xhi(0, 0, 0), zeta(0, 0, 0), theta(0, 0, 0), com(0, 0, 0);
    for (auto k : dir) {
        double q = 2 * pc::pi / spc.geo.getLength()[k];
        size_t N = 0;
        for (auto i : groups)
            for (auto &particle : spc.groups[i]) {
                theta[k] = particle.pos[k] * q;
                zeta[k] += std::sin(theta[k]);
                xhi[k] += std::cos(theta[k]);
                N++;
            }
        theta[k] = std::atan2(-zeta[k] / N, -xhi[k] / N) + pc::pi;
        com[k] = spc.geo.getLength()[k] * theta[k] / (2 * pc::pi);
    }
    spc.geo.boundary(com); // is this really needed?
    return com;
}

/**
 * @brief Calculates a gyration tensor of a range of particles
 *
 * The gyration tensor is computed from the atomic position vectors with respect to the reference point
 * which is always a center of mass,
 * \f$ t_{i} = r_{i} - r_\mathrm{cm} \f$:
 * \f$ S = (1 / \sum_{i=1}^{N} m_{i}) \sum_{i=1}^{N} m_{i} t_{i} t_{i}^{T} \f$
 *
 * Before the calculation, the molecule is made whole to moving it to the center or the
 * simulation box (0,0,0), then apply the given boundary function.
 *
 * @tparam iterator Iterator to `Particle` range
 * @param begin Iterator to first particle
 * @param end Iterator to end
 * @param mass_center The mass center used as reference and to remove PBC
 * @param boundary Function to apply periodic boundary functions (default: none)
 * @return gyration tensor; or zero tensor if empty particle range
 * @throws If total mass is non-positive
 */
template <typename iterator>
Tensor gyration(
    iterator begin, iterator end, const Point &mass_center, const BoundaryFunction boundary = [](auto &) {}) {
    Tensor S = Tensor::Zero();
    double total_mass = 0.0;
    std::for_each(begin, end, [&](const auto &particle) {
        const auto mass = particle.traits().mw;
        Point r = particle.pos - mass_center; // get rid...
        boundary(r);                          // ...of PBC (if any)
        S += mass * r * r.transpose();
        total_mass += mass;
    });
    if (total_mass > 0.0) {
        return S / total_mass;
    } else {
        throw std::runtime_error("gyration tensor: total mass must be positive");
    }
}

/**
 * @brief Calculates a gyration tensor of a range of particles
 *
 * The gyration tensor is computed from the atomic position vectors with respect to the reference point
 * which is always a center of mass,
 * \f$ t_{i} = r_{i} - r_\mathrm{cm} \f$:
 * \f$ S = (1 / \sum_{i=1}^{N} m_{i}) \sum_{i=1}^{N} m_{i} t_{i} t_{i}^{T} \f$
 *
 * Before the calculation, the molecule is made whole to moving it to the center or the
 * simulation box (0,0,0), then apply the given boundary function.
 *
 * @param begin Iterator to first position
 * @param end Iterator to past last position
 * @param mass Iterator to first mass
 * @param mass_center The mass center used as reference and to remove PBC
 * @param boundary Function to apply periodic boundary functions (default: none)
 * @return gyration tensor; or zero tensor if empty particle range
 * @throws If total mass is non-positive
 */
template <typename position_iterator, typename mass_iterator>
Tensor gyration(
    position_iterator begin, position_iterator end, mass_iterator mass, const Point &mass_center,
    const BoundaryFunction boundary = [](auto &) {}) {
    Tensor S = Tensor::Zero();
    double total_mass = 0.0;
    std::for_each(begin, end, [&](const auto &position) {
        Point r = position - mass_center; // get rid...
        boundary(r);                      // ...of PBC (if any)
        S += (*mass) * r * r.transpose();
        total_mass += (*mass);
        std::advance(mass, 1);
    });
    if (total_mass > 0.0) {
        return S / total_mass;
    } else {
        throw std::runtime_error("gyration tensor: total mass must be positive");
    }
}

/**
 * @brief Shape descriptors derived from gyration tensor
 *
 * The class is prepared with operator overloads to work with `AverageObj`
 * for averaging over multiple tensors
 */
struct ShapeDescriptors {
    double gyration_radius_squared = 0.0;
    double asphericity = 0.0;
    double acylindricity = 0.0;
    double relative_shape_anisotropy = 0.0; //!< relative shape anisotropy, kappa^2 (0=rod, 1=spherical)
    ShapeDescriptors() = default;
    ShapeDescriptors(const Tensor &gyration_tensor);             //!< Construct using an initial gyration tensor
    ShapeDescriptors &operator+=(const ShapeDescriptors &other); //!< Add another gyration tensor; req. for averaging
    ShapeDescriptors operator*(const double scale) const;        //!< Scale data; req. for averaging
};

void to_json(json &j, const ShapeDescriptors &shape); //!< Store Shape as json object

/**
 * @brief Calculates an inertia tensor of a molecular group
 *
 * The inertia tensor is computed from the atomic position vectors with respect to a reference point,
 * \f$ t_{i} = r_{i} - r_\mathrm{origin} \f$:
 * \f$ S = \sum_{i=1}^{N} m_{i} ( t_{i} \cdot t_{i} I  - t_{i} t_{i}^{T} ) \f$
 *
 * @param origin a reference point
 * @return inertia tensor (a zero tensor for an empty group)
 */
template <typename iter>
Tensor inertia(
    iter begin, iter end, const Point origin = {0, 0, 0}, const BoundaryFunction boundary = [](const Point &) {}) {
    Tensor I = Tensor::Zero();
    for (auto it = begin; it != end; ++it) {
        Point t = it->pos - origin;
        boundary(t);
        I += atoms.at(it->id).mw * (t.squaredNorm() * Tensor::Identity() - t * t.transpose());
    }
    return I;
}

/**
 * @brief Root-mean-square deviation of two data sets represented by iterators
 *
 * A binary function must be given that returns the difference between data points
 * in the two sets, for example `[](int a, int b){return a-b;}`.
 */
template <typename InputIt1, typename InputIt2, typename BinaryOperation>
double rootMeanSquareDeviation(InputIt1 begin, InputIt1 end, InputIt2 d_begin, BinaryOperation diff_squared_func) {
    assert(std::distance(begin, end) > 0);
    double sq_sum = 0;
    for (InputIt1 i = begin; i != end; ++i) {
        sq_sum += diff_squared_func(*i, *d_begin);
        ++d_begin;
    }
    return std::sqrt(sq_sum / std::distance(begin, end));
}

/**
 * @brief Scale particles to the surface of a sphere
 * @param particles Vector of particles
 *
 * The sphere radius is taken as the average radial distance
 * of all particles with respect to the mass center.
 * The _first_ particle of the given particles is excluded
 * from the COM calculation and re-positioned at the center
 * of the sphere. Therefore, make sure to add a dummy particle
 * to the beginning of the particle vector; its initial positions
 * are ignored and will be overwritten.
 * Similar to routine described in doi:10.1021/jp010360o
 */
ParticleVector mapParticlesOnSphere(const ParticleVector &);

/**
 * @brief Convert particles in hexagonal prism to space-filled cuboid
 * @param hexagon Input hexagonal prism
 * @param particles Particles in hexagonal prism
 * @return Cuboid and particle vector w. positions in cuboidal space
 *
 * The generated Cuboid has twice the volume of the hexagonal prism with
 * side-lengths [2 * inner_radius, 3 * outer_radius, height] and also twice
 * the number of particles
 */
template <typename Particles>
std::pair<Cuboid, ParticleVector> HexagonalPrismToCuboid(const HexagonalPrism &hexagon, const Particles &particles) {
    Cuboid cuboid({2.0 * hexagon.innerRadius(), 3.0 * hexagon.outerRadius(), hexagon.height()});
    ParticleVector cuboid_particles;
    cuboid_particles.reserve(2 * std::distance(particles.begin(), particles.end()));
    std::copy(particles.begin(), particles.end(), std::back_inserter(cuboid_particles)); // add central hexagon

    std::transform(particles.begin(), particles.end(), std::back_inserter(cuboid_particles), [&](auto particle) {
        particle.pos.x() += hexagon.innerRadius() * (particle.pos.x() > 0.0 ? -1.0 : 1.0);
        particle.pos.y() += hexagon.outerRadius() * (particle.pos.y() > 0.0 ? -1.5 : 1.5);
        assert(cuboid.collision(particle.pos) == false);
        return particle;
    }); // add the four corners; i.e. one extra, split hexagon
    assert(std::fabs(cuboid.getVolume() - 2.0 * hexagon.getVolume()) <= pc::epsilon_dbl);
    return {cuboid, cuboid_particles};
}

} // namespace Geometry
} // namespace Faunus
