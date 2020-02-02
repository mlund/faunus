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

struct Random;

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
typedef std::function<void(Point &)> BoundaryFunction;
typedef std::function<Point(const Point &, const Point &)> DistanceFunction;

//! Geometry variant used for Chameleon.
enum Variant { CUBOID = 0, SPHERE, CYLINDER, SLIT, HEXAGONAL, OCTAHEDRON, HYPERSPHERE2D };

//! Various methods of volume scaling, @see GeometryBase::setVolume.
enum VolumeMethod { ISOTROPIC, ISOCHORIC, XY, Z };

enum Coordinates { ORTHOGONAL, ORTHOHEXAGONAL, TRUNC_OCTAHEDRAL, NON3D };
enum Boundary { FIXED, PERIODIC };

/**
 * @brief A structure containing a type of boundary condition in each direction.
 *
 * A stub. It can be extended to fully json-configurable boundary conditions.
 */
struct BoundaryCondition {
    typedef Eigen::Matrix<Boundary, 3, 1> BoundaryXYZ;
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
    } //!< returns lambda to boundary()

    inline DistanceFunction getDistanceFunc() const {
        return [this](const Point &i, const Point &j) { return vdist(i, j); };
    } //!< returns lambda to vdist()

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

    //!< A unique pointer to a copy of self. To be used in copy constructors.
    virtual std::unique_ptr<GeometryImplementation> clone() const = 0;

    //!< Cereal serialisation
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
    Cuboid(const Point &p);
    Cuboid(double x, double y, double z);
    Cuboid(double x = 0.0);

    std::unique_ptr<GeometryImplementation> clone() const override {
        return std::make_unique<Cuboid>(*this);
    }; //!< A unique pointer to a copy of self.

    //!< Cereal serialisation
    template <class Archive> void serialize(Archive &archive) {
        archive(cereal::base_class<GeometryImplementation>(this), box, box_half, box_inv);
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
    void boundary(Point &a) const override;
    bool collision(const Point &a) const override;
    void randompos(Point &m, Random &rand) const override;
    void from_json(const json &j) override;
    void to_json(json &j) const override;
    Sphere(double radius = 0.0);

    std::unique_ptr<GeometryImplementation> clone() const override {
        return std::make_unique<Sphere>(*this);
    }; //!< A unique pointer to a copy of self.

    //!< Cereal serialisation
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

    //!< Cereal serialisation
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

    //!< Cereal serialisation
    template <class Archive> void serialize(Archive &archive) {
        archive(cereal::base_class<GeometryImplementation>(this), rhombic2cartesian, cartesian2rhombic, box);
    }
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
    void randompos(Point &m, Random &rand) const override;
    void from_json(const json &j) override;
    void to_json(json &j) const override;
    TruncatedOctahedron(double side = 0.0);

    std::unique_ptr<GeometryImplementation> clone() const override {
        return std::make_unique<TruncatedOctahedron>(*this);
    }; //!< A unique pointer to a copy of self.

    //!< Cereal serialisation
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
    Point len, len_half, len_inv; //!< Cached box dimensions, their half-values, and reciprocal values.
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
    // setLength() needed only for IO::FormatXTC::loadnextframe().
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
            if (std::fabs(a.z()) > len_half.y())
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

/**
 * More readable alternative, nearly same performance:
 *
 * ~~~ cpp
 * Point d(a - b);
 * for (size_t i = 0; i < 3; i++) {
 *     d[i] = std::fabs(d[i]);
 *     d[i] = d[i] - len_or_zero[i] * (d[i] > len_half[i]); // casting faster than branching
 * }
 * return d.squaredNorm();
 * ~~~
 */
inline double Chameleon::sqdist(const Point &a, const Point &b) const {
    if (geometry->boundary_conditions.coordinates == ORTHOGONAL) {
        Point d((a - b).cwiseAbs());
        return (d - (d.array() > len_half.array()).cast<double>().matrix().cwiseProduct(len_or_zero)).squaredNorm();
    } else
        return geometry->vdist(a, b).squaredNorm();
}

void to_json(json &, const Chameleon &);
void from_json(const json &, Chameleon &);

/*
   void unwrap( Point &a, const Point &ref ) const {
   a = vdist(a, ref) + ref;
   } //!< Remove PBC with respect to a reference point
 */
enum class weight { MASS, CHARGE, GEOMETRIC };

template <typename Titer, typename Tparticle = typename Titer::value_type, typename weightFunc>
Point anyCenter(Titer begin, Titer end, BoundaryFunction boundary, const weightFunc &weight,
                const Point &shift = {0, 0, 0}) {
    double sum = 0;
    Point c(0, 0, 0), t;
    for (auto &i = begin; i != end; ++i) {
        t = i->pos + shift; // translate to origo
        boundary(t);
        double w = weight(*i);
        c += w * t;
        sum += w;
    }
    if (sum <= pc::epsilon_dbl)
        return {0, 0, 0};
    else
        c = c / sum - shift;
    boundary(c);
    return c;
} //!< Mass, charge, or geometric center of a collection of particles

template <typename Titer, typename Tparticle = typename Titer::value_type>
Point massCenter(
    Titer begin, Titer end, BoundaryFunction boundary = [](Point &) {}, const Point &shift = {0, 0, 0}) {
    return anyCenter(
        begin, end, boundary, [](const Tparticle &p) { return atoms.at(p.id).mw; }, shift);
} // Mass center

template <class Titer = typename std::vector<T>::iterator>
void translate(
    Titer begin, Titer end, const Point &d, BoundaryFunction boundary = [](Point &) {}) {
    for (auto i = begin; i != end; ++i) {
        i->pos += d;
        boundary(i->pos);
    }
} //!< Vector displacement of a range of particles

template <typename Titer>
void cm2origo(
    Titer begin, Titer end, BoundaryFunction boundary = [](Point &) {}) {
    Point cm = massCenter(begin, end, boundary);
    translate(begin, end, -cm, boundary);
} //!< Translate to that mass center is in (0,0,0)

template <typename Titer>
void rotate(
    Titer begin, Titer end, const Eigen::Quaterniond &q, BoundaryFunction boundary = [](Point &) {},
    const Point &shift = Point(0, 0, 0)) {
    auto m = q.toRotationMatrix(); // rotation matrix
    for (auto i = begin; i != end; ++i) {
        i->rotate(q, m); // rotate internal coordinates
        i->pos += shift;
        boundary(i->pos);
        i->pos = q * i->pos;
        boundary(i->pos);
        i->pos -= shift;
        boundary(i->pos);
    }
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
 * @brief Calculates a gyration tensor of a molecular group
 *
 * The gyration tensor is computed from the atomic position vectors with respect to the reference point
 * which is always a center of mass,
 * \f$ t_{i} = r_{i} - r_\mathrm{cm} \f$:
 * \f$ S = (1 / \sum_{i=1}^{N} m_{i}) \sum_{i=1}^{N} m_{i} t_{i} t_{i}^{T} \f$
 *
 * @param origin center of mass of the molecular group
 * @return gyration tensor (a zero tensor for an empty group)
 */
template <typename iter>
Tensor gyration(iter begin, iter end, const Point origin = {0,0,0},
        const BoundaryFunction boundary = [](const Point &) {}) {
    Tensor S = Tensor::Zero();
    double mw_sum = 0;
    for (auto it = begin; it != end; ++it) {
        const auto mw = atoms.at(it->id).mw;
        Point t = it->pos - origin;
        boundary(t);
        mw_sum += mw;
        S += mw * t * t.transpose();
    }
    if (mw_sum != 0) {
        S /= mw_sum;
    } else {
        assert(S == Tensor::Zero()); // otherwise we have negative atom weights
    }
    return S;
}

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
Tensor inertia(iter begin, iter end, const Point origin = {0,0,0},
        const BoundaryFunction boundary = [](const Point &) {}) {
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

} // namespace Geometry
} // namespace Faunus
