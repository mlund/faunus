#pragma once
#include "core.h"
#include "atomdata.h"
#include "tensor.h"
#include <cereal/types/memory.hpp>

namespace Eigen {
using Matrix3d = Matrix<double, 3, 3>;
using Quaterniond = Quaternion<double>;
} // namespace Eigen

namespace Faunus {

/**
 * @brief Base class for particle properties
 * @todo Is this really needed?
 */
struct ParticlePropertyBase {
    virtual void to_json(json &j) const = 0;                           //!< Convert to JSON object
    virtual void from_json(const json &j) = 0;                         //!< Convert from JSON object
    void rotate(const Eigen::Quaterniond &q, const Eigen::Matrix3d &); //!< Internal rotation
    virtual ~ParticlePropertyBase() = default;
    template <class Archive> void serialize(Archive &) {} //!< Cereal serialisation
};

template <typename... Ts>
auto to_json(json &) -> typename std::enable_if<sizeof...(Ts) == 0>::type {} //!< Particle to JSON
template <typename T, typename... Ts> void to_json(json &j, const ParticlePropertyBase &a, const Ts &... rest) {
    a.to_json(j);
    to_json<Ts...>(j, rest...);
} //!< Particle to JSON

// JSON --> Particle
template <typename... Ts> auto from_json(const json &) -> typename std::enable_if<sizeof...(Ts) == 0>::type {}
template <typename T, typename... Ts> void from_json(const json &j, ParticlePropertyBase &a, Ts &... rest) {
    a.from_json(j);
    from_json<Ts...>(j, rest...);
}

struct Radius : public ParticlePropertyBase {
    double radius = 0.0; //!< Particle radius
    void to_json(json &j) const override;
    void from_json(const json &j) override;
    template <class Archive> void serialize(Archive &archive) { archive(radius); }
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
}; //!< Radius property

struct Charge : public ParticlePropertyBase {
    double charge = 0.0; //!< Particle radius
    void to_json(json &j) const override;
    void from_json(const json &j) override;
    template <class Archive> void serialize(Archive &archive) { archive(charge); }
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
}; //!< Charge (monopole) property

/** @brief Dipole properties
 *
 * Json (de)serialization:
 *
 * ```{.cpp}
 *     Dipole d = R"( "mu":[0,0,1], "mulen":10 )"_json
 * ```
 */
struct Dipole : public ParticlePropertyBase {
    Point mu = {0.0, 0.0, 0.0};                                        //!< dipole moment unit vector
    double mulen = 0.0;                                                //!< dipole moment scalar
    void rotate(const Eigen::Quaterniond &q, const Eigen::Matrix3d &); //!< Rotate dipole moment
    void to_json(json &j) const override;
    void from_json(const json &j) override;
    template <class Archive> void serialize(Archive &archive) { archive(mu, mulen); }
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

struct Polarizable : public ParticlePropertyBase {
    Point mui = {1.0, 0.0, 0.0};                                        //!< induced dipole moment unit vector
    double muilen = 0.0;                                                //!< induced dipole moment scalar
    Tensor alpha;                                                       //!< polarizability tensor
    void rotate(const Eigen::Quaterniond &q, const Eigen::Matrix3d &m); //!< Rotate polarizability tensor
    void to_json(json &j) const override;
    void from_json(const json &j) override;
    template <class Archive> void serialize(Archive &archive) { archive(mui, muilen, alpha); }
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

struct Quadrupole : public ParticlePropertyBase {
    Tensor Q;                                                           //!< Quadrupole
    void rotate(const Eigen::Quaterniond &q, const Eigen::Matrix3d &m); //!< Rotate quadrupole moment
    void to_json(json &j) const override;
    void from_json(const json &j) override;
    template <class Archive> void serialize(Archive &archive) { archive(Q); }
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
}; // Quadrupole property

struct Cigar : public ParticlePropertyBase {
    Point scdir = {1, 0, 0};                                           //!< Sphero-cylinder direction unit vector
    double sclen = 0;                                                  //!< Length
    void rotate(const Eigen::Quaterniond &q, const Eigen::Matrix3d &); //!< Rotate sphero-cylinder
    void to_json(json &j) const override;
    void from_json(const json &j) override;
    template <class Archive> void serialize(Archive &archive) { archive(scdir, sclen); }
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
}; //!< Sphero-cylinder properties

/**
 * @brief Particle template
 *
 * This is the Particle class used to store information about
 * particles. In addition to charge, positionn and ID, the particle
 * can store auxiliary information via the `shape` pointer. This
 * can is typically anisotropic properties which in this way
 * is stored in an alternative memory location.
 * */
template <typename... Properties> class ParticleTemplate : public Properties... {
  private:
    // rotate internal coordinates
    template <typename... Ts>
    auto _rotate(const Eigen::Quaterniond &, const Eigen::Matrix3d &) ->
        typename std::enable_if<sizeof...(Ts) == 0>::type {}
    template <typename T, typename... Ts>
    void _rotate(const Eigen::Quaterniond &q, const Eigen::Matrix3d &m, T &a, Ts &... rest) {
        a.rotate(q, m);
        _rotate<Ts...>(q, m, rest...);
    }

    // Cereal serialisation

    template <typename... Ts, class Archive>
    auto _serialize(Archive &) -> typename std::enable_if<sizeof...(Ts) == 0>::type {}

    template <typename T, typename... Ts, class Archive> void _serialize(Archive &archive, T &a, Ts &... rest) {
        a.serialize(archive);
        _serialize<Ts...>(archive, rest...);
    }

  public:
    ParticleTemplate() : Properties()... {};

    explicit ParticleTemplate(const AtomData &a) : Properties()... { *this = json(a).front(); }

    void rotate(const Eigen::Quaterniond &q, const Eigen::Matrix3d &m) {
        _rotate<Properties...>(q, m, dynamic_cast<Properties &>(*this)...);
    } //!< Rotate all internal coordinates if needed

    template <class Archive> void serialize(Archive &archive) {
        _serialize<Properties...>(archive, dynamic_cast<Properties &>(*this)...);
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

template <typename... Properties> void to_json(json &j, const ParticleTemplate<Properties...> &a) {
    to_json<Properties...>(j, Properties(a)...);
}

template <typename... Properties> void from_json(const json &j, ParticleTemplate<Properties...> &a) {
    from_json<Properties...>(j, dynamic_cast<Properties &>(a)...);
}

/**
 * @brief Particle class for storing positions, id, and other properties
 *
 * Particles carry `id`, `pos`, `charge` by default but can have additional
 * or _extended_ data stored using a different memory model. When serializing
 * from a json object, extended properties are automatically detected and
 * memory is automatically allocated
 *
 * @todo: memory model for extended properties not optimal
 */
class Particle {
  public:
    using ParticleExtension = ParticleTemplate<Dipole, Quadrupole, Cigar>;
    std::shared_ptr<ParticleExtension> ext = nullptr; //!< Point to extended properties
    int id = -1;                                      //!< Particle id/type
    double charge = 0.0;                              //!< Particle charge
    Point pos = {0.0, 0.0, 0.0};                      //!< Particle position vector

    Particle() = default;
    Particle(const AtomData &a, const Point &pos);
    Particle(const AtomData &a);           //!< construct from AtomData
    Particle(const Particle &);            //!< copy constructor
    Particle &operator=(const Particle &); //!< assignment operator
    const AtomData &traits() const;        //!< get properties from AtomData
    void rotate(const Eigen::Quaterniond &quaternion, const Eigen::Matrix3d &rotation_matrix); //!< internal rotation
    bool hasExtension() const;            //!< check if particle has extensions (dipole etc.)
    ParticleExtension &createExtension(); //!< Create extension
    inline ParticleExtension &getExt() { return ext ? *ext : createExtension(); } //!< get/create extension
    inline const ParticleExtension &getExt() const {
        assert(ext);
        return *ext;
    } //!< Get extended particle properties;

    /**
     * @brief Cereal serialisation
     * @param archive Archive to serialize to/from
     * @warning Still under construction
     */
    template <class Archive> void serialize(Archive &archive) {
        archive(ext, id, charge, pos);
        // if (ext != nullptr)
        //    ext->serialize(archive);
    } //!<
};

//! Storage type for collections of particles
using ParticleVector = std::vector<Particle>;

void from_json(const json &, Particle &);
void to_json(json &, const Particle &);

} // namespace Faunus
