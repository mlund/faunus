#pragma once
#include "core.h"
#include "atomdata.h"

namespace Faunus {

/** @brief Base class for particle properties */
struct ParticlePropertyBase {
    virtual void to_json(json &j) const = 0;   //!< Convert to JSON object
    virtual void from_json(const json &j) = 0; //!< Convert from JSON object
    void rotate(const Eigen::Quaterniond &q, const Eigen::Matrix3d &);
    inline virtual ~ParticlePropertyBase(){};
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

struct PositionAndID : public ParticlePropertyBase {
    int id = -1;           //!< Particle id/type
    Point pos = {0, 0, 0}; //!< Particle position vector
    const AtomData &traits();
    void to_json(json &j) const;
    void from_json(const json &j);
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
}; //!< Position and ID property

struct Radius : public ParticlePropertyBase {
    double radius = 0; //!< Particle radius
    void to_json(json &j) const override;
    void from_json(const json &j) override;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
}; //!< Radius property

struct Charge : public ParticlePropertyBase {
    double charge = 0; //!< Particle radius
    void to_json(json &j) const override;
    void from_json(const json &j) override;
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
    Point mu = {1, 0, 0};                                              //!< dipole moment unit vector
    double mulen = 0;                                                  //!< dipole moment scalar
    void rotate(const Eigen::Quaterniond &q, const Eigen::Matrix3d &); //!< Rotate dipole moment
    void to_json(json &j) const override;
    void from_json(const json &j) override;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

struct Polarizable : public ParticlePropertyBase {
    Point mui = {1, 0, 0};                                              //!< induced dipole moment unit vector
    double muilen = 0;                                                  //!< induced dipole moment scalar
    Tensor alpha;                                                       //!< polarizability tensor
    void rotate(const Eigen::Quaterniond &q, const Eigen::Matrix3d &m); //!< Rotate polarizability tensor
    void to_json(json &j) const override;
    void from_json(const json &j) override;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

struct Quadrupole : public ParticlePropertyBase {
    Tensor Q;                                                           //!< Quadrupole
    void rotate(const Eigen::Quaterniond &q, const Eigen::Matrix3d &m); //!< Rotate quadrupole moment
    void to_json(json &j) const override;
    void from_json(const json &j) override;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
}; // Quadrupole property

struct Cigar : public ParticlePropertyBase {
    Point scdir = {1, 0, 0};                                           //!< Sphero-cylinder direction unit vector
    double sclen = 0;                                                  //!< Length
    void rotate(const Eigen::Quaterniond &q, const Eigen::Matrix3d &); //!< Rotate sphero-cylinder
    void to_json(json &j) const override;
    void from_json(const json &j) override {
        scdir = j.value("scdir", scdir);
        sclen = j.value("sclen", sclen);
    }
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
}; //!< Sphero-cylinder properties

/**
 * @brief Particle
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

  public:
    ParticleTemplate() : Properties()... {}

    ParticleTemplate(const AtomData &a) : Properties()... { *this = json(a).front(); }

    void rotate(const Eigen::Quaterniond &q, const Eigen::Matrix3d &m) {
        _rotate<Properties...>(q, m, dynamic_cast<Properties &>(*this)...);
    } //!< Rotate all internal coordinates if needed

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

template <typename... Properties> void to_json(json &j, const ParticleTemplate<Properties...> &a) {
    to_json<Properties...>(j, Properties(a)...);
}

template <typename... Properties> void from_json(const json &j, ParticleTemplate<Properties...> &a) {
    from_json<Properties...>(j, dynamic_cast<Properties &>(a)...);
}

class Particle {
  public:
    int id = -1;           //!< Particle id/type
    double charge = 0;     //!< Particle charge
    Point pos = {0, 0, 0}; //!< Particle position vector
    std::shared_ptr<ParticleTemplate<Dipole, Quadrupole, Cigar>> shape = nullptr;
    const AtomData &traits();
    Particle() = default;
    Particle(const AtomData &a);
    void rotate(const Eigen::Quaterniond &q, const Eigen::Matrix3d &m);
};

void from_json(const json &j, Particle &p);

void to_json(json &j, const Particle &p);

// using Particle = ParticleTemplate<PositionAndID, Charge>;
using ParticleAllProperties = ParticleTemplate<PositionAndID, Radius, Dipole, Charge, Quadrupole, Cigar>;

#ifdef DOCTEST_LIBRARY_INCLUDED
// convert test to use `Particle::shape`
TEST_CASE("[Faunus] Particle") {
    using doctest::Approx;
    ParticleAllProperties p1, p2;
    p1.id = 100;
    p1.pos = {1, 2, 3};
    p1.charge = -0.8;
    p1.radius = 7.1;
    p1.mu = {0, 0, 1};
    p1.mulen = 2.8;
    p1.scdir = {-0.1, 0.3, 1.9};
    p1.sclen = 0.5;
    p1.Q = Tensor(1, 2, 3, 4, 5, 6);

    p2 = json(p1); // p1 --> json --> p2

    CHECK(json(p1) == json(p2)); // p1 --> json == json <-- p2 ?

    CHECK(p2.id == 100);
    CHECK(p2.pos == Point(1, 2, 3));
    CHECK(p2.charge == -0.8);
    CHECK(p2.radius == 7.1);
    CHECK(p2.mu == Point(0, 0, 1));
    CHECK(p2.mulen == 2.8);
    CHECK(p2.scdir == Point(-0.1, 0.3, 1.9));
    CHECK(p2.sclen == 0.5);
    CHECK(p2.Q == Tensor(1, 2, 3, 4, 5, 6));

    // check of all properties are rotated
    QuaternionRotate qrot(pc::pi / 2, {0, 1, 0});
    p1.mu = p1.scdir = {1, 0, 0};
    p1.rotate(qrot.first, qrot.second);

    CHECK(p1.mu.x() == Approx(0));
    CHECK(p1.mu.z() == Approx(-1));
    CHECK(p1.scdir.x() == Approx(0));
    CHECK(p1.scdir.z() == Approx(-1));

    CHECK(p1.Q(0, 0) == Approx(6));
    CHECK(p1.Q(0, 1) == Approx(5));
    CHECK(p1.Q(0, 2) == Approx(-3));
    CHECK(p1.Q(1, 1) == Approx(4));
    CHECK(p1.Q(1, 2) == Approx(-2));
    CHECK(p1.Q(2, 2) == Approx(1));
}
#endif

} // namespace Faunus
