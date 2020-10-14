#pragma once
#include <doctest/doctest.h>
#include "units.h"
#include "aux/pow_function.h"

namespace Faunus {
/**
 * @brief Returns ion-dipole interaction.
 * @param QBxMuA Product of ion B:s charge and dipole A:s scalar
 * @param muA Unit dipole moment vector of particel A
 * @param QAxMuB Product of ion A:s charge and dipole B:s scalar
 * @param muB Unit dipole moment vector of particel B
 * @param r Direction \f$ r_A - r_B \f$
 */
template <class Tvec> double q2mu(double QBxMuA, const Tvec &muA, double QAxMuB, const Tvec &muB, const Tvec &r) {
    double r2i = 1 / r.squaredNorm(); // B = sol(dip), A = ch(charge)
    double r1i = sqrt(r2i);
    double r3i = r1i * r2i;
    double W1 = QBxMuA * muA.dot(r) * r3i;
    double W2 = QAxMuB * muB.dot(-r) * r3i;
    return (W1 + W2);
}

/**
 * @brief Returns dipole-dipole interaction
 *
 * @param muA Unit dipole moment vector of particle A
 * @param muB Unit dipole moment vector of particle B
 * @param muAxmuB Product of dipole scalars
 * @param r Direction \f$ r_A - r_B \f$
 */
template <class Tvec>
double mu2mu(const Tvec &muA, const Tvec &muB, double muAxmuB, const Tvec &r, double a = 1.0, double b = 0.0) {
#ifdef FAU_APPROXMATH
    double r1i = invsqrtQuake(r.squaredNorm());
    double r2i = r1i * r1i;
#else
    double r2i = 1 / r.squaredNorm();
    double r1i = sqrt(r2i);
#endif
    double dot = muA.dot(muB);
    double T = (3 * muA.dot(r) * muB.dot(r) * r2i - dot) * a + dot * b;
    return -muAxmuB * T * r1i * r2i;
}

/**
 * @brief Returns ion-quadrupole interaction
 * @param qA Charge of particle A
 * @param quadB Quadrupole moment of particle B
 * @param qB Charge of particle B
 * @param quadA Quadrupole moment of particle A
 * @param r Direction \f$ r_A - r_B \f$
 */
template <class Tvec, class Tmat>
double q2quad(double qA, const Tmat &quadB, double qB, const Tmat &quadA, const Tvec &r) {
    double r2i = 1 / r.squaredNorm();
    double r1i = sqrt(r2i);
    double r3i = r1i * r2i;
    double r5i = r3i * r2i;
    double WAB = r.transpose() * quadB * r;
    WAB = 3 * WAB * r5i - quadB.trace() * r3i;
    double WBA = r.transpose() * quadA * r;
    WBA = 3 * WBA * r5i - quadA.trace() * r3i;
    return (qA * WAB + qB * WBA);
}

namespace Potential {

/**
 * @brief Returns the factorial of 'n'. Note that 'n' must be positive semidefinite.
 * @note Calculated at compile time and thus have no run-time overhead.
 */
constexpr unsigned int factorial(unsigned int n) { return n <= 1 ? 1 : n * factorial(n - 1); }

TEST_CASE("[Faunus] Factorial") {
    CHECK(factorial(0) == 1);
    CHECK(factorial(1) == 1);
    CHECK(factorial(2) == 2);
    CHECK(factorial(3) == 6);
    CHECK(factorial(10) == 3628800);
}

/**
 * @brief Help-function for `sfQpotential` using order 3.
*/
inline void dipoleDipoleQ2Help_3(double &a3, double &b3, double q) {
    double q2 = q*q;
    constexpr double one_third = 1.0 / 3.0, two_third = 2.0 / 3.0, eight_third = 8.0 / 3.0;
    a3 = (((-5.0 * q2 + eight_third * q) * q + q) * q + one_third) * q2 + 1.0;
    b3 = -two_third * ((15.0 * q2 - 10.0 * q - 5.0) * q2 + 1.0) * q2;
}

/**
 * @brief Help-function for `sfQpotential` using order 4.
*/
inline void dipoleDipoleQ2Help_4(double &a4, double &b4, double q) {
    double q2 = q*q;
    double q3 = q2*q;
    a4 = ( ( (10.0 * q2 - 9.0 * q - 8.0 )*q3 + 10.0 )*q3 - 2.0 )*q - 1.0;
    b4 = ( ( ( 21.0 * q2 - 16.0 * q - 35.0 / 3.0 )*q3 + 16.0 / 3.0 )*q3 + 1.0 / 3.0 )*q2 + 1.0;
}

/**
* @brief Help-function for `sfQ0potential`.
*/
inline double dipoleDipoleQ2Help(double q, int l=0, int P=300, bool a=true) {
    if(q >= 1.0 - (1.0/2400.0))
        return 0.0;
    if(q <= (1.0/2400.0) && a)
        return 1.0;
    if(q <= (1.0/2400.0) && !a)
        return 0.0;

    double qP = 1.0; // Will end as q-Pochhammer Symbol, (q^l;q)_P
    double fac = Faunus::powi(q, l);
    double sum1 = 0.0, sum2 = 0.0, sum4 = 0.0;

    for( int n = 1; n <= P; n++) {
        fac *= q; // q^(l+n)
        qP *= (1.0 - fac);

        double tmp0 = (l+n)*fac/(1.0 - fac);
        sum2 += tmp0;
        sum1 += tmp0*(fac + l + n - 1.0)/(1.0 - fac);
        sum4 += tmp0*(4.0*fac + l + n - 4.0)/(1.0 - fac);
    }
    sum2 = sum2*sum2;
    double ac = qP/3.0*(3.0 - sum4 + sum2);
    double bc = qP/3.0*(sum2 - sum1);

    if(a)
        return ac;
    return bc;
}

/**
 * @brief Help-function for the q-potential in class CoulombGalore
 *
 * More information here: http://mathworld.wolfram.com/q-PochhammerSymbol.html
 * P = 300 gives an error of about 10^-17 for k < 4
 */
inline double qPochhammerSymbol(double q, int k = 1, int P = 300) {
    double value = 1.0;
    double temp = Faunus::powi(q, k);
    for (int i = 0; i < P; i++) {
        value *= (1.0 - temp);
        temp *= q;
    }
    return value;
}

TEST_CASE("[Faunus] qPochhammerSymbol") {
    double q = 0.5;
    CHECK(qPochhammerSymbol(q, 0, 0) == 1);
    CHECK(qPochhammerSymbol(0, 0, 1) == 0);
    CHECK(qPochhammerSymbol(1, 0, 1) == 0);
    CHECK(qPochhammerSymbol(1, 1, 2) == 0);
    // add tests...
}

} // namespace Potential

/**
 * @brief Returns the total charge for a set of particles
 * @param begin First particle
 * @param end Last particle
 */
template <class Titer> double monopoleMoment(Titer begin, Titer end) {
    double z = 0;
    for (auto it = begin; it != end; ++it)
        z += it->charge;
    return z;
} //!< Calculates dipole moment vector for a set of particles

/**
 * @brief Returns the total dipole-moment for a set of particles
 * @param begin First particle
 * @param end Last particle
 * @param boundary Function to use for boundary
 * @param origin Origin, default (0,0,0)
 * @param cutoff Cut-off for included particles with regard to origin, default value is infinite
 */
template <class Titer, class BoundaryFunction>
Point dipoleMoment(
    Titer begin, Titer end, BoundaryFunction boundary = [](const Point &) {}, const Point origin = {0, 0, 0},
    double cutoff = pc::infty) {
    Point mu(0, 0, 0);
    std::for_each(begin, end, [&](const auto &particle) {
        Point r = particle.pos - origin;
        boundary(r);
        if (r.squaredNorm() < cutoff * cutoff) {
            mu += r * particle.charge;
        }
    });
    return mu;
} //!< Calculates dipole moment vector

TEST_CASE("[Faunus] dipoleMoment") {
    using doctest::Approx;
    ParticleVector p(2);
    p[0].pos = {10, 20, 30};
    p[1].pos = {-10, 0, -30};
    p[0].charge = -0.5;
    p[1].charge = 0.5;

    SUBCASE("Neutral molecule") {
        auto mu = dipoleMoment(p.begin(), p.end(), [](auto &) {}, {2, 3, 4}); // some origin
        CHECK(mu.squaredNorm() == Approx(10 * 10 + 10 * 10 + 30 * 30));
        mu = dipoleMoment(p.begin(), p.end(), [](auto &) {}, {20, 30, 40}); // another origin
        CHECK(mu.squaredNorm() == Approx(10 * 10 + 10 * 10 + 30 * 30));
    }

    SUBCASE("Charged molecule") {
        p[0].charge *= -1.0; // give molecule a net charge
        auto mu = dipoleMoment(p.begin(), p.end(), [](auto &) {}, {2, 3, 4});
        CHECK(mu.x() == Approx(-2));
        CHECK(mu.y() == Approx(7));
        CHECK(mu.z() == Approx(-4));
    }
}

/**
 * @brief Returns the total quadrupole-moment for a set of particles, note with trace!
 * @param begin First particle
 * @param end Last particle
 * @param boundary Function to use for boundary
 * @param origin Origin for quadrupole-moment, default (0,0,0)
 * @param cutoff Cut-off for included particles with regard to origin, default value is infinite
 */
template <class Titer, class BoundaryFunction>
Tensor quadrupoleMoment(Titer begin, Titer end, BoundaryFunction boundary = [](const Point &) {}, Point origin={0,0,0},
                        double cutoff = pc::infty) {
    Tensor theta;
    theta.setZero();
    for (auto it = begin; it != end; ++it) {
        Point t = it->pos - origin;
        boundary(t);
        if (t.squaredNorm() < cutoff * cutoff) {
            theta += t * t.transpose() * it->charge;
        }
    }
    return 0.5 * theta;
} //!< Calculates quadrupole moment tensor (with trace)

TEST_CASE("[Faunus] quadrupoleMoment") {
    using doctest::Approx;
    ParticleVector p(2);
    p[0].pos = {10, 20, 30};
    p[1].pos = {-10, 0, -30};
    p[0].charge = -0.5;
    p[1].charge = 0.3;
    auto mu = quadrupoleMoment(p.begin(), p.end(), [](auto &) {}, {2, 3, 4});
    CHECK(mu.trace() == Approx(-60.9));
    p[0].charge *= -1.0;
    mu = quadrupoleMoment(p.begin(), p.end(), [](auto &) {}, {2, 3, 4});
    CHECK(mu.trace() == Approx(453.6));
}

/**
 * @brief Converts a group to a multipole-particle
 * @param g Group
 * @param boundary Function to use for boundary
 * @param cutoff Cut-off for included particles with regard to origin, default value is infinite
 */
template <class Tgroup, class BoundaryFunction>
auto toMultipole(const Tgroup &g, BoundaryFunction boundary = [](const Point &) {}, double cutoff = pc::infty) {
    Particle m;
    m.pos = g.cm;
    m.charge = Faunus::monopoleMoment(g.begin(), g.end());                                // monopole
    m.getExt().mu = Faunus::dipoleMoment(g.begin(), g.end(), boundary, m.pos, cutoff);    // dipole
    m.getExt().Q = Faunus::quadrupoleMoment(g.begin(), g.end(), boundary, m.pos, cutoff); // quadrupole
    m.getExt().mulen = m.getExt().mu.norm();
    if (m.getExt().mulen > 1e-9)
        m.getExt().mu.normalize();
    return m;
} //<! Group --> Multipole

} // namespace Faunus
