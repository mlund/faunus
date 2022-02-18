#pragma once

#include "potentials_base.h"
#include "units.h"

namespace Faunus::SpheroCylinder {

/**
 * @brief Calculate perpendicular projection of the first vector versus the second vector
 *
 * Calculate projection of a vector to a plane define by second vector (normal of a plan)
 * @param a the first vector
 * @param b the second vector
 */
inline Point vec_perpproject(const Point& a, const Point& b) { return a - b * a.dot(b); }

/**
 * @brief Calculate minimum distance between two line segments
 *
 * Find closest distance between line segments and return its vector
 * gets orientations and lengths of line segments and the vector connecting
 * their center os masses (from vec1 to vec2)
 * Copyright 2001, softSurfer (www.softsurfer.com)
 * This code may be freely used and modified for any purpose
 * providing that this copyright notice is included with it.
 * SoftSurfer makes no warranty for this code, and cannot be held
 * liable for any real or imagined damage resulting from its use.
 * Users of this code must verify correctness for their application.
 *
 * @param dir1 Direction of first segment
 * @param halfl1 Half length of first segment
 * @param dir2 Direction of second segment
 * @param halfl2 Half length of second segment
 * @param r_cm Distance vector between the middle of the two segments
 */
Point mindist_segment2segment(const Point& dir1, double halfl1, const Point& dir2, double halfl2, const Point& r_cm);

/**
 * @param dir Direction of segment
 * @param halfl Half length of segment
 * @param r_cm Distance vector between the middle segment to point
 */
inline Point mindist_segment2point(const Point& dir, double halfl, const Point& r_cm) {
    double d;
    double c = dir.dot(r_cm);
    if (c > halfl) {
        d = halfl;
    } else {
        if (c > -halfl) {
            d = c;
        } else {
            d = -halfl;
        }
    }
    return -r_cm + (dir * d);
}

/**
 * Finds intersections of spherocylinder and plane defined by vector
 * "w_vec" and if they are in all-way patch then returns number of them (PSC)
 */
int find_intersect_plane(const Cigar& part1, const Cigar& part2, const Point& r_cm, const Point& w_vec,
                         double cutoff_squared, double cospatch, std::array<double, 5>& intersections);

/**
 * @brief Finds if vector "vec" has angular intersection w. patch of part1
 */
int test_intrpatch(const Cigar& part1, Point& vec, double cospatch, double ti, std::array<double, 5>& intersections);

/**
 * @brief Intersect of plane
 *
 * Finds intersections of plane defined by vector "w_vec"
 * and if they are in cylindrical patch then returns number of them (CPSC)
 */
int find_intersect_planec(const Cigar& part1, const Cigar& part2, const Point& r_cm, const Point& w_vec, double rcut2,
                          double cospatch, std::array<double, 5>& intersections);

/**
 * @brief Intersections of spherocylinder2 with a all-way patch of spherocylinder1 and return them (PSC)
 */
int psc_intersect(const Cigar& part1, const Cigar& part2, const Point& r_cm, std::array<double, 5>& intersections,
                  double rcut2);

/**
 * @brief Intersection of PSC2 with cylindrical patch of PSC1 and return them (CPSC)
 */
int cpsc_intersect(const Cigar& part1, const Cigar& part2, const Point& r_cm, std::array<double, 5>& intersections,
                   double rcut2);
} // namespace Faunus::SpheroCylinder

// ---------------------------------------------------------------------------------------------

namespace Faunus::Potential {

inline double fanglscale(const double a, const Cigar& cigar) {
    // a = r_ij * n_i
    if (a <= cigar.pcanglsw) {
        return 0.0;
    }
    if (a >= cigar.pcangl) {
        return 1.0;
    }
    return 0.5 - ((cigar.pcanglsw + cigar.pcangl) * 0.5 - a) / (cigar.pcangl - cigar.pcanglsw);
}

/** @brief Hard-sphere pair potential for spherocylinders */
class HardSpheroCylinder : public PairPotentialBase {
  public:
    double operator()(const Particle& particle1, const Particle& particle2, [[maybe_unused]] double d,
                      const Point& center_to_center_distance) const override {
        assert(particle1.hasExtension() && particle2.hasExtension());
        auto minimum_distance_squared = SpheroCylinder::mindist_segment2segment(
                                            particle1.ext->scdir, particle1.ext->half_length, particle2.ext->scdir,
                                            particle1.ext->half_length, center_to_center_distance)
                                            .squaredNorm();
        const auto contact_distance = 0.5 * (particle1.traits().sigma + particle2.traits().sigma);
        if (minimum_distance_squared < contact_distance * contact_distance) {
            return pc::infty;
        }
        return 0.0;
    } //!< Pair energy in units of kT

    HardSpheroCylinder();
    void to_json(json& j) const override;
    void from_json(const json& j) override;
};

/**
 * @brief Pair potential between a patchy sphero-cylinder (first particle) and a sphere (second particle)
 * @tparam PatchPotential Pair potential between sphere and point on patch (isotropic)
 * @tparam CylinderPotential Pair potential between sphere and closest cylinder part (isotropic)
 */
template <typename PatchPotential, typename CylinderPotential> class CigarWithSphere : public PairPotentialBase {
  private:
    PatchPotential patch_potential;
    CylinderPotential cylinder_potential;

  public:
    CigarWithSphere()
        : PairPotentialBase("cigar-sphere", ""s, false) {}

    inline double operator()(const Particle& cigar, const Particle& sphere, [[maybe_unused]] double distance_squared,
                             const Point& center_separation) const override {
        assert(cigar.hasExtension());
        const auto c = cigar.getExt().scdir.dot(center_separation);
        double contt = 0;
        if (c > cigar.ext->half_length) {
            contt = cigar.ext->half_length;
        } else {
            if (c > -cigar.ext->half_length) {
                contt = c;
            } else {
                contt = -cigar.ext->half_length;
            }
        }
        Point distvec = -center_separation + (cigar.ext->scdir * contt);
        if (cigar.traits().sphero_cylinder.type == SpheroCylinderData::PatchType::None &&
            sphere.traits().sphero_cylinder.type == SpheroCylinderData::PatchType::None) {
            return cylinder_potential(cigar, sphere, distvec.squaredNorm(), Point::Zero());
        }

        // patchy interaction
        const auto cutoff_squared = patch_potential.cutOffSquared(cigar.id, sphere.id);
        // scaling function: angular dependence of patch1
        Point vec1 = SpheroCylinder::vec_perpproject(distvec, cigar.ext->scdir).normalized();
        const auto s = vec1.dot(cigar.ext->patchdir);
        const auto f1 = fanglscale(s, cigar.getExt());

        // scaling function for the length of spherocylinder within cutoff
        auto ndist_squared = distvec.dot(distvec);
        auto t = sqrt(cutoff_squared - ndist_squared); // TODO cutoff
        double f0;
        if (contt + t > cigar.ext->half_length) {
            f0 = cigar.ext->half_length;
        } else {
            f0 = contt + t;
        }
        if (contt - t < -cigar.ext->half_length) {
            f0 -= -cigar.ext->half_length;
        } else {
            f0 -= contt - t;
        }
        return f1 * (f0 + 1.0) * patch_potential(cigar, sphere, ndist_squared, Point::Zero()) +
               cylinder_potential(cigar, sphere, ndist_squared, Point::Zero());
    }

    void to_json(json& j) const override {
        j["patch"] = static_cast<json>(patch_potential);
        j["cylinder"] = static_cast<json>(cylinder_potential);
    }

    void from_json(const json& j) override {
        patch_potential = j;
        cylinder_potential = j;
    }
};

/**
 * @brief Template for cigar-cigar interactions including patchy cigars
 *
 * This template takes care of cigar-cigar interactions. If there are patches
 * here we calculate scaling factors based on the overlapping segments of the two
 * cigars (spherocylinders) withn their patches. f0 is for size of overlapping segment
 * whicle f1 anf f2 are scaling fators for orientation of pacthes.
 * `PatchPotential` is an isotropic potential, typically CosAttract, used for the patch, while `IsotropicPairPotential`
 * is used for the remaining "isotropic" SC, typically WeeksChandlerAndersen.
 *
 * @tparam PatchPotential Pair potential between two points on the patches (isotropic, e.g. CosAttract)
 * @tparam CylinderPotential Pair potential between closest cylinder parts (isotropic, e.g. WCA)
 * @todo Energy calculation badly needs refactoring!
 */
template <typename PatchPotential, typename CylinderPotential> class CigarWithCigar : public PairPotentialBase {
  private:
    PatchPotential patch_potential;
    CylinderPotential cylinder_potential;

    double patchyPatchyEnergy(const Particle& particle1, const Particle& particle2,
                              const Point& center_separation) const { // patchy sc with patchy sc
        const auto cutoff_squared = patch_potential.cutOffSquared(particle1.id, particle1.id);
        std::array<double, 5> intersections;

        // distance for repulsion
        const auto rclose_squared =
            SpheroCylinder::mindist_segment2segment(particle1.ext->scdir, particle1.ext->half_length,
                                                    particle2.ext->scdir, particle2.ext->half_length, center_separation)
                .squaredNorm();
        // 1- do intersections of spherocylinder2 with patch of spherocylinder1 at.
        //  cut distance C
        int intrs = 0;
        intersections.fill(0.0);
        if (particle1.traits().sphero_cylinder.type == SpheroCylinderData::PatchType::Full) {
            intrs = SpheroCylinder::psc_intersect(particle1.getExt(), particle2.getExt(), center_separation,
                                                  intersections, cutoff_squared);
        } else {
            if (particle1.traits().sphero_cylinder.type == SpheroCylinderData::PatchType::Capped) {
                intrs = SpheroCylinder::cpsc_intersect(particle1.getExt(), particle2.getExt(), center_separation,
                                                       intersections, cutoff_squared);
            } else {
                throw std::runtime_error("unimplemented");
            }
        }
        if (intrs < 2) {
            // sc is all outside patch, attractive energy is 0
            return cylinder_potential(particle1, particle2, rclose_squared, Point::Zero());
        }
        auto T1 = intersections[0]; // points on sc2
        auto T2 = intersections[1];
        // 2- now do the same oposite way psc1 in patch of psc2
        intersections.fill(0.0);
        if (particle1.traits().sphero_cylinder.type ==
            SpheroCylinderData::PatchType::Full) { //!< @warning should this not be b.traits()?
            intrs = SpheroCylinder::psc_intersect(particle2.getExt(), particle1.getExt(), -center_separation,
                                                  intersections, cutoff_squared);
        } else {
            if (particle1.traits().sphero_cylinder.type ==
                SpheroCylinderData::PatchType::Capped) { //!< @warning should this not be particle2.traits()?
                intrs = SpheroCylinder::cpsc_intersect(particle2.getExt(), particle1.getExt(), -center_separation,
                                                       intersections, cutoff_squared);
            } else {
                throw std::runtime_error("unimplemented");
            }
        }
        if (intrs < 2) { // sc is all outside patch, attractive energy is 0
            return cylinder_potential(particle1, particle2, rclose_squared, Point::Zero());
        }
        auto S1 = intersections[0]; // points on sc1
        auto S2 = intersections[1];

        // 3- scaling function1: dependence on the length of intersetions
        auto v1 = fabs(S1 - S2) * 0.5;
        auto v2 = fabs(T1 - T2) * 0.5;
        auto f0 = v1 + v2;
        // 4a- with two intersection pices calculate vector between their CM
        //-this is for angular orientation
        Point vec1 = particle1.ext->scdir * (S1 + S2) * 0.5;
        Point vec2 = particle2.ext->scdir * (T1 + T2) * 0.5;
        Point vec_intrs = vec2 - vec1 - center_separation; // vec_intrs should be from sc1 t sc2

        // 5- scaling function2: angular dependence of patch1
        vec1 = SpheroCylinder::vec_perpproject(vec_intrs, particle1.ext->scdir);
        vec1.normalize();
        auto s = vec1.dot(particle1.ext->patchdir);
        auto f1 = fanglscale(s, particle1.getExt());

        // 6- scaling function3: angular dependence of patch2
        vec1 = SpheroCylinder::vec_perpproject(-vec_intrs, particle2.ext->scdir).normalized();
        s = vec1.dot(particle2.ext->patchdir);
        auto f2 = fanglscale(s, particle2.getExt());

        // 7 - calculate closest distance attractive energy from it
        auto vec_mindist =
            SpheroCylinder::mindist_segment2segment(particle1.ext->scdir, v1, particle2.ext->scdir, v2, vec_intrs);
        auto ndistsq = vec_mindist.dot(vec_mindist);

        // 8- put it all together and output scale
        return f0 * f1 * f2 * patch_potential(particle1, particle2, ndistsq, Point::Zero()) +
               cylinder_potential(particle1, particle2, rclose_squared, Point::Zero());
    }

    double isotropicIsotropicEnergy(const Particle& particle1, const Particle& particle2,
                                    const Point& center_separation) const { // isotropic sc with isotropic sc
        const auto mindist =
            SpheroCylinder::mindist_segment2segment(particle1.ext->scdir, particle1.ext->half_length,
                                                    particle2.ext->scdir, particle2.ext->half_length, center_separation)
                .squaredNorm();
        return patch_potential(particle1, particle2, mindist, Point::Zero()) +
               cylinder_potential(particle1, particle2, mindist, Point::Zero());
    }

  public:
    CigarWithCigar()
        : PairPotentialBase("cigar-cigar", ""s, false) {}

    double operator()(const Particle& particle1, const Particle& particle2,
                      [[maybe_unused]] double center_separation_squared,
                      const Point& center_separation) const override {

        assert(particle1.hasExtension() && particle2.hasExtension());

        if (particle1.traits().sphero_cylinder.type != SpheroCylinderData::PatchType::None &&
            particle2.traits().sphero_cylinder.type != SpheroCylinderData::PatchType::None) {
            return patchyPatchyEnergy(particle1, particle2, center_separation);
        }
        return isotropicIsotropicEnergy(particle1, particle2, center_separation);
    }

    void to_json(json& j) const override {
        j["patch"] = static_cast<json>(patch_potential);
        j["cylinder"] = static_cast<json>(cylinder_potential);
    }

    void from_json(const json& j) override {
        patch_potential = j;
        cylinder_potential = j;
    }
};

/**
 * @brief Set pair potential between cigar-cigar, sphere-sphere, cigar-sphere
 *
 * This takes three pair potentials that will be called dependent on the
 * nature of the two particles. If the sphero-cylinder has zero length it
 * is assumed to be a spherical, isotropic particle.
 *
 * @tparam PatchPotential Pair potential used for patch part (isotropic, e.g. CosAttract)
 * @tparam CylinderPotential Pair potential used for cylindrical path (isotropic, e.g. WCA)
 * @tparam CylinderPotential Pair potential used sphere-sphere interaction (isotropic, e.g. WCA)
 */
template <typename PatchPotential, typename CylinderPotential, typename SphereWithSphere = CylinderPotential>
class CompleteCigarPotential : public PairPotentialBase {
  private:
    SphereWithSphere sphere_sphere;                                  // pair potential between spheres
    CigarWithCigar<PatchPotential, CylinderPotential> cigar_cigar;   // pair potential between cigars
    CigarWithSphere<PatchPotential, CylinderPotential> cigar_sphere; // pair potential cigar <-> sphere

  public:
    inline double operator()(const Particle& a, const Particle& b, [[maybe_unused]] double distance_squared,
                             const Point& distance) const override {
        const double small_number = 1e-6;
        const auto a_is_sphere = !a.hasExtension() || a.ext->half_length < small_number;
        const auto b_is_sphere = !a.hasExtension() || b.ext->half_length < small_number;

        if (!a_is_sphere && !b_is_sphere) {
            return cigar_cigar(a, b, distance_squared, distance);
        }
        if (a_is_sphere && b_is_sphere) {
            return sphere_sphere(a, b, distance_squared, distance);
        }
        if (b_is_sphere) {
            return cigar_sphere(a, b, distance_squared, distance);
        }
        return cigar_sphere(b, a, distance_squared, distance);
    }

    void to_json([[maybe_unused]] json& j) const override { j = {sphere_sphere, cigar_cigar, cigar_sphere}; }

    void from_json(const json& j) override {
        sphere_sphere = j;
        cigar_cigar = j;
        cigar_sphere = j;
    }

    CompleteCigarPotential()
        : PairPotentialBase("complete cigar", ""s, false) {}
};

} // namespace Faunus::Potential