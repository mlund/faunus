#pragma once

#include "potentials_base.h"
#include "units.h"

/**
 * The algorithms found here are mainly direct conversions from
 * Robert Vacha's spherocylinder C code (~2008-2010).
 *
 * @todo cleanup and split into smaller functions.
 */
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
 * @param segment_direction Direction of segment
 * @param half_length Half length of segment
 * @param separation Distance vector between the middle segment to point
 */
inline Point mindist_segment2point(const Point& segment_direction, const double half_length, const Point& separation) {
    const auto c = segment_direction.dot(separation);
    double d;
    if (c > half_length) {
        d = half_length;
    } else {
        d = (c > -half_length) ? d = c : -half_length;
    }
    return -separation + (segment_direction * d);
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
int psc_intersect(const Cigar& particle1, const Cigar& particle2, const Point& r_cm,
                  std::array<double, 5>& intersections, double cutoff_squared);

/**
 * @brief Intersection of PSC2 with cylindrical patch of PSC1 and return them (CPSC)
 */
int cpsc_intersect(const Cigar& part1, const Cigar& part2, const Point& r_cm, std::array<double, 5>& intersections,
                   double rcut2);

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

} // namespace Faunus::SpheroCylinder

// --------------------

namespace Faunus::Potential {

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
    PatchPotential patch_potential;       //!< Isotropic pair-potential between patches
    CylinderPotential cylinder_potential; //!< Isotropic pair-potential between non-patchy parts

  public:
    CigarWithSphere();
    void to_json(json& j) const override;
    void from_json(const json& j) override;

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
        const Point distvec = -center_separation + (cigar.ext->scdir * contt);
        if (cigar.traits().sphero_cylinder.type == SpheroCylinderData::PatchType::None &&
            sphere.traits().sphero_cylinder.type == SpheroCylinderData::PatchType::None) {
            return cylinder_potential(cigar, sphere, distvec.squaredNorm(), Point::Zero());
        }

        // patchy interaction
        const auto cutoff_squared = patch_potential.cutOffSquared(cigar.id, sphere.id);
        // scaling function: angular dependence of patch1
        const Point vec1 = SpheroCylinder::vec_perpproject(distvec, cigar.ext->scdir).normalized();
        const auto s = vec1.dot(cigar.ext->patchdir);
        const auto f1 = SpheroCylinder::fanglscale(s, cigar.getExt());

        // scaling function for the length of spherocylinder within cutoff
        const auto ndist_squared = distvec.dot(distvec);
        const auto t = sqrt(cutoff_squared - ndist_squared); // TODO cutoff
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
};

/**
 * @brief Template for interactions between patchy sphero-cylinders (PSCs) interactions
 *
 * This template takes care of psc-psc interactions. If there are patches
 * we calculate scaling factors based on the overlapping segments of the two
 * PSCs withn their patches. f0 is for size of overlapping segment
 * whicle f1 anf f2 are scaling fators for orientation of pacthes.
 *
 * @tparam PatchPotential Pair potential between two points on the patches (isotropic, e.g. CosAttract)
 * @tparam CylinderPotential Pair potential between closest cylinder parts (isotropic, e.g. WCA)
 * @todo Energy calculation badly needs refactoring!
 */
template <typename PatchPotential, typename CylinderPotential> class CigarWithCigar : public PairPotentialBase {
  private:
    PatchPotential patch_potential;       //!< Isotropic pair-potential for patchy parts
    CylinderPotential cylinder_potential; //!< Isotropic pair-potential for cylindrical parts

    double patchyPatchyEnergy(const Particle& particle1, const Particle& particle2,
                              const Point& center_separation) const;

    double isotropicIsotropicEnergy(const Particle& particle1, const Particle& particle2,
                                    const Point& center_separation) const;

  public:
    CigarWithCigar();
    void to_json(json& j) const override;
    void from_json(const json& j) override;

    inline double operator()(const Particle& particle1, const Particle& particle2,
                             [[maybe_unused]] double center_separation_squared,
                             const Point& center_separation) const override {
        assert(particle1.hasExtension() && particle2.hasExtension());
        if (particle1.traits().sphero_cylinder.type != SpheroCylinderData::PatchType::None &&
            particle2.traits().sphero_cylinder.type != SpheroCylinderData::PatchType::None) {
            return patchyPatchyEnergy(particle1, particle2, center_separation);
        }
        return isotropicIsotropicEnergy(particle1, particle2, center_separation);
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

    void to_json([[maybe_unused]] json& j) const override;
    void from_json(const json& j) override;
    CompleteCigarPotential();
};

} // namespace Faunus::Potential