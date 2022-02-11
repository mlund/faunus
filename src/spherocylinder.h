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
 * @brief Pair potential between a sphere and a patchy sphero-cylinder
 */
template <typename SpherePatchPotential, typename SphereCylinderPotential>
class PatchyCigarSphere : public PairPotentialBase {
  private:
    SpherePatchPotential first;
    SphereCylinderPotential second;
  public:

    inline double operator()(const Particle& a, const Particle& b, [[maybe_unused]] double distance_squared,
                      const Point& center_separation) const override {
        // 0- isotropic, 1-PSC all-way patch,2 -CPSC cylindrical patch
        // b is sphere, a is spherocylinder
        double contt = 0;
        assert(a.ext->half_length < 1e-6); // First (a) should be cigar then sphere, not opposite!
        const auto c = a.ext->scdir.dot(center_separation);
        if (c > a.ext->half_length) {
            contt = a.ext->half_length;
        } else {
            if (c > -a.ext->half_length) {
                contt = c;
            } else {
                contt = -a.ext->half_length;
            }
        }
        Point distvec = -center_separation + (a.ext->scdir * contt);
        if (a.traits().sphero_cylinder.type == SpheroCylinderData::PatchType::None &&
            b.traits().sphero_cylinder.type == SpheroCylinderData::PatchType::None) {
            return second(a, b, distvec.squaredNorm(), Point::Zero());
        }

        // patchy interaction
        const auto cutoff_squared = first.cutoff_squared(a.id, b.id);
        // scaling function: angular dependence of patch1
        Point vec1 = SpheroCylinder::vec_perpproject(distvec, a.ext->scdir).normalized();
        auto s = vec1.dot(a.ext->patchdir);
        auto f1 = fanglscale(s, a.getExt());

        // scaling function for the length of spherocylinder within cutoff
        auto ndist_squared = distvec.dot(distvec);
        auto t = sqrt(cutoff_squared - ndist_squared); // TODO cutoff
        double f0;
        if (contt + t > a.ext->half_length) {
            f0 = a.ext->half_length;
        } else {
            f0 = contt + t;
        }
        if (contt - t < -a.ext->half_length) {
            f0 -= -a.ext->half_length;
        } else {
            f0 -= contt - t;
        }
        return first(a, b, ndist_squared, Point::Zero()) * f1 * (f0 + 1.0) + second(a, b, ndist_squared, Point::Zero());
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
 * There are two types of patches evaluated here:
 * patchtype 1 (PSC) runs along the whole axis including the ends
 * patchtype 2 (CPSC) is limited to the cylindrical part
 *
 * @todo Energy calculation badly needs refactoring!
 */
template <typename PatchPotential, typename SpheroCylinderPotential> class PatchyCigarCigar : public PairPotentialBase {
  private:
    PatchPotential patch_pairpotential;             // typically `CosAttact`
    SpheroCylinderPotential spherocylinder_pairpot; // typically `WeeksChandlerAndersen`

    double patchyPatchyEnergy(const Particle& particle1, const Particle& particle2,
                              const Point& center_separation) const { // patchy sc with patchy sc
        const auto cutoff_squared = patch_pairpotential.cutOffSquared();
        std::array<double, 5> intersections;

        // distance for repulsion
        const auto rclose = SpheroCylinder::mindist_segment2segment(particle1.ext->scdir, particle1.ext->half_length,
                                                                    particle2.ext->scdir, particle2.ext->half_length,
                                                                    center_separation);
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
            return spherocylinder_pairpot(particle1, particle2, rclose.squaredNorm(), Point::Zero());
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
                SpheroCylinderData::PatchType::Capped) { //!< @warning should this not be b.traits()?
                intrs = SpheroCylinder::cpsc_intersect(particle2.getExt(), particle1.getExt(), -center_separation,
                                                       intersections, cutoff_squared);
            } else {
                throw std::runtime_error("unimplemented");
            }
        }
        if (intrs < 2) {
            return spherocylinder_pairpot(particle1, particle2, rclose.squaredNorm(),
                                          Point::Zero()); // sc is all outside patch, attractive energy is 0
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
        return f0 * f1 * f2 * patch_pairpotential(particle1, particle2, ndistsq, Point::Zero()) +
               spherocylinder_pairpot(particle1, particle2, rclose.squaredNorm(), Point::Zero());
    }

    double isotropicIsotropicEnergy(const Particle& particle1, const Particle& particle2,
                                    const Point& center_separation) const { // isotropic sc with isotropic sc
        const auto mindist =
            SpheroCylinder::mindist_segment2segment(particle1.ext->scdir, particle1.ext->half_length,
                                                    particle2.ext->scdir, particle2.ext->half_length, center_separation)
                .squaredNorm();
        return patch_pairpotential(particle1, particle2, mindist, Point::Zero()) +
               spherocylinder_pairpot(particle1, particle2, mindist, Point::Zero());
    }

  public:
    PatchyCigarCigar()
        : PairPotentialBase("patchy-cigar-cigar", ""s, false) {}

    double operator()(const Particle& particle1, const Particle& particle2,
                      [[maybe_unused]] double center_separation_squared,
                      const Point& center_separation) const override {
        if (particle1.traits().sphero_cylinder.type != SpheroCylinderData::PatchType::None &&
            particle2.traits().sphero_cylinder.type != SpheroCylinderData::PatchType::None) {
            return patchyPatchyEnergy(particle1, particle2, center_separation);
        }
        if (particle1.traits().sphero_cylinder.type == SpheroCylinderData::PatchType::None &&
            particle2.traits().sphero_cylinder.type == SpheroCylinderData::PatchType::None) {
            return isotropicIsotropicEnergy(particle1, particle2, center_separation);
        }
        throw std::runtime_error("PSC w. isotropic cigar not implemented");
    }

    void to_json(json& j) const override {
        j["patch"] = static_cast<json>(patch_pairpotential);
        j["cylinder"] = static_cast<json>(spherocylinder_pairpot);
    }
    void from_json(const json& j) override {
        if (!j.contains("patch") || !j.contains("cylinder")) {
            throw ConfigurationError("patch and/or cylinder undefined");
        }
        patch_pairpotential = j["patch"];
        spherocylinder_pairpot = j["cylinder"];
    }
};

/**
 * @brief Set pair potential between cigar-cigar, sphere-sphere, cigar-sphere
 *
 * This takes three pair potentials that will be called dependent on the
 * nature of the two particles. If the sphero-cylinder has zero length it
 * is assumed to be a spherical, isotropic particle.
 */
template <typename Tcigarcigar, typename Tspheresphere, typename Tcigarsphere> class CigarSphereSplit {
  public:
    Tspheresphere pairpot_ss;
    PatchyCigarCigar<Tcigarcigar, Tcigarcigar> pairpot_cc;
    PatchyCigarSphere<Tcigarsphere, Tcigarsphere> pairpot_cs;

    inline double operator()(const Particle& a, const Particle& b, const Point& r_cm) {
        if (a.ext->half_length < 1e-6) {
            // a sphere - b sphere
            if (b.ext->half_length < 1e-6) {
                return pairpot_ss(a, b, r_cm.squaredNorm());
            }
            // a sphere - b cigar
            // PatchyCigarSphere(b,a)
            assert(b.ext->scdir.squaredNorm() > 1e-6 && "Direction vector of patchy spherocylinder has zero size.");
            assert(b.ext->patchdir.squaredNorm() > 1e-6 && "Patch direction vector of patchy spherocylinder has zero size.");
            assert(b.ext->patchsides[0].squaredNorm() > 1e-6 && "Vector associated with patch side has zero size. "
                                                           "Patchy spherocylinder were probably not initialized.");
            return pairpot_cs(b, a, r_cm);
        }
        // a cigar - b sphere
        if (b.ext->half_length < 1e-6) {
            // PatchyCigarSphere(a,b)
            assert(a.ext->scdir.squaredNorm() > 1e-6 && "Direction vector of patchy spherocylinder has zero size.");
            assert(a.ext->patchdir.squaredNorm() > 1e-6 && "Patch direction vector of patchy spherocylinder has zero size.");
            assert(a.ext->patchsides[0].squaredNorm() > 1e-6 && "Vector associated with patch side has zero size. "
                                                           "Patchy spherocylinder were probably not initialized.");
            return pairpot_cs(a, b, r_cm);
        }
        // a cigar - b cigar
        // PatchyCigarCigar
        assert(a.ext->scdir.squaredNorm() > 1e-6 && "Direction vector of patchy spherocylinder has zero size.");
        assert(b.ext->scdir.squaredNorm() > 1e-6 && "Direction vector of patchy spherocylinder has zero size.");
        return pairpot_cc(a, b, r_cm);
    }
};

} // namespace Faunus::Potential