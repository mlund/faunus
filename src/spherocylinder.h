#pragma once

#include "core.h"
#include "particle.h"
#include "potentials.h"

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
int find_intersect_plane(const Cigar& part1, const Cigar& part2, const Point& r_cm, const Point& w_vec, double rcut2,
                         double cospatch, double intersections[5]);

/**
 * @brief Finds if vector "vec" has angular intersection w. patch of part1
 */
int test_intrpatch(const Cigar& part1, Point& vec, double cospatch, double ti, double intersections[5]);

/**
 * @brief Intersect of plane
 *
 * Finds intersections of plane defined by vector "w_vec"
 * and if they are in cylindrical patch then returns number of them (CPSC)
 */
int find_intersect_planec(const Cigar& part1, const Cigar& part2, const Point& r_cm, const Point& w_vec, double rcut2,
                          double cospatch, double intersections[5]);

/**
 * @brief Intersections of spherocylinder2 with a all-way patch of spherocylinder1 and return them (PSC)
 */
int psc_intersect(const Cigar& part1, const Cigar& part2, const Point& r_cm, double intersections[5], double rcut2);

/**
 * @brief Intersection of PSC2 with cylindrical patch of PSC1 and return them (CPSC)
 */
int cpsc_intersect(const Cigar& part1, const Cigar& part2, const Point& r_cm, double intersections[5], double rcut2);
} // namespace Faunus::SpheroCylinder

// ---------------------------------------------------------------------------------------------

namespace Faunus::Potential {

inline double fanglscale(const double a, const Cigar& particle) {
    // a = r_ij * n_i
    if (a <= particle.pcanglsw) {
        return 0.0;
    }
    if (a >= particle.pcangl) {
        return 1.0;
    }
    return 0.5 - ((particle.pcanglsw + particle.pcangl) * 0.5 - a) / (particle.pcangl - particle.pcanglsw);
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
 * @brief Brief description here (one line)
 *
 * Detailed description here...
 * For patchy spherocylinder `Tcigarshere` should be a combined pair potential,
 * where `first` accounts for patchy interaction and `second` is isotropic, only.
 */
template <typename Tcigarsphere> class PatchyCigarSphere {
  public:
    Tcigarsphere pairpot;

    inline double operator()(const Particle& a, const Particle& b, const Point& r_cm) {
        // 0- isotropic, 1-PSC all-way patch,2 -CPSC cylindrical patch
        // b is sphere, a is spherocylinder
        double s, t, f0, f1, contt;

        assert(a.half_length < 1e-6 && "First (a) should be cigar then sphere, not opposite!");
        double c = a.ext->scdir.dot(r_cm);
        if (c > a.ext->half_length) {
            contt = a.ext->half_length;
        } else {
            if (c > -a.ext->half_length) {
                contt = c;
            } else {
                contt = -a.ext->half_length;
            }
        }
        Point distvec = -r_cm + (a.ext->scdir * contt);

        if (a.traits().sphero_cylinder.patch_type == 0 && b.traits().sphero_cylinder.patch_type == 0) {
            return pairpot(a, b, distvec.dot(distvec));
        }

        // patchy interaction
        auto rcut2 = pairpot.first.rcut2(a.id, b.id);
        // scaling function: angular dependence of patch1
        Point vec1 = SpheroCylinder::vec_perpproject(distvec, a.ext->scdir).normalized();
        s = vec1.dot(a.ext->patchdir);
        f1 = fanglscale(s, a.getExt());

        // scaling function for the length of spherocylinder within cutoff
        auto ndistsq = distvec.dot(distvec);
        t = sqrt(rcut2 - ndistsq); // TODO cutoff
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
        return pairpot.first(a, b, ndistsq) * f1 * (f0 + 1.0) + pairpot.second(a, b, ndistsq);
    }
};

/**
 * @brief Template for cigar-cigar interactions including patchy cigars
 *
 * This template takes care of cigar-cigar interactions. If there are patches
 * here we calculate scaling factors based on the overlapping segments of the two
 * cigars (spherocylinders) withn their patches. f0 is for size of overlapping segment
 * whicle f1 anf f2 are scaling fators for orientation of pacthes.
 * For patchy spherocylinder `Tcigarcigar` potential should be a combined pair potential,
 * where `first` accounts for patchy interaction and `second` is isotropic, only.
 * There are two types of patches evaluated here:
 * patchtype 1 (PSC) runs along the whole axis including the ends
 * patchtype 2 (CPSC) is limited to the cylindrical part
 */
template <typename Tcigarcigar> class PatchyCigarCigar {
  public:
    Tcigarcigar pairpot;

    inline double operator()(const Particle& a, const Particle& b, const Point& r_cm) {
        // 0- isotropic, 1-PSC all-way patch,2 -CPSC cylindrical patch
        if (a.traits().sphero_cylinder.patch_type > 0) {
            if (b.traits().sphero_cylinder.patch_type > 0) {
                // patchy sc with patchy sc
                int intrs;
                double rcut2 = pairpot.first.rcut2(a.id, b.id);
                double ndistsq;
                double v1, v2, f0, f1, f2, T1, T2, S1, S2, s;
                double intersections[5];
                Point vec1, vec2, vec_intrs, vec_mindist;

                assert(rcut2 > 1e-6 && "Cutoff for patchy cigar interaction has zero size.");
                // distance for repulsion
                Point rclose = SpheroCylinder::mindist_segment2segment(a.ext->scdir, a.ext->half_length, b.ext->scdir,
                                                                       b.ext->half_length, r_cm);
                intrs = 0;
                for (int i = 0; i < 5; i++) {
                    intersections[i] = 0;
                }
                // 1- do intersections of spherocylinder2 with patch of spherocylinder1 at.
                //  cut distance C
                if (a.traits().sphero_cylinder.patch_type == 1) {
                    intrs = SpheroCylinder::psc_intersect(a.getExt(), b.getExt(), r_cm, intersections, rcut2);
                } else {
                    if (a.traits().sphero_cylinder.patch_type == 2) {
                        intrs = SpheroCylinder::cpsc_intersect(a.getExt(), b.getExt(), r_cm, intersections, rcut2);
                    } else {
                        assert(!"Patchtype not implemented!");
                    }
                }
                if (intrs < 2) {
                    return pairpot.second(a, b, rclose.squaredNorm());
                }                      // sc is all outside patch, attractive energy is 0
                T1 = intersections[0]; // points on sc2
                T2 = intersections[1];
                // 2- now do the same oposite way psc1 in patch of psc2
                for (int i = 0; i < 5; i++) {
                    intersections[i] = 0;
                }
                if (a.traits().sphero_cylinder.patch_type == 1) {
                    intrs = SpheroCylinder::psc_intersect(b.getExt(), a.getExt(), -r_cm, intersections, rcut2);
                } else {
                    if (a.traits().sphero_cylinder.patch_type == 2) {
                        intrs = SpheroCylinder::cpsc_intersect(b.getExt(), a.getExt(), -r_cm, intersections, rcut2);
                    } else {
                        assert(!"Patchtype not implemented!");
                    }
                }
                if (intrs < 2) {
                    return pairpot.second(a, b,
                                          rclose.squaredNorm()); // sc is all outside patch, attractive energy is 0
                }
                S1 = intersections[0]; // points on sc1
                S2 = intersections[1];

                // 3- scaling function1: dependence on the length of intersetions
                v1 = fabs(S1 - S2) * 0.5;
                v2 = fabs(T1 - T2) * 0.5;
                f0 = v1 + v2;
                // 4a- with two intersection pices calculate vector between their CM
                //-this is for angular orientation
                vec1 = a.ext->scdir * (S1 + S2) * 0.5;
                vec2 = b.ext->scdir * (T1 + T2) * 0.5;
                vec_intrs = vec2 - vec1 - r_cm;
                // vec_intrs should be from sc1 t sc2
                // 4b - calculate closest distance attractive energy from it
                vec_mindist = SpheroCylinder::mindist_segment2segment(a.ext->scdir, v1, b.ext->scdir, v2, vec_intrs);
                ndistsq = vec_mindist.dot(vec_mindist);

                // 5- scaling function2: angular dependence of patch1
                vec1 = SpheroCylinder::vec_perpproject(vec_intrs, a.ext->scdir);
                vec1.normalize();
                s = vec1.dot(a.ext->patchdir);
                f1 = fanglscale(s, a.getExt());

                // 6- scaling function3: angular dependence of patch2
                vec1 = SpheroCylinder::vec_perpproject(-vec_intrs, b.ext->scdir).normalized();
                s = vec1.dot(b.ext->patchdir);
                f2 = fanglscale(s, b.getExt());

                // 8- put it all together and output scale
                return f0 * f1 * f2 * pairpot.first(a, b, ndistsq) + pairpot.second(a, b, rclose.squaredNorm());
            }
            assert(!"PSC w. isotropic cigar not implemented!");
        } else {
            if (b.traits().sphero_cylinder.patch_type > 0) {
                assert(!"PSC w. isotropic cigar not implemented!");
                // isotropic sc with patchy sc - we dont have at the moment
            } else {
                // isotropic sc with isotropic sc
                Point rclose = SpheroCylinder::mindist_segment2segment(a.ext->scdir, a.ext->half_length, b.ext->scdir,
                                                                       b.ext->half_length, r_cm);
                return pairpot(a, b, rclose.squaredNorm());
            }
        }
        assert(!"Something we have not implemented");
        return 0.0;
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
    PatchyCigarCigar<Tcigarcigar> pairpot_cc;
    PatchyCigarSphere<Tcigarsphere> pairpot_cs;

    inline double operator()(const Particle& a, const Particle& b, const Point& r_cm) {
        if (a.ext->half_length < 1e-6) {
            // a sphere - b sphere
            if (b.ext->half_length < 1e-6) {
                return pairpot_ss(a, b, r_cm.squaredNorm());
            }
            // a sphere - b cigar
            // PatchyCigarSphere(b,a)
            assert(b.dir.squaredNorm() > 1e-6 && "Direction vector of patchy spherocylinder has zero size.");
            assert(b.patchdir.squaredNorm() > 1e-6 && "Patch direction vector of patchy spherocylinder has zero size.");
            assert(b.patchsides[0].squaredNorm() > 1e-6 && "Vector associated with patch side has zero size. "
                                                           "Patchy spherocylinder were probably not initialized.");
            return pairpot_cs(b, a, r_cm);
        }
        // a cigar - b sphere
        if (b.ext->half_length < 1e-6) {
            // PatchyCigarSphere(a,b)
            assert(a.dir.squaredNorm() > 1e-6 && "Direction vector of patchy spherocylinder has zero size.");
            assert(a.patchdir.squaredNorm() > 1e-6 && "Patch direction vector of patchy spherocylinder has zero size.");
            assert(a.patchsides[0].squaredNorm() > 1e-6 && "Vector associated with patch side has zero size. "
                                                           "Patchy spherocylinder were probably not initialized.");
            return pairpot_cs(a, b, r_cm);
        }
        // a cigar - b cigar
        // PatchyCigarCigar
        assert(a.dir.squaredNorm() > 1e-6 && "Direction vector of patchy spherocylinder has zero size.");
        assert(b.dir.squaredNorm() > 1e-6 && "Direction vector of patchy spherocylinder has zero size.");
        return pairpot_cc(a, b, r_cm);
    }
};

} // namespace Faunus::Potential