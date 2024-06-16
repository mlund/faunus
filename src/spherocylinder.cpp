#include "spherocylinder.h"
#include <Eigen/Geometry>
#include "potentials.h"

namespace Faunus::SpheroCylinder {

/**
 * @param dir1
 * @param half_length1
 * @param dir2
 * @param half_length2
 * @param r_cm
 * @return
 */
Point mindist_segment2segment(const Point& dir1, const double half_length1, const Point& dir2,
                              const double half_length2, const Point& r_cm) {
    constexpr auto very_small_number = 0.00000001;
    Point u = dir1 * (half_length1 * 2);                        // S1.P1 - S1.P0;
    Point v = dir2 * (half_length2 * 2);                        // S2.P1 - S2.P0;
    Point w = dir2 * half_length2 - dir1 * half_length1 - r_cm; // S1.P0-S2.P0
    auto a = u.dot(u);                                          // always >= 0
    auto b = u.dot(v);
    auto c = v.dot(v); // always >= 0
    auto d = u.dot(w);
    auto e = v.dot(w);
    auto D = a * c - b * b; // always >= 0
    auto sc = D;
    auto sN = D;
    auto sD = D; // sc = sN / sD, default sD = D >= 0
    auto tc = D;
    auto tN = D;
    auto tD = D; // tc = tN / tD, default tD = D >= 0

    // compute the line parameters of the two closest points
    if (D < very_small_number) { // the lines are almost parallel
        sN = 0.0;                // force using point P0 on segment S1
        sD = 1.0;                // to prevent possible division by 0.0 later
        tN = e;
        tD = c;
    } else { // get the closest points on the infinite lines
        sN = (b * e - c * d);
        tN = (a * e - b * d);
        if (sN < 0.0) { // sc < 0 => the s=0 edge is visible
            sN = 0.0;
            tN = e;
            tD = c;
        } else if (sN > sD) { // sc > 1 => the s=1 edge is visible
            sN = sD;
            tN = e + b;
            tD = c;
        }
    }
    if (tN < 0.0) { // tc < 0 => the t=0 edge is visible
        tN = 0.0;
        // recompute sc for this edge
        if (-d < 0.0) {
            sN = 0.0;
        } else if (-d > a) {
            sN = sD;
        } else {
            sN = -d;
            sD = a;
        }
    } else if (tN > tD) { // tc > 1 => the t=1 edge is visible
        tN = tD;
        // recompute sc for this edge
        if ((-d + b) < 0.0) {
            sN = 0;
        } else if ((-d + b) > a) {
            sN = sD;
        } else {
            sN = (-d + b);
            sD = a;
        }
    }
    // finally do the division to get sc and tc
    sc = (std::fabs(sN) < very_small_number) ? 0.0 : sN / sD;
    tc = (std::fabs(tN) < very_small_number) ? 0.0 : tN / tD;

    // get the difference of the two closest points
    return u * sc + w - v * tc;
}

/**
 * @param part1 Particle 1
 * @param part2 Particle 2
 * @param r_cm
 * @param w_vec
 * @param cutoff_squared
 * @param cospatch
 * @param intersections
 * @return
 */
int find_intersect_plane(const Cigar& part1, const Cigar& part2, const Point& r_cm, const Point& w_vec,
                         const double cutoff_squared, const double cospatch, std::array<double, 5>& intersections) {
    Point nplane = part1.scdir.cross(w_vec).normalized();
    const auto a = nplane.dot(part2.scdir);
    if (std::fabs(a) <= pc::epsilon_dbl) {
        return 0; /* there is no intersection plane and sc are paralel*/
    }
    const auto ti = nplane.dot(r_cm) / a;
    if ((ti > part2.half_length) || (ti < -part2.half_length)) {
        return 0; /* there is no intersection plane sc is too short*/
    }
    Point intersection_to_cm = ti * part2.scdir - r_cm; /*vector from intersection point to CM*/
    const auto c = intersection_to_cm.dot(w_vec);
    if (c * cospatch < 0) {
        return 0; /* the intersection in plane is on other side of patch */
    }
    const auto d = std::fabs(intersection_to_cm.dot(part1.scdir)) - part2.half_length;
    const auto disti = (d > 0.0) ? (d * d + c * c) : (c * c); // inside either cylinder or patch?
    if (disti > cutoff_squared) {
        return 0; /* the intersection is outside sc */
    }
    int intrs = 1;
    int i = 0;
    while (intersections[i] != 0) {
        if (ti == intersections[i]) {
            intrs = 0;
        } /* found intersection we already have -it is at boundary*/
        i++;
    }
    if (intrs > 0) {
        intersections[i] = ti;
    }
    return intrs;
}

/**
 * @param part1
 * @param vec
 * @param cospatch
 * @param ti
 * @param intersections
 * @return
 */
int test_intrpatch(const Cigar& part1, Point& vec, double cospatch, double ti, std::array<double, 5>& intersections) {
    /*test if we have intersection*/
    /* do projection to patch plane*/
    vec = vec_perpproject(vec, part1.scdir).normalized();
    /* test angle distance from patch*/
    const auto a = part1.patchdir.dot(vec);
    int intrs = 0;
    if (a >= cospatch) {
        intrs = 1;
        int i = 0;
        while (intersections[i] != 0) {
            if (ti == intersections[i]) {
                intrs = 0;
            } /* found intersection we already have -it is at boundary*/
            i++;
        }
        if (intrs > 0) {
            intersections[i] = ti;
        }
    }
    return intrs;
}

/**
 * @param part1
 * @param part2
 * @param r_cm
 * @param w_vec
 * @param rcut2
 * @param cospatch
 * @param intersections
 * @return
 */
int find_intersect_planec(const Cigar& part1, const Cigar& part2, const Point& r_cm, const Point& w_vec, double rcut2,
                          double cospatch, std::array<double, 5>& intersections) {
    Point nplane = part1.scdir.cross(w_vec).normalized();
    auto a = nplane.dot(part2.scdir);
    if (std::fabs(a) <= pc::epsilon_dbl) {
        return 0; // there is no intersection plane and sc are paralel
    }
    const auto ti = nplane.dot(r_cm) / a;
    if ((ti > part2.half_length) || (ti < -part2.half_length)) {
        return 0; // there is no intersection plane sc is too short
    }
    Point d_vec = ti * part2.scdir - r_cm; /*vector from intersection point to CM*/
    const auto c = d_vec.dot(w_vec);
    if (c * cospatch < 0) {
        return 0; /* the intersection in plane is on other side of patch */
    }
    int intrs = 0;
    const auto d = fabs(d_vec.dot(part1.scdir)) - part2.half_length;
    if (d <= 0) {
        auto disti = c * c; /*is inside cylinder*/
        if (disti > rcut2) {
            return 0; /* the intersection is outside sc */
        }
        intrs = 1;
        int i = 0;
        while (intersections[i] != 0) {
            if (ti == intersections[i]) {
                intrs = 0; /* found intersection we already have -it is at boundary*/
            }
            i++;
        }
        if (intrs > 0) {
            intersections[i] = ti;
        }
    }
    return intrs;
}

/**
 * @param particle1
 * @param particle2
 * @param r_cm
 * @param intersections
 * @param cutoff_squared
 * @return
 * @todo Add documentation and split into smaller parts(!)
 */
int psc_intersect(const Cigar& particle1, const Cigar& particle2, const Point& r_cm,
                  std::array<double, 5>& intersections, const double cutoff_squared) {
    double a, b, c, d, e, x1, x2;
    Point cm21, vec1, vec2, vec3, vec4;

    /*1- do intersections of spherocylinder2 with patch of spherocylinder1 at
      cut distance C*/
    /*1a- test intersection with half planes of patch and look how far they are
      from spherocylinder. If closer than C  we got itersection*/

    /* plane1 */
    /* find intersections of part2 with plane by par1 and patchsides[0] */
    int intrs = 0;
    intrs += find_intersect_plane(particle1, particle2, r_cm, particle1.patchsides[0], cutoff_squared,
                                  particle1.pcanglsw, intersections);
    /* plane2 */
    /* find intersections of part2 with plane by par1 and patchsides[1] */
    intrs += find_intersect_plane(particle1, particle2, r_cm, particle1.patchsides[1], cutoff_squared,
                                  particle1.pcanglsw, intersections);

    if ((intrs == 2) && (particle1.pcanglsw < 0)) {
        assert("Patch>180 -> two segments.");
    }

    /*1b- test intersection with cylinder - it is at distance C*/
    if (intrs < 2) {
        cm21 = -r_cm;
        vec1 = cm21.cross(particle1.scdir);
        vec2 = particle2.scdir.cross(particle1.scdir);
        a = vec2.dot(vec2);
        b = 2 * vec1.dot(vec2);
        c = -cutoff_squared + vec1.dot(vec1);
        d = b * b - 4 * a * c;
        if (d >= 0) {                      /*there is intersection with infinite cylinder */
            x1 = (-b + sqrt(d)) * 0.5 / a; /*parameter on line of SC2 determining intersection*/
            if ((x1 >= particle2.half_length) || (x1 <= -particle2.half_length)) {
                intrs += 0;
            } /*intersection is outside sc2*/
            else {
                /* vectors from center os sc1 to intersection with infinite cylinder*/
                vec1 = particle2.scdir * x1 - r_cm;
                e = particle1.scdir.dot(vec1);
                if ((e >= particle1.half_length) || (e <= -particle1.half_length)) {
                    intrs += 0;
                } /*intersection is outside sc1*/
                else {
                    intrs += test_intrpatch(particle1, vec1, particle1.pcanglsw, x1, intersections);
                }
            }
            if (d > 0) {
                x2 = (-b - sqrt(d)) * 0.5 / a; /*parameter on line of SC2 determining intersection*/
                if ((x2 >= particle2.half_length) || (x2 <= -particle2.half_length)) {
                    intrs += 0;
                } /*intersection is outside sc2*/
                else {
                    vec2 = particle2.scdir * x2 - r_cm;
                    e = particle1.scdir.dot(vec2);
                    if ((e >= particle1.half_length) || (e <= -particle1.half_length))
                        intrs += 0; /*intersection is outside sc1*/
                    else {
                        intrs += test_intrpatch(particle1, vec2, particle1.pcanglsw, x2, intersections);
                    }
                }
            }
        }
    }
    /*1c- test intersection with spheres at the end - it is at distace C*/
    if (intrs < 2) {
        /*centers of spheres*/
        /*relative to the CM of sc2*/
        vec1 = particle1.scdir * particle1.half_length - r_cm;
        vec2 = -particle1.scdir * particle1.half_length - r_cm;

        /*sphere1*/
        a = particle2.scdir.dot(particle2.scdir);
        b = 2.0 * vec1.dot(particle2.scdir);
        c = vec1.dot(vec1) - cutoff_squared;
        d = b * b - 4 * a * c;
        if (d >= 0) {                      /*if d<0 there are no intersections*/
            x1 = (-b + sqrt(d)) * 0.5 / a; /*parameter on line of SC2 determining intersection*/
            if ((x1 >= particle2.half_length) || (x1 <= -particle2.half_length)) {
                intrs += 0;
            } /*intersection is outside sc2*/
            else {
                vec3 = particle2.scdir * x1 - r_cm;
                e = particle1.scdir.dot(vec3);
                if ((e >= particle1.half_length) ||
                    (e <= -particle1.half_length)) { /*if not intersection is inside sc1*/
                    intrs += test_intrpatch(particle1, vec3, particle1.pcanglsw, x1, intersections);
                }
            }
            if (d > 0) {
                x2 = (-b - sqrt(d)) * 0.5 / a; /*parameter on line of SC2 determining intersection*/
                if ((x2 >= particle2.half_length) || (x2 <= -particle2.half_length)) {
                    intrs += 0;
                } /*intersection is outside sc2*/
                else {
                    vec4 = particle2.scdir * x2 - r_cm;
                    e = particle1.scdir.dot(vec4);
                    if ((e >= particle1.half_length) ||
                        (e <= -particle1.half_length)) { /*if not intersection is inside sc1*/
                        intrs += test_intrpatch(particle1, vec4, particle1.pcanglsw, x2, intersections);
                    }
                }
            }
        }
        /*sphere2*/
        a = particle2.scdir.dot(particle2.scdir);
        b = 2.0 * vec2.dot(particle2.scdir);
        c = vec2.dot(vec2) - cutoff_squared;
        d = b * b - 4 * a * c;
        if (d >= 0) {                      /*if d<0 there are no intersections*/
            x1 = (-b + sqrt(d)) * 0.5 / a; /*parameter on line of SC2 determining intersection*/
            if ((x1 >= particle2.half_length) || (x1 <= -particle2.half_length)) {
                intrs += 0;
            } /*intersection is outside sc2*/
            else {
                vec3 = particle2.scdir * x1 - r_cm;
                e = particle1.scdir.dot(vec3);
                if ((e >= particle1.half_length) ||
                    (e <= -particle1.half_length)) { /*if not intersection is inside sc1*/
                    intrs += test_intrpatch(particle1, vec3, particle1.pcanglsw, x1, intersections);
                }
            }
            if (d > 0) {
                x2 = (-b - sqrt(d)) * 0.5 / a; /*parameter on line of SC2 determining intersection*/
                if ((x2 >= particle2.half_length) || (x2 <= -particle2.half_length)) {
                    intrs += 0;
                } /*intersection is outside sc2*/
                else {
                    vec4 = particle2.scdir * x2 - r_cm;
                    e = particle1.scdir.dot(vec4);
                    if ((e >= particle1.half_length) ||
                        (e <= -particle1.half_length)) { /*if not intersection is inside sc1*/
                        intrs += test_intrpatch(particle1, vec4, particle1.pcanglsw, x2, intersections);
                    }
                }
            }
        }
    }

    /*1d- if there is only one itersection shperocylinder ends within patch wedge
      set as second intersection end inside patch*/
    if (intrs < 2) {
        /*whole spherocylinder is in or all out if intrs ==0*/
        vec1 = particle2.scdir * particle2.half_length - r_cm;
        /*vector from CM of sc1 to end of sc2*/
        /*check is is inside sc1*/
        a = vec1.dot(particle1.scdir);
        vec3 = vec1 - particle1.scdir * a;
        b = vec3.dot(vec3);
        d = fabs(a) - particle1.half_length;
        if (d <= 0) {
            c = b;
        } /*is inside cylindrical part*/
        else {
            c = d * d + b;
        } /*is inside caps*/
        /*c is distance squared from line or end to test if is inside sc*/
        if (c < cutoff_squared) {
            intrs += test_intrpatch(particle1, vec1, particle1.pcanglsw, particle2.half_length, intersections);
        }
        if (intrs < 2) {
            vec2 = -particle2.scdir * particle2.half_length - r_cm;
            /*check is is inside sc1*/
            a = vec2.dot(particle1.scdir);
            vec4 = vec2 - particle1.scdir * a;
            b = vec4.dot(vec4);
            d = fabs(a) - particle1.half_length;
            if (d <= 0) {
                c = b;
            } /*is inside cylindrical part*/
            else {
                c = d * d + b;
            } /*is inside caps*/
            /*c is distance squared from line or end to test if is inside sc*/
            if (c < cutoff_squared) {
                intrs +=
                    test_intrpatch(particle1, vec2, particle1.pcanglsw, -1.0 * particle2.half_length, intersections);
            }
        }
    }
    return intrs;
}

/**
 * @param cigar1
 * @param cigar2
 * @param r_cm
 * @param intersections
 * @param cutoff_squared
 * @return
 */
int cpsc_intersect(const Cigar& cigar1, const Cigar& cigar2, const Point& r_cm, std::array<double, 5>& intersections,
                   const double cutoff_squared) {
    int intrs;
    double a, b, c, d, e, x1, x2;
    Point cm21, vec1, vec2, vec3, vec4;

    intrs = 0;
    /*1- do intersections of spherocylinder2 with patch of spherocylinder1 at
      cut distance C*/
    /*1a- test intersection with half planes of patch and look how far they are
      from spherocylinder. If closer than C  we got itersection*/

    /* plane1 */
    /* find intersections of part2 with plane by par1 and part1.patchsides[0] */
    intrs += find_intersect_planec(cigar1, cigar2, r_cm, cigar1.patchsides[0], cutoff_squared, cigar1.pcanglsw,
                                   intersections);
    /* plane2 */
    /* find intersections of part2 with plane by par1 and part1.patchsides[1] */
    intrs += find_intersect_planec(cigar1, cigar2, r_cm, cigar1.patchsides[1], cutoff_squared, cigar1.pcanglsw,
                                   intersections);

    if ((intrs == 2) && (cigar1.pcanglsw < 0)) {
        throw std::runtime_error("Patch larger than 180 deg -> two segments - this is not yet implemented");
    }

    /*1b- test intersection with cylinder - it is at distance C*/
    if (intrs < 2) {
        cm21 = -r_cm;
        vec1 = cm21.cross(cigar1.scdir);
        vec2 = cigar2.scdir.cross(cigar1.scdir);
        a = vec2.dot(vec2);
        b = 2 * vec1.dot(vec2);
        c = -cutoff_squared + vec1.dot(vec1);
        d = b * b - 4 * a * c;
        if (d >= 0) {                      /*there is intersection with infinite cylinder */
            x1 = (-b + sqrt(d)) * 0.5 / a; /*parameter on line of SC2 determining intersection*/
            if ((x1 >= cigar2.half_length) || (x1 <= -cigar2.half_length)) {
                intrs += 0;
            } /*intersection is outside sc2*/
            else {
                /* vectors from center os sc1 to intersection with infinite cylinder*/
                vec1 = cigar2.scdir * x1 - r_cm;
                e = cigar1.scdir.dot(vec1);
                if ((e >= cigar1.half_length) || (e <= -cigar1.half_length)) {
                    intrs += 0;
                } /*intersection is outside sc1*/
                else {
                    intrs += test_intrpatch(cigar1, vec1, cigar1.pcanglsw, x1, intersections);
                }
            }
            if (d > 0) {
                x2 = (-b - sqrt(d)) * 0.5 / a; /*parameter on line of SC2 determining intersection*/
                if ((x2 >= cigar2.half_length) || (x2 <= -cigar2.half_length)) {
                    intrs += 0; /*intersection is outside sc2*/
                } else {
                    vec2 = cigar2.scdir * x2 - r_cm;
                    e = cigar1.scdir.dot(vec2);
                    if ((e >= cigar1.half_length) || (e <= -cigar1.half_length)) {
                        intrs += 0;
                    } /*intersection is outside sc1*/
                    else {
                        intrs += test_intrpatch(cigar1, vec2, cigar1.pcanglsw, x2, intersections);
                    }
                }
            }
        }
    }

    /*1c- test intersection with plates at the end - it is at distace C and in wedge*/
    if (intrs < 2) {
        a = cigar1.scdir.dot(cigar2.scdir);
        if (a == 0.0) {
            intrs = 0;
        } /* there is no intersection plane and sc are paralel*/
        else {
            /*plane cap1*/
            vec1 = r_cm + cigar1.half_length * cigar1.scdir;
            x1 = cigar1.scdir.dot(vec1) / a; /*parameter on line of SC2 determining intersection*/
            if ((x1 > cigar2.half_length) || (x1 < -cigar2.half_length)) {
                intrs += 0;
            } /* there is no intersection plane sc is too short*/
            else {
                vec2 = x1 * cigar2.scdir - vec1; /*vector from ENDPOINT to intersection point */
                b = vec2.dot(vec2);
                if (b > cutoff_squared) {
                    intrs += 0;
                } /* the intersection is outside sc */
                else {
                    intrs += test_intrpatch(cigar1, vec2, cigar1.pcanglsw, x1, intersections);
                }
            }
            /*plane cap2*/
            vec1 = r_cm - cigar1.half_length * cigar1.scdir;
            x2 = cigar1.scdir.dot(vec1) / a; /*parameter on line of SC2 determining intersection*/
            if ((x2 > cigar2.half_length) || (x2 < -cigar2.half_length)) {
                intrs += 0;
            } /* there is no intersection plane sc is too short*/
            else {
                vec2 = x2 * cigar2.scdir - vec1; /*vector from ENDPOINT to intersection point */
                b = vec2.dot(vec2);
                if (b > cutoff_squared) {
                    intrs += 0;
                } /* the intersection is outside sc */
                else {
                    intrs += test_intrpatch(cigar1, vec2, cigar1.pcanglsw, x2, intersections);
                }
            }
        }
    }
    /*1d- if there is only one itersection shperocylinder ends within patch wedge
      set as second intersection end inside patch*/
    if (intrs < 2) {
        /*whole spherocylinder is in or all out if intrs ==0*/
        vec1 = cigar2.scdir * cigar2.half_length - r_cm;
        /*vector from CM of sc1 to end of sc2*/
        /*check is is inside sc1*/
        a = vec1.dot(cigar1.scdir);
        vec3 = vec1 - cigar1.scdir * a;
        b = vec3.dot(vec3);
        d = fabs(a) - cigar1.half_length;
        if (d <= 0) { /*is in cylindrical part*/
            /*c is distance squared from line or end to test if is inside sc*/
            if (b < cutoff_squared) {
                intrs += test_intrpatch(cigar1, vec1, cigar1.pcanglsw, cigar2.half_length, intersections);
            }
        }
        if (intrs < 2) {
            vec2 = -cigar2.scdir * cigar2.half_length - r_cm;
            /*check is is inside sc1*/
            a = vec2.dot(cigar1.scdir);
            vec4 = vec2 - cigar1.scdir * a;
            b = vec4.dot(vec4);
            d = fabs(a) - cigar1.half_length;
            if (d <= 0) {
                /*c is distance squared from line or end to test if is inside sc*/
                if (b < cutoff_squared) {
                    intrs += test_intrpatch(cigar1, vec2, cigar1.pcanglsw, -1.0 * cigar2.half_length, intersections);
                }
            }
        }
    }
    return intrs;
}
} // namespace Faunus::SpheroCylinder

namespace Faunus::pairpotential {

HardSpheroCylinder::HardSpheroCylinder()
    : PairPotential("hardspherocylinder", "", false) {}

void HardSpheroCylinder::to_json([[maybe_unused]] json& j) const {}

void HardSpheroCylinder::from_json([[maybe_unused]] const json& j) {}

/**
 * @param particle1 First PSC particle
 * @param particle2 Second PSC particle
 * @param center_separation
 * @return
 */
template <pairpotential::RequirePairPotential PatchPotential, pairpotential::RequirePairPotential CylinderPotential>
double CigarWithCigar<PatchPotential, CylinderPotential>::patchyPatchyEnergy(
    const Particle& particle1, const Particle& particle2,
    const Point& center_separation) const { // patchy sc with patchy sc
    const auto cutoff_squared = patch_potential.cutOffSquared(particle1.id, particle1.id);
    std::array<double, 5> intersections;

    // distance for repulsion
    const auto rclose_squared =
        SpheroCylinder::mindist_segment2segment(particle1.ext->scdir, particle1.ext->half_length, particle2.ext->scdir,
                                                particle2.ext->half_length, center_separation)
            .squaredNorm();
    // 1- do intersections of spherocylinder2 with patch of spherocylinder1 at.
    //  cut distance C
    int intrs = 0;
    intersections.fill(0.0);
    if (particle1.traits().sphero_cylinder.type == SpheroCylinderData::PatchType::Full) {
        intrs = SpheroCylinder::psc_intersect(particle1.getExt(), particle2.getExt(), center_separation, intersections,
                                              cutoff_squared);
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
        intrs = SpheroCylinder::psc_intersect(particle2.getExt(), particle1.getExt(), -center_separation, intersections,
                                              cutoff_squared);
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
    const auto S1 = intersections[0]; // points on sc1
    const auto S2 = intersections[1];

    // 3- scaling function1: dependence on the length of intersetions
    const auto v1 = fabs(S1 - S2) * 0.5;
    const auto v2 = fabs(T1 - T2) * 0.5;
    const auto f0 = v1 + v2;
    // 4a- with two intersection pices calculate vector between their CM
    //-this is for angular orientation
    Point vec1 = particle1.ext->scdir * (S1 + S2) * 0.5;
    Point vec2 = particle2.ext->scdir * (T1 + T2) * 0.5;
    Point vec_intrs = vec2 - vec1 - center_separation; // vec_intrs should be from sc1 t sc2

    // 5- scaling function2: angular dependence of patch1
    vec1 = SpheroCylinder::vec_perpproject(vec_intrs, particle1.ext->scdir);
    vec1.normalize();
    auto s = vec1.dot(particle1.ext->patchdir);
    const auto f1 = SpheroCylinder::fanglscale(s, particle1.getExt());

    // 6- scaling function3: angular dependence of patch2
    vec1 = SpheroCylinder::vec_perpproject(-vec_intrs, particle2.ext->scdir).normalized();
    s = vec1.dot(particle2.ext->patchdir);
    const auto f2 = SpheroCylinder::fanglscale(s, particle2.getExt());

    // 7 - calculate closest distance attractive energy from it
    auto vec_mindist =
        SpheroCylinder::mindist_segment2segment(particle1.ext->scdir, v1, particle2.ext->scdir, v2, vec_intrs);
    auto ndistsq = vec_mindist.dot(vec_mindist);

    // 8- put it all together and output scale
    return f0 * f1 * f2 * patch_potential(particle1, particle2, ndistsq, Point::Zero()) +
           cylinder_potential(particle1, particle2, rclose_squared, Point::Zero());
}

template <pairpotential::RequirePairPotential PatchPotential, pairpotential::RequirePairPotential CylinderPotential>
double CigarWithCigar<PatchPotential, CylinderPotential>::isotropicIsotropicEnergy(
    const Particle& particle1, const Particle& particle2,
    const Point& center_separation) const { // isotropic sc with isotropic sc
    const auto mindist =
        SpheroCylinder::mindist_segment2segment(particle1.ext->scdir, particle1.ext->half_length, particle2.ext->scdir,
                                                particle2.ext->half_length, center_separation)
            .squaredNorm();
    return patch_potential(particle1, particle2, mindist, Point::Zero()) +
           cylinder_potential(particle1, particle2, mindist, Point::Zero());
}

template <RequirePairPotential PatchPotential, RequirePairPotential CylinderPotential>
void CigarWithCigar<PatchPotential, CylinderPotential>::to_json(json& j) const {
    j["patch"] = static_cast<json>(patch_potential);
    j["cylinder"] = static_cast<json>(cylinder_potential);
}

template <RequirePairPotential PatchPotential, RequirePairPotential CylinderPotential>
void CigarWithCigar<PatchPotential, CylinderPotential>::from_json(const json& j) {
    pairpotential::from_json(j, patch_potential);
    pairpotential::from_json(j, cylinder_potential);
}

template <RequirePairPotential PatchPotential, RequirePairPotential CylinderPotential>
CigarWithCigar<PatchPotential, CylinderPotential>::CigarWithCigar()
    : PairPotential("cigar-cigar", ""s, false) {}

template class CigarWithCigar<CosAttractMixed, WeeksChandlerAndersen>; // explicit initialization

TEST_CASE("[Faunus] CigarWithCigar") {
    using CigarCigar = CigarWithCigar<CosAttractMixed, WeeksChandlerAndersen>;

    // Check that we use same temperature as R. Vacha. Note also that the epsilon value
    // in his code is in kT per nanometer.
    CHECK_EQ(1.4_kT / 1.0_kJmol, doctest::Approx(3.47056));

    // the Faunus epsilon value is in kJ/mol/Å
    Faunus::atoms = R"([
             { "A": { "sigma": 10.0, "eps": 0.347056, "rc": 11.2246205, "wc": 10.0, "psc": {"length": 60, "type": "capped", "patch_angle": 80, "patch_angle_switch": 5.0} } },
             { "B": { "sigma": 10.0, "eps": 0.347056 } }
             ])"_json.get<decltype(atoms)>();

    Particle a;
    Particle b;
    a.createExtension();
    b.createExtension();
    a = Faunus::atoms.at(0);
    b = Faunus::atoms.at(0);

    SUBCASE("maximum contact attraction") {
        // place two parallel CPSC in contact and with patches facing each other
        auto pairpot = CigarCigar(R"( {"cos2": {"mixing": "LB"}, "wca": {"mixing": "LB"}} )"_json);
        const auto& atom_data = Faunus::atoms.at(0); // atom type "A"
        Point psc_dir = {0.0, 1.0, 0.0};             // ↑
        Point patch_dir = {1.0, 0.0, 0.0};           // →
        a.getExt().setDirections(atom_data.sphero_cylinder, psc_dir, patch_dir);
        b.getExt().setDirections(atom_data.sphero_cylinder, psc_dir, -patch_dir);
        a.pos = {0.0, 0.0, 0.0};
        b.pos = {atom_data.sigma, 0.0, 0.0}; // place a and b at contact
        Point separation = a.pos - b.pos;
        CHECK_EQ(pairpot(a, b, separation.squaredNorm(), separation), doctest::Approx(-8.26)); // kT
    }
}

// -----------------------------

template <RequirePairPotential PatchPotential, RequirePairPotential CylinderPotential,
          RequirePairPotential SphereWithSphere>
void CompleteCigarPotential<PatchPotential, CylinderPotential, SphereWithSphere>::to_json(json& j) const {
    j = {sphere_sphere, cigar_cigar, cigar_sphere};
}

template <RequirePairPotential PatchPotential, RequirePairPotential CylinderPotential,
          RequirePairPotential SphereWithSphere>
void CompleteCigarPotential<PatchPotential, CylinderPotential, SphereWithSphere>::from_json(const json& j) {
    pairpotential::from_json(j, sphere_sphere);
    pairpotential::from_json(j, cigar_cigar);
    pairpotential::from_json(j, cigar_sphere);
}

template <RequirePairPotential PatchPotential, RequirePairPotential CylinderPotential,
          RequirePairPotential SphereWithSphere>
CompleteCigarPotential<PatchPotential, CylinderPotential, SphereWithSphere>::CompleteCigarPotential()
    : PairPotential("complete cigar", ""s, false) {}

template class CompleteCigarPotential<CosAttractMixed, WeeksChandlerAndersen>; // explicit initialization

// -------------------------------

template <RequirePairPotential PatchPotential, RequirePairPotential CylinderPotential>
CigarWithSphere<PatchPotential, CylinderPotential>::CigarWithSphere()
    : PairPotential("cigar-sphere", ""s, false) {}

template <RequirePairPotential PatchPotential, RequirePairPotential CylinderPotential>
void CigarWithSphere<PatchPotential, CylinderPotential>::to_json(json& j) const {
    j["patch"] = static_cast<json>(patch_potential);
    j["cylinder"] = static_cast<json>(cylinder_potential);
}

template <RequirePairPotential PatchPotential, RequirePairPotential CylinderPotential>
void CigarWithSphere<PatchPotential, CylinderPotential>::from_json(const json& j) {
    pairpotential::from_json(j, patch_potential);
    pairpotential::from_json(j, cylinder_potential);
}

template class CigarWithSphere<CosAttractMixed, WeeksChandlerAndersen>; // explicit initialization

} // namespace Faunus::Potential