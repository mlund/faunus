#include "spherocylinder.h"

Point Faunus::SpheroCylinder::mindist_segment2segment(const Point& dir1, double halfl1, const Point& dir2,
                                                      double halfl2, const Point& r_cm) {
    Point u = dir1 * (halfl1 * 2);                  // S1.P1 - S1.P0;
    Point v = dir2 * (halfl2 * 2);                  // S2.P1 - S2.P0;
    Point w = dir2 * halfl2 - dir1 * halfl1 - r_cm; // S1.P0-S2.P0
    auto a = u.dot(u);                              // always >= 0
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
    if (D < 0.00000001) { // the lines are almost parallel
        sN = 0.0;         // force using point P0 on segment S1
        sD = 1.0;         // to prevent possible division by 0.0 later
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
        if (-d < 0.0)
            sN = 0.0;
        else if (-d > a)
            sN = sD;
        else {
            sN = -d;
            sD = a;
        }
    } else if (tN > tD) { // tc > 1 => the t=1 edge is visible
        tN = tD;
        // recompute sc for this edge
        if ((-d + b) < 0.0)
            sN = 0;
        else if ((-d + b) > a)
            sN = sD;
        else {
            sN = (-d + b);
            sD = a;
        }
    }
    // finally do the division to get sc and tc
    if (std::fabs(sN) < 0.00000001)
        sc = 0.0;
    else
        sc = sN / sD;
    if (std::fabs(tN) < 0.00000001)
        tc = 0.0;
    else
        tc = tN / tD;

    // get the difference of the two closest points
    return u * sc + w - v * tc;
}
int Faunus::SpheroCylinder::find_intersect_plane(const Faunus::Cigar& part1, const Faunus::Cigar& part2,
                                                 const Faunus::Point& r_cm, const Faunus::Point& w_vec, double rcut2,
                                                 double cospatch, double* intersections) {
    int intrs;
    double disti;

    Point nplane = part1.scdir.cross(w_vec);
    nplane.normalize();
    auto a = nplane.dot(part2.scdir);
    if (a == 0.0)
        intrs = 0; /* there is no intersection plane and sc are paralel*/
    else {
        auto ti = nplane.dot(r_cm) / a;
        if ((ti > part2.halfl) || (ti < -part2.halfl))
            intrs = 0; /* there is no intersection plane sc is too short*/
        else {
            Point d_vec = ti * part2.scdir - r_cm; /*vector from intersection point to CM*/
            auto c = d_vec.dot(w_vec);
            if (c * cospatch < 0)
                intrs = 0; /* the intersection in plane is on other side of patch */
            else {
                double d = fabs(d_vec.dot(part1.scdir)) - part2.halfl;
                if (d <= 0)
                    disti = c * c; /*is inside cylinder*/
                else
                    disti = d * d + c * c; /*is inside patch*/
                if (disti > rcut2)
                    intrs = 0; /* the intersection is outside sc */
                else {
                    intrs = 1;
                    int i = 0;
                    while (intersections[i] != 0) {
                        if (ti == intersections[i])
                            intrs = 0; /* found intersection we already have -it is at boundary*/
                        i++;
                    }
                    if (intrs > 0) {
                        intersections[i] = ti;
                    }
                }
            }
        }
    }
    return intrs;
}
int Faunus::SpheroCylinder::test_intrpatch(const Faunus::Cigar& part1, Faunus::Point& vec, double cospatch, double ti,
                                           double* intersections) {
    /*test if we have intersection*/
    /* do projection to patch plane*/
    vec = vec_perpproject(vec, part1.scdir);
    vec.normalize();
    /* test angle distance from patch*/
    auto a = part1.patchdir.dot(vec);
    int intrs = 0;
    if (a >= cospatch) {
        intrs = 1;
        int i = 0;
        while (intersections[i] != 0) {
            if (ti == intersections[i])
                intrs = 0; /* found intersection we already have -it is at boundary*/
            i++;
        }
        if (intrs > 0)
            intersections[i] = ti;
    }
    return intrs;
}
int Faunus::SpheroCylinder::find_intersect_planec(const Faunus::Cigar& part1, const Faunus::Cigar& part2,
                                                  const Faunus::Point& r_cm, const Faunus::Point& w_vec, double rcut2,
                                                  double cospatch, double* intersections) {
    Point nplane = part1.scdir.cross(w_vec).normalized();
    auto a = nplane.dot(part2.scdir);
    if (a == 0.0) { // @todo unstable comparison
        return 0;   // there is no intersection plane and sc are paralel
    }
    auto ti = nplane.dot(r_cm) / a;
    if ((ti > part2.halfl) || (ti < -part2.halfl)) {
        return 0; // there is no intersection plane sc is too short
    }
    Point d_vec = ti * part2.scdir - r_cm; /*vector from intersection point to CM*/
    auto c = d_vec.dot(w_vec);
    if (c * cospatch < 0) {
        return 0; /* the intersection in plane is on other side of patch */
    }
    int intrs = 0;
    auto d = fabs(d_vec.dot(part1.scdir)) - part2.halfl;
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
int Faunus::SpheroCylinder::psc_intersect(const Faunus::Cigar& part1, const Faunus::Cigar& part2,
                                          const Faunus::Point& r_cm, double* intersections, double rcut2) {
    double a, b, c, d, e, x1, x2;
    Point cm21, vec1, vec2, vec3, vec4;

    /*1- do intersections of spherocylinder2 with patch of spherocylinder1 at
      cut distance C*/
    /*1a- test intersection with half planes of patch and look how far they are
      from spherocylinder. If closer then C  we got itersection*/

    /* plane1 */
    /* find intersections of part2 with plane by par1 and patchsides[0] */
    int intrs = 0;
    intrs += find_intersect_plane(part1, part2, r_cm, part1.patchsides[0], rcut2, part1.pcanglsw, intersections);
    /* plane2 */
    /* find intersections of part2 with plane by par1 and patchsides[1] */
    intrs += find_intersect_plane(part1, part2, r_cm, part1.patchsides[1], rcut2, part1.pcanglsw, intersections);

    if ((intrs == 2) && (part1.pcanglsw < 0)) {
        assert("Patch>180 -> two segments.");
    }

    /*1b- test intersection with cylinder - it is at distance C*/
    if (intrs < 2) {
        cm21 = -r_cm;
        vec1 = cm21.cross(part1.scdir);
        vec2 = part2.scdir.cross(part1.scdir);
        a = vec2.dot(vec2);
        b = 2 * vec1.dot(vec2);
        c = -rcut2 + vec1.dot(vec1);
        d = b * b - 4 * a * c;
        if (d >= 0) {                      /*there is intersection with infinite cylinder */
            x1 = (-b + sqrt(d)) * 0.5 / a; /*parameter on line of SC2 determining intersection*/
            if ((x1 >= part2.halfl) || (x1 <= -part2.halfl)) {
                intrs += 0;
            } /*intersection is outside sc2*/
            else {
                /* vectors from center os sc1 to intersection with infinite cylinder*/
                vec1 = part2.scdir * x1 - r_cm;
                e = part1.scdir.dot(vec1);
                if ((e >= part1.halfl) || (e <= -part1.halfl)) {
                    intrs += 0;
                } /*intersection is outside sc1*/
                else {
                    intrs += test_intrpatch(part1, vec1, part1.pcanglsw, x1, intersections);
                }
            }
            if (d > 0) {
                x2 = (-b - sqrt(d)) * 0.5 / a; /*parameter on line of SC2 determining intersection*/
                if ((x2 >= part2.halfl) || (x2 <= -part2.halfl)) {
                    intrs += 0;
                } /*intersection is outside sc2*/
                else {
                    vec2 = part2.scdir * x2 - r_cm;
                    e = part1.scdir.dot(vec2);
                    if ((e >= part1.halfl) || (e <= -part1.halfl))
                        intrs += 0; /*intersection is outside sc1*/
                    else {
                        intrs += test_intrpatch(part1, vec2, part1.pcanglsw, x2, intersections);
                    }
                }
            }
        }
    }
    /*1c- test intersection with spheres at the end - it is at distace C*/
    if (intrs < 2) {
        /*centers of spheres*/
        /*relative to the CM of sc2*/
        vec1 = part1.scdir * part1.halfl - r_cm;
        vec2 = -part1.scdir * part1.halfl - r_cm;

        /*sphere1*/
        a = part2.scdir.dot(part2.scdir);
        b = 2.0 * vec1.dot(part2.scdir);
        c = vec1.dot(vec1) - rcut2;
        d = b * b - 4 * a * c;
        if (d >= 0) {                      /*if d<0 there are no intersections*/
            x1 = (-b + sqrt(d)) * 0.5 / a; /*parameter on line of SC2 determining intersection*/
            if ((x1 >= part2.halfl) || (x1 <= -part2.halfl)) {
                intrs += 0;
            } /*intersection is outside sc2*/
            else {
                vec3 = part2.scdir * x1 - r_cm;
                e = part1.scdir.dot(vec3);
                if ((e >= part1.halfl) || (e <= -part1.halfl)) { /*if not intersection is inside sc1*/
                    intrs += test_intrpatch(part1, vec3, part1.pcanglsw, x1, intersections);
                }
            }
            if (d > 0) {
                x2 = (-b - sqrt(d)) * 0.5 / a; /*parameter on line of SC2 determining intersection*/
                if ((x2 >= part2.halfl) || (x2 <= -part2.halfl)) {
                    intrs += 0;
                } /*intersection is outside sc2*/
                else {
                    vec4 = part2.scdir * x2 - r_cm;
                    e = part1.dir.dot(vec4);
                    if ((e >= part1.halfl) || (e <= -part1.halfl)) { /*if not intersection is inside sc1*/
                        intrs += test_intrpatch(part1, vec4, part1.pcanglsw, x2, intersections);
                    }
                }
            }
        }
        /*sphere2*/
        a = part2.scdir.dot(part2.scdir);
        b = 2.0 * vec2.dot(part2.scdir);
        c = vec2.dot(vec2) - rcut2;
        d = b * b - 4 * a * c;
        if (d >= 0) {                      /*if d<0 there are no intersections*/
            x1 = (-b + sqrt(d)) * 0.5 / a; /*parameter on line of SC2 determining intersection*/
            if ((x1 >= part2.halfl) || (x1 <= -part2.halfl)) {
                intrs += 0;
            } /*intersection is outside sc2*/
            else {
                vec3 = part2.scdir * x1 - r_cm;
                e = part1.scdir.dot(vec3);
                if ((e >= part1.halfl) || (e <= -part1.halfl)) { /*if not intersection is inside sc1*/
                    intrs += test_intrpatch(part1, vec3, part1.pcanglsw, x1, intersections);
                }
            }
            if (d > 0) {
                x2 = (-b - sqrt(d)) * 0.5 / a; /*parameter on line of SC2 determining intersection*/
                if ((x2 >= part2.halfl) || (x2 <= -part2.halfl)) {
                    intrs += 0;
                } /*intersection is outside sc2*/
                else {
                    vec4 = part2.scdir * x2 - r_cm;
                    e = part1.scdir.dot(vec4);
                    if ((e >= part1.halfl) || (e <= -part1.halfl)) { /*if not intersection is inside sc1*/
                        intrs += test_intrpatch(part1, vec4, part1.pcanglsw, x2, intersections);
                    }
                }
            }
        }
    }

    /*1d- if there is only one itersection shperocylinder ends within patch wedge
      set as second intersection end inside patch*/
    if (intrs < 2) {
        /*whole spherocylinder is in or all out if intrs ==0*/
        vec1 = part2.scdir * part2.halfl - r_cm;
        /*vector from CM of sc1 to end of sc2*/
        /*check is is inside sc1*/
        a = vec1.dot(part1.scdir);
        vec3 = vec1 - part1.scdir * a;
        b = vec3.dot(vec3);
        d = fabs(a) - part1.halfl;
        if (d <= 0)
            c = b; /*is inside cylindrical part*/
        else
            c = d * d + b; /*is inside caps*/
        /*c is distance squared from line or end to test if is inside sc*/
        if (c < rcut2)
            intrs += test_intrpatch(part1, vec1, part1.pcanglsw, part2.halfl, intersections);
        if (intrs < 2) {
            vec2 = -part2.scdir * part2.halfl - r_cm;
            /*check is is inside sc1*/
            a = vec2.dot(part1.scdir);
            vec4 = vec2 - part1.scdir * a;
            b = vec4.dot(vec4);
            d = fabs(a) - part1.halfl;
            if (d <= 0)
                c = b; /*is inside cylindrical part*/
            else
                c = d * d + b; /*is inside caps*/
            /*c is distance squared from line or end to test if is inside sc*/
            if (c < rcut2)
                intrs += test_intrpatch(part1, vec2, part1.pcanglsw, -1.0 * part2.halfl, intersections);
        }
    }
    return intrs;
}
int Faunus::SpheroCylinder::cpsc_intersect(const Faunus::Cigar& part1, const Faunus::Cigar& part2,
                                           const Faunus::Point& r_cm, double* intersections, double rcut2) {
    int intrs;
    double a, b, c, d, e, x1, x2;
    Point cm21, vec1, vec2, vec3, vec4;

    intrs = 0;
    /*1- do intersections of spherocylinder2 with patch of spherocylinder1 at
      cut distance C*/
    /*1a- test intersection with half planes of patch and look how far they are
      from spherocylinder. If closer then C  we got itersection*/

    /* plane1 */
    /* find intersections of part2 with plane by par1 and part1.patchsides[0] */
    intrs += find_intersect_planec(part1, part2, r_cm, part1.patchsides[0], rcut2, part1.pcanglsw, intersections);
    /* plane2 */
    /* find intersections of part2 with plane by par1 and part1.patchsides[1] */
    intrs += find_intersect_planec(part1, part2, r_cm, part1.patchsides[1], rcut2, part1.pcanglsw, intersections);

    if ((intrs == 2) && (part1.pcanglsw < 0))
        throw std::runtime_error("Patch larger than 180 deg -> two segments - this is not yet implemented");

    /*1b- test intersection with cylinder - it is at distance C*/
    if (intrs < 2) {
        cm21 = -r_cm;
        vec1 = cm21.cross(part1.scdir);
        vec2 = part2.scdir.cross(part1.scdir);
        a = vec2.dot(vec2);
        b = 2 * vec1.dot(vec2);
        c = -rcut2 + vec1.dot(vec1);
        d = b * b - 4 * a * c;
        if (d >= 0) {                      /*there is intersection with infinite cylinder */
            x1 = (-b + sqrt(d)) * 0.5 / a; /*parameter on line of SC2 determining intersection*/
            if ((x1 >= part2.halfl) || (x1 <= -part2.halfl))
                intrs += 0; /*intersection is outside sc2*/
            else {
                /* vectors from center os sc1 to intersection with infinite cylinder*/
                vec1 = part2.scdir * x1 - r_cm;
                e = part1.scdir.dot(vec1);
                if ((e >= part1.halfl) || (e <= -part1.halfl))
                    intrs += 0; /*intersection is outside sc1*/
                else {
                    intrs += test_intrpatch(part1, vec1, part1.pcanglsw, x1, intersections);
                }
            }
            if (d > 0) {
                x2 = (-b - sqrt(d)) * 0.5 / a; /*parameter on line of SC2 determining intersection*/
                if ((x2 >= part2.halfl) || (x2 <= -part2.halfl))
                    intrs += 0; /*intersection is outside sc2*/
                else {
                    vec2 = part2.scdir * x2 - r_cm;
                    e = part1.scdir.dot(vec2);
                    if ((e >= part1.halfl) || (e <= -part1.halfl))
                        intrs += 0; /*intersection is outside sc1*/
                    else {
                        intrs += test_intrpatch(part1, vec2, part1.pcanglsw, x2, intersections);
                    }
                }
            }
        }
    }

    /*1c- test intersection with plates at the end - it is at distace C and in wedge*/
    if (intrs < 2) {
        a = part1.scdir.dot(part2.scdir);
        if (a == 0.0)
            intrs = 0; /* there is no intersection plane and sc are paralel*/
        else {
            /*plane cap1*/
            vec1 = r_cm + part1.halfl * part1.scdir;
            x1 = part1.scdir.dot(vec1) / a; /*parameter on line of SC2 determining intersection*/
            if ((x1 > part2.halfl) || (x1 < -part2.halfl))
                intrs += 0; /* there is no intersection plane sc is too short*/
            else {
                vec2 = x1 * part2.scdir - vec1; /*vector from ENDPOINT to intersection point */
                b = vec2.dot(vec2);
                if (b > rcut2)
                    intrs += 0; /* the intersection is outside sc */
                else {
                    intrs += test_intrpatch(part1, vec2, part1.pcanglsw, x1, intersections);
                }
            }
            /*plane cap2*/
            vec1 = r_cm - part1.halfl * part1.scdir;
            x2 = part1.scdir.dot(vec1) / a; /*parameter on line of SC2 determining intersection*/
            if ((x2 > part2.halfl) || (x2 < -part2.halfl))
                intrs += 0; /* there is no intersection plane sc is too short*/
            else {
                vec2 = x2 * part2.scdir - vec1; /*vector from ENDPOINT to intersection point */
                b = vec2.dot(vec2);
                if (b > rcut2)
                    intrs += 0; /* the intersection is outside sc */
                else {
                    intrs += test_intrpatch(part1, vec2, part1.pcanglsw, x2, intersections);
                }
            }
        }
    }
    /*1d- if there is only one itersection shperocylinder ends within patch wedge
      set as second intersection end inside patch*/
    if (intrs < 2) {
        /*whole spherocylinder is in or all out if intrs ==0*/
        vec1 = part2.scdir * part2.halfl - r_cm;
        /*vector from CM of sc1 to end of sc2*/
        /*check is is inside sc1*/
        a = vec1.dot(part1.scdir);
        vec3 = vec1 - part1.scdir * a;
        b = vec3.dot(vec3);
        d = fabs(a) - part1.halfl;
        if (d <= 0) { /*is in cylindrical part*/
            /*c is distance squared from line or end to test if is inside sc*/
            if (b < rcut2)
                intrs += test_intrpatch(part1, vec1, part1.pcanglsw, part2.halfl, intersections);
        }
        if (intrs < 2) {
            vec2 = -part2.scdir * part2.halfl - r_cm;
            /*check is is inside sc1*/
            a = vec2.dot(part1.scdir);
            vec4 = vec2 - part1.scdir * a;
            b = vec4.dot(vec4);
            d = fabs(a) - part1.halfl;
            if (d <= 0) {
                /*c is distance squared from line or end to test if is inside sc*/
                if (b < rcut2)
                    intrs += test_intrpatch(part1, vec2, part1.pcanglsw, -1.0 * part2.halfl, intersections);
            }
        }
    }
    return intrs;
}
