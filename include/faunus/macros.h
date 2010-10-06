#ifndef FAU_MACROS_H
#define FAU_MACROS_H

namespace Faunus {

/*!
 * \page macros Macros
 * Macros are provided to loop over particles in vectors and groups:
 *
 * \code
 *    double Ztot=0, Zgrp=0, Utot=0;
 *    vector<particle> pvec;
 *    group grp;
 *
 *    ...
 *
 *    FOR_PARTICLES_IN_VECTOR(p, pvec)
 *      Ztot += p->charge;
 *
 *    FOR_PARTICLES_IN_GROUP(p, pvec, grp)
 *      Zgrp += p->charge;
 *
 *    FOR_PAIRS_IN_VECTOR(pi, pj, pvec)
 *      Utot += energy(pi,pj);
 * \endcode
 */

#define FOR_PARTICLES_IN_VECTOR(p,vec) \
  for (vector<particle>::iterator p=vec.begin(); p<vec.end(); ++p)

#define FOR_PAIRS_IN_VECTOR(pi,pj,vec) \
  vector<particle>::iterator pi=vec.begin(); \
  vector<particle>::iterator pj=pi+1; \
  for (; pi<pj; ++pi) \
    for (pj=pi+1; pj<vec.end(); ++pj)

#define FOR_PARTICLES_IN_GROUP(p,vec,group) \
  for (vector<particle>::iterator p=vec.begin()+group.beg; p<=vec.begin()+group.end; ++p)

}//namespace
#endif
