#pragma once

#include "core.h"

namespace Faunus {
    namespace Potential {

        struct PairPotentialBase {}; //!< Base for pair-potentials

        template<class T1, class T2>
        struct CombinedPairPotential : public PairPotentialBase {
            T1 first;  //!< First pair potential of type T1
            T2 second; //!< Second pair potential of type T2

            CombinedPairPotential(const T1 &a, const T2 &b) : first(a), second(b) {}

            template<class Tparticle>
            inline double operator()(const Tparticle &a, const Tparticle &b, const Point &r) const {
                return first(a, b, r) + second(a, b, r);
            }
        };

        template<class T1, class T2,
                class = typename std::enable_if<std::is_base_of<PairPotentialBase, T1>::value>::type,
                class = typename std::enable_if<std::is_base_of<PairPotentialBase, T2>::value>::type>
        CombinedPairPotential<T1, T2> &operator+(const T1 &pot1, const T2 &pot2) {
            return *(new CombinedPairPotential<T1, T2>(pot1, pot2));
        } //!< Add two pair potentials

        struct Coulomb : public PairPotentialBase {
            double lB; //!< Bjerrum length
            template<typename... T>
            double operator()(const Particle<T...> &a, const Particle<T...> &b, const Point &r) const {
                return lB * a.charge * b.charge / r.norm();
            }
        };

        struct HardSphere : public PairPotentialBase {
            template<typename... T>
            double operator()(const Particle<T...> &a, const Particle<T...> &b, const Point &r2) const {
                auto d=a.radius+b.radius;
                return (d*d<r2.squaredNorm()) ? pc::infty : 0;
            }
        };
    }//end of namespace Potential
}//end of namespace Faunus
