#ifndef FAU_POT_DEBYEHUCKELP3VEC_H
#define FAU_POT_DEBYEHUCKELP3VEC_H
#include "faunus/potentials/pot_debyehuckelP3.h"
#include "faunus/legendre.h"
#include <Eigen/Dense>

namespace Eigen {
  template<typename Scalar>
    struct invsqrt_op {
      EIGEN_EMPTY_STRUCT_CTOR( invsqrt_op )
      inline const Scalar operator() (const Scalar& a) const { return Scalar(1)/ei_sqrt(a); }
      typedef typename ei_packet_traits<Scalar>::type Packet;
      inline Packet packetOp(const Packet& a) const { return _mm_rsqrt_ps(a); }
    };
  template<typename Scalar>
    struct ei_functor_traits<invsqrt_op<Scalar> > {
      enum {
        Cost = 2 * NumTraits<Scalar>::MulCost,
        PacketAccess = 
#ifdef EIGEN_VECTORIZE_SSE
          true
#else
          false
#endif
      };
    };
}

namespace Faunus {
  /*! \brief Debye-Huckel potential for periodic boundry 
   *         conditions in 3D, it is extended to preform 
   *         under conditions of constant pressure.
   *         See class isobaric->markovmove.h
   *  \author Mikael Lund/Bjoern Persson
   *  \date Lund/Prag 2008
   */
  class pot_debyehuckelP3vec : public pot_debyehuckelP3 {
    using pot_debyehuckelP3::k;
    private:
      Eigen::ArrayXf r2,qq,sigma,u,vdw6;
    
    public:
      pot_debyehuckelP3vec( inputfile &in ) : pot_debyehuckelP3(in) {
        name+=" - vectorized for Eigen";
      }

      double energy(const vector<particle> &p, const vector<int> &pairs) {
        float kappa=k;
        int l=0, n=pairs.size()/2;
        r2.resize(n);
        qq.resize(n);
        sigma.resize(n);
        for (int k=0; k<n; ++k) {
          int i=pairs[l++];
          int j=pairs[l++];
          r2[k] = sqdist( p[i], p[j] );
          qq[k] = p[i].charge*p[j].charge;
          sigma[k]  = p[i].radius+p[j].radius;
        }
        u    = r2.unaryExpr( Eigen::invsqrt_op<float>() ); // -> 1/r (approximate)
        vdw6 = (sigma*u).square();  // -> (s/r)^2
        vdw6 = vdw6*vdw6*vdw6;      // -> (s/r)^6
        u = qq * u * ((-kappa*r2*u).exp());
        return f*eps*(vdw6*vdw6 - vdw6).sum() + f*u.sum();
      }

      string info() {
        std::ostringstream o;
        o << pot_debyehuckelP3::info();
        return o.str();
      }
  };
} // namespace
#endif
