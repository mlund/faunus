#ifndef FAU_SCATTER_H
#define FAU_SCATTER_H

#include <faunus/common.h>
#include <faunus/inputfile.h>
#include <faunus/point.h>
#include <faunus/geometry.h>
#include <faunus/analysis.h>

namespace Faunus {
  /*!
   * \brief Routines related to scattering
   */
  namespace Scatter {

    /*!
     * \brief Form factor, \c F(q), for a hard sphere of radius \c R.
     */
    class FormFactorSphere {
      private:
        typedef float Tfloat;
        inline Tfloat j1(Tfloat x) const { // spherical Bessel function
          register Tfloat xinv=1/x;
          return xinv*( sin(x)*xinv - cos(x) );
        }
      public:
        /*!
         * \param q q value in inverse angstroms
         * \param a particle to take radius, \c R from.
         * \returns
         * \f$ I(q) = \left [\frac{3}{(qR)^3} \left ( \sin{qR} - qR\cos{qR} \right ) \right ]^2\f$
         */
        inline Tfloat operator()(Tfloat q, const particle &a) const {
          Tfloat qR=q*a.radius;
          qR = 3. / (qR*qR*qR) * ( sin(qR) - qR*cos(qR) );
          return qR*qR;
        }
    };

    /*!
     * \brief Tabulated particle form factors loaded from disk
     * \note doi:10.1186/1471-2105-11-429
     */
    class FormFactorTable {
      private:
        typedef Analysis::Table2D<float,float> Ttable;
        std::map<particle::Tid, Ttable > F;
        void addVariants();//!< Add for atomic variants for already loaded data
      public:
        bool load(string, bool=true); //!< Load atomic F(q) tables from disk

        inline float operator()(float q, const particle &a) {
          assert( ~F.empty() && "Did you forget to load F(q) from disk?");
          assert( F.find(a.id) != F.end() && "F(q) for particle not known!");
          // or should we return largest F(q) if out of table?
          return F[a.id](q);
        }
    };

    template<typename Tgeometry, typename Tformfactor> class DebyeFormula {
      private:
        Tformfactor F; // scattering from a single particle
        Tgeometry geo; // geometry to use for distance calculations
      public:
        std::map<float,Average<float> > I; // Average I(q)

        DebyeFormula(InputMap &in) : geo(in) {}

        void sample(const p_vec &p, float qmin, float qmax, float dq) {
          std::map<float,float> _I;
          int n=(int)p.size();
#pragma omp parallel for reduction (+:I) schedule (dynamic)
          for (int i=0; i<n-1; i++) {
            for (int j=i+1; j<n; j++) {
              float r = geo.dist(p[i],p[j]);
              for (float q=qmin; q<=qmax; q+=dq) {
                float _iq = F(q,p[i]) * F(q,p[j]) * sin(q*r) / (q*r);
#pragma omp critical
                _I[q] += _iq;
              }
            }
          }
          for (auto &i : _I)
            I(i.first)+=2*i.second; // add to average I(q)
        }
    };

  } // end of namespace
} //end of namespace

#endif
