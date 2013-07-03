#ifndef FAU_SCATTER_H
#define FAU_SCATTER_H

#include <faunus/common.h>
#include <faunus/inputfile.h>
#include <faunus/analysis.h>

namespace Faunus {

  /** @brief Routines related to scattering */
  namespace Scatter {

    /** @brief Form factor, `F(q)`, for a hard sphere of radius `R`.
    */
    template<class T=float>
      class FormFactorSphere {
        private:
          T j1(T x) const { // spherical Bessel function
            register T xinv=1/x;
            return xinv*( sin(x)*xinv - cos(x) );
          }
        public:
          /**
           * @param q q value in inverse angstroms
           * @param a particle to take radius, \c R from.
           * @returns
           * @f$I(q)=\left [\frac{3}{(qR)^3}\left (\sin{qR}-qR\cos{qR}\right )\right ]^2@f$
           */
          template<class Tparticle>
            T operator()(T q, const Tparticle &a) const {
              assert(q>0 && a.radius>0 && "Particle radius and q must be positive");
              T qR=q*a.radius;
              qR = 3. / (qR*qR*qR) * ( sin(qR) - qR*cos(qR) );
              return qR*qR;
            }
      };

    /**
     * @brief Unity form factor (q independent)
     */
    template<class T=float>
      struct FormFactorUnity {
        template<class Tparticle>
          T operator()(T q, const Tparticle &a) const {
            return 1;
          }
      };

    /**
     * @brief Tabulated particle form factors loaded from disk
     * @note doi:10.1186/1471-2105-11-429
     */
    template<typename T=float>
      class FormFactorTable {
        private:
          typedef Analysis::Table2D<T,T> Ttable;
          std::map<int, Ttable > F;

          /**
           * @brief Add for atomic variants for already loaded data
           *
           * This will top up the tables with related form factors
           * not loaded from disk. I.e. if you loaded data for "HIS",
           * this will be duplicated for "HHIS" if not already loaded
           * and if it is defined in the Atoms class.
           */
          void addVariants() {
            for (auto &m : F) {            // loop over loaded F(q)
              string mname = atom[m.first].name;
              for (auto &a : atom.list)    // loop over species
                if (a.id!=m.first)         // if non-tabulated species
                  if (a.name.find(mname)!=std::string::npos) // and it is a mutant
                    F[a.id] = F[m.first];  // duplicate!
            }
          }

        public:
          template<class Tparticle>
            T operator()(T q, const Tparticle &a) {
              assert( ~F.empty() && "Did you forget to load F(q) from disk?");
              assert( F.find(a.id) != F.end() && "F(q) for particle not known!");
              // or should we return largest F(q) if out of table?
              return F[a.id](q);
            }

          /**
           * @brief Load atomic F(q) tables from disk
           * Example of the file format, where the first line gives
           * the q-resolution in inverse angstroms:
           * @verbatim
           * l1: 0.015
           * l2: q     ALA   ARG
           * l3: 0.000 8.983 23.527
           * l4: 0.015 8.980 23.516
           * l5: ...
           * @endverbatim
           * @param filename Multi-column file with F(q) for different species
           * @param variants True (default) if atomic variants to loaded
           * date should be generated
           */
          bool load(string filename, bool variants=true) {
            std::ifstream f(filename.c_str());
            if (f) {
              T dq;
              f >> dq; // read resolution from line 1

              // read atom names from line 2
              std::vector<string> namevec;
              string line, name;
              std::getline(f, line);
              std::stringstream ss(line);
              while (ss>>name)
                if (name!="q") {
                  namevec.push_back(name);
                  F[ atom[name].id ] = Ttable(dq,Ttable::XYDATA);
                }

              // read F(q) starting from line 3 until eof
              while (!f.eof()) {
                int rescnt=0;
                T _f, _q;
                std::getline(f, line);
                std::stringstream ss(line);
                ss >> _q;
                while (ss >> _f) {
                  auto id = atom[ namevec.at(rescnt++) ].id;
                  F[id](_q) = _f;
                }
              }
              if (variants)
                addVariants();
              return true;
            }
            return false;
          }

      };

    /**
     * @brief Calculates scattering intensity, I(q) using the Debye formula
     */
    template<class Tformfactor,class Tgeometry=Geometry::Sphere,class T=float>
      class DebyeFormula {
        protected:
          Tformfactor F; // scattering from a single particle
          Tgeometry geo; // geometry to use for distance calculations
        public:
          std::map<T,Average<T> > I; //!< Sampled, average I(q)

          DebyeFormula(InputMap &in) : geo(in) {}

          /**
           * @brief Sample I(q) and add to average
           * @param p Particle vector
           * @param qmin Minimum q-value to sample (1/A)
           * @param qmax Maximum q-value to sample (1/A)
           * @param dq q spacing (1/A)
           */
          template<class Tpvec>
            void sample(const Tpvec &p, T qmin, T qmax, T dq) {
              if (qmin<1e-6)
                qmin=dq;  // assure q>0
              std::map<T,T> _I;
              int n=(int)p.size();
#pragma omp parallel for reduction (+:I) schedule (dynamic)
              for (int i=0; i<n-1; i++) {
                for (int j=i+1; j<n; j++) {
                  T r = geo.dist(p[i],p[j]);
                  for (T q=qmin; q<=qmax; q+=dq) {
                    T _iq = F(q,p[i]) * F(q,p[j]) * sin(q*r) / (q*r);
#pragma omp critical
                    _I[q] += _iq;
                  }
                }
              }
              for (auto &i : _I)
                I[i.first]+=2*i.second/n+1; // add to average I(q)
            }

          void save(string filename) {
            if (!I.empty()) {
              std::ofstream f(filename.c_str());
              if (f) {
                for (auto &i : I)
                  f << i.first << " " << i.second << "\n";
              }
            }
          }
      };

    /**
     * @brief Structor factor calculation
     *
     * This will take a vector of points or particles and calculate
     * the structure factor using different methods described below.
     */
    template<typename T=double>
      class StructureFactor : public DebyeFormula<FormFactorSphere<T>> {
        public:

          typedef DebyeFormula<FormFactorSphere<T>> base;

          StructureFactor(InputMap &in) : base(in) {}

          /**
           * @brief Sample S(q) using a double sum.
           *
           * This will directly sample
           *
           * @f[
           *   S(\mathbf{q})= \frac{1}{N}\left < \sum_i^N\sum_j^N
           *   \exp(-i\mathbf{qr}_{ij}) \right >
           *   =1+\frac{1}{N} \left <
           *   \sum_{i\neq j}\sum_{j\neq i}\cos(\mathbf{qr}_{ij})
           *   \right >
           * @f]
           * The direction of the q vector is randomly generated on a unit
           * sphere, i.e. during a simulation there will be an angular
           * averaging, producing the same result as the Debye formula. The
           * number of directions for each sample event can be set using
           * the `Nq` parameter which defaults to 1. 
           *
           * @param p Particle/point vector
           * @param qmin Minimum q-value to sample (1/A)
           * @param qmax Maximum q-value to sample (1/A)
           * @param dq q spacing (1/A)
           * @param Nq Number of random q direction (default: 1)
           */
          template<class Tpvec>
            void sample(const Tpvec &p, T qmin, T qmax, T dq, int Nq=10) {
              if (qmin<1e-6)
                qmin=dq;  // assure q>0
              int n=(int)p.size();
              for (int k=0; k<Nq; k++) { // random q directions
                std::map<T,T> _I;
                Point qdir;
                qdir.ranunit(slp_global);
                for (int i=0; i<n-1; i++) {
                  for (int j=i+1; j<n; j++) {
                    Point r = base::geo.vdist(p[i],p[j]);
                    for (T q=qmin; q<=qmax; q+=dq)
                      _I[q] += cos( (q*qdir).dot(r) );
                  }
                } // end of particle loop
                for (auto &i : _I)
                  base::I[i.first]+=2*i.second/n+1; // add to average I(q)
              } // end of q averaging 
            }

          /**
           * @brief Single sum evaluation of S(q)
           *
           * @f[ S(\mathbf{q}) = \frac{1}{N} \left <
           *    \left ( \sum_i^N \sin(\mathbf{qr}_i) \right )^2 +
           *    \left ( \sum_j^N \cos(\mathbf{qr}_j) \right )^2
           *   \right >
           * @f]
           *
           * Angulalar averaging is done as in the `sample()` function.
           * For more information, see <http://doi.org/d8zgw5>
           *
           * @todo Swap loop order
           */
          template<class Tpvec>
            void sample2(const Tpvec &p, T qmin, T qmax, T dq, int Nq=1) {
              if (qmin<1e-6)
                qmin=dq;  // assure q>0
              int n=(int)p.size();
              for (int k=0; k<Nq; k++) { // random q directions
                Point qdir(1,0,0);
                qdir.ranunit(slp_global);
                std::map<T,T> _cos, _sin;
                for (T q=qmin; q<=qmax; q+=dq) {
                  for (int i=0; i<n; i++) {
                    T qr = (q*qdir).dot(p[i]);
                    _sin[q] += sin(qr);
                    _cos[q] += cos(qr);
                  }
                }
                for (auto &i : _sin) {
                  T q=i.first;
                  base::I[q] += (pow(_sin[q],2) + pow(_cos[q],2))/n;
                }
              } // end of q averaging 
            }

          template<class Tpvec>
            void sample3(const Tpvec &p, T qmin, T qmax, T dq, int Nq=1) {
              if (qmin<1e-6)
                qmin=dq;  // assure q>0
              int n=(int)p.size();
              for (int k=0; k<Nq; k++) { // random q directions
                Point qdir(1,0,0);
                //qdir.ranunit(slp_global);
                std::map<T,T> _cos, _sin;
                for (T q=qmin; q<=qmax; q+=dq) {
                  for (int i=0; i<n-1; i++) {
                    for (int j=i+1; j<n; j++) {
                      T qr = (q*qdir).dot(p[i]);
                      _sin[q] += sin(qr);
                      _cos[q] += cos(qr);
                    }
                  }
                }
                for (auto &i : _sin) {
                  T q=i.first;
                  base::I[q] += pow(_sin[q],2) + pow(_cos[q],2);
                }
              } // end of q averaging 
            }


      };

  } // end of namespace
} //end of namespace

#endif
