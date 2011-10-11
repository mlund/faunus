#ifndef FAU_BONDED
#define FAU_BONDED

#include <faunus/common.h>
#include <faunus/textio.h>
#include <faunus/geometry.h>

namespace Faunus {
  namespace Energy {

    /*!
     * \brief Harmonic bond energy
     */
    class HarmonicBond {
      public:
        double req; //!< Equilibrium distance (AA)
        double k;   //!< Force constant (kT/AA^2)
        HarmonicBond(double=0, double=0);
        double energy(double) const;
        string info() const;
    };

    /*!
     * \brief General class for handling pairs of particles
     * \date Lund, October 2011
     */
    template<class Tpairtype>
      class ParticlePairs {
        protected:
          std::map<int, std::map<int, Tpairtype> > list;
          string name;
        public:
          void add(int i, int j, const Tpairtype &p) {
            assert(i!=j);   //debug only
            if (i!=j) {
              list[i][j]=p;
              list[j][i]=p; //speed over space!
            }
          } 
          //!< \brief Retrieve bond via (i,j) operator
          //!< \warning Avoid using this to modify bond properties as both (i,j) and (j,i) need to be set.
          //!< \todo Should we rather return copy instead of reference (safer, slower)
          Tpairtype& operator() (int i, int j) {
            return list[i][j];
          }
          string info() {
            using namespace Faunus::textio;
            std::ostringstream o;
            o << header(name) << indent(SUB) << std::left
              << setw(7) << "i" << setw(7) << "j" << endl;
            for (auto i : list)
              for (auto j : list[i.first])
                if (i.first < j.first)
                  o << indent(SUB) << std::left << setw(7) << i.first
                    << setw(7) << j.first << j.second.info() << endl;
            return o.str();
          }
      };

    /*!
     * \brief Class for handling bond pairs
     * \date Lund, October 2011
     */
    template<class Tbondtype=HarmonicBond>
      class ParticleBonds : public ParticlePairs<Tbondtype> {
        typedef ParticlePairs<Tbondtype> pairs;
        public:
        ParticleBonds() { pairs::name+="Particle Bonds"; }
        //!< Returns total bond energy of i'th particle (kT)
        double bondEnergy(Geometry::Geometrybase &geo, const p_vec &p, int i) {
          double u=0;
          for (auto &j : pairs::list[i])
            u+=j.second.energy( geo.dist( p[i], p[j.first] ) );
          return u;
        }
      };

  } // namespace
} // namespace
#endif
