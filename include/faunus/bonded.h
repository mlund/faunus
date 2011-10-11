#ifndef FAU_BONDED
#define FAU_BONDED

#include <faunus/common.h>
#include <faunus/textio.h>

namespace Faunus {
  namespace Energy {

    /*!
     * \brief General class for handling pairs of particles
     * \date Lund, October 2011
     */
    template<class Tpairtype>
      class ParticlePairs {
        protected:
          std::map<int, std::map<int, Tpairtype*> > list;
          string name;
        public:
          void add(int i, int j, Tpairtype *p) {
            assert(i!=j);   //debug only
            assert(p!=NULL);
            if (i!=j) {
              list[i][j]=p;
              list[j][i]=p;
            }
          } 
          //!< \brief Retrieve bond via (i,j) operator
          //!< \warning Avoid using this to modify bond properties as both (i,j) and (j,i) need to be set.
          //!< \todo Should we rather return copy instead of reference (safer, slower)
          Tpairtype& operator() (int i, int j) {
            assert( list[i][j]!=NULL );
            return *list[i][j];
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
                    << setw(7) << j.first << j.second->info() << endl;
            return o.str();
          }
      };

    class Bondbase {
      public:
        virtual double energy(double) const=0;
        virtual string info() const=0;
    };

    /*!
     * \brief Harmonic bond energy
     */
    class HarmonicBond : public Bondbase {
      public:
        double req; //!< Equilibrium distance (AA)
        double k;   //!< Force constant (kT/AA^2)
        HarmonicBond(double=0, double=0);
        double energy(double) const;
        string info() const;
    };

    /*!
     * \brief Class for handling bond pairs
     * \date Lund, October 2011
     * \author Mikael Lund
     *
     * Example:
     * \code
     *    int i=10, j=11;                      // particle index
     *    Energy::ParticleBonds bonds;
     *    bonds.add(i, j, new Energy::HarmonicBond(0.1,5.0) );
     *    std::cout << bonds.info();
     *    double rij = ... ;                   // distance between i and j
     *    double u = bonds(i,j).energy( rij ); // i j bond energy
     * \endcode
     */
    class ParticleBonds : public ParticlePairs<Bondbase> {
      typedef ParticlePairs<Bondbase> pairs;
      public:
      ParticleBonds();
      //!< Returns total bond energy of i'th particle (kT)
      double bondEnergy(Geometry::Geometrybase&, const p_vec&, int);
    };

  } // namespace
} // namespace
#endif
