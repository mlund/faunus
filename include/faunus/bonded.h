#ifndef FAU_BONDED_H
#define FAU_BONDED_H

#include <faunus/common.h>
#include <faunus/textio.h>

namespace Faunus {
  namespace Energy {

    /*!
     * \brief General class for handling pairs of particles
     * \date Lund, October 2011
     *
     * This is a general class for handling properties for pairs. One example is bonds between
     * particles, identified through their particle index. Properties are added through to
     * add() template function which can also handle derived classes of Tpairprop. Upon adding
     * new properties, space is dynamically allocated inside the class for each property object.
     * This memory is freed upon destruction or when calling clear().
     */
    template<class Tpairprop, typename Tij=int>
      class ParticlePairs {
        private:
          vector<Tpairprop*> newPtr;
        protected:
          std::map<Tij, std::map<Tij, Tpairprop*> > list;
          string name;
        public:
          ~ParticlePairs() {
            clear();
          }

          void clear() {
            for (auto &ptr : newPtr)
              delete(ptr);
            newPtr.clear();
            list.clear();
          }

          template<typename Tderived>
            void add(Tij i, Tij j, Tderived p) {
              assert(i!=j); //debug
              if (i!=j) {
                newPtr.push_back( new Tderived(p) ); 
                list[i][j]=newPtr.back();
                list[j][i]=newPtr.back();
              }
            }

          //!< \brief Retrieve reference to bond via (i,j) operator
          Tpairprop& operator() (Tij i, Tij j) {
            assert( list[i][j]!=NULL ); //debug
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
        virtual ~Bondbase();
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
    class ParticleBonds : public ParticlePairs<Bondbase,int> {
      typedef ParticlePairs<Bondbase> pairs;
      public:
      ParticleBonds();

      //!< Returns total bond energy of i'th particle (kT)
      double totalEnergy(Geometry::Geometrybase&, const p_vec&, int);

      double totalEnergy(Geometry::Geometrybase&, const p_vec&, const Group&);

      //!< Returns total bond energy of all bonds (kT)
      double totalEnergy(Geometry::Geometrybase&, const p_vec&);
    };

  } // namespace
} // namespace
#endif
