#ifndef FAU_BONDED_H
#define FAU_BONDED_H

#include <faunus/json.h>
#include "faunus/textio.h"

namespace Faunus {

  /**
   * @brief Namespace for bonded interactions
   */
  namespace Bonded {

    typedef picojson::object::value_type Tjson;

    /**
     * @brief Base class for two- and manybody potentials
     */
    struct BondDataBase {

      vector<int> index; //!< Particle index

      inline BondDataBase() {};

      bool find(int i) const { //!< Determine if particle exists in `index`
        return (std::find(index.begin(), index.end(), i)!=index.end()) ? true : false;
      }
    };

    /**
     * @brief Data for harmonic bond
     *
     * The constructor accepts a json object with the following format:
     *
     * ~~~~{.js}
     *   "0 2" : { "k":0.5, "req":4.0 }
     * ~~~~
     *
     * where the first filed specifies the two particle index; `k` (kT) and `req` (angstrom)
     * are the force constant and equilibrium distance, respectively.
     */
    struct HarmonicBondData : public BondDataBase {
      double k;   //!< Force constant
      double req; //!< Equilibrium distance

      inline HarmonicBondData( int i, int j, double k, double req ) : k(k), req(req) { index = {0,1}; }

      inline HarmonicBondData( const Tjson &d ) {
        index = textio::words2vec<int>( d.first );
        assert( index.size() == 2);
        k =   json::value<double>( d.second, "k", 0 );
        req = json::value<double>( d.second, "req", 0 );
      }

      /** @brief Write to stream */
      friend std::ostream &operator<<(std::ostream &o, const HarmonicBondData &b) {
        o << b.index[0] << " " << b.index[1] << " " << b.k << " " << b.req;
        return o;
      }

      /** @brief Read from stream */
      HarmonicBondData& operator<<(std::istream &in) {
        in >> index[0] >> index[1] >> k >> req;
        return *this;
      }

    };

    /**
     * @brief Example structure for a four body potential
     *
     * The constructor accepts a json object with the following format:
     *
     * ~~~~{.js}
     *   "0 2 3 4" : { "k1":0.5, "k2":4.0 }
     * ~~~~
     *
     * where the first filed specifies the two particle index; `k` (kT) and `req` (angstrom)
     * are the force constant and equilibrium distance, respectively.
     * 
     * @todo Under construction
     */
    struct DihedralData : public BondDataBase {
      double k1;  //!< Some constant
      double k2;  //!< ...andother constant

      inline DihedralData() { index={-1,-1,-1,-1}; };

      inline DihedralData( const Tjson &d ) {
        index = textio::words2vec<int>( d.first );
        assert( index.size() == 4);
        k1 =   json::value<double>( d.second, "k1", 0 );
        k2 = json::value<double>( d.second, "k2", 0 );
      }

      /** @brief Write to stream */
      friend std::ostream &operator<<(std::ostream &o, const DihedralData &b) {
        for (auto i : b.index)
          o << i << " ";
        o << b.k1 << " " << b.k2;
        return o;
      }

      /** @brief Read from stream */
      DihedralData& operator<<(std::istream &in) {
        index.resize(4);
        for (auto &i : index)
          in >> i;
        in >> k1 >> k2;
        return *this;
      }
    };

  }//namespace
}//namespace

#endif
