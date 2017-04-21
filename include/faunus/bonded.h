#ifndef FAU_BONDED_H
#define FAU_BONDED_H

#include <faunus/json.h>
#include "faunus/textio.h"
#include <faunus/potentials.h>

namespace Faunus
{

  /**
   * @brief Namespace for bonded interactions
   */
  namespace Bonded
  {

    /**
     * @brief Base class for two- and manybody potentials
     */
    struct BondedBase
    {
        vector<int> index; //!< Particle index

        inline BondedBase() {};

        bool find( int i ) const
        { //!< Determine if particle exists in `index`
            return (std::find(index.begin(), index.end(), i) != index.end()) ? true : false;
        }

        /** @brief Add offset to each index */
        void shift( int offset )
        {
            assert(offset >= 0);
            for ( auto &i : index )
                i += offset;
        }
    };

    /**
     * @brief Data for linear bonds
     *
     * The constructor accepts a json object with the following format:
     *
     * ~~~~{.js}
     *   "0 2" : { "k":0.5, "req":4.0, "type":"harmonic" }
     * ~~~~
     *
     * where the first filed specifies the two particle index; `k` (kT) and `req` (angstrom)
     * are the force constant and equilibrium distance, respectively. By default, the bond
     * type is set to `HARMONIC`.
     */
    struct BondData : public BondedBase
    {

        enum class Type : char { HARMONIC, FENE };

        double k;   //!< Force constant
        double req; //!< Equilibrium distance

        Type type;  //!< Bond type (default: HARMONIC)

        inline BondData( int i, int j, double k, double req, Type type = Type::HARMONIC )
            : k(k), req(req), type(type)
        {
            index = {i, j};
        }

        inline BondData( Tmjson::iterator &it )
        {
            index = textio::words2vec<int>(it.key());
            assert(index.size() == 2);
            k = it.value()["k"] | 0.0;
            req = it.value()["req"] | 0.0;
            string t = it.value()["type"] | string("harmonic");
            if ( t == "harmonic" )
                type = Type::HARMONIC;
        }

        /** @brief Write to stream */
        friend std::ostream &operator<<( std::ostream &o, const BondData &b )
        {
            o << b.index[0] << " " << b.index[1] << " " << b.k << " " << b.req;
            return o;
        }

        /** @brief Read from stream */
        BondData &operator<<( std::istream &in )
        {
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
    struct DihedralData : public BondedBase
    {
        double k1;  //!< Some constant
        double k2;  //!< ...andother constant

        inline DihedralData() { index = {-1, -1, -1, -1}; };

        inline DihedralData( Tmjson::iterator &it )
        {
            index = textio::words2vec<int>(it.key());
            assert(index.size() == 4);
            k1 = it.value()["k1"] | 0.0;
            k2 = it.value()["k2"] | 0.0;
        }

        /** @brief Write to stream */
        friend std::ostream &operator<<( std::ostream &o, const DihedralData &b )
        {
            for ( auto i : b.index )
                o << i << " ";
            o << b.k1 << " " << b.k2;
            return o;
        }

        /** @brief Read from stream */
        DihedralData &operator<<( std::istream &in )
        {
            index.resize(4);
            for ( auto &i : index )
                in >> i;
            in >> k1 >> k2;
            return *this;
        }
    };

  }//namespace
}//namespace

#endif
