#pragma once
#include "space.h"

namespace Faunus {

    namespace Analysis {

        /** @brief Example analysis */
        template<class T, class Enable = void>
            struct analyse {
                void sample(T &p) {
                    std::cout << "not a dipole!" << std::endl;
                } //!< Sample
            }; // primary template

        /** @brief Example analysis */
        template<class T>
            struct analyse<T, typename std::enable_if<std::is_base_of<Dipole, T>::value>::type> {
                void sample(T &p) {
                    std::cout << "dipole!" << std::endl;
                } //!< Sample
            }; // specialized template

    }//namespace

}//namespace
