#pragma once

#include <random>
#include <iterator>
#include "json.hpp"

namespace Faunus {

    /**
     * Example code:
     *
     * ```{.cpp}
     *     Random r1;                                     // default deterministic seed
     *     Random r2 = json(r1);                          // copy engine state
     *     Random r3 = R"( {"seed" : "hardware"} )"_json; // non-deterministic seed
     *     Random r1.seed();                              // non-deterministic seed
     * ```
     */
    struct Random {
        std::mt19937 engine; //!< Random number engine used for all operations
        std::uniform_real_distribution<double> dist01; //!< Uniform real distribution [0,1)

        Random();
        void seed();
        double operator()(); //!< Double in uniform range [0,1)
        int range(int, int); //!< Integer in uniform range [min:max]

        template<class Titer>
            Titer sample(const Titer &beg, const Titer &end)
            {
                auto i = beg;
                std::advance(i, range(0, std::distance(beg, end) - 1));
                return i;
            } //!< Iterator to random element in container (std::vector, std::map, etc)
    }; //!< Class for handling random number generation

    void to_json(nlohmann::json&, const Random&);   //!< Random to json conversion
    void from_json(const nlohmann::json&, Random&); //!< json to Random conversion

    static Random random; // global instance of Random

#ifdef DOCTEST_LIBRARY_INCLUDED
    TEST_CASE("[Faunus] Random")
    {
        Random slump; // local instance

        CHECK( slump() == random() );

        int min=10, max=0, N=1e6;
        double x=0;
        for (int i=0; i<N; i++) {
            int j = slump.range(0,9);
            if (j<min) min=j;
            if (j>max) max=j;
            x+=j;
        }
        CHECK( min==0 );
        CHECK( max==9 );
        CHECK( std::fabs(x/N) == doctest::Approx(4.5).epsilon(0.01) );

        Random r1 = R"( {"seed" : "hardware"} )"_json; // non-deterministic seed
        Random r2; // default is a deterministic seed
        CHECK( r1() != r2() );
        Random r3 = nlohmann::json(r1); // r1 --> json --> r3
        CHECK( r1() == r3() );

        // check if random_device works
        Random a, b;
        CHECK( a() == b() );
        a.seed();
        b.seed();
        CHECK( a() != b() );
    }
#endif

}
