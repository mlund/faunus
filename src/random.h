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

    extern Random random; // global instance of Random

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

    /**
     * @brief Stores a series of elements with given weight
     *
     * Elements is accessed with `get()` that will
     * randomly pick from the weighted distribution.
     * Add elements with `push_back()`
     * where the default weight is _unity_.
     */
    template<typename T>
        class WeightedDistribution {

            private:
                std::discrete_distribution<> dist;
                std::vector<double> weights;

            public:
                std::vector<T> vec; //!< raw vector of T
                size_t index; //!< index from latest element access via push_back or get
                auto size() const { return vec.size(); }
                bool empty() const { return vec.empty(); }

                void clear() {
                    vec.clear();
                    weights.clear();
                }

                template<typename Tnumber>
                    void setWeight(const std::vector<Tnumber> &w) {
                        if (w.size() not_eq vec.size())
                            throw std::runtime_error("number of weights must match data");
                        weights.resize(w.size());
                        std::copy(w.begin(), w.end(), weights.begin());
                        dist = std::discrete_distribution<>(weights.begin(), weights.end());
                        assert(size_t(dist.max()) == vec.size() - 1);
                }

                void push_back(const T &value, double weight=1) {
                    vec.push_back(value);
                    weights.push_back(weight);
                    setWeight(weights);
                    index = vec.size()-1;
                } //!< add data and it's weight (default = 1)

                const T& get() {
                    assert( not empty() && "no data!");
                    index = dist(random.engine);
                    return vec.at(index);
                } //!< retrieve data with given weight
        };

#ifdef DOCTEST_LIBRARY_INCLUDED
    TEST_CASE("[Faunus] WeightedDistribution")
    {
        WeightedDistribution<double> v;

        v.push_back(0.5);
        CHECK( v.index==0 );
        CHECK( v.size()==1 );

        v.push_back(0.1, 4);
        CHECK( v.index==1 );
        CHECK( v.size()==2 );
        CHECK( not v.empty() );

        int N=1e4;
        double sum=0;
        for (int i=0; i<N; i++)
            sum += v.get();
        CHECK( sum/N == doctest::Approx( (0.5*1+0.1*4) / (1+4) ).epsilon(0.05) );

        v.setWeight<int>( {2,1} );
        sum=0;
        for (int i=0; i<N; i++)
            sum += v.get();
        CHECK( sum/N == doctest::Approx( (0.5*2+0.1*1) / (2+1) ).epsilon(0.05) );

        CHECK_THROWS(v.setWeight<float>({2,1,1}));

        v.clear();
        CHECK( v.empty() );
    }
#endif
}
