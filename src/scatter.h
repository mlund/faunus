#pragma once

namespace Faunus {

    /** @brief Routines related to scattering */
    namespace Scatter {

        /** @brief Form factor, `F(q)`, for a hard sphere of radius `R`.
        */
        template<class T=float>
            class FormFactorSphere {
                private:
                    T j1( T x ) const { // spherical Bessel function
                        T xinv = 1 / x;
                        return xinv * (sin(x) * xinv - cos(x));
                    }

                public:
                    /**
                     * @param q q value in inverse angstroms
                     * @param a particle to take radius, \c R from.
                     * @returns
                     * @f$I(q)=\left [\frac{3}{(qR)^3}\left (\sin{qR}-qR\cos{qR}\right )\right ]^2@f$
                     */
                    template<class Tparticle>
                        T operator()( T q, const Tparticle &a ) const {
                            assert(q > 0 && a.radius > 0 && "Particle radius and q must be positive");
                            T qR = q * a.radius;
                            qR = 3. / (qR * qR * qR) * (sin(qR) - qR * cos(qR));
                            return qR * qR;
                        }
            };

        /**
         * @brief Unity form factor (q independent)
         */
        template<class T=float>
            struct FormFactorUnity {
                template<class Tparticle>
                    T operator()( T q, const Tparticle &a ) const {
                        return 1;
                    }
            };

        /**
         * @brief Calculates scattering intensity, I(q) using the Debye formula
         *
         * It is important to note that distances should be calculated without
         * periodicity and if molecules cross periodic boundaries, these
         * must be made whole before performing the analysis.
         * The JSON object is scanned for the following keywords:
         *
         * - `qmin` Minimum q value (1/angstrom)
         * - `qmax` Maximum q value (1/angstrom)
         * - `dr` q spacing (1/angstrom)
         * - `cutoff` Cutoff distance (angstrom). *Experimental!*
         *
         * See also <http://dx.doi.org/10.1016/S0022-2860(96)09302-7>
         */
        template<class Tformfactor, class Tgeometry=Geometry::Sphere, class T=float>
            class DebyeFormula {
                private:
                    T qmin, qmax, dq, rc;
                protected:
                    Tformfactor F; // scattering from a single particle
                    Tgeometry geo; // geometry to use for distance calculations
                public:
                    std::map<T, T> I; //!< Sampled, average I(q)
                    std::map<T, T> S; //!< Weighted number of samplings

                    DebyeFormula(const json &j) {
                        geo = R"( { "radius" : 1e9 } )"_json;
                        dq = j.at("dq").get<double>();
                        qmin = j.at("qmin").get<double>();
                        qmax = j.at("qmax").get<double>();
                        rc = j.value("cutoff", 1.0e9);

                        if (dq<=0 || qmin<=0 || qmax<=0 || qmin>qmax)
                            throw std::runtime_error("DebyeFormula: invalid q parameters");
                    }

                    /**
                     * @brief Sample I(q) and add to average
                     *
                     * The q range is read from input as `qmin`, `qmax`, `dq` in units of
                     * inverse angstrom.
                     * It is also possible to specify an isotropic correction beyond
                     * a given cut-off -- see for example
                     * <https://debyer.readthedocs.org/en/latest/>.
                     * The cut-off distance - which should be smaller than half the
                     * box length for a cubic system is read specific by the
                     * keyword `sofq_cutoff`.
                     */
                    template<class Tpvec>
                        void sample( const Tpvec &p, T f = 1, T V = -1 ) {
                            assert(qmin > 0 && qmax > 0 && dq > 0 && "q range invalid.");
                            sample(p, qmin, qmax, dq, f, V);
                        }

                    /**
                     * @brief Sample I(q) and add to average
                     * @param p Particle vector
                     * @param qmin Minimum q-value to sample (1/A)
                     * @param qmax Maximum q-value to sample (1/A)
                     * @param dq q spacing (1/A)
                     * @param f weight of sampled configuration in biased simulations
                     * @param V Simulation volume (angstrom^3) used only for cut-off correction
                     */
                    template<class Tpvec>
                        void sample( const Tpvec &p, T qmin, T qmax, T dq, T f = 1, T V = -1 ) {
                            if ( qmin < 1e-6 )
                                qmin = dq;              // ensure that q>0

                            // Temporary f(q) functions - initialized to
                            // enable O(N) complexity iteration in inner loop.
                            std::map<T, T> _I, _ff;
                            for ( T q = qmin; q <= qmax; q += dq )
                                _I[q] = _ff[q] = 0;

                            int N = (int) p.size();
                            for ( int i = 0; i < N - 1; ++i )
                            {
                                for ( int j = i + 1; j < N; ++j )
                                {
                                    T r = geo.sqdist(p[i], p[j]);
                                    if ( r < rc * rc )
                                    {
                                        r = sqrt(r);
                                        for ( auto &m : _I )
                                        { // O(N) complexity
                                            T q = m.first;
                                            m.second += F(q, p[i]) * F(q, p[j]) * sin(q * r) / (q * r);
                                        }
                                    }
                                }
                            }
                            for ( int i = 0; i < N; i++ )
                                for ( T q = qmin; q <= qmax; q += dq )
                                    _ff[q] += pow(F(q, p[i]), 2);

                            for ( auto &i : _I )
                            {
                                T q = i.first, Icorr = 0;
                                if ( rc < 1e9 && V > 0 )
                                    Icorr = 4 * pc::pi * N / (V * pow(q, 3)) *
                                        (q * rc * cos(q * rc) - sin(q * rc));
                                S[q] += f;
                                I[q] += ((2 * i.second + _ff[q]) / N + Icorr) * f; // add to average I(q)
                            }
                        }

                    /**
                     * @brief Sample between all groups
                     *
                     * Instead of looping over all particles, this will ignore all internal
                     * group distances.
                     *
                     * @param p Particle vector
                     * @param groupList Vector of group pointers - i.e. as returned from `Space::groupList()`.
                     * @note Particle form factors are always set to unity, i.e. F(q) is ignored.
                     *
                     * @warning Untested
                     */
                    template<class Tpvec, class Tg>
                        void sampleg2g( const Tpvec &p, Tg groupList )
                        {
                            assert(qmin > 0 && qmax > 0 && dq > 0 && "q range invalid.");
                            sampleg2g(p, qmin, qmax, dq, groupList);
                        }

                    template<class Tpvec, class Tg>
                        void
                        sampleg2g( const Tpvec &p, T qmin, T qmax, T dq, Tg groupList )
                        {
                            std::vector<T> _I(int((qmax - qmin) / dq + 0.5));
                            // loop over all pairs of groups, then over particles
                            for ( int k = 0; k < (int) groupList.size() - 1; k++ )
                                for ( int l = k + 1; l < (int) groupList.size(); l++ )
                                    for ( auto i : *groupList.at(k))
                                        for ( auto j : *groupList.at(l))
                                        {
                                            T r = geo.vdist(p[i].pos, p[j].pos).norm();
                                            int cnt = 0;
                                            for ( T q = qmin; q <= qmax; q += dq )
                                                _I.at(cnt++) += sin(q * r) / (q * r);
                                        }
                            int cnt = 0, N = p.size();
                            for ( T q = qmin; q <= qmax; q += dq )
                                I[q] += 2. * _I.at(cnt++) / N + 1; // add to average I(q)
                        }

                    /**
                     * @brief Save I(q) to disk
                     */
                    void
                        save( const std::string &filename )
                        {
                            if ( !I.empty())
                            {
                                std::ofstream f(filename.c_str());
                                if ( f )
                                {
                                    for ( auto &i : I )
                                        f << i.first << " " << i.second / S[i.first] << "\n";
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
                private:
                    T qmin, qmax, dq;
                    Point qdir={0,0,0};
                public:
                    typedef DebyeFormula<FormFactorSphere<T>> base;

                    StructureFactor(const json &j) : base(j) {
                        qmin = j.at("qmin").get<double>();
                        qmax = j.at("qmax").get<double>();
                        dq = j.at("dq").get<double>();
                    }

                    template<class Tpvec>
                        void sample( const Tpvec &p ) {
                            sample(p, qmin, qmax, dq);
                        }

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
                     *
                     * The direction of the q vector is randomly generated on a unit
                     * sphere, i.e. during a simulation there will be an angular
                     * averaging, producing the same result as the Debye formula. The
                     * number of directions for each sample event can be set using
                     * the `Nq` parameter which defaults to 1. Averaging can be disabled
                     * by specifically setting a q vector via the last argument.
                     *
                     * @param p Particle/point vector
                     * @param qmin Minimum q-value to sample (1/A)
                     * @param qmax Maximum q-value to sample (1/A)
                     * @param dq q spacing (1/A)
                     * @param Nq Number of random q direction (default: 1)
                     * @param qdir Specific q vector - overrides averaging (default: 0,0,0)
                     */
                    template<class Tpvec>
                        void
                        sample( const Tpvec &p, T qmin, T qmax, T dq, int Nq = 10, Point qdir = Point(0, 0, 0))
                        {
                            if ( qmin < 1e-6 )
                                qmin = dq;  // assure q>0

                            int n = (int) p.size();
                            for ( int k = 0; k < Nq; k++ )
                            { // random q directions
                                std::map<T, T> _I; // temp map for I(q) value

                                // Random q vector if none given in input
                                if ( Nq == 0 && qdir.squaredNorm() > 1e-6 )
                                    Nq = 1;
                                else
                                    qdir = ranunit(random);

                                // N^2 loop over all particles
                                for ( int i = 0; i < n - 1; i++ )
                                {
                                    for ( int j = i + 1; j < n; j++ )
                                    {
                                        auto r = base::geo.vdist(p[i].pos, p[j].pos);
                                        for ( T q = qmin; q <= qmax; q += dq )
                                            _I[q] += cos((q * qdir).dot(r));
                                    }
                                }

                                for ( auto &i : _I )
                                    base::I[i.first] += 2 * i.second / n + 1; // add to average I(q)

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
                     * Angulalar averaging is done as in `sample()` and is
                     * in fact completely equivalent.
                     * For more information, see <http://doi.org/d8zgw5>
                     *
                     * @todo Swap loop order
                     */
                    template<class Tpvec>
                        void
                        sample2( const Tpvec &p, T qmin, T qmax, T dq, int Nq = 1 )
                        {
                            if ( qmin < 1e-6 )
                                qmin = dq;  // assure q>0
                            int n = (int) p.size();
                            for ( int k = 0; k < Nq; k++ )
                            { // random q directions
                                Point qdir(1, 0, 0);
                                qdir = ranunit(random);
                                std::map<T, T> _cos, _sin;
                                for ( T q = qmin; q <= qmax; q += dq )
                                {
                                    double ssum = 0, csum = 0; // tmp to avoid map lookup
                                    for ( auto &i : p )
                                    {
                                        T qr = (q * qdir).dot(i);
                                        ssum += sin(qr);
                                        csum += cos(qr);
                                    }
                                    _cos[q] = csum;
                                    _sin[q] = ssum;
                                }
                                for ( T q = qmin; q <= qmax; q += dq )
                                    base::I[q] += (pow(_sin[q], 2) + pow(_cos[q], 2)) / n;
                            } // end of q averaging
                        }

                    template<class Tpvec>
                        void
                        sample3( const Tpvec &p, T qmin, T qmax, T dq, int Nq = 1 )
                        {
                            if ( qmin < 1e-6 )
                                qmin = dq;  // assure q>0
                            int n = (int) p.size();
                            for ( int k = 0; k < Nq; k++ )
                            { // random q directions
                                Point qdir(1, 0, 0);
                                //qdir.ranunit(slump);
                                std::map<T, T> _cos, _sin;
                                for ( T q = qmin; q <= qmax; q += dq )
                                {
                                    for ( int i = 0; i < n - 1; i++ )
                                    {
                                        for ( int j = i + 1; j < n; j++ )
                                        {
                                            T qr = (q * qdir).dot(p[i]);
                                            _sin[q] += sin(qr);
                                            _cos[q] += cos(qr);
                                        }
                                    }
                                }
                                for ( auto &i : _sin )
                                {
                                    T q = i.first;
                                    base::I[q] += pow(_sin[q], 2) + pow(_cos[q], 2);
                                }
                            } // end of q averaging
                        }
            };

    } // end of namespace
} //end of namespace
