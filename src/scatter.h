#pragma once

#include <fstream>
#include <algorithm>
#include <cmath>

namespace Faunus {

/** @brief Routines related to scattering.
 */
namespace Scatter {

enum Algorithm { SIMD, EIGEN, GENERIC }; //!< Selections for math algorithms

/** @brief Form factor, `F(q)`, for a hard sphere of radius `R`.
 */
template <class T = float> class FormFactorSphere {
  private:
    T j1(T x) const { // spherical Bessel function
        T xinv = 1 / x;
        return xinv * (sin(x) * xinv - cos(x));
    }

  public:
    /**
     * @param q q value in inverse angstroms
     * @param a particle to take radius, \c R from
     * @returns
     * @f$I(q)=\left [\frac{3}{(qR)^3}\left (\sin{qR}-qR\cos{qR}\right )\right ]^2@f$
     */
    template <class Tparticle> T operator()(T q, const Tparticle &a) const {
        assert(q > 0 && a.radius > 0 && "Particle radius and q must be positive");
        T qR = q * a.radius;
        qR = 3. / (qR * qR * qR) * (sin(qR) - qR * cos(qR));
        return qR * qR;
    }
};

/**
 * @brief Unity form factor (q independent).
 */
template <class T = float> struct FormFactorUnity {
    template <class Tparticle> T operator()(T, const Tparticle &) const { return 1; }
};

/**
 * @brief Calculate scattering intensity, I(q), on a mesh using the Debye formula.
 *
 * It is important to note that distances should be calculated without periodicity and if molecules cross
 * periodic boundaries, these must be made whole before performing the analysis.
 *
 * The JSON object is scanned for the following keywords:
 *
 * - `qmin` minimum q value (1/angstrom)
 * - `qmax` maximum q value (1/angstrom)
 * - `dq` q mesh spacing (1/angstrom)
 * - `cutoff` cutoff distance (angstrom); *Experimental!*
 *
 * @see http://dx.doi.org/10.1016/S0022-2860(96)09302-7
 */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdouble-promotion"

template <class Tformfactor, class T = float> class DebyeFormula {
    static constexpr T r_cutoff_infty = 1e9; //<! a cutoff distance in angstrom considered to be infinity
    T q_mesh_min, q_mesh_max, q_mesh_step; //<! q_mesh parameters in inverse angstrom; used for inline lambda-functions

    /**
     * @param m mesh point index
     * @return the scattering vector magnitude q at the mesh point m
     */
    #pragma omp declare simd uniform(this) linear(m:1)
    inline T q_mesh(int m) {
        return q_mesh_min + m * q_mesh_step;
    }

    /**
     * @brief Initialize mesh for intensity and sampling.
     * @param q_min Minimum q-value to sample (1/A)
     * @param q_max Maximum q-value to sample (1/A)
     * @param q_step Spacing between mesh points (1/A)
     */
    void init_mesh(T q_min, T q_max, T q_step) {
        if (q_step <= 0 || q_min <= 0 || q_max <= 0 || q_min > q_max ||
            q_step / q_max < 4 * std::numeric_limits<T>::epsilon()) {
            throw std::range_error("DebyeFormula: Invalid mesh parameters for q");
        }
        q_mesh_min = q_min < T(1e-6) ? q_step : q_min; // ensure that q > 0
        q_mesh_max = q_max;
        q_mesh_step = q_step;
        try {
            // resolution of the 1D mesh approximation of the scattering vector magnitude q
            const int q_resolution = numeric_cast<int>(1.0 + std::floor((q_max - q_min) / q_step));
            intensity.resize(q_resolution, 0.0);
            sampling.resize(q_resolution, 0.0);
        } catch (std::overflow_error &e) {
            throw std::range_error("DebyeFormula: Too many samples");
        }
    }

    Geometry::Sphere geo = Geometry::Sphere(r_cutoff_infty / 2); //!< geometry to use for distance calculations
    T r_cutoff;                   //!< cut-off distance for scattering contributions (angstrom)
    Tformfactor form_factor;      //!< scattering from a single particle
    std::vector<T> intensity;     //!< sampled average I(q)
    std::vector<T> sampling;      //!< weighted number of samplings

  public:
    DebyeFormula(T q_min, T q_max, T q_step, T r_cutoff) : r_cutoff(r_cutoff) { init_mesh(q_min, q_max, q_step); };

    DebyeFormula(T q_min, T q_max, T q_step) : DebyeFormula(r_cutoff_infty, q_min, q_max, q_step) {};

    explicit DebyeFormula(const json &j)
        : DebyeFormula(j.at("qmin").get<double>(), j.at("qmax").get<double>(), j.at("dq").get<double>(),
                       j.value("cutoff", r_cutoff_infty)){};

    /**
     * @brief Sample I(q) and add to average.
     * @param p particle vector
     * @param weight weight of sampled configuration in biased simulations
     * @param volume simulation volume (angstrom cubed) used only for cut-off correction
     *
     * An isotropic correction is added beyond a given cut-off distance. For physics details see for example
     * @see https://debyer.readthedocs.org/en/latest/.
     *
     * O(N^2) * O(M) complexity where N is the number of particles and M the number of mesh points. The quadratic
     * complexity in N comes from the fact that the radial distribution function has to be computed.
     * The current implementation supports OpenMP parallelization. Roughly half of the execution time is spend
     * on computing sin values, e.g., in sinf_avx2.
     */
    template <class Tpvec> void sample(const Tpvec &p, const T weight = 1, const T volume = -1) {
        const int N = (int) p.size(); // number of particles
        const int M = (int) intensity.size(); // number of mesh points
        std::vector<T> intensity_sum(M, 0.0);

        // Allow parallelization with a hand written reduction of intensity_sum at the end.
        // https://gcc.gnu.org/gcc-9/porting_to.html#ompdatasharing
        // #pragma omp parallel default(none) shared(N, M) shared(geo, r_cutoff, p) shared(intensity_sum)
        #pragma omp parallel default(shared) shared(intensity_sum)
        {
            std::vector<T> intensity_sum_private(M, 0.0); // a temporal private intensity_sum
            #pragma omp for schedule(dynamic)
            for (int i = 0; i < N - 1; ++i) {
                for (int j = i + 1; j < N; ++j) {
                    T r = T(geo.sqdist(p[i], p[j])); // the square root follows
                    if (r < r_cutoff * r_cutoff) {
                        r = std::sqrt(r);
                        // Black magic: The q_mesh function must be inlineable otherwise the loop cannot be unrolled
                        // using advanced SIMD instructions leading to a huge performance penalty (a factor of 4).
                        // The unrolled loop uses a different sin implementation, which may be spotted when profiling.
                        // TODO: Optimize also for other compilers than GCC by using a vector math library, e.g.,
                        // TODO: https://github.com/vectorclass/version2
                        // #pragma GCC unroll 16 // for diagnostics, GCC issues warning when cannot unroll
                        for (int m = 0; m < M; ++m) {
                            const T q = q_mesh(m);
                            intensity_sum_private[m] +=
                                form_factor(q, p[i]) * form_factor(q, p[j]) * std::sin(q * r) / (q * r);
                        }
                    }
                }
            }
            // reduce intensity_sum_private into intensity_sum
            #pragma omp critical
            std::transform(intensity_sum.begin(), intensity_sum.end(), intensity_sum_private.begin(),
                           intensity_sum.begin(), std::plus<T>());
        }

        // https://gcc.gnu.org/gcc-9/porting_to.html#ompdatasharing
        // #pragma omp parallel for default(none) shared(N, M, weight, volume) shared(p, r_cutoff, intensity_sum) shared(sampling, intensity)
        #pragma omp parallel for shared(sampling, intensity)
        for (int m = 0; m < M; ++m) {
            const T q = q_mesh(m);
            T intensity_self_sum = 0;
            for (int i = 0; i < N; ++i) {
                intensity_self_sum += std::pow(form_factor(q, p[i]), 2);
            }
            T intensity_corr = 0;
            if (r_cutoff < r_cutoff_infty && volume > 0) {
                intensity_corr = 4 * pc::pi * N / (volume * std::pow(q, 3)) *
                         (q * r_cutoff * std::cos(q * r_cutoff) - std::sin(q * r_cutoff));
            }
            sampling[m] += weight;
            intensity[m] += ((2 * intensity_sum[m] + intensity_self_sum) / N + intensity_corr) * weight;
        }
    }

    /**
     * @return a tuple of min, max, and step parameters of a q-mash
     */
    auto getQMeshParameters() {
        return std::make_tuple(q_mesh_min, q_mesh_max, q_mesh_step);
    }

    /**
     * @return a map containing q (key) and average intensity (value)
     */
    auto getIntensity() {
        std::map<T, T> averaged_intensity;
        for (size_t m = 0; m < intensity.size(); ++m) {
            const T average = intensity[m] / (sampling[m] != T(0.0) ? sampling[m] : T(1.0));
            averaged_intensity.emplace(q_mesh(m), average);
        }
        return averaged_intensity;
    }
};
#pragma GCC diagnostic pop

/**
 * A policy for collecting samples. To be used together with StructureFactor class templates.
 *
 * @tparam T float or double
 */
template <typename T> class SamplingPolicy {
  public:
    struct sampled_value {
        T value;
        T weight;
    };
    typedef std::map<T, sampled_value> TSampledValueMap;
  private:
    TSampledValueMap samples;
    const T precision = 10000.0; //!< precision of the key for better binning

  public:
    std::map<T, T> getSampling() const {
        std::map<T, T> average;
        for (auto [key, sample] : samples) {
            average.emplace(key, sample.value / sample.weight);
        }
        return average;
    }

    /**
     * @param key_approx is subject of rounding for better binning
     * @param value
     * @param weight
     */
    void addSampling(T key_approx, T value, T weight = 1.0) {
        const T key = std::round(key_approx * precision) / precision; // round |q| for better binning
        samples[key].value += value * weight;
        samples[key].weight += weight;
    }
};


/**
 * @brief Calculate structure factor using explicit q averaging.
 *
 * This averages over the thirteen permutations of the Miller index [100], [110], [101] using:
 *
 * @f[ S(\mathbf{q}) = \frac{1}{N} \left <
 *    \left ( \sum_i^N \sin(\mathbf{qr}_i) \right )^2 +
 *    \left ( \sum_j^N \cos(\mathbf{qr}_j) \right )^2
 *   \right >
 * @f]
 *
 * For more information, see @see http://doi.org/d8zgw5 and @see http://doi.org/10.1063/1.449987.
 */
template <typename T = double, Algorithm method = SIMD, typename TSamplingPolicy = SamplingPolicy<T>>
class StructureFactorPBC : private TSamplingPolicy {
    //! sample directions (h,k,l)
    const std::vector<Point> directions = {
        {1, 0, 0}, {0, 1, 0},  {0, 0, 1},                                      // 3 permutations
        {1, 1, 0}, {0, 1, 1},  {1, 0, 1},  {-1, 1, 0}, {-1, 0, 1}, {0, -1, 1}, // 6 permutations
        {1, 1, 1}, {-1, 1, 1}, {1, -1, 1}, {1, 1, -1}                          // 4 permutations
    };

    const int p_max;  //!< multiples of q to be sampled
    using TSamplingPolicy::addSampling;

  public:
    StructureFactorPBC(int q_multiplier) : p_max(q_multiplier){}

    template <class Tpositions> void sample(const Tpositions &positions, const Point &boxlength) {
        // https://gcc.gnu.org/gcc-9/porting_to.html#ompdatasharing
        // #pragma omp parallel for collapse(2) default(none) shared(directions, p_max, boxlength) shared(positions)
        #pragma omp parallel for collapse(2) default(shared)
        for (int i = 0; i < directions.size(); ++i) {
            for (int p = 1; p <= p_max; ++p) {                                         // loop over multiples of q
                const Point q = 2.0 * pc::pi * p * directions[i].cwiseQuotient(boxlength); // scattering vector
                T sum_sin = 0.0;
                T sum_cos = 0.0;
                if constexpr (method == SIMD) {
                    // When sine and cosine is computed in separate loops, advanced sine and cosine implementation
                    // utilizing SIMD instructions may be used to get at least 4 times performance boost.
                    // As of January 2020, only GCC exploits this using libmvec library if --ffast-math is enabled.
                    std::vector<T> qr_std(positions.size());
                    std::transform(positions.begin(), positions.end(), qr_std.begin(),
                                   [&q](const auto &r) { return q.dot(r); });
                    // as of January 2020 the std::transform_reduce is not implemented in libc++
                    for (const auto &qr : qr_std) {
                        sum_sin += std::sin(qr);
                    }
                    // as of January 2020 the std::transform_reduce is not implemented in libc++
                    for (const auto &qr : qr_std) {
                        sum_cos += std::cos(qr);
                    }
                } else if constexpr (method == EIGEN) {
                    // Map is a Nx3 matrix facade into original std::vector. Eigen does not accept `float`,
                    // hence `double` must be used instead of the template parameter T here.
                    auto qr = Eigen::Map<Eigen::MatrixXd, 0, Eigen::Stride<1, 3>>((double *)positions.data(),
                        positions.size(), 3) * q;
                    sum_sin = qr.array().cast<T>().sin().sum();
                    sum_cos = qr.array().cast<T>().cos().sum();
                } else if constexpr (method == GENERIC) {
                    // TODO: Optimize also for other compilers than GCC by using a vector math library, e.g.,
                    // TODO: https://github.com/vectorclass/version2
                    for (const auto &r : positions) { // loop over positions
                        T qr = q.dot(r);        // scalar product q*r
                        sum_sin += sin(qr);
                        sum_cos += cos(qr);
                    }
                };
                // collect average, `norm()` gives the scattering vector length
                const T sf = (sum_sin * sum_sin + sum_cos * sum_cos) / (T)(positions.size());
                #pragma omp critical
                // avoid race conditions when updating the map
                addSampling(q.norm(), sf, 1.0);
            }
        }
    }

    int getQMultiplier() {
        return p_max;
    }

    using TSamplingPolicy::getSampling;
};

/**
 * @brief Calculate structure factor using explicit q averaging in isotropic periodic boundary conditions (IPBC).
 *
 * The sample directions reduce to 3 compared to 13 in regular periodic boundary conditions. Overall simplification
 * shall yield roughly 10 times faster computation.
 */
template <typename T = float, typename TSamplingPolicy = SamplingPolicy<T>> class StructureFactorIPBC : private TSamplingPolicy {
    //! Sample directions (h,k,l).
    //! Due to the symmetry in IPBC we need not consider permutations of directions.
    std::vector<Point> directions = {{1, 0, 0}, {1, 1, 0}, {1, 1, 1}};

    int p_max;  //!< multiples of q to be sampled
    using TSamplingPolicy::addSampling;

  public:
    explicit StructureFactorIPBC(int q_multiplier) : p_max(q_multiplier) {}

    template <class Tpositions> void sample(const Tpositions &positions, const Point &boxlength) {
        // https://gcc.gnu.org/gcc-9/porting_to.html#ompdatasharing
        // #pragma omp parallel for collapse(2) default(none) shared(directions, p_max, positions, boxlength)
        #pragma omp parallel for collapse(2) default(shared)
        for (size_t i = 0; i < directions.size(); ++i) {
            for (int p = 1; p <= p_max; ++p) {                                             // loop over multiples of q
                const Point q = 2.0 * pc::pi * p * directions[i].cwiseQuotient(boxlength); // scattering vector
                T sum_cos = 0;
                for (const auto &r : positions) { // loop over positions
                    // if q[i] == 0 then its cosine == 1 hence we can avoid cosine computation for performance reasons
                    T product = std::cos(T(q[0] * r[0]));
                    if (q[1] != 0)
                        product *= std::cos(T(q[1] * r[1]));
                    if (q[2] != 0)
                        product *= std::cos(T(q[2] * r[2]));
                    sum_cos += product;
                }
                // collect average, `norm()` gives the scattering vector length
                const T ipbc_factor = std::pow(2, directions[i].count()); // 2 ^ number of non-zero elements
                const T sf = (sum_cos * sum_cos) / (float)(positions.size()) * ipbc_factor;
                #pragma omp critical
                // avoid race conditions when updating the map
                addSampling(q.norm(), sf, 1.0);
            }
        }
    }

    int getQMultiplier() {
        return p_max;
    }

    using TSamplingPolicy::getSampling;
};

} // namespace Scatter
} // namespace Faunus
