#pragma once

#include <fstream>

namespace Faunus {

/** @brief Routines related to scattering */
namespace Scatter {

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
     * @param a particle to take radius, \c R from.
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
 * @brief Unity form factor (q independent)
 */
template <class T = float> struct FormFactorUnity {
    template <class Tparticle> T operator()(T, const Tparticle &) const { return 1; }
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
template <class Tformfactor, class Tgeometry = Geometry::Chameleon, class T = float> class DebyeFormula {
  protected:
    Tformfactor F;    // scattering from a single particle
    Tgeometry geo;    // geometry to use for distance calculations
    std::map<T, T> I; //!< Sampled, average I(q)
    std::map<T, T> S; //!< Weighted number of samplings
    T rc;

  public:
    T qmin, qmax, dq;

    DebyeFormula(const json &j) {
        geo = R"( { "type": "sphere", "radius": 1e9 } )"_json;
        dq = j.at("dq").get<double>();
        qmin = j.at("qmin").get<double>();
        qmax = j.at("qmax").get<double>();
        rc = j.value("cutoff", 1.0e9);

        if (dq <= 0 || qmin <= 0 || qmax <= 0 || qmin > qmax)
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
    template <class Tpvec> void sample(const Tpvec &p, T f = 1, T V = -1) {
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
    template <class Tpvec> void sample(const Tpvec &p, T qmin, T qmax, T dq, T f = 1, T V = -1) {
        if (qmin < 1e-6)
            qmin = dq; // ensure that q>0

        // Temporary f(q) functions - initialized to
        // enable O(N) complexity iteration in inner loop.
        std::map<T, T> _I, _ff;
        for (T q = qmin; q <= qmax; q += dq)
            _I[q] = _ff[q] = 0;

        int N = (int)p.size();
        for (int i = 0; i < N - 1; ++i) {
            for (int j = i + 1; j < N; ++j) {
                T r = geo.sqdist(p[i], p[j]);
                if (r < rc * rc) {
                    r = sqrt(r);
                    for (auto &m : _I) { // O(N) complexity
                        T q = m.first;
                        m.second += F(q, p[i]) * F(q, p[j]) * sin(q * r) / (q * r);
                    }
                }
            }
        }
        for (int i = 0; i < N; i++)
            for (T q = qmin; q <= qmax; q += dq)
                _ff[q] += pow(F(q, p[i]), 2);

        for (auto &i : _I) {
            T q = i.first, Icorr = 0;
            if (rc < 1e9 && V > 0)
                Icorr = 4 * pc::pi * N / (V * pow(q, 3)) * (q * rc * cos(q * rc) - sin(q * rc));
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
    template <class Tpvec, class Tg> void sampleg2g(const Tpvec &p, Tg groupList) {
        assert(qmin > 0 && qmax > 0 && dq > 0 && "q range invalid.");
        sampleg2g(p, qmin, qmax, dq, groupList);
    }

    template <class Tpvec, class Tg> void sampleg2g(const Tpvec &p, T qmin, T qmax, T dq, Tg groupList) {
        std::vector<T> _I(int((qmax - qmin) / dq + 0.5));
        // loop over all pairs of groups, then over particles
        for (int k = 0; k < (int)groupList.size() - 1; k++)
            for (int l = k + 1; l < (int)groupList.size(); l++)
                for (auto i : *groupList.at(k))
                    for (auto j : *groupList.at(l)) {
                        T r = geo.vdist(p[i].pos, p[j].pos).norm();
                        int cnt = 0;
                        for (T q = qmin; q <= qmax; q += dq)
                            _I.at(cnt++) += sin(q * r) / (q * r);
                    }
        int cnt = 0, N = p.size();
        for (T q = qmin; q <= qmax; q += dq)
            I[q] += 2. * _I.at(cnt++) / N + 1; // add to average I(q)
    }

    /**
     * @brief Save I(q) to disk
     */
    void save(const std::string &filename) {
        if (!I.empty()) {
            std::ofstream f(filename.c_str());
            if (f) {
                for (auto &i : I) {
                    if (not S.empty())
                        f << i.first << " " << i.second / S[i.first] << "\n";
                    else
                        f << i.first << " " << i.second << "\n";
                }
            }
        }
    }
};

/**
 * @brief Structor factor calculation using explicit q averaging
 *
 * This averages over the thirteen permutations of
 * the Miller index [100], [110], [101] using:
 *
 * @f[ S(\mathbf{q}) = \frac{1}{N} \left <
 *    \left ( \sum_i^N \sin(\mathbf{qr}_i) \right )^2 +
 *    \left ( \sum_j^N \cos(\mathbf{qr}_j) \right )^2
 *   \right >
 * @f]
 *
 * For more information, see http://doi.org/d8zgw5 and http://doi.org/10.1063/1.449987
 *
 * @todo Add OpenMP pragmas and particle formfactors
 */
template <typename T = double> class StructureFactor {
  private:
    std::map<T, T> S; //!< Average S(q)
    std::map<T, T> W; //!< Weighted number of samplings

    // Sample directions (h,k,l)
    std::vector<Point> directions = {
        {1, 0, 0}, {0, 1, 0},  {0, 0, 1},                                      // 3 permutations
        {1, 1, 0}, {0, 1, 1},  {1, 0, 1},  {-1, 1, 0}, {-1, 0, 1}, {0, -1, 1}, // 6 permutations
        {1, 1, 1}, {-1, 1, 1}, {1, -1, 1}, {1, 1, -1}                          // 4 permutations
    };

  public:
    StructureFactor(int number_of_scalings) : pmax(number_of_scalings) {}

    int pmax; //!< Multiples of q to be sampled

    template <class Tpvec> void sample(const Tpvec &positions, double boxlength) {
        for (auto &dir : directions) {                         // loop over 3+6+4=13 directions
            for (int p = 1; p <= pmax; p++) {                  // loop over multiples of q
                Point _q = (2 * pc::pi * p / boxlength) * dir; // scattering vector
                double sum_sin = 0, sum_cos = 0;               // temporary sums
                for (auto &r : positions) {                    // loop over positions
                    double qr = _q.dot(r);                     // scalar product q*r
                    sum_sin += sin(qr);
                    sum_cos += cos(qr);
                }
                // collect average, `norm()` gives the scattering vector length
                S[_q.norm()] += (sum_sin * sum_sin + sum_cos * sum_cos) / double(positions.size());
                W[_q.norm()] += 1.0;
            }
        }
    }

    void save(const std::string &filename) {
        if (not S.empty()) {
            std::ofstream f(filename.c_str());
            if (f)
                for (auto &i : S)
                    f << i.first << " " << i.second / W[i.first] << "\n";
        }
    }
};

} // namespace Scatter
} // namespace Faunus
