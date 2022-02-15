#pragma once

#include <map>
#include <fstream>
#include <range/v3/view/transform.hpp>
#if __cplusplus > 201703L
#include <concepts>
#endif

namespace Faunus {

/**
 * @brief Histogram for an arbitrary set of values using a sparse memory layout (map)
 *
 * Builds a histogram by binning given values to a specified resolution. Values are stored
 * in a memory efficient map-structure with log(N) lookup complexity.
 */
#if __cplusplus > 201703L
template <std::floating_point T = double> class SparseHistogram {
#else
template <typename T = double> class SparseHistogram {
#endif
    const T resolution;
    using map_type = std::map<int, unsigned int>;
    map_type data;

  public:
    explicit SparseHistogram(T resolution) : resolution(resolution) {}
    const T getResolution() const { return resolution; }
    void add(const T value) {
        if (std::isfinite(value)) {
            const auto index = static_cast<map_type::key_type>(std::round(value / resolution));
            data[index]++;
        } else {
            faunus_logger->warn("histogram: skipping inf/nan number");
        }
    }
    friend auto& operator<<(std::ostream& stream, const SparseHistogram& histogram) {
        std::for_each(histogram.data.begin(), histogram.data.end(), [&](const auto& sample) {
            stream << fmt::format("{:.6E} {}\n", static_cast<T>(sample.first) * histogram.resolution, sample.second);
        });
        return stream;
    }
};
} // namespace Faunus
