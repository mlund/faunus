#pragma once

#include <map>
#include <fstream>

namespace Faunus {

/**
 * @brief Histogram for an arbitrary set of values using a sparse memory layout (map)
 *
 * Builds a histogram by binning given values to a specified resolution. Values are stored
 * in a memory efficient map-structure with log(N) lookup complexity.
 */
template <typename T = double> class SparseHistogram {
    T resolution;
    std::map<int, unsigned int> data;

  public:
    explicit SparseHistogram(T resolution) : resolution(resolution) {}
    void add(const T value) {
        if (std::isfinite(value)) {
            data[static_cast<int>(std::round(value / resolution))]++;
        } else {
            faunus_logger->warn("histogram: skipping inf/nan number");
        }
    }
    friend auto& operator<<(std::ostream& stream, const SparseHistogram& histogram) {
        std::for_each(histogram.data.begin(), histogram.data.end(), [&](const auto& sample) {
            stream << fmt::format("{:.6E} {}\n", T(sample.first) * histogram.resolution, sample.second);
        });
        return stream;
    }
};
} // namespace Faunus
