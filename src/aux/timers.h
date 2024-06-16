#pragma once
#include <chrono>
#include <iostream>

namespace Faunus {
/**
 * @brief Timer for measuring relative time consumption
 *
 * Time t=0 is set upon construction whereafter combined `start()`/
 * `stop()` calls can be made multiple times. The result is
 * the fraction of total time, consumed in between start/stop calls.
 */
template <typename Tunit = std::chrono::microseconds> class TimeRelativeOfTotal {
  private:
    Tunit delta;
    std::chrono::steady_clock::time_point t0, tx;

  public:
    TimeRelativeOfTotal()
        : delta(0)
    {
        t0 = std::chrono::steady_clock::now();
    }

    operator bool() const { return delta.count() != 0 ? true : false; }

    void start() { tx = std::chrono::steady_clock::now(); }

    void stop() { delta += std::chrono::duration_cast<Tunit>(std::chrono::steady_clock::now() - tx); }

    double result() const
    {
        auto now = std::chrono::steady_clock::now();
        auto total = std::chrono::duration_cast<Tunit>(now - t0);
        return delta.count() / double(total.count());
    }
};

/**
 * @brief Simple struct to measure duration of code
 *
 * Example:
 *
 * ~~~ cpp
 * Stopwatch w;
 * w.start();
 * // do something expensive...
 * w.stop(); // --> stdio
 * ~~~
 *
 * `stop()` takes an optional template parameter to specify
 * the time resolution. Default is milliseconds.
 */
class Stopwatch {
    using clock = std::chrono::high_resolution_clock;
    std::chrono::time_point<clock> starting_time, ending_time;

  public:
    inline Stopwatch() { start(); }
    inline void start() { starting_time = clock::now(); }
    inline void stop() { ending_time = clock::now(); }
};

} // namespace Faunus