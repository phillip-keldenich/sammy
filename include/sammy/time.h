#ifndef SAMMY_TIME_H_INCLUDED_
#define SAMMY_TIME_H_INCLUDED_

#include <chrono>

namespace sammy {

using Clock = std::chrono::steady_clock;

template <typename TP1, typename TP2>
inline double seconds_between(const TP1& before, const TP2& after) noexcept {
    return std::chrono::duration_cast<std::chrono::duration<double>>(after -
                                                                     before)
        .count();
}

} // namespace sammy

#endif
