#ifndef HS_RNG_H_INCLUDED_
#define HS_RNG_H_INCLUDED_

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <random>

namespace sammy {

inline std::mt19937_64& rng() noexcept {
    static std::atomic<std::int32_t> counter(0);
    thread_local std::mt19937_64 res(1337 + counter++);
    return res;
}

template<typename ContainerIn, typename ContainerOut, typename RNG>
inline void sample_from_range(const ContainerIn& container_in,
                              ContainerOut& container_out,
                              std::size_t goal_size, RNG& rng) 
{
    using IndexDist = std::uniform_int_distribution<std::size_t>;
    using Value = typename ContainerIn::value_type;
    if (container_in.empty()) {
        return;
    }
    if (container_in.size() <= goal_size) {
        container_out.insert(container_out.end(), container_in.begin(),
                             container_in.end());
        return;
    }
    std::size_t expected_max_out(1.1 * goal_size + container_out.size());
    container_out.reserve(expected_max_out);
    double probability = double(goal_size) / container_in.size();
    if (probability > 0.1) {
        IndexDist cdist(0, container_in.size() - 1);
        for (const Value& v : container_in) {
            std::size_t x = cdist(rng);
            if (x < goal_size) {
                container_out.push_back(v);
            }
        }
    } else {
        std::geometric_distribution<std::size_t> tdist(probability);
        auto i = container_in.begin(), e = container_in.end();
        i += tdist(rng);
        while (i < e) {
            container_out.push_back(*i);
            i += tdist(rng) + 1;
        }
    }
}

template <typename Container, typename RNG>
inline std::vector<typename Container::value_type>
sample_from_range(const Container& container, std::size_t goal_size, RNG& rng) {
    using Value = typename Container::value_type;
    std::vector<Value> result;
    sample_from_range(container, result, goal_size, rng);
    return result;
}

} // namespace sammy

#endif
