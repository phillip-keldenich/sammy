#ifndef HS_RNG_H_INCLUDED_
#define HS_RNG_H_INCLUDED_

#include <atomic>
#include <cstdint>
#include <random>
#include <cstddef>

namespace sammy {

inline std::mt19937_64& rng() noexcept {
    static std::atomic<std::int32_t> counter(0);
    thread_local std::mt19937_64 res(1337 + counter++);
    return res;
}

template<typename Container, typename RNG>
inline std::vector<typename Container::value_type>
sample_from_range(const Container& container, std::size_t goal_size, RNG& rng)
{
    using Value = typename Container::value_type;
    std::vector<Value> result;
    if(container.begin() == container.end()) return result;
    if(container.size() <= goal_size) {
        result.insert(result.end(), container.begin(), container.end());
        return result;
    }
    double probability = double(goal_size) / container.size();
    result.reserve(std::size_t(std::round(1.1 * goal_size)));
    if(probability > 0.1) {
        std::uniform_int_distribution<std::size_t> cdist(0, container.size() - 1);
        for(const Value& v : container) {
            std::size_t x = cdist(rng);
            if(x < goal_size) result.push_back(v);
        }
    } else {
        std::geometric_distribution<std::size_t> tdist(probability);
        auto i = container.begin(), e = container.end();
        i += tdist(rng);
        while(i < e) {
            result.push_back(*i);
            i += tdist(rng) + 1;
        }
    }
    return result;
}

} // namespace hs

#endif
