#ifndef SAMMY_SHALLOW_MEM_ESTIMATE_H_INCLUDED_
#define SAMMY_SHALLOW_MEM_ESTIMATE_H_INCLUDED_

#include "dynamic_bitset.h"
#include <climits>
#include <cstddef>
#include <type_traits>
#include <utility>
#include <vector>

namespace sammy {

template <typename T, typename A,
          std::enable_if_t<!std::is_same_v<T, bool>, int> = 0>
inline std::size_t shallow_memory_estimate(const std::vector<T, A>& vec) {
    return vec.capacity() * sizeof(T);
}

template <typename A>
inline std::size_t shallow_memory_estimate(const std::vector<bool, A>& vec) {
    return vec.capacity() / CHAR_BIT;
}

inline std::size_t shallow_memory_estimate(const DynamicBitset& bitset) {
    return shallow_memory_estimate(bitset.blocks());
}

template <typename T, typename A>
inline std::size_t depth1_memory_estimate(const std::vector<T, A>& vec) {
    std::size_t result = shallow_memory_estimate(vec);
    if (!vec.empty()) {
        result += vec.size() * shallow_memory_estimate(vec[0]);
    }
    return result;
}

} // namespace sammy

#endif
