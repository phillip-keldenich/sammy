#ifndef SAMMY_STAMP_SET_H_INCLUDED_
#define SAMMY_STAMP_SET_H_INCLUDED_

#include <cstdint>
#include <limits>
#include <type_traits>
#include <vector>

namespace sammy {

/**
 * @brief A data structure for efficient management of a set of elements.
 * 
 * The StampSet allows for efficient insertion, deletion, reset, and membership
 * checks of elements. It is particularly well-suited for scenarios where the
 * universe of elements is small and consists of integers indexed from 0 to n-1.
 * The implementation uses an array (`m_stamps`) to track membership by comparing
 * `m_stamps[i]` with the current stamp value. While a simple `std::vector<bool>`
 * could be used, this variant allows to perform `clear` in O(1) time instead of
 * O(n) time, as it just increases the current stamp value, invalidating all
 * previously inserted elements. This limits the number of feasible `clear`
 * calls to `std::numeric_limits<StampType>::max()`, but this is usually not a
 * problem in practice.
 * 
 * @tparam ValueType An integral type representing the elements in the set.
 * @tparam StampType An unsigned integral type used for stamps. It must be large
 * enough to safely handle the number of expected `clear` calls.
 */
template <typename ValueType, typename StampType = std::uint32_t>
class StampSet {
  private:
    static_assert(std::is_integral_v<ValueType>, "Value type must be integral");
    static_assert(std::is_integral_v<StampType> &&
                      std::is_unsigned_v<StampType>,
                  "Stamp type must be unsigned");

    std::vector<StampType> m_stamps;
    StampType m_current_stamp;

  public:
    explicit StampSet(ValueType universe_size)
        : m_stamps(universe_size, StampType(0)), m_current_stamp(1) {}

    // Basic constructor, copy constructor, and assignment operator.
    StampSet(const StampSet&) = default;
    StampSet& operator=(const StampSet&) = default;
    StampSet(StampSet&&) noexcept = default;
    StampSet& operator=(StampSet&&) noexcept = default;

    std::size_t universe_size() const noexcept { return m_stamps.size(); }

    // Resets the set to an empty state.
    void clear() noexcept {
        if (++m_current_stamp == 0) {
            std::fill(m_stamps.begin(), m_stamps.end(), StampType(0));
            m_current_stamp = 1;
        }
    }

    // Clears the set and assigns the elements from the range [begin, end).
    template <typename ForwardIterator>
    void assign(ForwardIterator begin, ForwardIterator end) noexcept {
        clear();
        insert(begin, end);
    }

    // Inserts the elements from the range [begin, end) into the set.
    template <typename ForwardIterator>
    void insert(ForwardIterator begin, ForwardIterator end) noexcept {
        std::for_each(begin, end, [&](ValueType l) { insert(l); });
    }

    // Inserts a single element into the set.
    void insert(ValueType v) noexcept { m_stamps[v] = m_current_stamp; }

    // Erases a single element from the set.
    void erase(ValueType v) noexcept { m_stamps[v] = 0; }

    // inserts a single element into the set and returns true if it was not
    // already present.
    bool check_insert(ValueType v) noexcept {
        StampType& s = m_stamps[v];
        bool result = (s != m_current_stamp);
        s = m_current_stamp;
        return result;
    }

    // Erases a single element from the set and returns true if it was present.
    bool check_erase(ValueType v) noexcept {
        StampType& s = m_stamps[v];
        bool result = (s == m_current_stamp);
        s = 0;
        return result;
    }

    // Checks if the set contains a single element.
    bool count(ValueType v) const noexcept {
        return m_stamps[v] == m_current_stamp;
    }

    // Checks if the set contains a single element.
    bool contains(ValueType v) const noexcept { return count(v); }
};

} // namespace sammy

#endif
