#ifndef SAMMY_RANGE_H_INCLUDED_
#define SAMMY_RANGE_H_INCLUDED_

#include <boost/iterator/iterator_facade.hpp>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <random>
#include <type_traits>
#include <utility>

namespace sammy {

/**
 * @brief A range spanned by two iterators.
 *
 * @tparam Iterator
 */
template <typename Iterator> class IteratorRange {
  public:
    IteratorRange(Iterator b, Iterator e) : m_beg(b), m_end(e) {}

    explicit IteratorRange(const std::pair<Iterator, Iterator>& p)
        : m_beg(p.first), m_end(p.second) {}

    Iterator begin() const noexcept { return m_beg; }
    Iterator end() const noexcept { return m_end; }
    std::size_t size() const noexcept { return std::distance(m_beg, m_end); }

  private:
    Iterator m_beg, m_end;
};

/**
 * @brief A range of random integers of a certain length.
 *
 * @tparam Value
 * @tparam RngType
 */
template <typename Value, typename RngType> class RandomIntRange {
    RngType* m_rng;
    std::size_t m_length;
    std::uniform_int_distribution<Value> m_dist;

  public:
    class Iterator {
        friend class RandomIntRange;
        RngType* m_rng;
        std::uniform_int_distribution<Value> m_dist;
        std::size_t m_offs;
        Value m_curr;

        Iterator(RngType* rng, std::uniform_int_distribution<Value> dist,
                 std::size_t offs, Value curr) noexcept
            : m_rng(rng), m_dist(dist), m_offs(offs), m_curr(curr) {}

      public:
        using iterator_category = std::input_iterator_tag;
        using value_type = Value;
        using difference_type = std::ptrdiff_t;
        using pointer = const Value*;
        using reference = Value;

        Iterator() = default;

        bool operator==(const Iterator& o) const noexcept {
            return m_offs == o.m_offs;
        }

        bool operator!=(const Iterator& o) const noexcept {
            return m_offs != o.m_offs;
        }

        Iterator& operator++() noexcept {
            ++m_offs;
            m_curr = m_dist(*m_rng);
            return *this;
        }

        Iterator operator++(int) noexcept {
            Iterator result(*this);
            ++*this;
            return result;
        }

        Value operator*() const noexcept { return m_curr; }

        const Value* operator->() const noexcept { return &m_curr; }
    };

    Iterator begin() { return Iterator(m_rng, m_dist, 0, m_dist(*m_rng)); }

    Iterator end() { return Iterator(m_rng, m_dist, m_length, 0); }

    RandomIntRange(RngType& rng, std::size_t length,
                   std::uniform_int_distribution<Value> dist = {})
        : m_rng(&rng), m_length(length), m_dist(dist) {}
};

template <typename ValueType, typename RngType>
inline auto random_int_range(std::size_t length, RngType& rng) noexcept {
    return RandomIntRange<ValueType, RngType>(rng, length);
}

template <typename ValueType, typename RngType>
inline auto random_int_range(
    std::size_t length, RngType& rng,
    const std::uniform_int_distribution<ValueType>& dist) noexcept {
    return RandomIntRange<ValueType, RngType>(rng, length, dist);
}

// deduction guide
template <typename Iterator1, typename Iterator2>
explicit IteratorRange(Iterator1, Iterator2)
    -> IteratorRange<std::remove_cv_t<std::remove_reference_t<Iterator1>>>;

template <typename RandomAccessIterator, typename Callable>
static inline void
split_range(RandomAccessIterator begin, RandomAccessIterator end,
            std::size_t max_partitions, Callable&& callback) {
    auto count = std::size_t(end - begin);
    if (max_partitions > count)
        max_partitions = count;
    std::size_t elements_per_partition = count / max_partitions;
    std::size_t rem_mod = count % max_partitions;
    RandomAccessIterator current = begin;
    for (std::size_t i = 0; i < max_partitions; ++i) {
        std::size_t current_count = elements_per_partition;
        if (rem_mod > 0) {
            --rem_mod;
            ++current_count;
        }
        RandomAccessIterator current_end = current + current_count;
        callback(RandomAccessIterator(current),
                 RandomAccessIterator(current_end));
        current = current_end;
    }
}

template <typename RNGType, typename Iterator>
class RandomSkipIterator
    : boost::iterator_facade<RandomSkipIterator<RNGType, Iterator>,
                             typename std::iterator_traits<Iterator>::value,
                             std::forward_iterator_tag,
                             typename std::iterator_traits<Iterator>::reference,
                             std::ptrdiff_t> {
  public:
    RandomSkipIterator(RNGType* rng, Iterator iter, Iterator end,
                       double prob) noexcept
        : rng(rng), iter(iter), end(end), prob(prob) {}

    RandomSkipIterator() noexcept : rng(nullptr), iter(), end(), prob(0.0) {}

  private:
    friend class boost::iterator_core_access;

    typename std::iterator_traits<Iterator>::reference
    dereference() const noexcept {
        return *iter;
    }

    bool equal(const RandomSkipIterator& other) const noexcept {
        return iter == other.iter;
    }

    void increment() noexcept {
        std::geometric_distribution<std::size_t> jump_dist(prob);
        std::size_t dist = jump_dist(*rng) + 1;
        if (std::size_t(std::distance(iter, end)) <= dist) {
            iter = end;
        } else {
            std::advance(iter, dist);
        }
    }

    RNGType* rng;
    Iterator iter;
    Iterator end;
    double prob;
};

} // namespace sammy

#endif
