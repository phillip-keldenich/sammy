#ifndef SAMMY_ALGORITHM_EX_H_INCLUDED_
#define SAMMY_ALGORITHM_EX_H_INCLUDED_

#include "range.h"
#include <algorithm>
#include <boost/iterator/counting_iterator.hpp>

namespace sammy {

template <typename ForwardIterator, typename Callable>
static inline ForwardIterator
find_pair_if(ForwardIterator begin, ForwardIterator end, Callable&& callable) {
    ForwardIterator last = begin;
    if (begin != end) {
        for (++begin; begin != end; ++begin) {
            if (callable(*last, *begin)) {
                return last;
            }
            last = begin;
        }
    }
    return end;
}

template <typename InputIterator, typename OutputIterator,
          typename Transformation, typename Predicate>
static inline OutputIterator
copy_transformed_if(InputIterator begin, InputIterator end, OutputIterator out,
                    Transformation&& transform, Predicate&& predicate) {
    begin = std::find_if(begin, end, std::forward<Predicate>(predicate));
    while (begin != end) {
        *out = transform(*begin);
        ++begin;
        ++out;
        begin = std::find_if(begin, end, std::forward<Predicate>(predicate));
    }
    return out;
}

template <typename Container>
static inline void sort_unique(Container& container) {
    std::sort(container.begin(), container.end());
    container.erase(std::unique(container.begin(), container.end()),
                    container.end());
}

template <typename SetType, typename Container>
static inline void nosort_unique(Container& container) {
    using Reference = typename Container::const_reference;
    SetType set;
    auto new_end = std::remove_if(
        container.begin(), container.end(),
        [&set](Reference elem) { return !set.insert(elem).second; });
    container.erase(new_end, container.end());
}

template <typename Container, typename InputIterator>
static inline void add_and_make_unique(Container& container,
                                       InputIterator begin, InputIterator end) {
    container.insert(container.end(), begin, end);
    std::sort(container.begin(), container.end());
    container.erase(std::unique(container.begin(), container.end()),
                    container.end());
}

template <typename ForwardIterator, typename IndexForwardIterator>
static inline ForwardIterator remove_indices(ForwardIterator begin,
                                             ForwardIterator end,
                                             IndexForwardIterator remove_begin,
                                             IndexForwardIterator remove_end) {
    if (remove_begin == remove_end)
        return end;
    auto first_removed = *remove_begin++;
    ForwardIterator out = begin;
    std::advance(out, first_removed);
    ForwardIterator in = out;
    ++in;
    auto current_index = first_removed;
    ++current_index;
    while (remove_begin != remove_end) {
        for (auto next_removed = *remove_begin; current_index < next_removed;
             ++current_index, ++in, ++out)
        {
            *out = std::move(*in);
        }
        ++in;
        ++current_index;
        ++remove_begin;
    }
    for (; in != end; ++in, ++out) {
        *out = std::move(*in);
    }
    return out;
}

template <typename Range> static inline auto vector(const Range& range) {
    using std::begin;
    using std::end;
    using Iterator = decltype(begin(range));
    using ValueType = typename std::iterator_traits<Iterator>::value_type;
    return std::vector<ValueType>{begin(range), end(range)};
}

template <typename Iterator>
static inline auto vector(Iterator begin, Iterator end) {
    using ValueType = typename std::iterator_traits<Iterator>::value_type;
    return std::vector<ValueType>{begin, end};
}

template <typename IntType>
static inline auto range(IntType begin, IntType end) {
    end = (std::max)(begin, end);
    return IteratorRange{boost::counting_iterator<IntType>(begin),
                         boost::counting_iterator<IntType>(end)};
}

template <typename IntType> static inline auto range(IntType end) {
    return range(IntType(0), end);
}

// helpers for std::variant visitation
template <class... Ts> struct overloaded : Ts... {
    using Ts::operator()...;
};
template <class... Ts> overloaded(Ts...) -> overloaded<Ts...>;

} // namespace sammy

#endif
