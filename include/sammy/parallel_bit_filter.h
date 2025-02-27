#ifndef SAMMY_PARALLEL_BIT_FILTER_H_INCLUDED_
#define SAMMY_PARALLEL_BIT_FILTER_H_INCLUDED_

#include "dynamic_bitset.h"
#include "thread_group.h"
#include <vector>

namespace sammy {

/**
 * Buffer for parallel bitset operations.
 */
class BitsetOperationsBuffer {
  public:
    explicit BitsetOperationsBuffer(ThreadGroup<void>* pool)
        : thread_pool(pool), bitsets(pool->num_threads() + 1, Bitset{}),
          used_threads(pool->num_threads() + 1, 0) {}

    void prepare() { used_threads.assign(thread_pool->num_threads() + 1, 0); }

    template <typename Combiner>
    Bitset& initialize(std::size_t index, const Bitset& from_bs,
                       Combiner&& combiner) {
        Bitset& s = bitsets[index];
        if (used_threads[index]) {
            // re-initialization (if one thread got two tasks)
            combiner(s, from_bs);
        } else {
            // first time this thread is used this run
            used_threads[index] = 1;
            s = from_bs;
        }
        return s;
    }

    template <typename Reducer>
    void reduce(Bitset& initial_and_result, Reducer&& reducer) {
        for (std::size_t i = 0, nt = bitsets.size(); i < nt; ++i) {
            if (!used_threads[i])
                continue;
            std::forward<Reducer>(reducer)(initial_and_result, bitsets[i]);
        }
    }

    ThreadGroup<void>& thread_group() { return *thread_pool; }

  private:
    ThreadGroup<void>* thread_pool;
    std::vector<Bitset> bitsets;
    std::vector<std::size_t> used_threads;
};

/**
 * Given an existing thread group, a range of bitsets y_1, ..., y_n,
 * and an initial set x, compute x - y_1 - ... - y_n, possibly in parallel
 * using the thread group of op_buffer, and write the result to
 * initial_and_result.
 */
template <typename BitsetsIterator>
inline void bitwise_filter(BitsetOperationsBuffer& op_buffer,
                           DynamicBitset& initial_and_result,
                           BitsetsIterator begin, BitsetsIterator end) {
    auto& tgroup = op_buffer.thread_group();
    auto nt = tgroup.num_threads() + 1;
    std::size_t num_sets = std::distance(begin, end);
    std::size_t total_size = initial_and_result.size() * num_sets / 8;
    if (nt <= 1 || num_sets < 4 * nt || total_size <= 1024 * 1024) {
        std::for_each(begin, end, [&](const DynamicBitset& bs) {
            initial_and_result -= bs;
        });
    } else {
        op_buffer.prepare();
        tgroup.context_function_parallel_foreach_iterator(
            begin, end,
            [&](std::size_t thread_index,
                const BitsetsIterator& first_iter) -> Bitset& {
                return op_buffer.initialize(
                    thread_index, *first_iter,
                    [](Bitset& bout, const Bitset& bin) { bout |= bin; });
            },
            [&](Bitset& output, std::size_t, const BitsetsIterator& iter) {
                output |= *iter;
            },
            [](Bitset&, std::size_t) {});
        op_buffer.reduce(
            initial_and_result,
            [](Bitset& in_out, const Bitset& in) { in_out -= in; });
    }
}

/**
 * Given an existing thread group, a range of bitsets y_1, ..., y_n,
 * and an initial set x, compute x - y_1 - ... - y_n, possibly in parallel
 * using the thread group tgroup, and write the result to initial_and_result.
 */
template <typename BitsetsIterator>
inline void bitwise_and(BitsetOperationsBuffer& op_buffer,
                        DynamicBitset& initial_and_result,
                        BitsetsIterator begin, BitsetsIterator end) {
    auto& tgroup = op_buffer.thread_group();
    auto nt = tgroup.num_threads() + 1;
    std::size_t num_sets = std::distance(begin, end);
    std::size_t total_size = initial_and_result.size() * num_sets / 8;
    if (nt <= 1 || num_sets < 4 * nt || total_size <= 1024 * 1024) {
        std::for_each(begin, end, [&](const DynamicBitset& bs) {
            initial_and_result &= bs;
        });
    } else {
        op_buffer.prepare();
        tgroup.context_function_parallel_foreach_iterator(
            begin, end,
            [&](std::size_t thread_index,
                const BitsetsIterator& first_iter) -> Bitset& {
                return op_buffer.initialize(
                    thread_index, *first_iter,
                    [](Bitset& bout, const Bitset& bin) { bout &= bin; });
            },
            [&](Bitset& output, std::size_t, const BitsetsIterator& iter) {
                output &= *iter;
            },
            [](Bitset&, std::size_t) {});
        op_buffer.reduce(
            initial_and_result,
            [](Bitset& in_out, const Bitset& in) { in_out &= in; });
    }
}

} // namespace sammy

#endif
