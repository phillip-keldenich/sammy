#ifndef SAMMY_PAIR_INFEASIBILITY_MAP_H_
#define SAMMY_PAIR_INFEASIBILITY_MAP_H_

#include "clause_db.h"
#include "cuda_iteration.h"
#include "dynamic_bitset.h"
#include "literals.h"
#include "rng.h"
#include "shared_db_propagator.h"
#include "thread_group.h"
#include "shallow_mem_estimate.h"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <memory>
#include <new>

namespace sammy {

using Index = Lit;

inline Bitset make_bitset(std::size_t length, bool value = false) {
    return Bitset(length, value);
}

inline Bitset make_larger_bitset(const Bitset& from) {
    Bitset bs;
    bs.reserve(from.size() * 2);
    bs = from;
    return bs;
}

#ifdef SAMMY_CUDA_SUPPORTED
static inline void cuda_extract_feasibilities(
    Var n_concrete, std::vector<DynamicBitset>& definitely_feasible_matrix,
    const std::vector<DynamicBitset>& literals_in_class) {
    if (literals_in_class.empty())
        return;

    const Lit nclit = 2 * n_concrete;
    if (definitely_feasible_matrix.size() != nclit ||
        definitely_feasible_matrix[0].size() != nclit)
    {
        throw std::invalid_argument(
            "definitely_feasible_matrix has wrong size");
    }
    std::size_t u32_per_bitset = literals_in_class[0].blocks().size() *
                                 sizeof(DynamicBitset::Block) /
                                 sizeof(std::uint32_t);
    std::vector<std::uint32_t> host_prepare_buffer;
    host_prepare_buffer.reserve(literals_in_class.size() * u32_per_bitset);
    for (const DynamicBitset& bset : literals_in_class) {
        for (DynamicBitset::Block b : bset.blocks()) {
            detail::to_prepare_buffer(b, host_prepare_buffer);
        }
    }
    detail::CUDADevicePointer<const std::uint32_t> device_bit_data(
        host_prepare_buffer);
    detail::CUDADevicePointer<std::uint32_t> device_output_buffer(
        nclit * u32_per_bitset);
    Var v = 0;
    while (v < n_concrete) {
        Var num_rows = n_concrete - v;
        if (num_rows > 2 * detail::GOAL_ROWS_PER_CALL()) {
            num_rows = 2 * detail::GOAL_ROWS_PER_CALL();
        }
        detail::call_cuda_extract_kernel(device_bit_data.get(), u32_per_bitset,
                                         device_output_buffer.get(), nclit, v,
                                         num_rows, literals_in_class.size());
        v += num_rows;
    }
    std::vector<std::uint32_t> host_output =
        device_output_buffer.get_host_copy();
    for (Lit l = 0; l < nclit; ++l) {
        definitely_feasible_matrix[l].binary_or(
            &host_output[l * u32_per_bitset]);
    }
}
#endif

/**
 * Matrix that stores pairs that are
 * known-feasible or known-infeasible.
 */
class PairInfeasibilityMap {
  public:
    using Row = Bitset;
    using Matrix = std::vector<Row>;

  private:
    std::size_t num_vars;
    Matrix m_matrix;
    Matrix m_def_feasible;
    std::vector<Lit> m_incorporate_buffer;

    explicit PairInfeasibilityMap(std::size_t n_concrete,
                                  const std::vector<std::vector<bool>>& bits)
        : num_vars(n_concrete), m_matrix(), m_def_feasible() {
        m_matrix.reserve(bits.size());
        m_def_feasible.reserve(bits.size());
        for (std::size_t i = 0; i < bits.size(); ++i) {
            m_matrix.emplace_back(bits[i]);
            m_def_feasible.emplace_back(m_matrix.back());
            m_def_feasible.back().flip();
        }
    }

  public:
    std::size_t get_memory_size() const noexcept {
        std::size_t result = sizeof(PairInfeasibilityMap) +
            shallow_memory_estimate(m_matrix) +
            shallow_memory_estimate(m_def_feasible) +
            shallow_memory_estimate(m_incorporate_buffer);
        if(!m_matrix.empty()) {
            result += m_matrix.size() * shallow_memory_estimate(m_matrix[0]);
        }
        if(!m_def_feasible.empty()) {
            result += m_def_feasible.size() * shallow_memory_estimate(m_def_feasible[0]);
        }
        return result;
    }

    std::size_t get_n_concrete() const noexcept { return num_vars; }

    explicit PairInfeasibilityMap(std::size_t n_concrete)
        : num_vars(n_concrete),
          m_matrix(2 * num_vars, Bitset(2 * num_vars, false)),
          m_def_feasible(2 * num_vars, Bitset(2 * num_vars, false)) {
        Lit l = 0;
        for (auto& r : m_matrix) {
            r[lit::negate(l)] = true;
            ++l;
        }
    }

    Row& operator[](Lit l) noexcept { return m_matrix[l]; }

    const Row& operator[](Lit l) const noexcept { return m_matrix[l]; }

    void literal_infeasible(Lit l) noexcept {
        const auto nl = 2 * num_vars;
        m_matrix[l].set();
        for (Lit i = 0; i < nl; ++i) {
            m_matrix[i][l] = true;
        }
    }

    void literal_pair_infeasible(Lit l1, Lit l2) noexcept {
        m_matrix[l1][l2] = true;
        m_matrix[l2][l1] = true;
    }

    bool operator()(Lit lmin, Lit lmax) const noexcept {
        return m_matrix[lmin][lmax];
    }

    bool is_definitely_feasible(Lit lmin, Lit lmax) const noexcept {
        return m_def_feasible[lmin][lmax];
    }

    void set_definitely_feasible(Lit lmin, Lit lmax) noexcept {
        m_def_feasible[lmin][lmax] = true;
        m_def_feasible[lmax][lmin] = true;
    }

    std::vector<std::vector<bool>> export_bits() const {
        std::vector<std::vector<bool>> result;
        result.reserve(m_matrix.size());
        for (const auto& row : m_matrix) {
            result.emplace_back(static_cast<std::vector<bool>>(row));
        }
        return result;
    }

    static PairInfeasibilityMap
    import_bits(const std::vector<std::vector<bool>>& bits) {
        if (bits.size() == 0)
            throw std::runtime_error("Empty matrix in import_bits!");
        if (bits.size() % 2 == 1)
            throw std::runtime_error("Odd matrix in import_bits!");
        if (bits[0].size() != bits.size())
            throw std::runtime_error("Non-square matrix in import_bits!");
        std::size_t num_concrete = bits.size() / 2;
        return PairInfeasibilityMap{num_concrete, bits};
    }

    void incorporate_complete_class(const SharedDBPropagator& propagator) {
        const Lit nc_lit = 2 * num_vars;
        m_incorporate_buffer.clear();
        std::copy_if(propagator.get_trail().begin(),
                     propagator.get_trail().end(),
                     std::back_inserter(m_incorporate_buffer),
                     [&](Lit l) { return l < nc_lit; });
        std::sort(m_incorporate_buffer.begin(), m_incorporate_buffer.end());
        for (auto i = m_incorporate_buffer.begin(),
                  e = m_incorporate_buffer.end();
             i != e; ++i)
        {
            auto& row = m_def_feasible[*i];
            for (auto j = m_incorporate_buffer.begin(); j != i; ++j) {
                row[*j] = true;
            }
            for (auto j = std::next(i); j != e; ++j) {
                row[*j] = true;
            }
        }
    }

    void incorporate_complete_classes(
        const std::vector<DynamicBitset>& literals_in_class,
        ThreadGroup<void>& tpool) {
        auto handle_var = [&](Var v) {
            Lit plit = lit::positive_lit(v);
            DynamicBitset& positive_row = m_def_feasible[plit];
            DynamicBitset& negative_row = m_def_feasible[lit::negative_lit(v)];
            for (const DynamicBitset& bset : literals_in_class) {
                DynamicBitset& out_row =
                    (bset[plit] ? positive_row : negative_row);
                out_row |= bset;
            }
        };

        tpool.parallel_foreach_iterator(Var(0), Var(num_vars), handle_var);
    }

#ifdef SAMMY_CUDA_SUPPORTED
    void cuda_incorporate_complete_classes(
        const std::vector<DynamicBitset>& literals_in_class,
        ThreadGroup<void>& tpool) {
        try {
            cuda_extract_feasibilities(num_vars, m_def_feasible,
                                       literals_in_class);
            return;
        } catch (const CUDAError& err) {
            std::cerr << "CUDA error in cuda_incorporate_complete_classes: "
                      << err.what() << std::endl;
            had_cuda_error(err);
        }
        incorporate_complete_classes(literals_in_class, tpool);
    }
#endif

    /*
        template<typename PropIterator>
        void incorporate_complete_classes(PropIterator begin, PropIterator end,
       ThreadGroup<void>& tpool) { std::vector<Lit> prop_trails; std::size_t
       nclasses = std::size_t(end - begin); const Var ncvars = num_vars; const
       Lit nclits = 2 * ncvars; prop_trails.reserve(nclasses * num_vars);
            std::for_each(begin, end, [&] (const SharedDBPropagator& prop) {
                const auto& t = prop.get_trail();
                std::copy_if(t.begin(), t.end(),
       std::back_inserter(prop_trails),
                             [&] (Lit l) { return l < nclits; });
            });
            auto trail_range = [&] (std::size_t class_index) {
                auto tbegin = prop_trails.begin() + (class_index * ncvars);
                return std::make_pair(tbegin, tbegin + ncvars);
            };
            auto handle_var = [&] (Var v) {
                const Lit p = lit::positive_lit(v), n = lit::negative_lit(v);
                for(std::size_t iclass = 0; iclass < nclasses; ++iclass) {
                    Row& out = begin[iclass].is_true(p) ? m_def_feasible[p] :
       m_def_feasible[n]; auto [tbegin, tend] = trail_range(iclass);
                    std::for_each(tbegin, tend, [&] (Lit l) {out[l].set();});
                }
            };
            tpool.parallel_foreach_iterator(Var(0), ncvars, handle_var);
        } */

    std::size_t total_memory_usage() const noexcept {
        return 2 * (m_matrix.capacity() * sizeof(Row) +
                    m_matrix[0].bytes_used() * m_matrix.size()) +
               sizeof(PairInfeasibilityMap);
    }

    std::size_t count_vertices() const noexcept {
        auto res = std::transform_reduce(
            m_def_feasible.begin(), m_def_feasible.end(), std::size_t(0),
            std::plus<>{}, [](const Row& r) { return r.count(); });
        res /= 2;
        return res;
    }

    std::vector<Vertex>
    collect_vertices(std::size_t reserve_size = 0) const noexcept {
        std::vector<Vertex> result;
        if (reserve_size)
            result.reserve(reserve_size);
        for (Lit i = 0, nclit = 2 * num_vars; i < nclit - 2; ++i) {
            std::transform(
                m_def_feasible[i].ones_from_begin(i + 1),
                m_def_feasible[i].ones_end(), std::back_inserter(result),
                [&](std::size_t lmax) { return Vertex{i, Lit(lmax)}; });
        }
        return result;
    }

    std::vector<Vertex> sample_vertices(std::size_t target_size,
                                        std::size_t total_size = 0) {
        if (!total_size) {
            total_size = count_vertices();
        }
        double sample_prob = double(target_size) / total_size;
        sample_prob = std::min(1.0, sample_prob);
        std::geometric_distribution<std::size_t> dist{sample_prob};
        auto& rng = sammy::rng();
        std::size_t skip_count = dist(rng);
        std::vector<Vertex> result;
        result.reserve(std::size_t(1.1 * target_size));
        for (Lit i = 0, nclit = 2 * num_vars; i < nclit - 2; ++i) {
            for (std::size_t j : m_def_feasible[i].ones_from(i + 1)) {
                if (!skip_count) {
                    skip_count = dist(rng);
                    result.emplace_back(i, Lit(j));
                } else {
                    --skip_count;
                }
            }
        }
        return result;
    }
};

} // namespace sammy

#endif
