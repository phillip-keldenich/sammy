#ifndef SAMMY_UNIVERSE_SUBGRAPH_H_INCLUDED_
#define SAMMY_UNIVERSE_SUBGRAPH_H_INCLUDED_

#include "vertex_operations.h"
#include "shared_db_propagator.h"
#include "dynamic_bitset.h"
#include "thread_group.h"
#include "parallel_bit_filter.h"
#include "clique_or_indset_builder.h"
#include "pair_infeasibility_map.h"
#include "thread_interrupt.h"

namespace sammy {

/**
 * UniverseSubgraph represents a growable subgraph of
 * the universe of interactions with an edge between
 * two interactions if they are mutually exclusive.
 */
class UniverseSubgraph {
  public:
    UniverseSubgraph(const UniverseSubgraph& o) :
        propagator(o.propagator),
        infeasibility_map(o.infeasibility_map),
        vertices(o.vertices),
        matrix(o.matrix),
        vertices_with_literal(o.vertices_with_literal),
        vertex_index_map(o.vertex_index_map),
        degree(o.degree),
        parallel_bits(&o.parallel_bits.thread_group())
    {}

    UniverseSubgraph &operator=(const UniverseSubgraph& o) {
        propagator = o.propagator;
        infeasibility_map = o.infeasibility_map;
        vertices = o.vertices;
        matrix = o.matrix;
        vertices_with_literal = o.vertices_with_literal;
        vertex_index_map = o.vertex_index_map;
        degree = o.degree;
        return *this;
    }

    UniverseSubgraph(UniverseSubgraph&& o) noexcept :
        propagator(std::move(o.propagator)),
        infeasibility_map(o.infeasibility_map),
        vertices(std::move(o.vertices)),
        matrix(std::move(o.matrix)),
        vertices_with_literal(std::move(o.vertices_with_literal)),
        vertex_index_map(std::move(o.vertex_index_map)),
        degree(std::move(o.degree)),
        parallel_bits(&o.parallel_bits.thread_group())
    {}

    UniverseSubgraph &operator=(UniverseSubgraph&& o) noexcept {
        std::swap(propagator, o.propagator);
        std::swap(infeasibility_map, o.infeasibility_map);
        std::swap(vertices, o.vertices);
        std::swap(matrix, o.matrix);
        std::swap(vertices_with_literal, o.vertices_with_literal);
        std::swap(vertex_index_map, o.vertex_index_map);
        std::swap(degree, o.degree);
        return *this;
    }

    bool operator==(const UniverseSubgraph& o) const noexcept {
        return vertices == o.vertices &&
               matrix == o.matrix;
    }

    bool operator!=(const UniverseSubgraph& o) const noexcept {
        return !(*this == o);
    }

    UniverseSubgraph(ClauseDB* all_clauses, ThreadGroup<void>* thread_pool,
                     const PairInfeasibilityMap* infeasibility_map,
                     std::vector<Vertex> vertices)
        : propagator(all_clauses), infeasibility_map(infeasibility_map),
          vertices(std::move(vertices)),
          matrix(this->vertices.size(),
                 DynamicBitset(this->vertices.size(), false)),
          vertices_with_literal(2 * all_clauses->num_vars(),
                                DynamicBitset(this->vertices.size(), false)),
          parallel_bits(thread_pool)
    {
        p_build_vertex_index_map();
        p_build_vertices_with_literal();
        p_build_matrix();
    }

    SharedDBPropagator& get_propagator() noexcept { return propagator; }

    const SharedDBPropagator& get_propagator() const noexcept {
        return propagator;
    }

    DynamicBitset& matrix_row(std::size_t index) noexcept {
        return matrix[index];
    }

    const DynamicBitset& matrix_row(std::size_t index) const noexcept {
        return matrix[index];
    }

    std::size_t n() const noexcept { return vertices.size(); }

    const std::vector<std::size_t>& get_degrees() const noexcept {
        return degree;
    }

    const std::vector<Vertex>& vertex_set() const noexcept { return vertices; }

    Vertex vertex(std::size_t index) const noexcept { return vertices[index]; }

    bool has_vertex(Vertex v) const noexcept {
        return vertex_index_map.count(v);
    }

    std::size_t vertex_index(Vertex v) const { return vertex_index_map.at(v); }

    using CliqueBuilder = sammy::CliqueBuilder<UniverseSubgraph>;
    using IndependentSetBuilder = sammy::IndependentSetBuilder<UniverseSubgraph>;

    CliqueBuilder clique_builder() {
        return CliqueBuilder{this, &parallel_bits};
    }

    IndependentSetBuilder independent_set_builder() {
        return IndependentSetBuilder{this, &parallel_bits};
    }

    bool is_edge(Vertex v, Vertex w) const noexcept {
        return matrix[vertex_index(v)][vertex_index(w)];
    }

    bool is_edge(std::size_t i1, std::size_t i2) const noexcept {
        return matrix[i1][i2];
    }

    UniverseSubgraph restrict_to_degree_at_least(std::size_t min_degree) {
        std::vector<std::size_t> new_degrees = degree;
        std::vector<std::size_t> removal_queue;
        for (std::size_t i = 0, n = this->n(); i < n; ++i) {
            if (new_degrees[i] < min_degree) {
                removal_queue.push_back(i);
            }
        }
        std::size_t queue_pos = 0;
        while (queue_pos < removal_queue.size()) {
            std::size_t r = removal_queue[queue_pos++];
            for (std::size_t neighbor : matrix[r].ones()) {
                if (new_degrees[neighbor]-- == min_degree) {
                    removal_queue.push_back(neighbor);
                }
            }
        }
        std::sort(removal_queue.begin(), removal_queue.end());
        std::vector<Vertex> new_vertices;
        std::vector<DynamicBitset> new_matrix;
        std::vector<DynamicBitset> new_vertices_with_literal;
        removal_queue.push_back(this->n());
        p_walk_removals(removal_queue, new_vertices, new_matrix,
                        new_vertices_with_literal);
        return UniverseSubgraph(propagator, infeasibility_map,
                                std::move(new_vertices), std::move(new_matrix),
                                std::move(new_vertices_with_literal), parallel_bits.thread_group());
    }

    /**
     * @brief Add a bulk of new vertices.
     *
     * @param vertices
     */
    void add_vertices(const std::vector<Vertex>& new_vertices) {
        std::size_t old_n = vertices.size();
        std::size_t new_n = old_n + new_vertices.size();
        std::size_t res_size = 1;
        while (res_size <= new_n) {
            res_size <<= 1;
        }
        for (auto& row : matrix) {
            row.reserve(res_size);
            row.resize(new_n, false);
        }
        matrix.reserve(res_size);
        matrix.resize(new_n, DynamicBitset(new_n, false));
        for (auto& vwl : vertices_with_literal) {
            vwl.reserve(res_size);
            vwl.resize(new_n, false);
        }
        vertices.insert(vertices.end(), new_vertices.begin(),
                        new_vertices.end());
        for(std::size_t i = old_n; i < new_n; ++i) {
            vertex_index_map.try_emplace(vertices[i], i);
        }
        degree.resize(vertices.size(), 0);
        p_extend_vertices_with_literal(old_n, new_n);
        p_extend_matrix(old_n, new_n);
    }

    Var get_n_concrete() const noexcept {
        return infeasibility_map->get_n_concrete();
    }

    const PairInfeasibilityMap* get_infeasibility_map() const noexcept {
        return infeasibility_map;
    }

    void nonedge_to_edge(std::size_t index1, std::size_t index2) {
        matrix[index1][index2] = true;
        matrix[index2][index1] = true;
        ++degree[index1];
        ++degree[index2];
    }

    template<typename InputIterator>
    void add_new_vertices(InputIterator begin, InputIterator end)
    {
        std::vector<Vertex> vnew;
        std::copy_if(begin, end, std::back_inserter(vnew),
                     [&] (Vertex v) { return !has_vertex(v); });
        add_vertices(vnew);
    }

    template<typename VertexContainer>
    std::vector<std::size_t> to_indices(const VertexContainer& vertices) const noexcept {
        std::vector<std::size_t> result;
        result.reserve(vertices.size());
        std::transform(vertices.begin(), vertices.end(), std::back_inserter(result),
                       [&] (Vertex v) {return vertex_index(v);});
        return result;
    }

    /**
     * @brief Beside the edges implied by conflicts found in any case,
     *        we may also try to extend the information by
     *        checking, for each pair that is not already excluded,
     *        whether both vertices can be pushed simultaneously.
     *        This is MUCH more expensive than building the basic matrix.
     * @return The number of additional edges found.
     */
    std::size_t extend_matrix_by_propagation(bool interruptible = false) {
        if(is_extended) {
            return 0;
        }
        std::size_t count = 0, check_count = 0;
        for (std::size_t i = 1, nv = n(); i < nv; ++i) {
            Vertex v = vertices[i];
            reset_and_push_noresolve(propagator, v);
            DynamicBitset& row_i = matrix[i];
            for (std::size_t j = 0; j < i; ++j) {
                if (!row_i[j]) {
                    Vertex w = vertices[j];
                    if (!can_push(propagator, w)) {
                        row_i[j] = true;
                        matrix[j][i] = true;
                        ++degree[i];
                        ++degree[j];
                        ++count;
                    }
                }
            }
            if(interruptible && ++check_count % 1024 == 0) {
                throw_if_interrupted();
            }
        }
        is_extended = true;
        return count;
    }

    bool is_extended_by_propagation() const noexcept {
        return is_extended;
    }

    ThreadGroup<void>& thread_group() noexcept {
        return parallel_bits.thread_group();
    }

    std::size_t get_degree(std::size_t vertex_index) const noexcept {
        return degree[vertex_index];
    }

    std::size_t get_degree(Vertex v) const noexcept {
        return degree[vertex_index(v)];
    }

    BitsetOperationsBuffer* get_parallel_bits() noexcept {
        return &parallel_bits;
    }

  private:
    UniverseSubgraph(const SharedDBPropagator& propagator,
                     const PairInfeasibilityMap* infeasibility_map,
                     std::vector<Vertex> vertices,
                     std::vector<DynamicBitset> matrix,
                     std::vector<DynamicBitset> vertices_with_literal,
                     ThreadGroup<void>& thread_pool)
        : propagator(propagator), infeasibility_map(infeasibility_map),
          vertices(std::move(vertices)), matrix(std::move(matrix)),
          vertices_with_literal(std::move(vertices_with_literal)),
          degree(this->vertices.size(), 0),
          parallel_bits(&thread_pool)
    {
        this->propagator.reset_or_throw();
        p_build_vertex_index_map();
        for (std::size_t i = 0, n = this->n(); i < n; ++i) {
            degree[i] = this->matrix[i].count();
        }
    }

    void p_walk_removals_vertices(const std::vector<std::size_t>& removed,
                                  std::vector<Vertex>& nvertices) {
        std::size_t new_n = vertices.size() - removed.size();
        nvertices.reserve(new_n);
        auto rem_iter = removed.begin();
        for (std::size_t i = 0, nold = n(); i < nold; ++i) {
            if (i == *rem_iter) {
                ++rem_iter;
                continue;
            }
            nvertices.push_back(vertices[i]);
        }
    }

    void p_walk_removals_matrix(const std::vector<std::size_t>& removed,
                                std::vector<DynamicBitset>& nmatrix) {
        auto out_rem_iter = removed.begin();
        const std::size_t nold = n();
        std::size_t out_new_i = 0;
        for (std::size_t out_i = 0; out_i < nold; ++out_i) {
            if (out_i == *out_rem_iter) {
                ++out_rem_iter;
                continue;
            }
            auto& row = nmatrix[out_new_i];
            const auto& old_row = matrix[out_i];
            auto in_rem_iter = removed.begin();
            std::size_t in_new_i = 0;
            for (std::size_t in_i = 0; in_i < nold; ++in_i) {
                if (in_i == *in_rem_iter) {
                    ++in_rem_iter;
                    continue;
                }
                if (old_row[in_i]) {
                    row[in_new_i].set();
                }
                ++in_new_i;
            }
            ++out_new_i;
        }
    }

    void p_walk_removals_vwl(const std::vector<std::size_t>& removed,
                             std::vector<DynamicBitset>& nv_with_lit) {
        Lit nalit = 2 * propagator.db().num_vars();
        for (Lit l = 0; l < nalit; ++l) {
            auto rem_iter = removed.begin();
            std::size_t new_i = 0;
            const auto& old_row = vertices_with_literal[l];
            auto& new_row = nv_with_lit[l];
            for (std::size_t i = 0, n = this->n(); i < n; ++i) {
                if (i == *rem_iter) {
                    ++rem_iter;
                    continue;
                }
                if (old_row[i]) {
                    new_row[new_i].set();
                }
                ++new_i;
            }
        }
    }

    void p_walk_removals(const std::vector<std::size_t>& removed,
                         std::vector<Vertex>& nvertices,
                         std::vector<DynamicBitset>& nmatrix,
                         std::vector<DynamicBitset>& nv_with_lit) {
        p_walk_removals_vertices(removed, nvertices);
        auto nalit = 2 * propagator.db().num_vars();
        std::size_t new_n = nvertices.size();
        nmatrix.assign(new_n, DynamicBitset(new_n, false));
        nv_with_lit.assign(nalit, DynamicBitset(new_n, false));
        p_walk_removals_matrix(removed, nmatrix);
        p_walk_removals_vwl(removed, nv_with_lit);
    }

    void p_build_vertex_index_map() {
        vertex_index_map.reserve(vertices.size());
        for (std::size_t vi = 0, n = vertices.size(); vi < n; ++vi) {
            if (!vertex_index_map.try_emplace(vertices[vi], vi).second) {
                throw std::logic_error(
                    "Duplicate vertex given to UniverseSubgraph!");
            }
        }
    }

    void p_assert_symmetry() {
        for (std::size_t i = 1, n = this->n(); i < n; ++i) {
            const auto& row = matrix[i];
            for (std::size_t j = 0; j < i; ++j) {
                if (bool(row[j]) != bool(matrix[j][i])) {
                    std::cerr << "Asymmetry: matrix[" << i << "," << j
                              << "] != matrix[" << j << "," << i << "]!\n"
                              << std::flush;
                    std::abort();
                }
            }
        }
    }

    void p_build_matrix_row_range(std::size_t begin, std::size_t end, std::size_t old_n) {
        SharedDBPropagator local_prop = propagator;
        const std::size_t n = vertices.size();
        for (std::size_t vi = begin; vi != end; ++vi) {
            Vertex v = vertices[vi];
            auto& row = matrix[vi];
            reset_and_push_noresolve(local_prop, v);
            for (Lit lpos : local_prop.get_trail()) {
                Lit lneg = lit::negate(lpos);
                row |= vertices_with_literal[lneg];
            }
            if(is_extended) {
                // already checked this for old rows
                // during extension of existing rows
                for(std::size_t vj = 0; vj != old_n; ++vj) {
                    if(!row[vj]) {
                        row[vj] = matrix[vj][vi];
                    }
                }
                for(std::size_t vj = old_n; vj != vi; ++vj) {
                    if(row[vj]) continue;
                    Vertex w = vertices[vj];
                    if(!can_push(local_prop, w)) row[vj] = true;
                }
                for(std::size_t vj = vi + 1; vj < n; ++vj) {
                    if(row[vj]) continue;
                    Vertex w = vertices[vj];
                    if(!can_push(local_prop, w)) row[vj] = true;
                }
            }
            degree[vi] = row.count();
        }
    }

    void p_extend_matrix_row_range(std::size_t begin, std::size_t end,
                                   std::size_t start_from) {
        SharedDBPropagator local_prop = propagator;
        for (std::size_t vi = begin; vi != end; ++vi) {
            Vertex v = vertices[vi];
            reset_and_push_noresolve(local_prop, v);
            for (Lit lpos : local_prop.get_trail()) {
                Lit lneg = lit::negate(lpos);
                matrix[vi].binary_or(vertices_with_literal[lneg], start_from);
            }
            if(is_extended) {
                auto& row = matrix[vi];
                for(std::size_t vj = start_from, n = vertices.size(); vj != n; ++vj) {
                    if(row[vj]) continue;
                    Vertex w = vertices[vj];
                    if(!can_push(local_prop, w)) {
                        row[vj] = true;
                    }
                }
            }
            degree[vi] = matrix[vi].count();
        }
    }

    void p_build_matrix(std::size_t row_begin, std::size_t row_end) {
        propagator.reset_or_throw();
        std::size_t num_threads = std::thread::hardware_concurrency();
        std::size_t tp_threads = parallel_bits.thread_group().num_threads() + 1;
        std::size_t num_vertices = row_end - row_begin;
        num_threads = (std::min)(num_threads, tp_threads);
        num_threads = (std::min)(num_threads, num_vertices);
        if(num_threads == 1) {
            p_build_matrix_row_range(row_begin, row_end, row_begin);
            return;
        }
        std::unique_ptr<std::thread[]> builders =
            std::make_unique<std::thread[]>(num_threads);
        std::size_t vertices_per_thread = num_vertices / num_threads;
        for(std::size_t i = 0, cur = row_begin; i < num_threads; ++i, cur += vertices_per_thread) {
            std::size_t cend = cur + vertices_per_thread;
            if(i == num_threads - 1) cend = row_end;
            builders[i] = std::thread(
                [this, &row_begin](std::size_t b, std::size_t e) {
                    p_build_matrix_row_range(b, e, row_begin);
                },
                cur, cend);
        }
        for (std::size_t i = 0; i < num_threads; ++i) {
            builders[i].join();
        }
        static_cast<void>(&UniverseSubgraph::p_assert_symmetry);
        assert((p_assert_symmetry(),
                "Graph/adjacency matrix should be symmetric!"));
    }

    void p_build_matrix() {
        degree.resize(vertices.size(), 0);
        p_build_matrix(0, vertices.size());
    }

    void p_extend_existing_rows(std::size_t existing_rows, std::size_t begin_new) {
        std::size_t num_threads = std::thread::hardware_concurrency();
        std::size_t tp_threads = parallel_bits.thread_group().num_threads() + 1;
        num_threads = (std::min)(num_threads, tp_threads);
        num_threads = (std::min)(num_threads, existing_rows);
        if(num_threads == 1) {
            p_extend_matrix_row_range(0, existing_rows, begin_new);
            return;
        }
        std::unique_ptr<std::thread[]> builders =
            std::make_unique<std::thread[]>(num_threads);
        std::size_t vertices_per_thread = existing_rows / num_threads;
        for (std::size_t i = 0, cur = 0; i < num_threads;
             ++i, cur += vertices_per_thread)
        {
            std::size_t cend = cur + vertices_per_thread;
            if(i == num_threads - 1) cend = existing_rows;
            builders[i] = std::thread(
                [&] (std::size_t b, std::size_t e) {
                    p_extend_matrix_row_range(b, e, begin_new);
                }, cur, cend);
        }
        for(std::size_t i = 0; i < num_threads; ++i) {
            builders[i].join();
        }
    }

    void p_extend_matrix(std::size_t old_n, std::size_t new_n) {
        propagator.reset_or_throw();
        p_extend_existing_rows(old_n, old_n);
        p_build_matrix(old_n, new_n);
    }

    void p_extend_vertices_with_literal(std::size_t vindex_begin,
                                        std::size_t vindex_end) 
    {
        for (std::size_t vi = vindex_begin; vi != vindex_end; ++vi) {
            Vertex v = vertices[vi];
            reset_and_push_noresolve(propagator, v);
            for (Lit l : propagator.get_trail()) {
                vertices_with_literal[l][vi].set();
            }
        }
    }

    void p_build_vertices_with_literal() {
        p_extend_vertices_with_literal(0, vertices.size());
    }

    SharedDBPropagator propagator;
    const PairInfeasibilityMap* infeasibility_map;
    std::vector<Vertex> vertices;
    std::vector<DynamicBitset> matrix;
    std::vector<DynamicBitset> vertices_with_literal;
    VertexMapTo<std::size_t> vertex_index_map;
    std::vector<std::size_t> degree;
    mutable BitsetOperationsBuffer parallel_bits;
    bool is_extended = false;
};

}

#endif
