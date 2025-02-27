#ifndef SAMMY_CLIQUE_OR_INDSET_BUILDER_H_INCLUDED_
#define SAMMY_CLIQUE_OR_INDSET_BUILDER_H_INCLUDED_

#include "dynamic_bitset.h"
#include "literals.h"
#include "parallel_bit_filter.h"
#include <boost/iterator/transform_iterator.hpp>

namespace sammy {

template <typename GraphType, bool BuildingClique>
class CliqueOrIndependentSetBuilder {
  public:
    explicit CliqueOrIndependentSetBuilder(
        GraphType* graph, BitsetOperationsBuffer* parallel_bits)
        : graph(graph), parallel_bits(parallel_bits),
          possible_vertices(graph->n(), true) {}

    void reset_vertices() {
        auto g_n = graph->n();
        if (possible_vertices.size() != g_n) {
            possible_vertices.resize(g_n);
        }
        possible_vertices.set();
        clique_vertices.clear();
    }

    void reset_vertices(const std::vector<Vertex>& initial) {
        reset_vertices();
        if (initial.empty()) {
            return;
        }
        auto transform_to_index = [this](Vertex v) {
            return graph->vertex_index(v);
        };
        std::transform(initial.begin(), initial.end(),
                       std::back_inserter(clique_vertices), transform_to_index);
        for (std::size_t ci : clique_vertices) {
            possible_vertices[ci] = false;
        }
        auto index_to_bitset = [&](std::size_t i) -> const DynamicBitset& {
            return graph->matrix_row(i);
        };
        auto bitset_begin =
            boost::make_transform_iterator<decltype(index_to_bitset)>(
                clique_vertices.begin());
        auto bitset_end =
            boost::make_transform_iterator<decltype(index_to_bitset)>(
                clique_vertices.end());
        if constexpr (BuildingClique) {
            sammy::bitwise_and(*parallel_bits, possible_vertices, bitset_begin,
                               bitset_end);
        } else {
            sammy::bitwise_filter(*parallel_bits, possible_vertices,
                                  bitset_begin, bitset_end);
        }
    }

    template <typename RNG> void randomly_make_maximal(RNG& rng) {
        for (;;) {
            bool success = false;
            std::uniform_int_distribution<std::size_t> indices(0,
                                                               graph->n() - 1);
            for (std::size_t trial_count = 0; trial_count < 100; ++trial_count)
            {
                auto index = indices(rng);
                if (possible_vertices[index]) {
                    add_vertex(index);
                    success = true;
                    break;
                }
            }
            if (!success) {
                possible_buffer.clear();
                std::copy(possible_vertices.ones_begin(),
                          possible_vertices.ones_end(),
                          std::back_inserter(possible_buffer));
                while (!possible_buffer.empty()) {
                    std::uniform_int_distribution<std::size_t> iindex_dist(
                        0, possible_buffer.size() - 1);
                    auto iindex = iindex_dist(rng);
                    auto vindex = possible_buffer[iindex];
                    add_vertex(vindex);
                    possible_buffer.erase(
                        std::remove_if(possible_buffer.begin(),
                                       possible_buffer.end(),
                                       [&](std::size_t v) {
                                           return !possible_vertices[v];
                                       }),
                        possible_buffer.end());
                }
                return;
            }
        }
    }

    void add_vertex(std::size_t index) {
        clique_vertices.push_back(index);
        possible_vertices[index] = false;
        if (BuildingClique) {
            possible_vertices &= graph->matrix_row(index);
        } else {
            possible_vertices -= graph->matrix_row(index);
        }
    }

    void add_vertex(Vertex v) { add_vertex(graph->vertex_index(v)); }

    std::size_t size() const noexcept { return clique_vertices.size(); }

    const std::vector<std::size_t>& get_indices() const noexcept {
        return clique_vertices;
    }

    std::vector<Vertex> get_vertices() const {
        std::vector<Vertex> result;
        result.reserve(clique_vertices.size());
        std::transform(clique_vertices.begin(), clique_vertices.end(),
                       std::back_inserter(result),
                       [&](std::size_t i) { return graph->vertex(i); });
        return result;
    }

    template <typename VertexIterator>
    void greedily_extend_on(VertexIterator begin, VertexIterator end) {
        for (auto v : IteratorRange(begin, end)) {
            auto index = graph->vertex_index(v);
            if (possible_vertices[index]) {
                add_vertex(index);
            }
        }
    }

    template <typename VertexIndexIterator>
    void greedily_extend_on_indices(VertexIndexIterator begin,
                                    VertexIndexIterator end) {
        for (auto index : IteratorRange(begin, end)) {
            if (possible_vertices[index]) {
                add_vertex(index);
            }
        }
    }

    template <typename WeightedIndexIterator>
    double greedily_extend_weighted(WeightedIndexIterator begin,
                                    WeightedIndexIterator end) {
        double total = 0.0;
        for (auto [index, weight] : IteratorRange(begin, end)) {
            if (possible_vertices[index]) {
                add_vertex(index);
                total += weight;
            }
        }
        return total;
    }

  private:
    GraphType* graph;
    BitsetOperationsBuffer* parallel_bits;
    std::vector<std::size_t> clique_vertices;
    DynamicBitset possible_vertices;
    std::vector<std::size_t> possible_buffer;
};

template <typename GraphType>
using CliqueBuilder = CliqueOrIndependentSetBuilder<GraphType, true>;
template <typename GraphType>
using IndependentSetBuilder = CliqueOrIndependentSetBuilder<GraphType, false>;

} // namespace sammy

#endif
