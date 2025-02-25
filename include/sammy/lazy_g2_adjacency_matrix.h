#ifndef SAMMY_LAZY_G2_ADJACENCY_MATRIX_H_INCLUDED_
#define SAMMY_LAZY_G2_ADJACENCY_MATRIX_H_INCLUDED_

#include "literals.h"
#include "dynamic_bitset.h"
#include "pair_infeasibility_map.h"
#include "shared_db_propagator.h"
#include "vertex_operations.h"
#include "thread_interrupt.h"


namespace sammy {

/**
 * Adjacency matrix of a subgraph of G2 on a
 * fixed subset of all possible vertices; this
 * may either be a reduced vertex set of the 
 * full universe or a real subset, e.g., during LNS.
 */
class LazyG2AdjacencyMatrix {
  public:
    using Vertices = std::vector<Vertex>;
    using Indices = std::vector<std::size_t>;
    using VertexMap = VertexMapTo<std::size_t>;

    LazyG2AdjacencyMatrix(
        std::vector<Vertex> considered_vertices,
        ClauseDB& clause_db, 
        std::size_t n_concrete
    );

    /**
     * Get a list of vertex indices containing the given concrete literal l.
     */
    const Indices& vertices_with_concrete_literal(Lit l) const noexcept {
        assert(l < m_vertices_containing_concrete_literal.size());
        return m_vertices_containing_concrete_literal[l];
    }

    /**
     * A vector-based map of:
     *  - concrete literal l -> vertex indices implying l by UP.
     */
    const std::vector<DynamicBitset>& vertices_implying_literal() const noexcept {
        return m_definitive_nonedges;
    }

    /**
     * Get a bitset of vertices that imply the given literal l.
     */
    const DynamicBitset& vertices_implying_literal(Lit l) const noexcept {
        assert(l < m_vertices_implying_literals.size());
        return m_vertices_implying_literals[l];
    }

    /**
     * Get a vector-based map of:
     * concrete literal l -> vertex indices containing l.
     */
    const std::vector<Indices>& vertices_with_concrete_literal() const noexcept {
        return m_vertices_containing_concrete_literal;
    }

    /**
     * Get a reference to the list of all vertices.
     * Vertex indices point into this list.
     */
    const Vertices& all_vertices() const noexcept {
        return m_all_vertices;
    }

    /**
     * Get the number of vertices.
     */
    std::size_t n() const noexcept {
        return m_all_vertices.size();
    }

    /**
     * Get a vertex by its index.
     */
    Vertex vertex(std::size_t index) const noexcept {
        return m_all_vertices[index];
    }

    /**
     * Get a temporary propagator reference.
     * Clear after each use.
     */
    SharedDBPropagator& temp_propagator() noexcept {
        return m_empty_propagator;
    }

    /**
     * Get a reference to the given row of the adjacency matrix.
     */
    DynamicBitset &row(std::size_t index) noexcept {
        return m_adjacency_matrix[index];
    }

    /**
     * Get a const reference to the given row of the adjacency matrix.
     */
    const DynamicBitset &row(std::size_t index) const noexcept {
        return m_adjacency_matrix[index];
    }

    /**
     * Check if there definitely is no edge between 
     * vertices with index i and j; if necessary, does
     * the pairwise propagation to confirm this,
     * and cache the result.
     */
    bool is_definitive_nonedge(std::size_t i, std::size_t j) noexcept {
        if(m_adjacency_matrix[i][j]) return false;
        if(m_definitive_nonedges[i][j]) return true;
        if(push_vertex_pair(m_empty_propagator, m_all_vertices[i], m_all_vertices[j])) {
            m_empty_propagator.reset_to_zero();
            m_definitive_nonedges[i][j].set();
            m_definitive_nonedges[j][i].set();
            return true;
        } else {
            m_adjacency_matrix[i][j].set();
            m_adjacency_matrix[j][i].set();
            return false;
        }
    }

    std::size_t num_vars() const noexcept {
        return m_n_all;
    }

    std::size_t num_concrete() const noexcept {
        return m_n_concrete;
    }

    void nonedge_to_edge(std::size_t v, std::size_t w) noexcept {
        m_adjacency_matrix[v][w].set();
        m_adjacency_matrix[w][v].set();
        m_definitive_nonedges[v][w].reset();
        m_definitive_nonedges[w][v].reset();
    }

    std::size_t index_of(Vertex v) const noexcept {
        return m_vertex_indices.at(v);
    }

    template<typename Range>
    std::vector<std::size_t> indices_of(Range&& r) const {
        std::vector<std::size_t> result;
        for(Vertex v : r) {
            result.push_back(index_of(v));
        }
        return result;
    }

    template<typename Range>
    std::vector<Vertex> vertices_of(Range&& r) const {
        std::vector<Vertex> result;
        for(std::size_t i : r) {
            result.push_back(vertex(i));
        }
        return result;
    }

  private:
    void p_initialize_vertices_implying_literal();
    void p_initialize_vertices_containing_concrete_literal();
    void p_initialize_matrix_from_implied_literals();

    /**
     * The number of all literals.
     */
    std::size_t m_n_all;

    /**
     * The number of concrete literals.
     */
    std::size_t m_n_concrete;

    /**
     * A propagator we can reuse.
     * After use, it must be emptied again.
     */
    SharedDBPropagator m_empty_propagator;

    /**
     * The set of all vertices; in other places, they 
     * are referenced by their index in this vector.
     */
    std::vector<Vertex> m_all_vertices;

    /**
     * For each concrete literal, the list of vertices containing that literal.
     */
    std::vector<std::vector<std::size_t>> m_vertices_containing_concrete_literal;

    /**
     * For each (concrete or non-concrete) literal, the set of vertices
     * _implying_ that literal by propagation.
     */
    std::vector<DynamicBitset> m_vertices_implying_literals;

    /**
     * Adjacency matrix that contains ones for definite edges;
     * zeros might turn out to be edges in G2, but not in the
     * current subgraph.
     */
    std::vector<DynamicBitset> m_adjacency_matrix;

    /**
     * A 'true' means that there is definitely no edge,
     * i.e., propagation of interactions m_all_vertices[i]
     * and m_all_vertices[j] at the same time found no conflict.
     */
    std::vector<DynamicBitset> m_definitive_nonedges;

    /**
     * Get the index of a given vertex by hash-table lookup.
     */
    VertexMap m_vertex_indices;
};

LazyG2AdjacencyMatrix::LazyG2AdjacencyMatrix(
    std::vector<Vertex> considered_vertices,
    ClauseDB& clause_db, 
    std::size_t n_concrete
) :
    m_n_all(clause_db.num_vars()),
    m_n_concrete(n_concrete),
    m_empty_propagator(&clause_db),
    m_all_vertices(std::move(considered_vertices))
{
    m_vertex_indices.reserve(m_all_vertices.size());
    for(std::size_t i = 0, n = m_all_vertices.size(); i < n; ++i) {
        m_vertex_indices[m_all_vertices[i]] = i;
    }
    p_initialize_vertices_containing_concrete_literal();
    p_initialize_vertices_implying_literal();
    p_initialize_matrix_from_implied_literals();
}

/**
 * Initialize the list of vertices containing each concrete literal.
 */
void LazyG2AdjacencyMatrix::p_initialize_vertices_containing_concrete_literal() {
    m_vertices_containing_concrete_literal.assign(2 * m_n_concrete, Indices{});
    for(std::size_t vi = 0, vn = m_all_vertices.size(); vi < vn; ++vi) {
        Vertex v = m_all_vertices[vi];
        assert(lit::var(v.first) < m_n_concrete);
        assert(lit::var(v.second) < m_n_concrete);
        m_vertices_containing_concrete_literal[v.first].push_back(vi);
        m_vertices_containing_concrete_literal[v.second].push_back(vi);
    }
}

/**
 * Initialize the vertices_implying_literals member.
 */
void LazyG2AdjacencyMatrix::p_initialize_vertices_implying_literal() {
    m_vertices_implying_literals.reserve(2 * m_n_all);
    std::generate_n(std::back_inserter(m_vertices_implying_literals), 2 * m_n_all,
                    [&] () { return DynamicBitset(m_all_vertices.size(), false); });
    for(std::size_t vi = 0, n = m_all_vertices.size(); vi < n; ++vi) {
        assert(m_empty_propagator.get_current_level() == 0);
        Vertex v = m_all_vertices[vi];
        if(push_vertex(m_empty_propagator, v) < 0)
            throw std::logic_error("Infeasible interaction in m_all_vertices!");
        for(Lit l : m_empty_propagator.get_trail()) {
            m_vertices_implying_literals[l][vi].set();
        }
        m_empty_propagator.reset_to_zero();
    }
}

void LazyG2AdjacencyMatrix::p_initialize_matrix_from_implied_literals() {
    m_adjacency_matrix.assign(m_all_vertices.size(), DynamicBitset(m_all_vertices.size(), false));
    m_definitive_nonedges.assign(m_all_vertices.size(), DynamicBitset(m_all_vertices.size(), false));
    assert(m_adjacency_matrix.size() == m_all_vertices.size());
    assert(m_adjacency_matrix.empty() || 
           m_adjacency_matrix[0].size() == m_all_vertices.size());
    for(std::size_t vi = 0, n = m_all_vertices.size(); vi < n; ++vi) {
        assert(m_empty_propagator.get_current_level() == 0);
        Vertex v = m_all_vertices[vi];
        DynamicBitset& row = m_adjacency_matrix[vi];
        push_vertex(m_empty_propagator, v);
        for(Lit lpos : m_empty_propagator.get_trail()) {
            row |= m_vertices_implying_literals[lit::negate(lpos)];
        }
        m_empty_propagator.reset_to_zero();
        if(vi % 16 == 15) throw_if_interrupted();
    }
}

}

#endif
