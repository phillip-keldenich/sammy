#ifndef SAMMY_FAST_CLIQUE_H_INCLUDED_
#define SAMMY_FAST_CLIQUE_H_INCLUDED_

#include "literals.h"
#include "rng.h"
#include "thread_group.h"
#include "shared_db_propagator.h"

namespace sammy {

/**
 * @brief Quickly compute some maximal clique (i.e. a set of mutually exclusive
 * literal pairs) on a subgraph of vertices (i.e., literal pairs). The vertices
 * are assumed to be feasible if they can be pushed to a propagator; normally,
 * with vertices extracted from a sample, this should be a non-issue. Running
 * this heuristic should be feasible even if storing/computing the complete
 * subgraph is too expensive. Mainly useful to quickly get lower bounds.
 */
class FastCliqueBuilder {
  public:
    explicit FastCliqueBuilder(SharedDBPropagator prop)
        : m_prop(std::move(prop)) {}

    std::vector<Vertex> compute_clique(Vertex start_vertex,
                                       const std::vector<Vertex>& subgraph) 
    {
        p_check_prop();
        m_subg = subgraph;
        if (!p_can_push(start_vertex))
            return {};
        p_filter_invalid();
        return p_compute_clique(start_vertex);
    }

    /**
     * @brief Compute a clique starting from a set of vertices. All vertices
     *        (starting vertices and subgraph) are assumed to be feasible interactions.
     */
    template<typename StartVerticesIterator> std::vector<Vertex> 
        compute_clique_known_valid(StartVerticesIterator begin, StartVerticesIterator end,
                                   const std::vector<Vertex>& subgraph)
    {
        std::vector<Vertex> result;
        p_check_prop();
        m_subg = subgraph;
        for(auto v : IteratorRange{begin, end}) {
            result.push_back(v);
            p_filter_nonneighbors(v);
        }
        while(!m_subg.empty()) {
            Vertex v = p_select_next_vertex();
            result.push_back(v);
            p_filter_nonneighbors(v);
        }
        return result;
    }

    std::vector<Vertex>
    random_multistart_best_clique(std::size_t iterations,
                                  const std::vector<Vertex>& subgraph) 
    {
        p_check_prop();
        std::vector<Vertex> result;
        m_subg = subgraph;
        p_filter_invalid();
        std::vector<Vertex> filtered_whole = m_subg;
        for (std::size_t i = 0; i < iterations; ++i) {
            m_subg = filtered_whole;
            Vertex vbegin = p_select_next_vertex();
            std::vector<Vertex> current = p_compute_clique(vbegin);
            if (current.size() > result.size())
                result = current;
        }
        return result;
    }

    std::vector<Vertex>
    random_multistart_best_clique_known_valid(std::size_t iterations,
                                              const std::vector<Vertex>& subgraph)
    {
        p_check_prop();
        std::vector<Vertex> result;
        for(std::size_t i = 0; i < iterations; ++i) {
            m_subg = subgraph;
            Vertex vbegin = p_select_next_vertex();
            std::vector<Vertex> current = p_compute_clique(vbegin);
            if(current.size() > result.size())
                result = current;
        }
        return result;
    }

    // compute degree sequence (deg(vertices[i]) == result[i]);
    // much slower than random_multistart_best_clique!
    std::vector<std::size_t>
    subgraph_degrees(const std::vector<Vertex>& vertices) 
    {
        p_check_prop();
        std::vector<std::size_t> result;
        for (auto i = vertices.begin(), e = vertices.end(); i != e; ++i) {
            Vertex v_i = *i;
            if (!p_push(v_i)) {
                result.push_back(0);
            } else {
                auto count = static_cast<std::size_t>(std::count_if(
                    vertices.begin(), vertices.end(), [&](Vertex v_j) {
                        return v_i != v_j && !p_can_push(v_j);
                    }));
                result.push_back(count);
            }
            while (m_prop.get_current_level() > 0)
                m_prop.pop_level();
        }
        return result;
    }

    std::vector<double>
    subgraph_degree_estimates(const std::vector<Vertex>& vertices,
                              std::size_t num_edge_samples) 
    {
        p_check_prop();
        if (vertices.empty())
            return {};
        if (vertices.size() == 1)
            return std::initializer_list<double>{0.0};
        std::size_t total_edges = (vertices.size() * (vertices.size() - 1)) / 2;
        double sample_prob = double(num_edge_samples) / total_edges;
        std::geometric_distribution<std::size_t> skip_dist(sample_prob);
        double inv_prob = 1.0 / sample_prob;
        auto& rng = sammy::rng();
        std::vector<double> estimates(vertices.size(), 0.0);
        auto b = vertices.begin(), i = b + 1, j = b, e = vertices.end();
        p_push(*i);
        std::size_t skip = skip_dist(rng);
        for (;;) {
            if (p_advance_pair(i, j, b, e, skip)) {
                while (m_prop.get_current_level() > 0)
                    m_prop.pop_level();
                if (i == e)
                    return estimates;
                p_push(*i);
            }
            if (!p_can_push(*j)) {
                estimates[i - b] += inv_prob;
                estimates[j - b] += inv_prob;
            }
            skip = 1 + skip_dist(rng);
        }
    }

    SharedDBPropagator& propagator() noexcept { return m_prop; }

  private:
    void p_check_prop() {
        if (m_prop.get_current_level() != 0)
            throw std::logic_error("Propagator not at level 0");
        m_prop.incorporate_or_throw();
    }

    template <typename Iterator>
    bool p_advance_pair(Iterator& i, Iterator& j, Iterator begin, Iterator end,
                        std::size_t dist) {
        bool result = false;
        while (i != end && dist > 0) {
            std::size_t row_remaining(i - j);
            if (row_remaining > dist) {
                j += dist;
                break;
            } else {
                dist -= row_remaining;
                ++i;
                j = begin;
                result = true;
            }
        }
        return result;
    }

    std::vector<Vertex> p_compute_clique(Vertex start_vertex) {
        std::vector<Vertex> result;
        p_filter_nonneighbors(start_vertex);
        while (!m_subg.empty()) {
            Vertex v = p_select_next_vertex();
            p_filter_nonneighbors(v);
            result.push_back(v);
        }
        return result;
    }

    Vertex p_select_next_vertex() const {
        std::uniform_int_distribution<std::size_t> dist(0, m_subg.size() - 1);
        return m_subg[dist(sammy::rng())];
    }

    bool p_push(Vertex v) {
        if (m_prop.is_false(v.first) || m_prop.is_false(v.second))
            return false;
        bool pushed = false;
        if (m_prop.is_open(v.first)) {
            if (!m_prop.push_level(v.first) || m_prop.is_false(v.second)) {
                m_prop.pop_level();
                return false;
            }
            pushed = true;
        }
        if (m_prop.is_open(v.second)) {
            if (!m_prop.push_level(v.second)) {
                m_prop.pop_level();
                if (pushed)
                    m_prop.pop_level();
                return false;
            }
        }
        return true;
    }

    bool p_can_push(Vertex v) {
        auto old_level = m_prop.get_current_level();
        if (!p_push(v))
            return false;
        while (m_prop.get_current_level() > old_level)
            m_prop.pop_level();
        return true;
    }

    void p_filter_invalid() {
        m_subg.erase(std::remove_if(m_subg.begin(), m_subg.end(),
                                    [&](Vertex v) { return !p_can_push(v); }),
                     m_subg.end());
    }

    void p_filter_nonneighbors(Vertex vnew) {
        p_push(vnew);
        auto nend = std::remove_if(m_subg.begin(), m_subg.end(), [&](Vertex v) {
            return v == vnew || p_can_push(v);
        });
        m_subg.erase(nend, m_subg.end());
        while (m_prop.get_current_level() > 0)
            m_prop.pop_level();
    }

    SharedDBPropagator m_prop;
    std::vector<Vertex> m_subg;
};

class ParallelFastCliqueBuilder {
  public:
    explicit ParallelFastCliqueBuilder(SharedDBPropagator base_prop, ThreadGroup<void>* thread_pool) :
        m_base_prop(std::move(base_prop)),
        m_thread_pool(thread_pool)
    {}

    std::vector<Vertex> 
        random_multistart_best_clique(std::size_t iterations_per_thread,
                                      const std::vector<Vertex>& subgraph)
    {
        m_base_prop.incorporate_or_throw();
        std::mutex m_out_lock;
        std::vector<Vertex> result;
        m_thread_pool->run_n_copies(m_thread_pool->num_threads() + 1,
            [&] () {
                FastCliqueBuilder builder{m_base_prop};
                auto rclique = builder.random_multistart_best_clique(iterations_per_thread, subgraph);
                {
                    std::unique_lock<std::mutex> l{m_out_lock};
                    if(rclique.size() > result.size()) {
                        result = std::move(rclique);
                    }
                }
            }
        );
        return result;
    }

  private:
    SharedDBPropagator m_base_prop;
    ThreadGroup<void>* m_thread_pool;
};

} // namespace sammy

#endif
