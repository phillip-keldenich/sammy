#ifndef SAMMY_COLORING_H_INCLUDED_
#define SAMMY_COLORING_H_INCLUDED_

#include "universe_subgraph.h"
#include "kissat_solver.h"
#include "output.h"
#include "experiment_flags.h"
#include <boost/iterator/counting_iterator.hpp>

namespace sammy {

inline std::pair<
        std::vector<SharedDBPropagator>,
        std::vector<Vertex> 
> coloring_to_classes(
    const SharedDBPropagator& prop_template,
    const std::vector<Vertex>& vertices,
    const std::vector<std::size_t>& coloring,
    std::size_t num_colors) 
{
    std::cout << "To classes: " << num_colors << std::endl;
    std::vector<std::vector<std::size_t>> color_classes(num_colors);
    std::vector<SharedDBPropagator> classes(num_colors, prop_template);
    std::vector<Vertex> spawners;
    for(std::size_t i = 0, n = vertices.size(); i < n; ++i) {
        auto c = coloring[i];
        color_classes[c].push_back(i);
    }
    for(std::size_t i = 0, n = color_classes.size(); i < n; ++i) {
        std::cout << "Class " << i << ": " << color_classes[i].size() << std::endl;
    }
    for(std::size_t i = 0; i < num_colors; ++i) {
        auto& cc = color_classes[i];
        auto& ci = classes[i];
        spawners.push_back(vertices[cc.front()]);
        cc.erase(std::remove_if(cc.begin(), cc.end(),
                 [&] (std::size_t vi) { return push_vertex(ci, vertices[vi]) >= 0; }),
                 cc.end());
    }
    for(std::size_t i = 0; i < num_colors; ++i) {
        auto& cc = color_classes[i];
        for(std::size_t vi : cc) {
            bool found = false;
            for(std::size_t k = 0, cc = classes.size(); k < cc; ++k) {
                if(push_vertex(classes[k], vertices[vi]) >= 0) {
                    found = true;
                    break;
                }
            }
            if(!found) {
                classes.push_back(prop_template);
                push_vertex(classes.back(), vertices[vi]);
                spawners.push_back(vertices[vi]);
            }
        }
    }
    return {std::move(classes), std::move(spawners)};
}

struct ColoringExtraConstraints {
    using ExtraConstraint = boost::container::small_vector<std::size_t, 4>;
    using ExtraConstraints = std::vector<ExtraConstraint>;

    ExtraConstraints constraints;
    std::vector<std::vector<std::size_t>> extra_constraints_by_vertex;

    void update_vertex_count(const UniverseSubgraph& subgraph) {
        extra_constraints_by_vertex.resize(subgraph.n());
    }

    void add_constraint(const ExtraConstraint& extra_constraint) {
        constraints.push_back(extra_constraint);
        p_update_extra_constraints_by_vertex();
    }

    void add_constraint(ExtraConstraint&& extra_constraint) {
        constraints.emplace_back(std::move(extra_constraint));
        p_update_extra_constraints_by_vertex();
    }

    template<typename Iterator>
    void add_constraint(Iterator begin, Iterator end)
    {
        constraints.emplace_back(begin, end);
        p_update_extra_constraints_by_vertex();
    }

    std::vector<SharedDBPropagator> coloring_to_classes(
        SharedDBPropagator& template_propagator,
        const UniverseSubgraph& subgraph,
        const std::vector<std::size_t>& coloring,
        std::size_t num_colors) 
    {
        std::vector<std::vector<std::size_t>> color_classes(num_colors);
        for(std::size_t vi = 0; vi < coloring.size(); ++vi) {
            color_classes[coloring[vi]].push_back(vi);
        }
        std::vector<SharedDBPropagator> result;
        for(std::size_t cc = 0; cc < num_colors; ++cc) {
            auto& cclass = color_classes[cc];
            if(cclass.empty()) continue;
            template_propagator.reset_or_throw();
            for(auto vi_it = cclass.begin(), vi_end = cclass.end(); vi_it != vi_end; ++vi_it) {
                Vertex v = subgraph.vertex(*vi_it);
                if(push_vertex(template_propagator, v) < 0) {
                    p_constraint_from_prop(template_propagator, cclass.begin(), vi_it + 1, subgraph);
                    p_coloring_to_classes_failed(template_propagator, subgraph, color_classes, cc + 1);
                    result.clear();
                    return result;
                }
            }
            result.push_back(template_propagator);
        }
        return result;
    }

  private:
    void p_coloring_to_classes_failed(
        SharedDBPropagator& template_propagator,
        const UniverseSubgraph& subgraph,
        const std::vector<std::vector<std::size_t>>& color_classes,
        std::size_t start_with_class) 
    {
        for(std::size_t i = start_with_class; i < color_classes.size(); ++i) {
            auto& cclass = color_classes[i];
            if(cclass.size() <= 2) continue;
            template_propagator.reset_or_throw();
            for(auto vi_it = cclass.begin(), vi_end = cclass.end(); vi_it != vi_end; ++vi_it) {
                Vertex v = subgraph.vertex(*vi_it);
                if(push_vertex(template_propagator, v) < 0) {
                    p_constraint_from_prop(template_propagator, cclass.begin(),
                                           vi_it + 1, subgraph);
                    break;
                }
            }
        }
    }

    template<typename RandomAccessIterator>
    void p_constraint_from_prop(SharedDBPropagator& propagator,
                                RandomAccessIterator pushed_begin,
                                RandomAccessIterator pushed_end,
                                const UniverseSubgraph& subgraph)
    {
        // pushing this vertex failed;
        // the propagator has all previous vertices pushed
        std::size_t failed_index = *(pushed_end-1);
        Vertex failed = subgraph.vertex(failed_index);
        if(propagator.is_false(failed.first)) {
            // case 1a: the first literal of the vertex is false
            p_learn_constraint(propagator, pushed_begin, pushed_end-1, lit::negate(failed.first), failed_index, subgraph);
            return;
        }
        if(propagator.is_false(failed.second)) {
            // case 1b: the second literal of the vertex is false
            p_learn_constraint(propagator, pushed_begin, pushed_end-1, lit::negate(failed.second), failed_index, subgraph);
            return;
        }
        if(propagator.is_open(failed.first)) {
            if(!propagator.push_level(failed.first)) {
                // case 2a: pushing the first literal from the vertex causes a conflict
                p_learn_constraint_from_conflict(propagator, pushed_begin, pushed_end, subgraph);
                return;
            }
            if(propagator.is_false(failed.second)) {
                // case 3: pushing the first literal from the vertex causes no conflict,
                // but implies the negation of the second literal
                // TODO: this case needs more thought...
                add_constraint(pushed_begin, pushed_end);
                return;
            }
        }
        // case 2b: pushing the second literal from the vertex causes a conflict
        if(propagator.is_open(failed.second) && !propagator.push_level(failed.second)) {
            p_learn_constraint_from_conflict(propagator, pushed_begin, pushed_end, subgraph);
            return;
        }
        throw std::logic_error("Invalid case in p_constraint_from_prop!");
    }

    template<typename RandomAccessIterator>
    std::vector<DynamicBitset> p_compute_vertex_consequences(
        const SharedDBPropagator& propagator_template,
        RandomAccessIterator vbegin, RandomAccessIterator vend,
        const UniverseSubgraph& subgraph) 
    {
        SharedDBPropagator propagator{propagator_template};
        std::size_t n = propagator.db().num_vars();
        std::vector<DynamicBitset> vertex_consequences;
        for(std::size_t vindex : IteratorRange{vbegin, vend}) {
            propagator.reset_or_throw();
            Vertex v = subgraph.vertex(vindex);
            if(propagator.is_open(v.first)) propagator.push_level(v.first);
            if(propagator.is_open(v.second)) propagator.push_level(v.second);
            vertex_consequences.emplace_back(2 * n, false);
            auto& b = vertex_consequences.back();
            for(Lit l : propagator.get_trail()) {
                b[l] = true;
            }
        }
        return vertex_consequences;
    }

    template<typename RandomAccessIterator>
    void p_learn_constraint_from_conflict(
        SharedDBPropagator& propagator, 
        RandomAccessIterator pushed_begin, RandomAccessIterator pushed_end, 
        const UniverseSubgraph& subgraph)
    {
        auto vertex_consequences = p_compute_vertex_consequences(propagator, pushed_begin, pushed_end, subgraph);
        auto [conflict_literal, conflict_reason] = propagator.get_conflict();
        DynamicBitset active(2 * propagator.db().num_vars(), false);
        Lit coverage_block = lit::negate(conflict_literal);
        active[coverage_block] = true;
        for(Lit l : conflict_reason.lits(propagator.db())) {
            if(l != conflict_literal) {
                auto lneg = lit::negate(l);
                active[lneg] = true;
            }
        }
        auto best_cover = p_walk_trail(propagator, propagator.get_trail().size() - 1, coverage_block,
                                       active, vertex_consequences, pushed_begin, pushed_end);
        add_constraint(best_cover);
    }

    template<typename RandomAccessIterator>
    ExtraConstraint
         p_walk_trail(SharedDBPropagator& propagator,
                      std::size_t start_index, Lit initial_blocker,
                      DynamicBitset& active,
                      const std::vector<DynamicBitset>& vertex_consequences,
                      RandomAccessIterator pushed_begin,
                      RandomAccessIterator pushed_end)
    {
        ExtraConstraint best_cover;
        const auto& trail = propagator.get_trail();
        const auto& reasons = propagator.get_reasons();
        Lit coverage_blocker = initial_blocker;
        DynamicBitset tmp;
        DynamicBitset coverable = vertex_consequences.front();
        std::for_each(vertex_consequences.begin() + 1, vertex_consequences.end(),
                      [&] (const DynamicBitset& b) {coverable |= b;});
        for(std::size_t tindex = start_index; tindex <= start_index; --tindex) {
            Lit tlit = trail[tindex];
            if(!active[tlit]) continue;
            Reason r = reasons[tindex];
            switch(r.reason_length) {
                case 0: break;
                case 1: active[tlit] = false; break;
                default:
                    for(Lit l : r.lits(propagator.db())) {
                        if(l != tlit) {
                            active[lit::negate(l)] = true;
                        }
                    }
                    active[tlit] = false;
                break;
            }
            if(!active[coverage_blocker]) {
                tmp = active;
                tmp -= coverable;
                auto blocker = tmp.ones_begin();
                if(blocker == tmp.ones_end()) {
                    p_compute_cover(active, tmp, vertex_consequences, pushed_begin, pushed_end, best_cover);
                    if(best_cover.size() <= 2) break;
                } else {
                    coverage_blocker = Lit(*blocker);
                }
            }
        }
        if(best_cover.empty()) 
            throw std::logic_error("Could not cover conflict literals with vertices!");
        return best_cover;
    }

    template<typename RandomAccessIterator>
    void p_learn_constraint(SharedDBPropagator& propagator,
                            RandomAccessIterator pushed_begin, RandomAccessIterator pushed_end,
                            Lit problem_literal, std::size_t extra_vertex, const UniverseSubgraph& subgraph)
    {
        auto vertex_consequences = p_compute_vertex_consequences(propagator, pushed_begin, pushed_end, subgraph);
        DynamicBitset active(2 * propagator.db().num_vars(), false);
        active[problem_literal] = true;
        std::size_t problem_index = propagator.get_trail_index(problem_literal);
        auto best_cover = p_walk_trail(propagator, problem_index, problem_literal, active, vertex_consequences,
                                       pushed_begin, pushed_end);
        best_cover.push_back(extra_vertex);
        add_constraint(best_cover);
    }

    /**
     * Cover the current conflict literals
     * with vertices in a greedy fashion
     * (their decisions & consequences).
     */
    template<typename RandomAccessIterator>
    void p_compute_cover(DynamicBitset& active, DynamicBitset& tmp, 
                         const std::vector<DynamicBitset>& vconsq,
                         RandomAccessIterator begin, RandomAccessIterator end,
                         ExtraConstraint& best_cover)
    {
        DynamicBitset still_active = active;
        ExtraConstraint cover;
        while(still_active.any()) {
            std::size_t best_num_covered = 0;
            std::size_t best_vertex = 0;
            std::size_t pindex = 0, best_pindex = 0;
            for(RandomAccessIterator p = begin; p != end; ++p, ++pindex) {
                tmp = still_active;
                tmp &= vconsq[pindex];
                auto count = tmp.count();
                if(count > best_num_covered) {
                    best_num_covered = count;
                    best_vertex = *p;
                    best_pindex = pindex;
                }
            }
            cover.push_back(best_vertex);
            still_active -= vconsq[best_pindex];
        }
        if(best_cover.empty() || cover.size() < best_cover.size()) {
            best_cover = std::move(cover);
        }
    }

    void p_update_extra_constraints_by_vertex() {
        ExtraConstraint& last = constraints.back();
        std::cout << "CONSTRAINT ADDED:";
        for(std::size_t vi : last) {
            std::cout << " " << vi;
        }
        std::cout << std::endl;
        std::size_t idx = constraints.size() - 1;
        for(std::size_t v : last) {
            extra_constraints_by_vertex[v].push_back(idx);
        }
    }
};

class SATKColoringSolver {
  public:
    SATKColoringSolver(UniverseSubgraph* subgraph,
                       EventRecorder* recorder,
                       std::vector<std::size_t> clique,
                       std::size_t k) :
        m_subgraph(subgraph),
        m_recorder(recorder),
        m_clique(std::move(clique)),
        m_k(k)
    {}

    void set_extra_constraints(ColoringExtraConstraints* extra) {
        m_extra_cons = extra;
    }

    std::optional<bool> solve(double time_limit = 
                              std::numeric_limits<double>::infinity()) 
    {
        m_recorder->store_event("BEGIN_SOLVE_K_COLORING", 
                                {{"k", m_k}, {"n", m_subgraph->n()}}, 
                                "k", "n");
        p_reduce_graph();
        std::vector<std::size_t> core_coloring;
        if(!m_new_to_old.empty()) {
            std::vector<bool> solution;
            auto optres = p_sat_solve(solution, time_limit);
            if(!optres) {
                m_recorder->store_event("DONE_SOLVE_K_COLORING", {{"result", "aborted"}}, "result");
                return std::nullopt;
            }
            m_recorder->store_event("DONE_SOLVE_K_COLORING",
                                    {{"result", *optres ? "coloring" : "bound"}}, 
                                    "result");
            if(*optres) {
                core_coloring = p_core_coloring_from_solution(solution);
            } else {
                return false;
            }
        }
        p_compact_core_coloring(core_coloring);
        p_extend_core_coloring(core_coloring);
        return true;
    }

    const std::vector<std::size_t>& get_coloring() const noexcept {
        return m_coloring;
    }

  private:
    std::vector<std::size_t> p_core_coloring_from_solution(const std::vector<bool>& sol) {
        std::size_t nnew = m_new_to_old.size();
        std::vector<std::size_t> result(nnew, std::numeric_limits<std::size_t>::max());
        std::size_t offs = 0;
        for(std::size_t vi = 0; vi < nnew; ++vi) {
            for(std::size_t c = 0; c < m_k; ++c, ++offs) {
                if(sol[offs]) {
                    result[vi] = c;
                    offs += m_k - c;
                    break;
                }
            }
        }
        return result;
    }

    auto p_var_range(std::size_t new_vertex_index) {
        using Lit = KissatSolver::Lit;
        return IteratorRange{
            boost::make_counting_iterator(Lit(new_vertex_index * m_k + 1)),
            boost::make_counting_iterator(Lit(new_vertex_index * m_k + m_k + 1))
        };
    }

    void p_add_clique(KissatSolver& solver, std::size_t nnew) {
        std::size_t clique_next_color = 0;
        for(std::size_t cvi : m_clique) {
            std::size_t new_i = m_old_to_new[cvi];
            if(new_i < nnew) {
                solver.add_short_clause(p_var_range(new_i).begin()[clique_next_color++]);
            }
        }
    }

    void p_add_edges(KissatSolver& solver, std::size_t nold, std::size_t nnew)
    {
        using Lit = KissatSolver::Lit;
        for(std::size_t i = 0; i < nold; ++i) {
            std::size_t new_i = m_old_to_new[i];
            if(new_i >= nnew) continue;
            auto vri = p_var_range(new_i);
            for(std::size_t j : m_subgraph->matrix_row(i).ones_from(i + 1)) {
                std::size_t new_j = m_old_to_new[j];
                if(new_j >= nnew) continue;
                auto vrj = p_var_range(new_j);
                auto vi = vri.begin();
                for(Lit lj : vrj) {
                    solver.add_short_clause(-lj, -*vi);
                    ++vi;
                }
            }
        }
    }

    void p_add_extra(KissatSolver& solver) 
    {
        ColoringExtraConstraints::ExtraConstraint extra;
        for(std::size_t i = 0, ncons = m_extra_cons->constraints.size(); i < ncons; ++i) {
            if(!m_extra_cons_active[i]) continue;
            extra.clear();
            const auto& cons = m_extra_cons->constraints[i];
            std::transform(cons.begin(), cons.end(), std::back_inserter(extra),
                           [&] (std::size_t vold) { return m_old_to_new[vold]; });
            for(std::size_t c = 0; c < m_k; ++c) {
                std::for_each(extra.begin(), extra.end(), [&] (std::size_t vnew) {
                    solver.add_literal(-p_var_range(vnew).begin()[c]);
                });
                solver.finish_clause();
            }
        }
    }

    std::optional<bool> p_sat_solve(std::vector<bool>& solution,
                                    double time_limit) 
    {
        std::size_t nnew = m_new_to_old.size();
        std::size_t nold = m_subgraph->n();
        KissatSolver solver;
        solver.reserve(m_k * nnew);
        for(std::size_t i = 0; i < nnew; ++i) {
            auto r = p_var_range(i);
            solver.add_clause(r.begin(), r.end());
            for(auto c1i = r.begin(), e = r.end(); c1i != e; ++c1i) {
                for(auto c2i = c1i + 1; c2i != e; ++c2i) {
                    solver.add_short_clause(-*c1i, -*c2i);
                }
            }
        }
        p_add_clique(solver, nnew);
        p_add_edges(solver, nold, nnew);
        if(m_extra_cons) p_add_extra(solver);
        auto r = solver.solve(time_limit);
        if(r && *r) solution = std::move(solver.get_model().raw());
        return r;
    }

    void p_compact_core_coloring(std::vector<std::size_t>& core_coloring) {
        DynamicBitset unused(m_k, true);
        for(std::size_t c : core_coloring) {
            unused[c].reset();
        }
        std::vector<std::size_t> uco;
        std::copy(unused.ones_begin(), unused.ones_end(), std::back_inserter(uco));
        if(uco.empty()) return;
        std::vector<std::size_t> remap;
        auto uco_iter = uco.begin();
        std::size_t next_color = 0;
        for(std::size_t c = 0; c < m_k; ++c) {
            if(c == *uco_iter) {
                remap.push_back(0);
                if(++uco_iter == uco.end()) {
                    for(++c; c < m_k; ++c) {
                        remap.push_back(next_color++);
                    }
                    break;
                }
            } else {
                remap.push_back(next_color++);
            }
        }
        for(std::size_t& c : core_coloring) {
            c = remap[c];
        }
    }

    void p_extend_core_coloring(const std::vector<std::size_t>& core_coloring) 
    {
        m_coloring.assign(m_subgraph->n(), 0);
        for(std::size_t i = 0, nnew = core_coloring.size(); i < nnew; ++i) {
            std::size_t old_index = m_new_to_old[i];
            m_coloring[old_index] = core_coloring[i];
        }
        DynamicBitset avail(m_k, true);
        for(std::size_t deleted : IteratorRange{m_deletion_order.rbegin(), 
                                                m_deletion_order.rend()}) 
        {
            avail.assign(m_k, true);
            for(std::size_t neigh : m_subgraph->matrix_row(deleted).ones()) {
                if(m_degrees[neigh] >= m_k) {
                    avail[m_coloring[neigh]].reset();
                }
            }
            if(m_extra_cons) {
                for(std::size_t neigh : m_extra_pseudo_edges[deleted]) {
                    if(m_degrees[neigh] >= m_k) {
                        avail[m_coloring[neigh]].reset();
                    }
                }
            }
            assert(avail.count() != 0);
            m_coloring[deleted] = *avail.ones_begin();
            m_degrees[deleted] = m_k; // make sure the new vertex now counts as present
        }
    }

    bool p_degree_too_low(std::size_t vindex) const noexcept {
        return m_degrees[vindex] < m_k;
    }

    void p_reduce_graph() {
        p_compute_degrees();
        p_init_deletion_order();
        if(m_deletion_order.empty()) {
            m_new_to_old.assign(boost::make_counting_iterator(std::size_t(0)),
                                boost::make_counting_iterator(m_subgraph->n()));
            m_old_to_new = m_new_to_old;
            m_extra_cons_active.assign(m_extra_cons->constraints.size(), true);
            m_extra_pseudo_edges.assign(m_subgraph->n(), boost::container::small_vector<std::size_t,4>{});
        } else {
            p_bfs_deletion_order();
            p_extract_new_to_old();
        }
    }

    void p_init_deletion_order() {
        m_deletion_order.clear();
        std::copy_if(
            boost::make_counting_iterator(std::size_t(0)),
            boost::make_counting_iterator(m_subgraph->n()),
            std::back_inserter(m_deletion_order),
            [&] (std::size_t vi) { return p_degree_too_low(vi); });
    }

    void p_bfs_deletion_order() {
        if(m_extra_cons) {
            m_extra_cons_active.assign(m_extra_cons->constraints.size(), true);
            m_extra_pseudo_edges.assign(m_subgraph->n(), boost::container::small_vector<std::size_t,4>{});
        }
        std::size_t queue_pos = 0;
        while(queue_pos < m_deletion_order.size()) {
            std::size_t deleted = m_deletion_order[queue_pos++];
            for(std::size_t neigh : m_subgraph->matrix_row(deleted).ones()) {
                if(--m_degrees[neigh] == m_k - 1) {
                    m_deletion_order.push_back(neigh);
                }
            }
            if(m_extra_cons) {
                for(std::size_t cons_idx : m_extra_cons->extra_constraints_by_vertex[deleted]) {
                    if(m_extra_cons_active[cons_idx]) {
                        p_deactivate_cons(deleted, cons_idx);
                    }
                }
            }
        }
    }

    void p_deactivate_cons(std::size_t deleted_vertex, std::size_t deactivated_cons) {
        m_extra_cons_active[deactivated_cons] = false;
        const auto& cons = m_extra_cons->constraints[deactivated_cons];
        std::size_t other_vertex = cons[0];
        if(other_vertex == deleted_vertex) {
            other_vertex = cons[1];
        }
        m_extra_pseudo_edges[deleted_vertex].push_back(other_vertex);
        m_extra_pseudo_edges[other_vertex].push_back(deleted_vertex);
        for(const auto& vertex : cons) {
            if(--m_degrees[vertex] == m_k - 1) {
                m_deletion_order.push_back(vertex);
            }
        }
    }

    void p_extract_new_to_old() {
        m_new_to_old.clear();
        m_old_to_new.clear();
        std::size_t num_used = 0;
        for(std::size_t i = 0, n = m_subgraph->n(); i < n; ++i) {
            if(m_degrees[i] >= m_k) {
                m_new_to_old.push_back(i);
                m_old_to_new.push_back(num_used++);
            } else {
                m_old_to_new.push_back(std::numeric_limits<std::size_t>::max());
            }
        }
    }

    void p_compute_degrees() {
        m_degrees = m_subgraph->get_degrees();
        if(m_extra_cons) {
            for(std::size_t i = 0, n = m_subgraph->n(); i < n; ++i) {
                m_degrees[i] += m_extra_cons->extra_constraints_by_vertex[i].size();
            }
        }
    }

    UniverseSubgraph* m_subgraph;
    EventRecorder* m_recorder;
    std::vector<std::size_t> m_clique;
    std::vector<std::size_t> m_deletion_order;
    std::vector<std::size_t> m_new_to_old;
    std::vector<std::size_t> m_old_to_new;
    std::vector<std::size_t> m_degrees;
    std::vector<std::size_t> m_coloring;
    std::size_t m_k;
    ColoringExtraConstraints *m_extra_cons = nullptr;
    std::vector<bool> m_extra_cons_active;
    std::vector<boost::container::small_vector<std::size_t,4>> m_extra_pseudo_edges;
};


class SATColoringSolver {
  public:
    SATColoringSolver(UniverseSubgraph* subgraph,
                      EventRecorder* recorder,
                      std::vector<Vertex> clique,
                      std::vector<std::size_t> best_known_coloring) :
        m_subgraph(subgraph),
        m_recorder(recorder),
        best_known_lower_bound(clique.size()),
        best_clique(std::move(clique)),
        best_coloring(std::move(best_known_coloring)),
        best_num_colors(*std::max_element(best_coloring.begin(), best_coloring.end()) + 1)
    {
        clique_indices.reserve(best_clique.size());
        std::transform(best_clique.begin(), best_clique.end(),
                       std::back_inserter(clique_indices),
                       [&] (Vertex v) { return m_subgraph->vertex_index(v); });
        best_num_colors = p_compact_coloring(best_coloring, best_num_colors);
    }

    void added_vertices() {
        const auto old_n = best_coloring.size();
        const auto new_n = m_subgraph->n();
        if(old_n >= new_n) return;
        for(std::size_t i = old_n; i < new_n; ++i) {
            
        }
    }

    std::optional<bool> run_with_num_colors(std::size_t num_colors, double time_limit =
                                            std::numeric_limits<double>::infinity()) 
    {
        if(num_colors >= best_num_colors) return true;
        if(num_colors < best_known_lower_bound) return false;
        const std::size_t n = m_subgraph->n();
        auto x_vc = [&] (std::size_t vindex, std::size_t color) -> Lit {
            return Lit(color * n + vindex + 1);
        };
        KissatSolver solver;
        solver.reserve(num_colors * n);
        for(std::size_t i = 0; i < clique_indices.size(); ++i) {
            solver.add_short_clause(x_vc(clique_indices[i], i));
        }
        std::vector<Lit> clause_buffer(num_colors, 0);
        for(std::size_t i = 0; i < n; ++i) {
            for(std::size_t j : m_subgraph->matrix_row(i).ones_from(i + 1)) {
                for(std::size_t k = 0; k < num_colors; ++k) {
                    solver.add_short_clause(-x_vc(j, k), -x_vc(i, k));
                }
            }
            for(std::size_t k = 0; k < num_colors; ++k) {
                clause_buffer[k] = x_vc(i, k);
            }
            solver.add_clause(clause_buffer.begin(), clause_buffer.end());
        }
        auto result = solver.solve(time_limit);
        if(!result) return std::nullopt;
        if(!*result) {
            m_recorder->store_event("NEW_LOWER_BOUND", {{"prev", best_known_lower_bound}, 
                                                        {"new", num_colors + 1}}, 
                                    "prev", "new");
            best_known_lower_bound = num_colors + 1;
            return false;
        }
        auto model = solver.get_model();
        const auto& solution = model.raw();
        m_recorder->store_event("NEW_COLORING", {{"prev", best_num_colors}, {"new", num_colors}}, 
                                "prev", "new");
        for(std::size_t k = 0, offs = 0; k < num_colors; ++k) {
            for(std::size_t v = 0; v < n; ++v, ++offs) {
                if(solution[offs]) {
                    best_coloring[v] = k;
                }
            }
        }
        best_num_colors = p_compact_coloring(best_coloring, num_colors);
        return true;
    }

    bool optimize_coloring(double time_limit = std::numeric_limits<double>::infinity()) {
        if(time_limit <= 0) return false;
        auto begin = Clock::now();
        while(best_known_lower_bound < best_num_colors) {
            std::size_t mid = (best_known_lower_bound + best_num_colors) / 2;
            double trem = time_limit - seconds_between(begin, Clock::now());
            m_recorder->store_event("BEGIN_SAT_COLORING", {{"best_lb", best_known_lower_bound},
                                                           {"best_coloring", best_num_colors},
                                                           {"query", mid},
                                                           {"n", m_subgraph->n()}},
                                    "best_lb", "best_coloring", "query", "n");
            auto result = run_with_num_colors(mid, trem);
            if(!result) {
                m_recorder->store_event("SAT_COLORING_TIMEOUT");
                return false;
            }
            m_recorder->store_event("DONE_SAT_COLORING",
                                    {{"best_lb", best_known_lower_bound},
                                     {"best_coloring", best_num_colors}}, 
                                    "best_lb", "best_coloring");
        }
        m_recorder->store_event("OPTIMUM_SUBGRAPH_COLORING", {{"best_coloring", best_num_colors}}, "best_coloring");
        return true;
    }

    std::pair<
        std::vector<SharedDBPropagator>,
        std::vector<Vertex> > coloring_to_classes() 
    {
        auto& prop = m_subgraph->get_propagator();
        prop.reset_or_throw();
        return sammy::coloring_to_classes(prop, m_subgraph->vertex_set(), best_coloring, best_num_colors);
    }

  private:
    std::size_t p_compact_coloring(std::vector<std::size_t>& coloring, std::size_t num_colors) {
        const auto n = coloring.size();
        std::vector<std::vector<std::size_t>> color_classes(num_colors);
        for(std::size_t i = 0; i < n; ++i) {
            color_classes[coloring[i]].push_back(i);
        }
        std::size_t num_used = 0;
        std::vector<DynamicBitset> vertex_compatibility(num_colors, DynamicBitset(n, true));
        for(std::size_t k = num_colors - 1; k < num_colors; --k) {
            const auto& cclass = color_classes[k];
            for(std::size_t v : cclass) {
                for(std::size_t cand = 0; cand < num_colors; ++cand) {
                    if(vertex_compatibility[cand][v]) {
                        if(cand >= num_used) num_used = cand + 1;
                        vertex_compatibility[cand] -= m_subgraph->matrix_row(v);
                        coloring[v] = cand;
                        break;
                    }
                }
            }
        }
        return num_used;
    }

    UniverseSubgraph* m_subgraph;
    EventRecorder* m_recorder;
    std::size_t best_known_lower_bound;
    std::vector<Vertex> best_clique;
    std::vector<std::size_t> clique_indices;
    std::vector<std::size_t> best_coloring;
    std::size_t best_num_colors;
};

}

#endif
