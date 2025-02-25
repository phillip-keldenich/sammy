#ifndef SAMMY_SATDSATUR_COLORER_H_INCLUDED_
#define SAMMY_SATDSATUR_COLORER_H_INCLUDED_

#include "shared_db_propagator.h"
#include "universe_subgraph.h"
#include "algorithm_ex.h"
#include "best_k.h"
#include "vertex_operations.h"
#include "class_completer.h"
#include "partial_solution.h"
#include <variant>
#include <atomic>


namespace sammy {

template<typename IncrementalSolver>
class SATDSaturCoverer {
    // Options for keeping open_classes updated?
    // - use only the graph datastructure
    //   * iterate neighbors on addition
    //   * no deletion/conflict learning
    //   * un 'unexpected' conflict, remove class
    //   * resets of classes on SAT coverage
    // - use graph datastructure and propagated literals
    //   * iterate neighbors on addition
    //   * have/iterate list of vertices with literal
    //   * no deletion/conflict learning
    //   * on unexpected conflict, remove class
    //   * the only way to remove vertices from a class
    //     is a successful SAT return (recoloring)
    //   * how to track covering classes?
    struct VertexInfo {
        DynamicBitset open_classes; //< a superset of the open classes for this vertex
        std::size_t num_open_classes; //< number of set bits in open_classes
        bool in_some_class; //< number of classes containing this vertex
        std::size_t degree; //< degree of the vertex in the graph

        /**
         * Mark the given class as unavailable for this vertex.
         */
        bool remove_from(std::size_t class_index) noexcept {
            if(open_classes[class_index]) {
                open_classes[class_index].reset();
                --num_open_classes;
                return true;
            }
            return false;
        }
    };

    struct CompareByInfo {
        bool operator()(std::size_t v1, std::size_t v2) const noexcept {
            const VertexInfo& i1 = that->m_vertex_info[v1];
            const VertexInfo& i2 = that->m_vertex_info[v2];
            if(i2.in_some_class) return !i1.in_some_class;
            if(i1.in_some_class) return false;
            if(i1.num_open_classes < i2.num_open_classes) return true;
            if(i2.num_open_classes < i1.num_open_classes) return false;
            return i2.degree < i1.degree;
        }

        const SATDSaturCoverer *that;
    };

    struct ColorClass {
        ColorClass(SATDSaturCoverer* that, const SharedDBPropagator& propagator, std::size_t index) :
            propagator(propagator),
            index(index)
        {
            p_init_info_from_trail(that);
        }

        ColorClass(SATDSaturCoverer* that,
                   const SharedDBPropagator& propagator,
                   std::size_t index,
                   std::size_t initial_vertex) :
            propagator(propagator),
            index(index)
        {
            push_vertex(this->propagator, that->m_subgraph->vertex(initial_vertex));
            p_init_info_from_trail(that);
        }

        void add_vertex(SATDSaturCoverer* that, std::size_t vertex) {
            UniverseSubgraph* g = that->m_subgraph;
            std::size_t tlength = propagator.get_trail().size();
            if(push_vertex(propagator, g->vertex(vertex)) < 0) {
                throw std::logic_error("Called add_vertex on incompatible class!");
            }
            for(std::size_t vo : g->matrix_row(vertex).ones()) {
                that->m_vertex_info[vo].remove_from(index);
            }
            for(Lit l : IteratorRange{propagator.get_trail().begin()+tlength,
                                      propagator.get_trail().end()}) 
            {
                Lit lneg = lit::negate(l);
                for(std::size_t v_pot : that->m_vertices_with_literal[l]) {
                    Vertex v = g->vertex(v_pot);
                    if(propagator.is_true(v.first) && propagator.is_true(v.second)) {
                        that->m_vertex_info[v_pot].in_some_class = true;
                    }
                }
                for(std::size_t v_out : that->m_vertices_with_literal[lneg]) {
                    that->m_vertex_info[v_out].remove_from(index);
                }
            }
        }

        bool can_add(SATDSaturCoverer* that, std::size_t vertex) {
            Vertex v = that->m_subgraph->vertex(vertex);
            return can_push(propagator, v);
        }

        void reset(SATDSaturCoverer* that, const DynamicBitset& expl_vertex_is_in) {
            propagator.reset_or_throw();
            p_init_info_from_trail(that);
            for(std::size_t vi_i : expl_vertex_is_in.ones()) {
                std::size_t vi = that->m_explicit_cover_order[vi_i];
                add_vertex(that, vi);
            }
        }

        SharedDBPropagator propagator;
        std::size_t index;

      private:
        void p_init_info_from_trail(SATDSaturCoverer* that) {
            UniverseSubgraph* g = that->m_subgraph;
            for(Lit l : propagator.get_trail()) {
                for(std::size_t v_pot : that->m_vertices_with_literal[l]) {
                    Vertex v = g->vertex(v_pot);
                    if(v.first == l) {
                        if(propagator.is_true(v.second)) {
                            that->m_vertex_info[v_pot].in_some_class = true;
                        }
                    }
                }
                Lit lneg = lit::negate(l);
                for(std::size_t v_out : that->m_vertices_with_literal[lneg]) {
                    that->m_vertex_info[v_out].remove_from(index);
                }
            }
        }
    };

    using SatLit = typename IncrementalSolver::Lit;
    using LitOrFixed = std::variant<SatLit, bool>;

    UniverseSubgraph* m_subgraph;
    std::vector<std::vector<std::size_t>> m_vertices_with_literal;
    std::vector<VertexInfo> m_vertex_info;
    std::vector<ColorClass> m_color_classes;
    BestK<std::size_t, CompareByInfo> m_best_infos;
    std::size_t m_clique_size, m_lower_bound, m_upper_bound;
    std::vector<std::size_t> m_explicit_cover_order;
    std::vector<std::size_t> m_true_zero_vertices;
    IncrementalSolver m_incremental_solver;
    // model variables: variable copies for each class;
    // entry [c][v] is the v-th variable of class c
    std::vector<std::vector<LitOrFixed>> m_class_literals;
    // model variables: vertex-in-class variables
    // entry [c][vi] encodes if vertex m_explicit_cover_order[vi] is in class c
    std::vector<std::vector<LitOrFixed>> m_vertex_in_class;
    // model variables: variables indicating that
    // more than some number of colors are needed;
    // entry [c] says that more than c + m_clique_size
    // colors are necessary
    std::vector<LitOrFixed> m_not_enough_colors;
    // buffer for building clauses
    std::vector<SatLit> m_clause_buffer;
    // flag to abort the search
    std::atomic<bool> m_aborted{false};
    // indicates whether we found an optimal coloring
    // or proved no improvement is possible, rather
    // than having been aborted
    bool m_success{false};

    void p_fill_vertices_with_literal() {
        m_vertices_with_literal.resize(2 * m_subgraph->get_propagator().db().num_vars());
        std::size_t index = 0;
        for(Vertex v : m_subgraph->vertex_set()) {
            m_vertices_with_literal[v.first].push_back(index);
            m_vertices_with_literal[v.second].push_back(index);
            ++index;
        }
    }

    void p_init_vertex_info() {
        const auto n = m_subgraph->n();
        VertexInfo info{DynamicBitset(m_lower_bound, true), 
                        m_lower_bound, false, 0};
        m_vertex_info.resize(n, info);
        for(std::size_t i = 0, n = m_subgraph->n(); i < n; ++i) {
            m_vertex_info[i].degree = m_subgraph->get_degree(i);
        }
    }

    void p_init_classes(std::vector<std::size_t> clique) {
        m_color_classes.reserve(clique.size());
        std::size_t index = 0;
        UniverseSubgraph* g = m_subgraph;
        SharedDBPropagator &propagator = g->get_propagator();
        propagator.reset_or_throw();
        for(std::size_t vi : clique) {
            m_color_classes.emplace_back(this, propagator, index++, vi);
        }
        m_explicit_cover_order = std::move(clique);
    }

    std::size_t p_find_true_next_vertex(const std::vector<std::size_t>& candidates) {
        m_true_zero_vertices.clear();
        for(std::size_t c : candidates) {
            VertexInfo& info = m_vertex_info[c];
            if(info.in_some_class) continue;
            std::size_t count = 0;
            for(std::size_t cindex : info.open_classes.ones()) {
                if(m_color_classes[cindex].can_add(this, c)) {
                    ++count;
                } else {
                    info.open_classes[cindex].reset();
                    --info.num_open_classes;
                }
            }
            if(count == 0) {
                m_true_zero_vertices.push_back(c);
            }
        }
        if(!m_true_zero_vertices.empty()) {
            return m_true_zero_vertices.front();
        }
        return *std::min_element(candidates.begin(), candidates.end(), CompareByInfo{this});
    }

    std::optional<std::size_t> p_find_next_vertex() {
        m_best_infos.clear();
        for(std::size_t i : range(m_subgraph->n())) {
            m_best_infos.push(i);
        }
        auto is_covered = [&] (std::size_t v) {
            return m_vertex_info[v].in_some_class;
        };
        const auto& best = m_best_infos.elements();
        if(std::all_of(best.begin(), best.end(), is_covered)) {
            return std::nullopt;
        }
        return p_find_true_next_vertex(best);
    }

    static bool p_is_level_vertex(const SharedDBPropagator& propagator,
                                  std::int32_t level, Vertex vertex)
    {
        Lit l = *propagator.level_begin(level);
        return vertex.first == l || vertex.second == l;
    }

    auto p_implied_by(const SharedDBPropagator& propagator, Vertex vertex) {
        if(propagator.get_current_level() == 0 || !p_is_level_vertex(propagator, 1, vertex)) {
            return IteratorRange{propagator.get_trail().begin(), propagator.level_end(0)};
        }
        if(propagator.get_current_level() < 2 || !p_is_level_vertex(propagator, 2, vertex)) {
            return IteratorRange{propagator.get_trail().begin(), propagator.level_end(1)};
        }
        return IteratorRange{propagator.get_trail().begin(), propagator.level_end(2)};
    }

    LitOrFixed p_value_in_class(Lit l, std::size_t class_index) {
        std::vector<LitOrFixed>& class_vars = m_class_literals[class_index];
        const bool negate = lit::negative(l);
        LitOrFixed& var = class_vars[lit::var(l)];
        return std::visit(overloaded{
            [&] (bool& fixed) { return LitOrFixed{std::in_place_index<1>, fixed ^ negate}; },
            [&] (SatLit& lit) { return LitOrFixed{std::in_place_index<0>, negate ? -lit : lit};}
        }, var);
    }

    template<typename LitIterator>
    void p_class_literals_embed_clause(LitIterator clause_begin, LitIterator clause_end,
                                       std::size_t class_index)
    {
        m_clause_buffer.clear();
        bool found_true = false;
        for(Lit l : IteratorRange{clause_begin, clause_end}) {
            LitOrFixed value = p_value_in_class(l, class_index);
            std::visit(overloaded{
                [&] (bool& fixed) { if(fixed) found_true = true; },
                [&] (SatLit& lit) { m_clause_buffer.push_back(lit); }
            }, value);
            if(found_true) return;
        }
        if(m_clause_buffer.empty()) {
            throw std::logic_error("UNSAT class (has empty clause)!");
        }
        m_incremental_solver.add_clause(m_clause_buffer.begin(), m_clause_buffer.end());
    }

    void p_class_literals_embed_formula(const ClauseDB& clauses, std::size_t class_index) {
        for(auto [l1, l2] : clauses.binary_clauses()) {
            const Lit l[2] = {l1,l2};
            p_class_literals_embed_clause(+l, l+2, class_index);
        }
        for(CRef c = 1, ndb = clauses.literal_db_size(); c < ndb; c = clauses.next_clause(c)) {
            auto lits = clauses.lits_of(c);
            p_class_literals_embed_clause(lits.begin(), lits.end(), class_index);
        }
    }

    template<typename FixedRange>
    void p_init_class_literals_with_fixed_range(std::size_t class_index, FixedRange fixed_range) {
        const SharedDBPropagator& prop = m_color_classes[class_index].propagator;
        const Var num_vars = prop.db().num_vars();
        DynamicBitset fixed{num_vars, false};
        for(Lit l : fixed_range) {
            fixed[lit::var(l)].set();
        }
        std::vector<LitOrFixed> class_literals;
        for(Var v : range(num_vars)) {
            if(fixed[v]) {
                bool value = prop.is_true(lit::positive_lit(v));
                class_literals.emplace_back(std::in_place_index<1>, value);
            } else {
                class_literals.emplace_back(std::in_place_index<0>, m_incremental_solver.new_var());
            }
        }
        m_class_literals.emplace_back(std::move(class_literals));
        p_class_literals_embed_formula(prop.db(), class_index);
    }

    void p_init_class_literals_clique(std::size_t class_index) {
        const SharedDBPropagator& prop = m_color_classes[class_index].propagator;
        std::size_t clique_vertex_index = m_explicit_cover_order[class_index];
        Vertex clique_vertex = m_subgraph->vertex(clique_vertex_index);
        p_init_class_literals_with_fixed_range(class_index, p_implied_by(prop, clique_vertex));
    }

    void p_init_class_literals_not_clique(std::size_t class_index) {
        const SharedDBPropagator& prop = m_color_classes[class_index].propagator;
        p_init_class_literals_with_fixed_range(
            class_index, IteratorRange{prop.level_begin(0), prop.level_end(0)});
    }

    void p_init_vertex_in_class(std::size_t class_index)
    {
        if(class_index > m_vertex_in_class.size())
            throw std::logic_error("Error: class_index out of range!");
        if(class_index == m_vertex_in_class.size()) {
            m_vertex_in_class.emplace_back();
        }
        std::vector<LitOrFixed>& v_in_c = m_vertex_in_class[class_index];
        auto one_fixed = [&] (bool f1, SatLit l2) {
            if(!f1) { v_in_c.emplace_back(std::in_place_index<1>, false); }
            else { v_in_c.emplace_back(std::in_place_index<0>, l2); }
        };
        for(std::size_t vi_i : range(v_in_c.size(), m_explicit_cover_order.size())) {
            std::size_t vi = m_explicit_cover_order[vi_i];
            Vertex v = m_subgraph->vertex(vi);
            LitOrFixed val1 = p_value_in_class(v.first, class_index);
            LitOrFixed val2 = p_value_in_class(v.second, class_index);
            std::visit(overloaded{
                [&] (bool& f1, bool& f2) {v_in_c.emplace_back(std::in_place_index<1>, f1 && f2);},
                [&] (SatLit& l1, SatLit& l2) {
                    SatLit var = m_incremental_solver.new_var();
                    m_incremental_solver.add_short_clause(-var, l1);
                    m_incremental_solver.add_short_clause(-var, l2);
                    m_incremental_solver.add_short_clause(var, -l1, -l2);
                    v_in_c.emplace_back(std::in_place_index<0>, var);
                },
                [&] (bool& f1, SatLit& l2) { one_fixed(f1, l2); },
                [&] (SatLit& l1, bool& f2) { one_fixed(f2, l1); }
            }, val1, val2);
        }
    }

    void p_init_class_clique(std::size_t class_index) {
        p_init_class_literals_clique(class_index);
        p_init_vertex_in_class(class_index);
    }

    void p_init_class_not_clique(std::size_t class_index) {
        p_init_class_literals_not_clique(class_index);
        p_init_vertex_in_class(class_index);
    }

    void p_init_class_literals() {
        for(std::size_t i : range(m_clique_size)) {
            p_init_class_clique(i);
        }
        for(std::size_t i : range(m_clique_size, m_lower_bound)) {
            p_init_class_not_clique(i);
        }
    }

    void p_setup_vertex_constraints(SatLit or_not_enough, std::size_t begin_vertex, 
                                    std::size_t end_vertex) 
    {
        auto is_fixed = [] (const LitOrFixed& f) { return f.index() == 1; };
        for(std::size_t vi_i : range(begin_vertex, end_vertex)) {
            bool found_true = false;
            m_clause_buffer.clear();
            for(std::size_t cindex : range(m_lower_bound)) {
                const LitOrFixed& f = m_vertex_in_class[cindex][vi_i];
                if(is_fixed(f)) {
                    if(std::get<1>(f)) {
                        found_true = true;
                        break;
                    }
                } else {
                    m_clause_buffer.push_back(std::get<0>(f));
                }
            }
            if(found_true) continue;
            m_clause_buffer.push_back(or_not_enough);
            m_incremental_solver.add_clause(m_clause_buffer.begin(), m_clause_buffer.end());
        }
    }

    void p_setup_vertex_constraints(std::size_t old_n) {
        std::size_t not_enough_idx = m_lower_bound - m_clique_size;
        std::size_t new_n = m_vertex_in_class[0].size();
        while(not_enough_idx > m_not_enough_colors.size()) {
            m_not_enough_colors.emplace_back(std::in_place_index<1>, true);
        }
        if(not_enough_idx == m_not_enough_colors.size()) {
            SatLit var = m_incremental_solver.new_var();
            m_not_enough_colors.emplace_back(std::in_place_index<0>, var);
            p_setup_vertex_constraints(var, 0, new_n);
        } else {
            p_setup_vertex_constraints(
                std::get<0>(m_not_enough_colors[not_enough_idx]),
                old_n, new_n);
        }     
    }

    void p_add_class() {
        for(VertexInfo& info : m_vertex_info) {
            info.open_classes.push_back(true);
            ++info.num_open_classes;
        }
        m_color_classes.emplace_back(this, m_subgraph->get_propagator(), m_lower_bound-1);
    }

    std::size_t p_update_classes() {
        std::size_t old_num_vertices = 0;
        if(m_class_literals.empty()) {
            p_init_class_literals();
        } else {
            old_num_vertices = m_vertex_in_class[0].size();
            for(std::size_t cindex : range(m_class_literals.size())) {
                p_init_vertex_in_class(cindex);
            }
            for(std::size_t cindex : range(m_class_literals.size(), m_lower_bound)) {
                p_init_class_not_clique(cindex);
            }
        }
        return old_num_vertices;
    }

    void p_no_class_available() {
        m_explicit_cover_order.insert(m_explicit_cover_order.end(),
                                      m_true_zero_vertices.begin(), m_true_zero_vertices.end());
        for(;;) {
            if(m_aborted.load()) return;
            std::size_t old_num_vertices = p_update_classes();
            p_setup_vertex_constraints(old_num_vertices);
            std::vector<SatLit> assumptions;
            assumptions.push_back(-std::get<0>(m_not_enough_colors.back()));
            if(m_aborted.load()) return;
            auto result = m_incremental_solver.solve(assumptions);
            if(!result) { m_aborted.store(true); }
            if(m_aborted.load()) return;
            if(*result) {
                // SAT
                const auto& model = m_incremental_solver.get_model();
                p_handle_sat_assignment(model);
                return;
            } else {
                // UNSAT - repeat with more colors
                m_incremental_solver.fix(std::get<0>(m_not_enough_colors.back()));
                if(++m_lower_bound == m_upper_bound) return;
                p_add_class();
            }
        }
    }

    template<typename ModelType>
    bool p_vertex_in_class(const ModelType& model, std::size_t class_index, std::size_t vi_i) {
        LitOrFixed& lf = m_vertex_in_class[class_index][vi_i];
        return std::visit(overloaded{
            [&] (bool& b) -> bool { return b; },
            [&] (SatLit& l) -> bool { return model[l]; }
        }, lf);
    }

    void p_reset_info() {
        for(VertexInfo& v : m_vertex_info) {
            v.open_classes.set();
            v.num_open_classes = m_lower_bound;
            v.in_some_class = false;
        }
    }

    template<typename ModelType>
    void p_handle_sat_assignment(const ModelType& model) {
        p_reset_info();
        DynamicBitset vertex_is_in{m_explicit_cover_order.size(), false};
        for(std::size_t cindex : range(m_lower_bound)) {
            vertex_is_in.reset();
            for(std::size_t vi_i : range(m_explicit_cover_order.size())) {
                if(p_vertex_in_class(model, cindex, vi_i)) {
                    vertex_is_in[vi_i].set();
                }
            }
            m_color_classes[cindex].reset(this, vertex_is_in);
        }
    }

    void p_reset() {
        p_reset_info();
        DynamicBitset explicit_is_in{m_explicit_cover_order.size(), false};
        for(std::size_t cindex : range(m_clique_size)) {
            explicit_is_in[cindex].set();
            m_color_classes[cindex].reset(this, explicit_is_in);
            explicit_is_in[cindex].reset();
        }
        for(std::size_t cindex : range(m_clique_size, m_lower_bound)) {
            m_color_classes[cindex].reset(this, explicit_is_in);
        }
        m_explicit_cover_order.resize(m_clique_size);
        m_incremental_solver = IncrementalSolver{};
        m_class_literals.clear();
        m_not_enough_colors.clear();
        m_vertex_in_class.clear();
    }

    bool p_finalize() {
        struct EmptyHandler {};
        EmptyHandler handler;
        Var n_all = m_subgraph->get_propagator().db().num_vars();
        Var n_concrete = m_subgraph->get_n_concrete();
        ClassCompleter<EmptyHandler> completer{n_concrete, n_all, &handler};
        bool result = true;
        for(std::size_t cindex : range(m_color_classes.size())) {
            SharedDBPropagator& prop = m_color_classes[cindex].propagator;
            if(!completer.complete_class(prop)) {
                result = false;
            }
        }
        if(result) return true;
        p_reset();
        return false;
    }

  public:
    SATDSaturCoverer(UniverseSubgraph* subgraph, std::vector<std::size_t> clique, std::size_t upper_bound) :
        m_subgraph(subgraph),
        m_best_infos(20, CompareByInfo{this}),
        m_clique_size(clique.size()),
        m_lower_bound(clique.size()),
        m_upper_bound(upper_bound)
    {
        p_fill_vertices_with_literal();
        p_init_vertex_info();
        p_init_classes(std::move(clique));
    }

    void abort() {
        m_aborted.store(true);
        m_incremental_solver.terminate();
    }

    void run_coloring(bool until_complete = true) {
        for(;;) {
            if(m_aborted.load()) return;
            auto n = p_find_next_vertex();
            if(!n) {
                if(!until_complete || p_finalize()) {
                    m_success = true;
                    return;
                }
                continue;
            }
            const VertexInfo& next = m_vertex_info[*n];
            if(next.num_open_classes == 0) {
                p_no_class_available();
                if(m_lower_bound == m_upper_bound) {
                    m_success = true;
                    return;
                }
            } else {
                for(std::size_t cindex : next.open_classes.ones()) {
                    m_color_classes[cindex].add_vertex(this, *n);
                    m_explicit_cover_order.push_back(*n);
                    break;
                }
            }
        }
    }

    bool was_successful() const {
        return m_success;
    }

    bool did_improve() const {
        return m_lower_bound < m_upper_bound;
    }

    std::size_t get_lower_bound() const noexcept {
        return m_lower_bound;
    }

    std::vector<SharedDBPropagator> get_incomplete_solution() const {
        std::vector<SharedDBPropagator> props;
        props.reserve(m_color_classes.size());
        std::transform(m_color_classes.begin(), m_color_classes.end(),
                       std::back_inserter(props), [] (const ColorClass& cc) {
                           return cc.propagator;
                       });
        return props;
    }

    PartialSolution get_partial_solution() {
        if(!was_successful() || !did_improve()) {
            throw std::logic_error("Tried to obtain partial solution from unsucessful or non-improving coverer!");
        }
        PartialSolution solution{m_subgraph->get_propagator().db().num_vars(), m_subgraph->get_infeasibility_map()};
        for(const auto& cc : m_color_classes) {
            solution.add_class(cc.propagator);
        }
        return solution;
    }
};

}

#endif
