#ifndef SAMMY_SAT_DSATUR_H_INCLUDED_
#define SAMMY_SAT_DSATUR_H_INCLUDED_

#include "barrage.h"
#include "barrage_lns_subproblem.h"
#include "clause_db.h"
#include "shared_db_propagator.h"
#include <cassert>
#include <bitset>

namespace sammy {

/**
 * In most cases, the number of classes is small enough
 * to fit into a small static bitset, say 64 or 128.
 */
template<std::size_t NClassesMax>
class StaticOpenClassTracker {
  public:
    explicit StaticOpenClassTracker(std::size_t num_classes)  :
        m_open_classes(),
        m_num_open_classes(num_classes)
    {
        assert(num_classes <= NClassesMax);
        m_open_classes.set();
        if(num_classes < NClassesMax) {
            m_open_classes >>= (NClassesMax - num_classes);
        }
        assert(m_open_classes.count() == num_classes);
    }

    bool close_class(std::size_t class_index) {
        assert(class_index < NClassesMax);
        if(m_open_classes[class_index]) {
            m_open_classes[class_index] = false;
            --m_num_open_classes;
            return true;
        }
        return false;
    }

    void reset(std::size_t num_classes) {
        assert(num_classes <= NClassesMax);
        m_open_classes.set();
        m_num_open_classes = num_classes;
        if(num_classes < NClassesMax) {
            m_open_classes >>= (NClassesMax - num_classes);
        }
        assert(m_open_classes.count() == num_classes);
    }

    std::size_t num_open() const noexcept {
        return m_num_open_classes;
    }

    template<typename Func/*(std::size_t one_index) -> bool*/>
    void iterate(Func&& func) {
        for(std::size_t i = 0, r = m_num_open_classes; r != 0; ++i) {
            if(m_open_classes[i]) {
                if(!func(i)) {
                    break;
                }
                --r;
            }
        }
    }

  private:
    std::bitset<NClassesMax> m_open_classes;
    std::size_t m_num_open_classes;
};

/**
 * In some cases, the number of classes is large enough
 * to require a dynamic bitset.
 */
class DynamicOpenClassTracker {
  public:
    explicit DynamicOpenClassTracker(std::size_t num_classes) :
        m_open_classes(num_classes, true),
        m_num_open_classes(num_classes)
    {}

    std::size_t num_open() const noexcept {
        return m_num_open_classes;
    }

    bool close_class(std::size_t class_index) {
        assert(class_index < m_open_classes.size());
        if (m_open_classes[class_index]) {
            m_open_classes[class_index] = false;
            --m_num_open_classes;
            return true;
        }
        return false;
    }

    void reset(std::size_t num_classes) {
        assert(num_classes == m_open_classes.size());
        m_open_classes.set();
        m_num_open_classes = num_classes;
    }

    template<typename Func/*(std::size_t one_index) -> bool*/>
    void iterate(Func&& func) {
        for(std::size_t i : m_open_classes.ones()) {
            if(!func(i)) {
                break;
            }
        }
    }

  private:
    DynamicBitset m_open_classes;
    std::size_t m_num_open_classes;
};

/**
 * Fixed mutually exclusive set solver
 * based on the DSatur-style heuristic
 * for determining which interactions to
 * cover explicitly next.
 */
template<typename IncrementalSatSolver>
class FixedMESSatDSaturSolver {
    using SatSolver = IncrementalSatSolver;
    using Lit = typename SatSolver::Lit;
    using SLit = sammy::Lit;
    using LitOrVal = std::variant<Lit, bool>;
    
    std::atomic<bool> m_abort{false};
    PortfolioSolver* m_portfolio;
    LNSSubproblem m_subproblem;
    SharedDBPropagator m_propagator;
    EventRecorder* m_recorder;
    std::size_t m_worker_id;
    std::size_t m_n_all;
    std::size_t m_n_concrete;
    bool m_infeasible_by_construction{false};
    std::size_t m_num_configs_allowed;
    std::vector<std::vector<LitOrVal>> m_config_vars;
    SatSolver m_solver;
    std::vector<Lit> m_buffer;
    std::vector<DynamicBitset> m_vertices_implying_literal;
    std::vector<std::vector<std::size_t>> m_vertices_containing_concrete;

    template<typename Tracker>
    struct VertexInfo {
        VertexInfo(std::size_t num_classes) :
            class_tracker(num_classes)
        {}

        std::size_t num_open() const noexcept {
            return class_tracker.num_open();
        }

        Tracker class_tracker;
        bool in_some_class{false}, explicitized{false};
    };
    std::vector<VertexInfo<StaticOpenClassTracker<64>>> m_vertex_info_64;
    std::vector<VertexInfo<StaticOpenClassTracker<128>>> m_vertex_info_128;
    std::vector<VertexInfo<DynamicOpenClassTracker>> m_vertex_info_d;
    std::vector<SharedDBPropagator> m_color_classes;
    std::vector<std::size_t> m_explicitized;
    // m_covered_in_class[i][ci]: coverage of vertex m_explicitized[i]
    // in class m_color_classes[ci]
    std::vector<std::vector<LitOrVal>> m_covered_in_class;
    std::vector<std::size_t> m_low_class_tiers[5];
    std::vector<std::int32_t> m_symmetry_breaking_levels;
    std::optional<std::vector<DynamicBitset>> m_final_solution;

  public:
    FixedMESSatDSaturSolver(
        PortfolioSolver* portfolio, LNSSubproblem&& subproblem,
        SharedDBPropagator prop, EventRecorder* recorder,
        std::size_t worker_id,
        std::size_t max_configs_allowed = 
            std::numeric_limits<std::size_t>::max()
    ) :
        m_portfolio(portfolio),
        m_subproblem(std::move(subproblem)),
        m_propagator(std::move(prop)),
        m_recorder(recorder),
        m_worker_id(worker_id),
        m_n_all(m_propagator.db().num_vars()),
        m_n_concrete(m_subproblem.num_concrete),
        m_num_configs_allowed(p_num_classes(max_configs_allowed))
    {
        try {
            if(m_subproblem.mutually_exclusive_set.size() > 
               m_num_configs_allowed) {
                m_infeasible_by_construction = true;
                return;
            }
            p_init_vertices_implying_literal();
            p_init_vertices_containing_concrete();
            p_init_vertex_info();
            p_init_color_classes();
            p_assign_mes_to_color_classes();
            p_init_tiers();
        } catch(...) {
            subproblem = std::move(m_subproblem);
            throw;
        }
    }

    LNSSubproblem move_out_subproblem() {
        return std::move(m_subproblem);
    }

    const std::vector<Vertex>& mes_vertices() const {
        return m_subproblem.mutually_exclusive_set;
    }

    static std::string name() {
        return std::string("FixedMESSatDSatur<") + SatSolver::name() + ">";
    }

    std::string strategy_name() const {
        return name();
    }

    void abort() {
        m_abort.store(true, std::memory_order_relaxed);
        m_solver.terminate();
    }

    const std::vector<DynamicBitset>& get_solution() const {
        return *m_final_solution;
    }

    std::optional<bool> solve() {
        if(m_infeasible_by_construction) {
            return false;
        }
        if(m_final_solution) {
            return true;
        }
        try {
            for(;;) {
                int res = p_continue();
                if(res == 0) {
                    return false;
                } else if(res == 1) {
                    return true;
                } else if(res == -1) {
                    return std::nullopt;
                }
            }
        } catch(InterruptError&) {
            return std::nullopt;
        }
    }

    bool check_solution_coverage(const std::vector<DynamicBitset>& solution) {
        std::vector<Vertex> unchecked = m_subproblem.uncovered_universe;
        for(const DynamicBitset& config : solution) {
            auto new_end = std::remove_if(
                unchecked.begin(), unchecked.end(),
                [&config](const Vertex& v) {
                    sammy::Var v1 = lit::var(v.first), v2 = lit::var(v.second);
                    bool neg1 = lit::negative(v.first);
                    bool neg2 = lit::negative(v.second);
                    return config[v1] != neg1 && config[v2] != neg2;
                }
            );
            unchecked.erase(new_end, unchecked.end());
        }
        return unchecked.empty();
    }

  private:
    std::size_t p_num_classes(std::size_t max_allowed) {
        if(max_allowed == std::numeric_limits<std::size_t>::max()) {
            return m_subproblem.removed_configurations.size() - 1;
        }
        return max_allowed;
    }

    void p_init_vertex_info() {
        if(m_num_configs_allowed <= 64) {
            p_init_vertex_info(&FixedMESSatDSaturSolver::m_vertex_info_64);
        } else if(m_num_configs_allowed <= 128) {
            p_init_vertex_info(&FixedMESSatDSaturSolver::m_vertex_info_128);
        } else {
            p_init_vertex_info(&FixedMESSatDSaturSolver::m_vertex_info_d);
        }
    }

    template<typename VertexInfosType>
    void p_init_tiers(VertexInfosType& vertex_info) {
        constexpr std::size_t num_tiers = sizeof(m_low_class_tiers) / 
                                          sizeof(m_low_class_tiers[0]);
        constexpr std::size_t ulim_tier = num_tiers - 1;
        std::size_t n = m_subproblem.uncovered_universe.size();
        for(std::size_t vi : range(n)) {
            auto& vinfo = vertex_info[vi];
            if(vinfo.in_some_class) {
                continue;
            }
            std::size_t num_open = vinfo.num_open();
            std::size_t tier = (std::min)(num_open, ulim_tier);
            m_low_class_tiers[tier].push_back(vi);
        }
        if(!m_low_class_tiers[1].empty()) {
            // additional symmetry breaking possible
            p_continue_tier1(vertex_info);
        }
        if(!m_low_class_tiers[0].empty()) {
            // just adding MES interactions gave infeasible interactions
            m_infeasible_by_construction = true;
        }
    }

    void p_init_tiers() {
        if(m_num_configs_allowed <= 64) {
            p_init_tiers(m_vertex_info_64);
        } else if(m_num_configs_allowed <= 128) {
            p_init_tiers(m_vertex_info_128);
        } else {
            p_init_tiers(m_vertex_info_d);
        }
        if(!m_infeasible_by_construction) {
            m_symmetry_breaking_levels.reserve(m_num_configs_allowed);
            for(const auto& p : m_color_classes) {
                m_symmetry_breaking_levels.push_back(p.get_current_level());
            }
        }
    }

    template<typename MemberPtrType>
    void p_init_vertex_info(MemberPtrType member_ptr) {
        auto& vertex_info = this->*member_ptr;
        std::size_t n = m_subproblem.uncovered_universe.size();
        vertex_info.reserve(n);
        for(std::size_t unused : range(n)) {
            static_cast<void>(unused);
            vertex_info.emplace_back(m_num_configs_allowed);
        }
    }

    /**
     * Sort the vertices in the subproblem universe
     * and initialize the m_vertices_implying_literal matrix.
     */
    void p_init_vertices_implying_literal() {
        std::transform(m_subproblem.uncovered_universe.begin(),
                       m_subproblem.uncovered_universe.end(),
                       m_subproblem.uncovered_universe.begin(),
                       [](Vertex v) {
                           return std::make_pair((std::min)(v.first, v.second), 
                                                 (std::max)(v.first, v.second));
                       });
        std::sort(m_subproblem.uncovered_universe.begin(),
                  m_subproblem.uncovered_universe.end());
        m_vertices_implying_literal.reserve(2 * m_n_all);
        std::size_t usize = m_subproblem.uncovered_universe.size();
        for(std::size_t i = 0; i < 2 * m_n_all; ++i) {
            m_vertices_implying_literal.emplace_back(usize, false);
        }
        m_propagator.reset_to_zero();
        SLit last_first = sammy::NIL;
        for(std::size_t vi = 0; vi < usize; ++vi) {
            Vertex v = m_subproblem.uncovered_universe[vi];
            SLit first = v.first, second = v.second;
            if(last_first != first) {
                m_propagator.reset_to_zero();
                if(!m_propagator.is_true(first)) {
                    if(m_propagator.is_false(first) || 
                       !m_propagator.push_level(first)) {
                        throw std::logic_error("Infeasible interaction "
                                               "in subproblem universe!");
                    }
                }
                last_first = first;
            }
            bool second_pushed = false;
            if(!m_propagator.is_true(second)) {
                if(m_propagator.is_false(second) || 
                   !(second_pushed = m_propagator.push_level(second))) {
                    throw std::logic_error("Infeasible interaction "
                                           "in subproblem universe!");
                }
            }
            for(SLit l : m_propagator.get_trail()) {
                m_vertices_implying_literal[l][vi].set();
            }
            if(second_pushed) {
                m_propagator.pop_level();
            }
        }
        m_propagator.reset_to_zero();
    }

    void p_init_vertices_containing_concrete() {
        m_vertices_containing_concrete.resize(2 * m_n_concrete);
        for(std::size_t vi = 0, n = m_subproblem.uncovered_universe.size();
            vi < n; ++vi) {
            Vertex v = m_subproblem.uncovered_universe[vi];
            SLit first = v.first, second = v.second;
            m_vertices_containing_concrete[first].push_back(vi);
            m_vertices_containing_concrete[second].push_back(vi);
        }
    }

    void p_init_color_classes() {
        m_color_classes.assign(m_num_configs_allowed, m_propagator);
    }

    void p_assign_mes_to_color_classes() {
        if(m_num_configs_allowed <= 64) {
            p_assign_mes_to_color_classes(
                &FixedMESSatDSaturSolver::m_vertex_info_64);
        } else if(m_num_configs_allowed <= 128) {
            p_assign_mes_to_color_classes(
                &FixedMESSatDSaturSolver::m_vertex_info_128);
        } else {
            p_assign_mes_to_color_classes(
                &FixedMESSatDSaturSolver::m_vertex_info_d);
        }
    }

    std::size_t p_find_vertex(Vertex v) const {
        v = std::make_pair((std::min)(v.first, v.second), 
                       (std::max)(v.first, v.second));
        auto it = std::lower_bound(m_subproblem.uncovered_universe.begin(),
                                   m_subproblem.uncovered_universe.end(), v);
        if(it == m_subproblem.uncovered_universe.end() || *it != v) {
            throw std::logic_error("Vertex not found in uncovered universe!");
        }
        return std::distance(m_subproblem.uncovered_universe.begin(), it);
    }

    template<typename MemberPtrType>
    void p_assign_mes_to_color_classes(MemberPtrType member_ptr) {
        auto& vertex_info = this->*member_ptr;
        std::size_t ci = 0;
        for(Vertex v : m_subproblem.mutually_exclusive_set) {
            std::size_t vi = p_find_vertex(v);
            assert(!vertex_info[vi].in_some_class);
            assert(m_subproblem.uncovered_universe.at(vi) == v);
            p_add_to_class(ci++, vi, v, vertex_info);
            m_explicitized.push_back(vi);
            vertex_info[vi].explicitized = true;
        }
    }

    template<typename VertexInfoType>
    void p_add_to_class(std::size_t class_index, std::size_t vertex_index,
                        VertexInfoType& vertex_info)
    {
        p_add_to_class(class_index, vertex_index, 
                       m_subproblem.uncovered_universe[vertex_index],
                       vertex_info);
    }

    template<typename VertexInfosType>
    void p_add_to_class(std::size_t class_index, std::size_t vertex_index,
                        Vertex vertex, VertexInfosType& vertex_info)
    {
        auto& info = vertex_info[vertex_index];
        auto& prop = m_color_classes[class_index];
        std::size_t old_trail_length = prop.get_trail().size();
        int pres = push_vertex(prop, vertex);
        if(pres < 0) {
            throw std::logic_error("Infeasible interaction in "
                                   "subproblem universe!");
        } else {
            info.in_some_class = true;
            p_added_to_class(class_index, vertex_index, prop,
                             vertex_info, old_trail_length);
        }
    }

    template<typename VertexInfosType>
    bool p_try_add_to_class(std::size_t class_index, std::size_t vertex_index,
                            VertexInfosType& vertex_info)
    {
        auto& prop = m_color_classes[class_index];
        std::size_t old_trail_length = prop.get_trail().size();
        Vertex vertex = m_subproblem.uncovered_universe[vertex_index];
        int pres = push_vertex(prop, vertex);
        if(pres < 0) {
            return false;
        }
        auto& info = vertex_info[vertex_index];
        info.in_some_class = true;
        info.explicitized = true;
        m_explicitized.push_back(vertex_index);
        p_added_to_class(class_index, vertex_index, prop,
                         vertex_info, old_trail_length);
        return true;
    }

    template<typename VertexInfosType>
    void p_added_to_class_new_concrete(SharedDBPropagator& class_prop,
                                       SLit new_lit,
                                       VertexInfosType& vertex_info)
    {
        for(std::size_t vpot : m_vertices_containing_concrete[new_lit]) {
            Vertex vo = m_subproblem.uncovered_universe[vpot];
            SLit other = vo.first;
            if(other == new_lit) {
                other = vo.second;
            }
            if(class_prop.is_true(other)) {
                vertex_info[vpot].in_some_class = true;
            }
        }
    }

    template<typename VertexInfosType>
    void p_added_to_class(std::size_t class_index, std::size_t vertex_index,
                          SharedDBPropagator& class_prop, 
                          VertexInfosType& vertex_info,
                          std::size_t old_trail_length)
    {
        const auto& trail = class_prop.get_trail();
        auto nclit = 2 * m_n_concrete;
        constexpr std::size_t lim_tier = 
            sizeof(m_low_class_tiers) / sizeof(m_low_class_tiers[0]) - 2;
        for(std::size_t i = old_trail_length; i < trail.size(); ++i) {
            SLit new_lit = trail[i];
            if (new_lit < nclit) {
                p_added_to_class_new_concrete(class_prop, new_lit, vertex_info);
            }
            SLit new_neg = lit::negate(new_lit);
            for(std::size_t vi : m_vertices_implying_literal[new_neg].ones()) {
                auto &vinf = vertex_info[vi];
                if(vinf.in_some_class) {
                    continue;
                }
                if(vinf.class_tracker.close_class(class_index)) {
                    std::size_t new_open = vinf.num_open();
                    if(new_open <= lim_tier) {
                        m_low_class_tiers[new_open].push_back(vi);
                    }
                }
            }
        }
    }

    static bool p_vertex_true_in(const SharedDBPropagator& prop, Vertex v) {
        return prop.is_true(v.first) && prop.is_true(v.second);
    }

    template<typename VertexInfosType>
    bool p_check_in_some_class(VertexInfosType& vertex_info) const {
        // check the sanity of .in_some_class, and .explicitized;
        // also check tiers
        HashSet<std::size_t> uncovered;
        std::size_t n = m_subproblem.uncovered_universe.size();
        for(std::size_t vi : range(n)) {
            const auto& vinfo = vertex_info[vi];
            Vertex v = m_subproblem.uncovered_universe[vi];
            if(vinfo.explicitized && !vinfo.in_some_class) {
                return false;
            }
            if(vinfo.in_some_class) {
                bool found = false;
                for(std::size_t ci : range(m_num_configs_allowed)) {
                    if(p_vertex_true_in(m_color_classes[ci], v)) {
                        found = true;
                        break;
                    }
                }
                if(!found) {
                    return false;
                }
            } else {
                uncovered.insert(vi);
            }
        }
        for(const auto& tier : m_low_class_tiers) {
            for(std::size_t vi : tier) {
                uncovered.erase(vi);
            }
        }
        if(!uncovered.empty()) {
            return false;
        }
        return true;
    }

    template<typename VertexInfosType>
    int p_continue(VertexInfosType& vertex_info) {
        constexpr std::size_t lim_tier = 
            sizeof(m_low_class_tiers) / sizeof(m_low_class_tiers[0]) - 2;
        if(m_abort.load(std::memory_order_relaxed)) {
            return -1;
        }
        assert(p_check_in_some_class(vertex_info));
        static_cast<void>(&FixedMESSatDSaturSolver::
                          p_check_in_some_class<VertexInfosType>);
        if(!m_low_class_tiers[0].empty()) {
            int res = p_continue_using_sat(vertex_info);
            return res;
        }
        if(!m_low_class_tiers[1].empty()) {
            // handle interactions with only one open class
            auto res = p_continue_tier1(vertex_info);
            // return if we actually changed the state
            if(res) {
                return *res;
            }
        }
        for(std::size_t t = 2; t <= lim_tier; ++t) {
            if(!m_low_class_tiers[t].empty()) {
                // select interactions from tier t
                auto res = p_continue_lim_tier(t, vertex_info);
                // return if we actually changed the state
                if(res) {
                    return *res;
                }
            }
        }
        if(!m_low_class_tiers[lim_tier + 1].empty()) {
            // select interactions from the unlimited #options tier
            auto res = p_continue_ulim_tier(lim_tier + 1, vertex_info);
            // return if we actually changed the state
            if(res) {
                return *res;
            }
        }
        // no more uncovered interactions - complete all partial configs
        int cc = p_continue_complete_all_classes(vertex_info);
        return cc;
    }

    template<typename VertexInfosType>
    bool p_find_open_class(VertexInfosType& vertex_info, std::size_t vi) {
        auto& vinfo = vertex_info[vi];
        Vertex v = m_subproblem.uncovered_universe[vi];
        bool found = false;
        vinfo.class_tracker.iterate([&] (std::size_t class_index) {
            auto& prop = m_color_classes[class_index];
            if(prop.is_true(v.first) || prop.is_true(v.second)) {
                if(p_try_add_to_class(class_index, vi, vertex_info)) {
                    found = true;
                    return false;
                } else {
                    vinfo.class_tracker.close_class(class_index);
                }
            }
            return true;
        });
        if(found) return true;
        vinfo.class_tracker.iterate([&] (std::size_t class_index) {
            if(p_try_add_to_class(class_index, vi, vertex_info)) {
                found = true;
                return false;
            } else {
                vinfo.class_tracker.close_class(class_index);
            }
            return true;
        });
        return found;
    }

    template<typename VertexInfosType>
    std::optional<int> p_continue_lim_tier(std::size_t tier,
                                           VertexInfosType& vertex_info)
    {
        auto& t1 = m_low_class_tiers[tier];
        auto& tl = m_low_class_tiers[tier - 1];
        std::size_t added = 0;
        while(!t1.empty() && tl.empty()) {
            std::size_t vi = t1.back();
            t1.pop_back();
            auto& vinfo = vertex_info[vi];
            assert(!vinfo.explicitized || vinfo.in_some_class);
            if(vinfo.in_some_class) {
                continue;
            }
            if(p_find_open_class(vertex_info, vi)) {
                ++added;
            } else {
                m_low_class_tiers[0].push_back(vi);
                return 2;
            }
        }
        if(added)
            return 2;
        return std::nullopt;
    }

    template<typename VertexInfosType>
    std::optional<int> p_continue_ulim_tier(std::size_t tier, 
                                            VertexInfosType& vertex_info) 
    {
        auto& t = m_low_class_tiers[tier];
        auto& tl = m_low_class_tiers[tier - 1];
        std::size_t min_open = std::numeric_limits<std::size_t>::max();
        std::vector<std::size_t> at_min_open;
        auto new_end = std::remove_if(t.begin(), t.end(),
            [&] (std::size_t vi) {
                auto& vinfo = vertex_info[vi];
                if(vinfo.in_some_class) {
                    return true;
                }
                std::size_t num_open = vinfo.num_open();
                if(num_open < min_open) {
                    min_open = num_open;
                    at_min_open.clear();
                    at_min_open.push_back(vi);
                } else if(num_open == min_open) {
                    at_min_open.push_back(vi);
                }
                return false;
            }
        );
        t.erase(new_end, t.end());
        if(t.empty()) {
            return std::nullopt;
        }
        while(!at_min_open.empty() && tl.empty()) {
            std::size_t vi = at_min_open.back();
            at_min_open.pop_back();
            auto& vinfo = vertex_info[vi];
            if(vinfo.in_some_class) {
                continue;
            }
            if(!p_find_open_class(vertex_info, vi)) {
                m_low_class_tiers[0].push_back(vi);
                return 2;
            }
        }
        return 2;
    }

    int p_continue() {
        if(m_num_configs_allowed <= 64) {
            return p_continue(m_vertex_info_64);
        } else if(m_num_configs_allowed <= 128) {
            return p_continue(m_vertex_info_128);
        } else {
            return p_continue(m_vertex_info_d);
        }
    }

    template<typename VertexInfosType>
    std::optional<int> p_continue_tier1(VertexInfosType& vertex_info) {
        std::size_t added = 0;
        auto& t = m_low_class_tiers[1];
        auto& t0 = m_low_class_tiers[0];
        while(!t.empty()) {
            std::size_t vi = t.back();
            t.pop_back();
            auto& vinfo = vertex_info[vi];
            assert(!vinfo.explicitized || vinfo.in_some_class);
            if(vinfo.in_some_class) {
                continue;
            }
            assert(vinfo.num_open() == 1);
            vinfo.class_tracker.iterate([&](std::size_t class_index) {
                if(p_try_add_to_class(class_index, vi, vertex_info)) {
                    ++added;
                } else {
                    vinfo.class_tracker.close_class(class_index);
                    t0.push_back(vi);
                }
                return false;
            });
            if(!t0.empty()) {
                return 2;
            }
        }
        if(added) {
            return 2;
        }
        return std::nullopt;
    }

    template<typename VertexInfosType>
    int p_continue_using_sat(VertexInfosType& vertex_info) {
        assert(std::none_of(m_low_class_tiers[0].begin(),
                            m_low_class_tiers[0].end(),
                            [&](std::size_t vi) {
                                return vertex_info[vi].in_some_class;
                            }));
        constexpr std::size_t tiers = 
            sizeof(m_low_class_tiers) / sizeof(m_low_class_tiers[0]);
        std::size_t prev_sat_explicitized = m_covered_in_class.size();
        if(prev_sat_explicitized == 0) {
            assert(m_config_vars.empty());
            p_explicitize_tier(0, vertex_info);
            p_explicitize_tier(1, vertex_info);
            for(std::size_t i = 2; i < tiers; ++i) {
                if(m_explicitized.size() >= 100) {
                    break;
                }
                std::size_t new_limit = 100 - m_explicitized.size();
                p_sample_explicitize_tier(i, vertex_info, new_limit);
            }
            if(!p_sat_first_setup()) {
                return 0;
            }
        } else {
            assert(m_config_vars.size() == m_num_configs_allowed);
            double universe = m_subproblem.uncovered_universe.size();
            double now_sat_explicitized = 1.33 * prev_sat_explicitized;
            if(now_sat_explicitized >= 0.75 * universe) {
                // explicitize all interactions
                p_explicitize_all(vertex_info);
            } else {
                std::size_t limit = std::size_t(now_sat_explicitized);
                p_explicitize_tier(0, vertex_info);
                for(std::size_t tier = 1; tier < tiers; ++tier) {
                    std::size_t cexp = m_explicitized.size();
                    if(cexp >= limit) {
                        break;
                    }
                    auto &t = m_low_class_tiers[tier];
                    if(cexp + t.size() <= limit) {
                        p_explicitize_tier(tier, vertex_info);
                    } else {
                        p_sample_explicitize_tier(tier, vertex_info, 
                                                  limit - cexp);
                    }
                }
            }
        }
        if(!p_sat_extend_covered_in_class()) {
            return 0;
        }
        return p_sat_solve(vertex_info);
    }

    template<typename VertexInfosType>
    int p_sat_solve(VertexInfosType& vertex_info) {
        auto res = m_solver.solve();
        if(!res) {
            // SAT solver was interrupted
            return -1;
        }
        if(!*res) {
            // known UNSAT; no possible improvement
            return 0;
        }
        // SAT solver found a solution
        if(m_explicitized.size() == m_subproblem.uncovered_universe.size()) {
            // all interactions were explicitized
            p_sat_solution_to_final_solution();
            return 1;
        }
        p_sat_solution_to_classes(vertex_info);
        return 2;
    }

    void p_sat_solution_to_final_solution() {
        std::vector<DynamicBitset> solution;
        auto model_map = m_solver.get_model();
        for(std::size_t ci : range(m_num_configs_allowed)) {
            DynamicBitset assignment(m_n_all, false);
            const auto& ci_vals = m_config_vars[ci];
            for(std::size_t li = 0; li < m_n_all; ++li) {
                bool v = std::visit(overloaded{
                        [&] (Lit l) -> bool {
                            return model_map[l];
                        },
                        [&] (bool b) -> bool {
                            return b;
                        }
                    }, ci_vals[li]);
                if(v) {
                    assignment[li].set();
                }
            }
            solution.push_back(std::move(assignment));
        }
        m_final_solution = std::move(solution);
    }

    /**
     * Validate that the propagator believes each
     * configuration in the solution from the
     * SAT solver to be feasible.
     * @param model_map The model map (solution) from the SAT solver.
     */
    template<typename ModelMap>
    bool p_check_sat_solution_vs_propagator(const ModelMap& model_map) {
        auto check_class = [&] (std::size_t ci) -> bool {
            auto checked_push = [&] (SLit l) {
                if(m_propagator.is_true(l)) return;
                if(m_propagator.is_false(l) || !m_propagator.push_level(l)) {
                    throw std::logic_error("SAT solution is infeasible!");
                }
            };
            m_propagator.reset_to_zero();
            const auto& ci_vals = m_config_vars[ci];
            assert(ci_vals.size() == m_n_all);
            for(std::size_t li = 0; li < m_n_all; ++li) {
                auto v = static_cast<sammy::Var>(li);
                std::visit(overloaded{
                    [&] (const Lit& l) {
                        if(model_map[l]) {
                            checked_push(lit::positive_lit(v));
                        } else {
                            checked_push(lit::negative_lit(v));
                        }
                    }, [&] (const bool& b) {
                        if(b) {
                            checked_push(lit::positive_lit(v));
                        } else {
                            checked_push(lit::negative_lit(v));
                        }
                    }
                }, ci_vals[li]);
            }
            return !m_propagator.is_conflicting();
        };
        auto r = range(m_num_configs_allowed);
        assert(std::all_of(std::begin(r), std::end(r), check_class));
        assert((m_propagator.reset_to_zero(), true));
        static_cast<void>(r);
        static_cast<void>(check_class);
        return true;
    }

    template<typename VertexInfosType>
    void p_sat_solution_to_classes(VertexInfosType& vertex_info) {
        auto model_map = m_solver.get_model();
        p_check_sat_solution_vs_propagator(model_map);
        auto bool_value = [&] (LitOrVal l) -> bool {
            return std::visit(overloaded{
                [&] (const Lit& l) -> bool {
                    return model_map[l];
                },
                [&] (const bool& b) -> bool {
                    return b;
                }
            }, l);
        };

        // reset our classes and tiers and in_some_class flags
        for(auto& cclass : m_color_classes) {
            cclass.reset_to_zero();
        }
        for(auto& t : m_low_class_tiers) {
            t.clear();
        }
        std::vector<std::size_t> still_to_check;
        for(std::size_t vi : range(vertex_info.size())) {
            auto& vinfo = vertex_info[vi];
            if(!vinfo.explicitized) {
                vinfo.in_some_class = false;
                vinfo.class_tracker.reset(m_num_configs_allowed);
                still_to_check.push_back(vi);
            }
        }

        // add explicitized
        for(std::size_t exi = 0, n_ex = m_explicitized.size(); 
            exi < n_ex; ++exi) 
        {
            std::size_t vi = m_explicitized[exi];
            Vertex v = m_subproblem.uncovered_universe[vi];
            const auto& covered_in = m_covered_in_class[exi];
            bool found = false;
            for(std::size_t ci : range(m_num_configs_allowed)) {
                if(bool_value(covered_in[ci])) {
                    auto& prop = m_color_classes[ci];
                    if(push_vertex(prop, v) < 0) {
                        throw std::logic_error("SAT solution is infeasible!");
                    }
                    found = true;
                    break;
                }
            }
            if(!found) {
                throw std::logic_error("SAT solution does not cover "
                                       "explicitized interaction!");
            }
        }

        // we need to recompute the flags and regenerate our tiers
        // begin by computing the in_some_class flags
        const auto& universe = m_subproblem.uncovered_universe;
        for(std::size_t ci = 0; ci < m_num_configs_allowed; ++ci) {
            SLit skip_first = sammy::NIL;
            SLit cov_first = sammy::NIL;
            const auto& prop = m_color_classes[ci];
            auto removal = [&] (std::size_t vi) -> bool {
                Vertex v = universe[vi];
                if(v.first == skip_first) {
                    return false;
                }
                if(v.first == cov_first) {
                    if(prop.is_true(v.second)) {
                        vertex_info[vi].in_some_class = true;
                        return true;
                    }
                    return false;
                }
                if(prop.is_true(v.first)) {
                    cov_first = v.first;
                    if(prop.is_true(v.second)) {
                        vertex_info[vi].in_some_class = true;
                        return true;
                    }
                } else {
                    skip_first = v.first;
                }
                return false;
            };
            auto new_end = std::remove_if(still_to_check.begin(),
                                          still_to_check.end(), removal);
            still_to_check.erase(new_end, still_to_check.end());
        }

        // now we compute the open classes where needed
        for(std::size_t ci = 0; ci < m_num_configs_allowed; ++ci) {
            const auto& prop = m_color_classes[ci];
            for(SLit l : prop.get_trail()) {
                SLit lneg = lit::negate(l);
                for(std::size_t vi : m_vertices_implying_literal[lneg].ones()) {
                    auto& vinfo = vertex_info[vi];
                    if(vinfo.in_some_class) {
                        continue;
                    }
                    vinfo.class_tracker.close_class(ci);
                }
            }
        }
        // finally, we can recompute the tiers
        constexpr std::size_t last_tier = 
            sizeof(m_low_class_tiers) / sizeof(m_low_class_tiers[0]) - 1;
        for(std::size_t vi : still_to_check) {
            auto& vinfo = vertex_info[vi];
            std::size_t tier = (std::min)(last_tier, vinfo.num_open());
            m_low_class_tiers[tier].push_back(vi);
        }
    }

    template<typename VertexInfosType>
    void p_explicitize_tier(std::size_t tier, VertexInfosType& vertex_info) 
    {
        auto& t = m_low_class_tiers[tier];
        for(std::size_t vi : t) {
            auto& vinfo = vertex_info[vi];
            if(vinfo.explicitized) {
                continue;
            }
            m_explicitized.push_back(vi);
            vinfo.explicitized = true;
            vinfo.in_some_class = true;
        }
    }

    template<typename VertexInfosType>
    void p_explicitize_all(VertexInfosType& vertex_info) {
        for(std::size_t vi = 0, n = m_subproblem.uncovered_universe.size();
            vi < n; ++vi) 
        {
            auto& vinfo = vertex_info[vi];
            if(vinfo.explicitized) {
                continue;
            }
            m_explicitized.push_back(vi);
            vinfo.explicitized = true;
            vinfo.in_some_class = true;
        }
    }

    template<typename VertexInfosType>
    void p_sample_explicitize_tier(std::size_t tier,
                                   VertexInfosType& vertex_info,
                                   std::size_t limit_new)
    {
        auto& t = m_low_class_tiers[tier];
        auto& rng = sammy::rng();
        if(limit_new * 5 >= t.size()) {
            // limit is relatively large compared to tier
            p_reduce_tier(t, vertex_info);
            if(t.size() <= limit_new) {
                p_explicitize_tier(tier, vertex_info);
                return;
            }
            std::shuffle(t.begin(), t.end(), rng);
            for(std::size_t i = 0; i < limit_new; ++i) {
                std::size_t vi = t[i];
                auto& vinfo = vertex_info[vi];
                m_explicitized.push_back(vi);
                vinfo.explicitized = true;
                vinfo.in_some_class = true;
            }
            return;
        }

        // limit seems relatively small compared to tier
        std::uniform_int_distribution<std::size_t> dist(0, t.size() - 1);
        std::size_t since_last_success = 0;
        std::size_t failure_limit = 30;
        while(limit_new > 0) {
            std::size_t ti = dist(rng);
            std::size_t vi = t[ti];
            auto& vinfo = vertex_info[vi];
            if(vinfo.explicitized) {
                if(++since_last_success >= failure_limit) {
                    break;
                }
            } else {
                since_last_success = 0;
                m_explicitized.push_back(vi);
                vinfo.explicitized = true;
                vinfo.in_some_class = true;
                --limit_new;
            }
        }
        if(!limit_new) {
            return;
        }
        p_reduce_tier(t, vertex_info);
        p_sample_explicitize_tier(tier, vertex_info, limit_new);
    }

    template<typename VertexInfosType>
    void p_reduce_tier(std::vector<std::size_t>& t,
                       VertexInfosType& vertex_info)
    {
        auto new_end = std::remove_if(
            t.begin(), t.end(),
            [&] (std::size_t vi) { return vertex_info[vi].explicitized; });
        t.erase(new_end, t.end());
    }

    static LitOrVal p_negate(LitOrVal v) {
        return std::visit(
            overloaded{
                [] (bool& b) -> LitOrVal {
                    return !b;
                },
                [] (Lit& l) -> LitOrVal {
                    return -l;
                }
            }, v
        );
    }

    LitOrVal p_get_config_var(std::size_t config_index, SLit lit) {
        auto var = lit::var(lit);
        bool neg = lit::negative(lit);
        LitOrVal val = m_config_vars[config_index][var];
        return neg ? p_negate(val) : val;
    }

    void p_sat_add_formula_on(const std::vector<LitOrVal>& config_vars) {
        // no need to add unary clauses; already fixed.
        const ClauseDB& clauses = m_propagator.db();

        // add all binary clauses that are not satisfied by fixed literals
        for(auto pair : clauses.binary_clauses()) {
            SLit l1 = pair.first, l2 = pair.second;
            auto lv1 = lit::var(l1), lv2 = lit::var(l2);
            bool neg1 = lit::negative(l1), neg2 = lit::negative(l2);
            auto val1 = config_vars[lv1], val2 = config_vars[lv2];
            // if one of the literals is fixed, we can ignore this clause:
            // if the fixed literal is false, UP has also fixed the other;
            // if it is true, the clause is satisfied.
            if(std::holds_alternative<bool>(val1) || 
               std::holds_alternative<bool>(val2)) 
            {
                continue;
            }
            Lit solver_lit1 = *std::get_if<Lit>(&val1);
            Lit solver_lit2 = *std::get_if<Lit>(&val2);
            if(neg1) solver_lit1 = -solver_lit1;
            if(neg2) solver_lit2 = -solver_lit2;
            m_solver.add_short_clause(solver_lit1, solver_lit2);
        }
        // add all longer clauses that are not satisfied by fixed literals,
        // shortening them by removing fixed dissatisfied literals
        std::vector<Lit> cbuffer;
        for(CRef cl = 1, cn = clauses.literal_db_size(); cl < cn; 
            cl = clauses.next_clause(cl))
        {
            cbuffer.clear();
            bool satisfied = false;
            for(SLit l : clauses.lits_of(cl)) {
                auto lv = lit::var(l);
                auto val = config_vars[lv];
                bool neg = lit::negative(l);
                if(std::holds_alternative<bool>(val)) {
                    bool bv = *std::get_if<bool>(&val);
                    if(bv != neg) {
                        satisfied = true;
                        break;
                    }
                } else {
                    Lit solver_lit = *std::get_if<Lit>(&val);
                    if(neg) solver_lit = -solver_lit;
                    cbuffer.push_back(solver_lit);
                }
            }
            if(satisfied) {
                continue;
            }
            if(cbuffer.empty()) {
                throw std::logic_error("Empty clause in SAT solver!");
            }
            m_solver.add_clause(cbuffer.begin(), cbuffer.end());
        }
    }

    void p_sat_setup_config_vars() {
        assert(m_config_vars.empty());
        assert(m_symmetry_breaking_levels.size() == m_num_configs_allowed);
        m_config_vars.reserve(m_num_configs_allowed);
        DynamicBitset fixed(m_n_all, false);
        for(std::size_t ci : range(m_num_configs_allowed)) {
            fixed.reset();
            std::vector<LitOrVal> config_vars;
            config_vars.reserve(m_n_all);
            const auto& prop = m_color_classes[ci];
            std::uint32_t sym_level = m_symmetry_breaking_levels[ci];
            for(SLit l : IteratorRange{prop.get_trail().begin(), 
                                       prop.level_end(sym_level)}) 
            {
                auto lv = lit::var(l);
                fixed[lv].set();
            }
            for(std::size_t v = 0; v < m_n_all; ++v) {
                if(fixed[v]) {
                    config_vars.emplace_back(
                        std::in_place_type<bool>,
                        prop.is_true(lit::positive_lit(v))
                    );
                } else {
                    config_vars.emplace_back(
                        std::in_place_type<Lit>,
                        m_solver.new_var()
                    );
                }
            }
            p_sat_add_formula_on(config_vars);
            m_config_vars.emplace_back(std::move(config_vars));
        }
        m_symmetry_breaking_levels.clear();
        assert(m_config_vars.size() == m_num_configs_allowed);
    }

    /**
     * Add clauses to ensure that
     * each concrete literal that occurs in an interaction in the subuniverse
     * has at least one configuration in which it occurs.
     */
    bool p_sat_add_vertical() {
        HashSet<SLit> occurs;
        SLit prev_first = sammy::NIL;
        for(Vertex v : m_subproblem.uncovered_universe) {
            if(v.first != prev_first) {
                prev_first = v.first;
                occurs.insert(v.first);
            }
            occurs.insert(v.second);
        }
        std::vector<Lit> cbuffer;
        for(SLit l : occurs) {
            cbuffer.clear();
            bool satisfied = false;
            for(std::size_t ci : range(m_num_configs_allowed)) {
                LitOrVal cval = p_get_config_var(ci, l);
                std::visit(overloaded{
                    [&] (bool& b) { if(b) satisfied = true; },
                    [&] (Lit& l) { cbuffer.push_back(l); }
                }, cval);
                if(satisfied) break;
            }
            if(satisfied) continue;
            if(cbuffer.empty()) return false;
            m_solver.add_clause(cbuffer.begin(), cbuffer.end());
        }
        return true;
    }

    /**
     * Called to complete all configurations to full configurations
     * after successfully covering all interactions, either by
     * explicitizing them or by implicit coverage.
     */
    template<typename VertexInfosType>
    int p_continue_complete_all_classes(VertexInfosType& vertex_info) {
        if(m_config_vars.empty()) {
            // have not previously used the SAT solver
            p_sat_first_setup();
        }
        std::vector<Lit> assumptions;
        for(std::size_t class_index : range(m_num_configs_allowed)) {
            const auto& prop = m_color_classes[class_index];
            auto& config_vars = m_config_vars[class_index];
            for(SLit l : prop.get_trail()) {
                auto lvar = lit::var(l);
                bool neg = lit::negative(l);
                LitOrVal config_var = config_vars[lvar];
                std::visit(overloaded{
                    [&] (bool b) {
                        // if lvar is fixed to true and l is negated,
                        // or vice versa, that is an error.
                        if(b == neg) {
                            throw std::logic_error("Broken assumptions "
                                                   "in propagator!");
                        }
                    },
                    [&] (Lit l) {
                        if(neg) {
                            assumptions.push_back(-l);
                        } else {
                            assumptions.push_back(l);
                        }
                    }
                }, config_var);
            }
        }
        auto res = m_solver.solve(assumptions);
        if(!res) {
            // SAT solver was interrupted
            return -1;
        }
        if(!*res) {
            // we could not improve, but that's no proof due to our assumptions;
            // invoke the SAT solver to find a working way of covering all
            // interactions that are currently explicitized
            p_sat_extend_covered_in_class();
            return p_sat_solve(vertex_info);
        }
        // SAT solver found a solution
        p_sat_solution_to_final_solution();
        return 1;
    }

    /**
     * Extend m_covered_in_class with the newly
     * explicitized interactions.
     */
    bool p_sat_extend_covered_in_class() {
        std::size_t prev_explicitized = m_covered_in_class.size();
        std::size_t now_explicitized = m_explicitized.size();
        std::vector<std::size_t> needs_extra_var;
        for(std::size_t ex_index = prev_explicitized;
            ex_index < now_explicitized; ++ex_index) 
        {
            std::size_t vi = m_explicitized[ex_index];
            Vertex v = m_subproblem.uncovered_universe[vi];
            SLit first = v.first, second = v.second;
            std::vector<LitOrVal> covered_in_class;
            covered_in_class.reserve(m_num_configs_allowed);
            bool satisfied = false;
            needs_extra_var.clear();
            for(std::size_t ci : range(m_num_configs_allowed)) {
                LitOrVal cval1 = p_get_config_var(ci, first);
                LitOrVal cval2 = p_get_config_var(ci, second);
                bool fixed1 = std::holds_alternative<bool>(cval1);
                bool fixed2 = std::holds_alternative<bool>(cval2);
                if((fixed1 && !*std::get_if<bool>(&cval1)) ||
                   (fixed2 && !*std::get_if<bool>(&cval2)))
                {
                    covered_in_class.emplace_back(
                        std::in_place_type<bool>, false);
                    continue;
                }
                if(fixed1 && fixed2) {
                    covered_in_class.emplace_back(
                        std::in_place_type<bool>, true);
                    satisfied = true;
                    continue;
                }
                if(fixed1) {
                    covered_in_class.emplace_back(
                        std::in_place_type<Lit>, *std::get_if<Lit>(&cval2)
                    );
                    continue;
                }
                if(fixed2) {
                    covered_in_class.emplace_back(
                        std::in_place_type<Lit>, *std::get_if<Lit>(&cval1)
                    );
                    continue;
                }
                // if we don't find this satisfied anyways,
                // we need to add a variable; defer creation
                needs_extra_var.push_back(ci);
                covered_in_class.emplace_back(
                    std::in_place_type<bool>, false
                );
            }
            if(satisfied) {
                m_covered_in_class.push_back(std::move(covered_in_class));
                continue;
            }
            for(std::size_t ci : needs_extra_var) {
                Lit new_var = m_solver.new_var();
                covered_in_class[ci].template emplace<Lit>(new_var);
                LitOrVal cvl1 = p_get_config_var(ci, first);
                LitOrVal cvl2 = p_get_config_var(ci, second);
                Lit cl1 = *std::get_if<Lit>(&cvl1);
                Lit cl2 = *std::get_if<Lit>(&cvl2);
                m_solver.add_short_clause(-cl1, -cl2, new_var);
                m_solver.add_short_clause(-new_var, cl1);
                m_solver.add_short_clause(-new_var, cl2);
            }
            std::size_t clause_size = 0;
            for(LitOrVal cval : covered_in_class) {
                Lit* l = std::get_if<Lit>(&cval);
                if(l) {
                    m_solver.add_literal(*l);
                    ++clause_size;
                }
            }
            m_covered_in_class.push_back(std::move(covered_in_class));
            if(!clause_size) {
                return false;
            }
            m_solver.finish_clause();
        }
        return true;
    }

    /**
     * Setup the SAT model the first time it is needed.
     * Clears m_symmetry_breaking_levels, sets up
     * m_config_vars and m_sat_covered_in_class.
     */
    bool p_sat_first_setup() {
        assert(m_config_vars.empty());
        p_sat_setup_config_vars();
        if(!p_sat_add_vertical()) {
            return false;
        }
        return true;
    }


};

}

#endif
