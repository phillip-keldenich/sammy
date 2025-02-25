#ifndef SAMMY_BARRAGE_WORKER_LNS_H_INCLUDED_
#define SAMMY_BARRAGE_WORKER_LNS_H_INCLUDED_

#include "barrage.h"
#include "gurobi_clique_solver_g2.h"
#include "satdsatur_colorer.h"
#include "cmsat5_solver.h"
#include "implied_vertices.h"
#include "algorithm_ex.h"
#include "clique_sat_dsatur.h"
#include "lns_destroy.h"
#include "thread_interrupt.h"
#include "barrage_worker_with_core.h"


namespace sammy {

static constexpr double RANDOM_DESTRUCTION_PROB  = 0.3;
static constexpr double RANDOM_WITH_TABLE_PROB = 0.2;
static constexpr double RANDOM_AND_LEAST_UNIQUE_DESTRUCTION_PROB = 0.5;

template<typename SubproblemSolverType> class SubproblemLNSSolverCore {
  public:
    using SubproblemSolver = SubproblemSolverType;

    static std::unique_ptr<SubproblemLNSSolverCore> factory(
        PortfolioSolver* solver, 
        PortfolioElementWithCore<SubproblemLNSSolverCore>* element
    ) {
        return std::make_unique<SubproblemLNSSolverCore>(solver, element);
    }

    SubproblemLNSSolverCore(PortfolioSolver* solver, PortfolioElementWithCore<SubproblemLNSSolverCore>* element) :
        m_portfolio(solver),
        m_element_mutex(&element->get_mutex()),
        m_termination_flag(get_interrupt_flag_ptr()),
        m_interrupt_flag(false),
        m_local_recorder(element->get_mutable_recorder()),
        m_source("LNSCore<" + SubproblemSolverType::name() + ">"),
        m_clauses(&solver->get_clauses()),
        m_propagator(m_clauses),
        m_destroy(m_local_recorder, m_worker_counter++, &solver->implied_cache(), 
                  solver->get_best_solution(), solver->get_best_mes(), m_propagator),
        m_destroy_seen_solution(m_destroy.total_num_configs())
    {
        set_interrupt_flag_ptr(&m_interrupt_flag);
        if(m_termination_flag->load()) {
            throw InterruptError();
        }
    }

    /**
     * Called by the PortfolioElementWithCore
     * if it finds the termination flag to be set.
     */
    void termination_flag_set() {
        assert(m_termination_flag->load());
        p_interrupt();
    }

    void interrupt_if_necessary(const InterruptionCheckInfo& info) {
        if(m_termination_flag->load()) {
            p_interrupt();
        } else if(info.best_upper_bound < m_destroy_seen_solution) {
            m_destroy_seen_solution = info.best_upper_bound;
            p_interrupt();
        }
    }

    void main() {
        while(!m_termination_flag->load()) {
            auto phase_result = p_destroy_phase();
            if(phase_result == DestroyPhaseResult::GLOBAL_OPT) {
                break;
            } else if(phase_result == DestroyPhaseResult::INTERRUPTED) {
                continue;
            }
            LNSSubproblem subproblem = m_destroy.move_out_subproblem();
            if(m_portfolio->subproblem_reporting_enabled()) {
                m_portfolio->report_subproblem(
                    subproblem, m_destroy.get_partial(),
                    m_destroy.best_global_mes().size(),
                    m_portfolio->get_best_lower_bound(),
                    m_source.c_str()
                );
            }
            auto before = std::chrono::steady_clock::now();
            try {
                p_create_subsolver(std::move(subproblem));
            } catch(InterruptError&) {
                m_destroy.return_subproblem_on_abort(std::move(subproblem));
                continue;
            }
            auto res = m_solver->solve();
            double time_taken = seconds_between(before, std::chrono::steady_clock::now());
            if(!res) {
                m_destroy.return_subproblem_on_abort(m_solver->move_out_subproblem());
                continue;
            }
            LNSSubproblem tmp = m_solver->move_out_subproblem();
            if(!*res) {
                auto old_lb = m_portfolio->get_best_lower_bound();
                if(old_lb < tmp.removed_configurations.size()) {
                    m_portfolio->report_lower_bound(tmp.removed_configurations.size(),
                                                    tmp.uncovered_universe, m_source.c_str());
                }
                m_portfolio->lns_report_failure(tmp.removed_configurations.size(), time_taken);
                m_destroy.improvement_impossible(std::move(tmp), m_solver->mes_vertices());
            } else {
                m_portfolio->lns_report_success(tmp.removed_configurations.size(), time_taken);
                const auto& improvement = m_solver->get_solution();
                PartialSolution improved = m_destroy.improve_destroyed(improvement);
                if(m_portfolio->report_solution(improved, m_source.c_str())) {
                    m_destroy.return_subproblem_on_success(std::move(tmp), std::move(improved));
                    std::size_t new_size = m_destroy.total_num_configs();
                    {
                        std::unique_lock l{*m_element_mutex};
                        if(new_size < m_destroy_seen_solution) {
                            m_destroy_seen_solution = new_size;
                        }
                    }
                } else {
                    m_destroy.return_subproblem_on_abort(std::move(tmp));
                }
            }
        }
        m_local_recorder->store_event("LNS_WORKER_EXITING", {{"worker_id", m_destroy.worker_index()}}, "worker_id");
    }

  private:
    enum class DestroyPhaseResult {
        GLOBAL_OPT,
        INTERRUPTED,
        SUBPROBLEM_FOUND
    };

    DestroyPhaseResult p_destroy_phase() {
        {
            std::unique_lock l{*m_element_mutex};
            if(m_solver) {
                m_solver.reset();
            }
        }
        try {
            p_update_destroy_if_needed();
            if(get_and_clear_interrupt_flag() || m_termination_flag->load()) {
                return DestroyPhaseResult::INTERRUPTED;
            }
            std::size_t num_removed = m_destroy.destroy(m_portfolio->lns_select_removed_class_count());
            if(num_removed == 0) {
                m_portfolio->report_lower_bound(
                    m_destroy.total_num_configs(), 
                    m_portfolio->get_infeasibility_map().collect_vertices(), 
                    "LNS whole solution destruction"
                );
                return DestroyPhaseResult::GLOBAL_OPT;
            }
            return DestroyPhaseResult::SUBPROBLEM_FOUND;
        } catch(InterruptError&) {
            return DestroyPhaseResult::INTERRUPTED;
        }
    }

    void p_update_destroy_if_needed() {
        if(p_mes_update_looks_necessary()) {
            m_destroy.update_global_mes_if_better(m_portfolio->get_best_mes());
        }
        if(!p_update_looks_necessary()) return;
        PartialSolution new_full = m_portfolio->get_best_solution();
        {
            std::unique_lock l{*m_element_mutex};
            m_destroy_seen_solution = new_full.size();
        }
        m_destroy.update_full_solution_if_better(std::move(new_full));
    }

    bool p_mes_update_looks_necessary() const {
        return m_portfolio->get_best_mes_size() > m_destroy.best_global_mes().size();
    }

    bool p_update_looks_necessary() const {
        std::unique_lock l{*m_element_mutex};
        return m_portfolio->get_best_solution_size() < m_destroy_seen_solution;
    }

    void p_interrupt() {
        m_interrupt_flag.store(true);
        if(m_solver) {
            m_solver->abort();
        }
    }

    void p_create_subsolver(LNSSubproblem&& subproblem) {
        SharedDBPropagator prop{m_clauses};
        // relies on the subsolver guarantee that, 
        // on interrupt exceptions, subproblem remains usable
        auto subsolver = std::make_unique<SubproblemSolver>(
            m_portfolio,
            std::move(subproblem),
            std::move(prop),
            m_local_recorder,
            m_destroy.worker_index()
        );
        std::unique_lock l{*m_element_mutex};
        m_solver = std::move(subsolver);
    }

    static std::atomic<std::size_t> m_worker_counter;
    PortfolioSolver* m_portfolio;
    std::mutex *m_element_mutex;
    std::atomic<bool> *m_termination_flag;
    std::atomic<bool> m_interrupt_flag;
    EventRecorder *m_local_recorder;
    std::string m_source;
    ClauseDB* m_clauses;
    SharedDBPropagator m_propagator;
    LNSDestroy m_destroy;
    std::size_t m_destroy_seen_solution;
    std::unique_ptr<SubproblemSolver> m_solver;
};

template<typename S> std::atomic<std::size_t> SubproblemLNSSolverCore<S>::m_worker_counter{0};

/**
 * Portfolio worker that does LNS using 
 * the CliqueSatDSaturSolver class template.
 * At any point in time, this is in one of
 * two stages: either we are looking for a
 * subproblem that we can hopefully solve
 * to improve the current solution (m_solver == std::nullopt),
 * or we are currently solving such a subproblem
 * (m_solver != std::nullopt).
 */
template<typename IncrementalSolverType>
class CliqueSatDSaturLNSElement : public PortfolioElement {
  public:
    using IncrementalSolver = IncrementalSolverType;
    using CSDSSolver = CliqueSatDSaturSolver<IncrementalSolver>;
    
    CliqueSatDSaturLNSElement(PortfolioSolver* solver, std::size_t upper_bound, 
                              std::size_t worker_index) :
        PortfolioElement(solver),
        m_worker_index(worker_index),
        m_clauses(nullptr),
        m_last_seen_upper_bound(upper_bound),
        m_our_old_solution(upper_bound),
        m_current_subproblem(solver->get_best_solution()),
        m_should_interrupt(false)
    {}

    void interrupt_if_necessary(const InterruptionCheckInfo& info) override {
        if(should_terminate.load()) {
            p_interrupt();
        }
        if(events & static_cast<EventMask>(PortfolioEvent::BETTER_UPPER_BOUND)) {
            if(info.best_upper_bound < m_last_seen_upper_bound) {
                m_last_seen_upper_bound = info.best_upper_bound;
                if(m_last_seen_upper_bound < m_our_old_solution) {
                    p_interrupt();
                }
            }
        }
        events = 0;
    }

    void main() override {
        p_init();
        set_interrupt_flag_ptr(&m_should_interrupt);
        m_last_start_time = Clock::now();
        bool begun_selection = false;
        while(!should_terminate.load()) {
            if(!begun_selection) m_local_recorder.store_event("BEGIN_SUBPROBLEM_SELECTION");
            else m_local_recorder.store_event("RETRY_SUBPROBLEM_SELECTION");
            begun_selection = true;
            p_switch_to_subproblem_selection_stage();
            if(!p_select_next_subproblem()) continue;
            if(!p_create_solver()) continue;
            p_switch_to_solving_stage();
            p_solve_cover();
            begun_selection = false;
        }
        m_local_recorder.store_event("LNS_WORKER_EXITING", 
                                     {{"description", get_description()}}, "description");
    }

    const EventRecorder* get_recorder() const override {
        return &m_local_recorder;
    }

    void synchronize_recorder(const EventRecorder& other) override {
        m_local_recorder.synchronize_with(other);
    }

    void set_recorder_quiet(bool quiet) override {
        m_local_recorder.set_print_events(!quiet);
    }

    std::string get_description() const override {
        return "C&P/SATDSatur LNS Worker #" + std::to_string(m_worker_index);
    }

  private:
    void p_init() {
        m_clauses = &solver->get_clauses();
        m_empty_propagator.emplace(m_clauses);
        m_fast_clique_builder.emplace(*m_empty_propagator);
    }

    /**
     * Select the next subproblem we try to solve.
     * If we found a subproblem, we also have computed
     * the filtered vertex set and a heuristic clique.
     * 
     * @return true if we found a subproblem, false if 
     *         we were told to terminate/interrupt or
     *         did not find a suitable subproblem.
     */
    bool p_select_next_subproblem() {
        m_current_subproblem = solver->get_best_solution();
        if(!p_updated_current_subproblem()) return false;
        if(!p_destroy()) return false;
        m_local_recorder.store_event("SELECTED_NEXT_SUBPROBLEM");
        return true;
    }

    /**
     * Destroy the current solution to find the next subproblem
     * we try to solve to improve the solution.
     * Must set the filtered vertex set and heuristic clique.
     * 
     * @return true if we should continue with this subproblem,
     *         false if we should retry.
     */
    bool p_destroy() {
        std::size_t goal_removed = solver->lns_select_removed_class_count();
        if(goal_removed < 3) {
            goal_removed = 3;
        }
        if(goal_removed > m_current_subproblem.size()) {
            m_local_recorder.store_event("DESTROY_WHOLE_SOLUTION");
            goal_removed = m_current_subproblem.size();
            return p_destroy_random(goal_removed);
        }
        double destroy_method_sel = std::uniform_real_distribution<double>{0.0,1.0}(sammy::rng());
        if(destroy_method_sel < RANDOM_DESTRUCTION_PROB) {
            m_local_recorder.store_event("DESTROY_RANDOM", {{"goal_removed", goal_removed}}, "goal_removed");
            return p_destroy_random(goal_removed);
        } else if(destroy_method_sel < RANDOM_DESTRUCTION_PROB + RANDOM_WITH_TABLE_PROB) {
            p_update_cached_mes_and_table();
            m_local_recorder.store_event("DESTROY_RANDOM_USING_BEST_MES", 
                                         {{"goal_removed", goal_removed}}, "goal_removed");
            return p_destroy_random_with_table(goal_removed);
        } else {
            m_local_recorder.store_event("DESTROY_RANDOM_AND_LEAST_UNIQUE", {{"goal_removed", goal_removed}}, "goal_removed");
            return p_destroy_random_and_least_unique(goal_removed);
        }
    }

    /**
     * Update the cached MES and the table
     * of classes covering MES vertices.
     */
    void p_update_cached_mes_and_table() {
        m_cached_best_mes = solver->get_best_mes();
        m_class_index_to_mes_index.clear();
        m_class_index_to_mes_index.resize(m_current_subproblem.size());
        m_mes_vertices_nonuniquely_covered.assign(m_cached_best_mes.size(), false);
        for(auto& v : m_mes_index_to_class_indices) v.clear();
        m_mes_index_to_class_indices.resize(m_cached_best_mes.size());
        std::size_t idx = 0;
        for(Vertex v : m_cached_best_mes) {
            m_current_subproblem.find_covering_classes(v.first, v.second, [&] (Index cls) {
                m_class_index_to_mes_index[cls] = idx;
                m_mes_index_to_class_indices[idx].push_back(cls);
            });
            if(m_mes_index_to_class_indices[idx].size() != 1) {
                m_mes_vertices_nonuniquely_covered[idx].set();
            }
            ++idx;
        }
    }

    /**
     * Check whether there is some hope for finding an improvement
     * by removing the classes with the given indices, or if
     * each of the indices covers a unique vertex in the MES.
     */
    bool p_index_set_could_improve(const std::vector<Index>& indices) {
        DynamicBitset mes_vertices_covered = m_mes_vertices_nonuniquely_covered;
        m_initial_clique.clear();
        for(Index i : indices) {
            auto v = m_class_index_to_mes_index[i];
            if(v) {
                if(!mes_vertices_covered[*v]) {
                    mes_vertices_covered[*v].set();
                    m_initial_clique.push_back(m_cached_best_mes[*v]);
                }
            }
        }
        if(m_initial_clique.size() == indices.size()) return false;
        return true;
    }

    /**
     * Remove goal_removed randomly selected classes.
     * Repeat the removal process up to max_iterations times
     * if the removal process as long as the removal process
     * results in a 'hopeless' subproblem according to the best MES.
     */
    bool p_remove_random_if_improvement_possible(std::size_t goal_removed, std::size_t max_iterations) {
        p_update_cached_mes_and_table();
        if(p_interrupt_check()) return false;
        auto& rng = sammy::rng();
        for(std::size_t iter = 0; iter < max_iterations; ++iter) {
            std::vector<Index> indices = vector(range(Index(m_current_subproblem.size())));
            std::shuffle(indices.begin(), indices.end(), rng);
            indices.resize(goal_removed);
            std::sort(indices.begin(), indices.end());
            if(p_index_set_could_improve(indices)) {
                m_currently_removed = m_current_subproblem.remove_assignments(indices);
                return true;
            }
        }
        return false;
    }

    /**
     * Destroy the current solution, making use of the
     * table of classes covering MES vertices to find
     * a deletion that may cause an improvement.
     */
    bool p_destroy_random_with_table(std::size_t goal_removed) {
        std::vector<Index> classes_without_mes_vertex;
        std::vector<std::size_t> multiply_covered_mes_vertices;
        for(Index i : range(Index(m_current_subproblem.size()))) {
            if(!m_class_index_to_mes_index[i]) {
                classes_without_mes_vertex.push_back(i);
            }
        }
        for(std::size_t vi : m_mes_vertices_nonuniquely_covered.ones()) {
            multiply_covered_mes_vertices.push_back(vi);
        }
        if(classes_without_mes_vertex.empty() && multiply_covered_mes_vertices.empty()) {
            // we actually have optimality here; may reach this depending on race condition
            should_terminate.store(true);
            return false;
        }
        auto& rng = sammy::rng();
        std::uniform_real_distribution<double> sel(0.0, 1.0);
        std::vector<Index> selected_classes;
        if(!classes_without_mes_vertex.empty() && (multiply_covered_mes_vertices.empty() || sel(rng) < 0.5)) {
            std::size_t select = (std::min)(goal_removed / 2 + 1, classes_without_mes_vertex.size());
            std::shuffle(classes_without_mes_vertex.begin(), classes_without_mes_vertex.end(), rng);
            selected_classes.assign(classes_without_mes_vertex.begin(), classes_without_mes_vertex.begin() + select);
        } else {
            std::size_t srnd = std::uniform_int_distribution<std::size_t>(0, multiply_covered_mes_vertices.size() - 1)(rng);
            std::size_t mvi = multiply_covered_mes_vertices[srnd];
            std::vector<Index> classes = m_mes_index_to_class_indices[mvi];
            selected_classes = m_mes_index_to_class_indices[mvi];
            if(selected_classes.size() > goal_removed) {
                std::shuffle(selected_classes.begin(), selected_classes.end(), rng);
                selected_classes.resize(goal_removed);
            }
        }
        m_initial_clique.clear();
        std::sort(selected_classes.begin(), selected_classes.end());
        auto is_in = [&] (Index idx) {
            return std::binary_search(selected_classes.begin(), selected_classes.end(), idx);
        };
        std::vector<Index> additional = vector(range(Index(m_current_subproblem.size())));
        additional.erase(std::remove_if(additional.begin(), additional.end(), is_in), additional.end());
        std::shuffle(additional.begin(), additional.end(), rng);
        additional.resize(goal_removed - selected_classes.size());
        auto inserted = selected_classes.insert(selected_classes.end(), additional.begin(), additional.end());
        std::sort(inserted, selected_classes.end());
        std::inplace_merge(selected_classes.begin(), inserted, selected_classes.end());
        m_currently_removed = m_current_subproblem.remove_assignments(selected_classes);
        return p_post_destroy();
    }

    /**
     * Destroy m_current_solution by removing 
     * goal_removed randomly selected elements.
     */
    bool p_destroy_random(std::size_t goal_removed) {
        bool res = p_remove_random_if_improvement_possible(goal_removed, 10);
        if(!res) {
            m_local_recorder.store_event("DESTROY_RANDOM_FAILED", {{"goal_removed", goal_removed}}, "goal_removed");
            if(p_interrupt_check()) return false;
            return p_destroy_random_with_table(goal_removed);
        }
        return p_post_destroy();
    }

    /**
     * Destroy m_current_solution by removing
     * roughly 1/3 * goal_removed randomly selected configurations
     * followed by roughly 2/3 * goal_removed configurations
     * that cover the least amount of uniquely covered elements.
     */
    bool p_destroy_random_and_least_unique(std::size_t goal_removed) {
        std::size_t goal_removed_random = (std::min)(goal_removed, goal_removed / 3 + 1);
        if(!p_remove_random_if_improvement_possible(goal_removed_random, 10)) {
            if(p_interrupt_check()) return false;
            return p_destroy_random_with_table(goal_removed);
        }
        std::size_t goal_removed_least_unique = goal_removed - goal_removed_random;
        std::vector<std::size_t> unique_per_class(m_current_subproblem.size(), 0);
        m_current_subproblem.iterate_all_uniquely_covered_with_class(
            [&] (Index i, Lit, Lit) {
                unique_per_class[i] += 1;
            });
        std::vector<Index> indices = vector(range(Index(m_current_subproblem.size())));
        std::nth_element(
            indices.begin(), indices.begin() + goal_removed_least_unique, indices.end(),
            [&] (Index a, Index b) { 
                return unique_per_class[a] < unique_per_class[b]; });
        indices.resize(goal_removed_least_unique);
        std::sort(indices.begin(), indices.end());
        auto further_removed = m_current_subproblem.remove_assignments(indices);
        for(auto& r : further_removed) m_currently_removed.emplace_back(std::move(r));
        return p_post_destroy();
    }

    /**
     * After destruction of m_current_solution,
     * compute the filtered set of uncovered vertices
     * and an initial clique/MES.
     */
    bool p_post_destroy() {
        if(p_interrupt_check()) return false;
        m_local_recorder.store_event("BEGIN_FILTER_UNCOVERED");
        m_filtered_uncovered.clear();
        m_current_subproblem.iterate_all_uncovered([&] (Lit lmin, Lit lmax) {
            m_filtered_uncovered.emplace_back(lmin, lmax);
        });
        solver->implied_cache().remove_implied(m_filtered_uncovered, *m_empty_propagator);
        m_local_recorder.store_event("DONE_FILTER_UNCOVERED");
        if(m_filtered_uncovered.empty()) {
            solver->report_solution(m_current_subproblem, "empty uncovered set");
            return false;
        }
        m_local_recorder.store_event("BEGIN_INITIAL_MES");
        m_initial_clique = p_find_initial_clique();
        m_local_recorder.store_event("DONE_INITIAL_MES", {{"size", m_initial_clique.size()}}, "size");
        if(m_initial_clique.size() >= m_currently_removed.size()) {
            m_local_recorder.store_event("INITIAL_MES_PRECLUDES_IMPROVEMENT");
            solver->report_mes(m_initial_clique, "LNS initial MES heuristic");
            return false;
        }
        return !p_interrupt_check();
    }

    /**
     * Notify that we updated the current subproblem from
     * the best solution from the PortfolioSolver.
     * @return true if we should continue with this subproblem, 
     *         false if we should retry or terminate.
     */
    bool p_updated_current_subproblem() {
        std::size_t new_size = m_current_subproblem.size();
        std::unique_lock l{mutex};
        m_our_old_solution = new_size;
        if(m_last_seen_upper_bound < m_our_old_solution) {
            return false;
        }
        if(m_our_old_solution < m_last_seen_upper_bound) {
            m_last_seen_upper_bound = m_our_old_solution;
        }
        return !should_terminate.load();
    }

    /**
     * Find an initial clique on m_filtered_uncovered
     * by heuristic means.
     */
    std::vector<Vertex> p_find_initial_clique() {
        assert(m_empty_propagator->get_current_level() == 0);
        std::vector<Vertex> best_clique = m_fast_clique_builder->random_multistart_best_clique_known_valid(10, m_filtered_uncovered);
        if(m_initial_clique.size() > best_clique.size()) {
            std::vector<Vertex> missing = m_initial_clique;
            auto rem = std::remove_if(missing.begin(), missing.end(), [&] (Vertex v) {
                return std::binary_search(m_filtered_uncovered.begin(), m_filtered_uncovered.end(), v);
            });
            missing.erase(rem, missing.end());
            if(!missing.empty()) {
                std::sort(missing.begin(), missing.end());
                auto ins = m_filtered_uncovered.insert(m_filtered_uncovered.end(), missing.begin(), missing.end());
                std::inplace_merge(m_filtered_uncovered.begin(), ins, m_filtered_uncovered.end());
            }
            best_clique = m_initial_clique;
            best_clique = m_fast_clique_builder->compute_clique_known_valid(best_clique.begin(), 
                                                                            best_clique.end(), 
                                                                            m_filtered_uncovered);
        }
        if(p_walk_cached_cliques(best_clique)) {
            std::size_t old_size = best_clique.size();
            best_clique = m_fast_clique_builder->compute_clique_known_valid(best_clique.begin(), 
                                                                            best_clique.end(), 
                                                                            m_filtered_uncovered);
            if(best_clique.size() > old_size) {
                solver->add_clique(best_clique);
            }
        }
        return best_clique;
    }
    
    bool p_walk_cached_cliques(std::vector<Vertex>& best_clique) {
        auto cview = solver->clique_cache_view();
        std::vector<Vertex> current;
        auto used = cview.end();
        for(auto i = cview.begin(), e = cview.end(); i != e; ++i) {
            auto clique = *i;
            if(clique.size() <= best_clique.size()) continue;
            current.clear();
            std::copy_if(clique.begin(), clique.end(), std::back_inserter(current),
                         [&] (Vertex v) {
                            return std::binary_search(m_filtered_uncovered.begin(), 
                                                      m_filtered_uncovered.end(), v);});
            if(current.size() > best_clique.size()) {
                best_clique = current;
                used = i;
            }
        }
        if(used != cview.end()) {
            solver->clique_was_used(cview, used);
            return true;
        } else {
            solver->add_clique(best_clique);
            return false;
        }
    }

    /**
     * Create the solver.
     * @return true if we should continue, false if we were interrupted.
     */
    bool p_create_solver() {
        m_local_recorder.store_event("BEGIN_CREATE_SUBPROBLEM_SOLVER");
        try {
            std::size_t best_mes_size = solver->get_best_mes_size();
            std::size_t best_lower_bound = solver->get_best_lower_bound();
            solver->report_subproblem(m_filtered_uncovered, m_initial_clique, 
                                      m_current_subproblem, m_currently_removed,
                                      best_mes_size, best_lower_bound, "satdsatur");
            m_solver.emplace(m_filtered_uncovered, &solver->get_infeasibility_map(), *m_clauses,
                         m_initial_clique, best_mes_size, 
                         best_lower_bound, m_currently_removed);
            m_solver->set_event_recorder(&m_local_recorder);
        } catch(const InterruptError&) {
            m_local_recorder.store_event("ABORTED_CREATE_SUBPROBLEM_SOLVER");
            return false;
        }
        m_local_recorder.store_event("DONE_CREATE_SUBPROBLEM_SOLVER");
        if(get_and_clear_interrupt_flag()) return false;
        return !p_interrupt_check();
    }

    /**
     * Start the solver and wait for it to finish or be aborted.
     * Automatically reports improved solutions to the
     * PortfolioSolver instance controlling this worker.
     */
    void p_solve_cover() {
        auto state = m_solver->solve();
        if(state == CSDSSolver::SolveResult::ABORTED) {
            m_local_recorder.store_event("ABORTED_SUBPROBLEM_SOLVE");
            return;
        }
        m_local_recorder.store_event("DONE_SUBPROBLEM_SOLVE");
        if(m_solver->improved_global_bound()) {
            solver->report_lower_bound(m_solver->get_best_bound(), 
                                       m_solver->get_best_bound_subgraph(), 
                                       "C&P/SATDSatur LNS");
        }
        if(m_solver->improved_mes()) {
            solver->report_mes(m_solver->get_best_mes(), "C&P/SATDSatur LNS");
        }
        double time_taken = seconds_between(m_last_start_time, Clock::now());
        if(state == CSDSSolver::SolveResult::SOLUTION_WAS_OPTIMAL) {
            solver->lns_report_failure(m_currently_removed.size(), time_taken);
        } else {
            PartialSolution subgraph_solution = m_solver->get_partial_solution();
            for(const auto& assignment : subgraph_solution.assignments()) {
                m_current_subproblem.add_assignment(assignment);
            }
            solver->report_solution(m_current_subproblem, "C&P/SATDSatur LNS");
            solver->lns_report_success(m_currently_removed.size(), time_taken);
        }
    }

    /**
     * Check whether we should interrupt what we are doing
     * (called mostly during subproblem selection stage).
     * This can be due to a better solution, optimality, 
     * or some termination criterion.
     * 
     * @return true if we should interrupt, false otherwise.
     */
    bool p_interrupt_check() {
        std::unique_lock l{mutex};
        if(should_terminate.load() || m_our_old_solution > m_last_seen_upper_bound) {
            return true;
        }
        return false;
    }

    /**
     * Cause the solver to terminate 
     * as soon as possible.
     * Called with lock held.
     */
    void p_interrupt() {
        if(m_have_solver) {
            m_solver->abort();
        } else {
            m_should_interrupt.store(true);
        }
    }

    /**
     * Switch to the subproblem selection stage.
     * This means taking the lock, resetting m_have_solver,
     * and destroying m_solver.
     */
    void p_switch_to_subproblem_selection_stage() {
        std::unique_lock l{mutex};
        m_have_solver = false;
        m_solver.reset();
    }

    /**
     * Switch to the subproblem solving stage.
     * This means taking the lock and setting m_have_solver.
     * m_solver MUST BE non-empty before calling this.
     */
    void p_switch_to_solving_stage() {
        std::unique_lock l{mutex};
        m_local_recorder.store_event("BEGIN_SUBPROBLEM_SOLVE", 
                                     {{"removed_configurations", m_currently_removed.size()}, 
                                      {"old_solution_value", m_our_old_solution},
                                      {"reduced_uncovered", m_solver->get_graph_size()}},
                                     "removed_configurations", "reduced_uncovered");
        assert(!!m_solver);
        m_have_solver = true;
        m_last_start_time = Clock::now();
    }

    /**
     * The stage indicator variable.
     * Invariant: m_solver is != std::nullopt if m_have_solver is true.
     * The worker thread will always populate
     * m_solver before setting it to true.
     * This variable is only changed with the mutex held.
     */
    bool m_have_solver{false};

    /**
     * The solver, if we currently have it.
     */
    std::optional<CSDSSolver> m_solver;

    /**
     * Index of this worker (mostly for debugging/logging/time measurement purposes).
     */
    std::size_t m_worker_index;

    /**
     * Pointer to our clauses.
     */
    ClauseDB* m_clauses;

    /**
     * Empty propagator we can reuse.
     */
    std::optional<SharedDBPropagator> m_empty_propagator;

    /**
     * Clique heuristic.
     */
    std::optional<FastCliqueBuilder> m_fast_clique_builder;

    /**
     * Our local recorder.
     */
    EventRecorder m_local_recorder;

    /**
     * The last seen upper bound.
     */
    std::size_t m_last_seen_upper_bound;

    /**
     * The value of the solution we destroyed to obtain
     * our current repair subproblem.
     */
    std::size_t m_our_old_solution;

    /**
     * The last time we started looking for an improvement.
     */
    Clock::time_point m_last_start_time;

    /**
     * The partial solution we currently try to extend.
     */
    PartialSolution m_current_subproblem;

    /**
     * The set of configurations we removed to
     * obtain the current subproblem.
     */
    std::vector<DynamicBitset> m_currently_removed;

    /**
     * A filtered set of uncovered vertices
     * that guarantees all uncovered vertices are covered.
     */
    std::vector<Vertex> m_filtered_uncovered;

    /**
     * Initial clique found on m_filtered_uncovered.
     */
    std::vector<Vertex> m_initial_clique;

    /**
     * The best global MES at the time of destruction of
     * the solution we currently work on.
     */
    std::vector<Vertex> m_cached_best_mes;

    /**
     * Index of the vertex in m_cached_best_mes
     * that is covered by class with index i, if any.
     */
    std::vector<std::optional<std::size_t>> m_class_index_to_mes_index;

    /**
     * Classes covering vertex i in m_cached_best_mes.
     */
    std::vector<std::vector<Index>> m_mes_index_to_class_indices;

    /**
     * A bitset of MES vertices that are not uniquely covered.
     */
    DynamicBitset m_mes_vertices_nonuniquely_covered;

    /**
     * Interruption flag for the worker thread.
     */
    std::atomic<bool> m_should_interrupt;
};

}

#endif
