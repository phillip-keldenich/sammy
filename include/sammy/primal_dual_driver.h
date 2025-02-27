#ifndef SAMMY_PRIMAL_DUAL_DRIVER_H_INCLUDED_
#define SAMMY_PRIMAL_DUAL_DRIVER_H_INCLUDED_

#include "algorithm_ex.h"
#include "coloring.h"
#include "experiment_flags.h"
#include "fast_clique.h"
#include "gurobi_clique_solver_g2.h"
#include "initial_coloring_heuristic.h"
#include "output.h"
#include "pair_infeasibility_map.h"
#include "simplification.h"
#include "universe_subgraph.h"

namespace sammy {

static constexpr std::size_t MAX_CUTS_FROM_PRIMAL = 5;
static constexpr double MIN_VERTICES_PER_PRICING = 10.0;
static constexpr double MIN_RELATIVE_PER_PRICING = 0.01;
static constexpr std::size_t MAX_CUT_ROUNDS_PER_PRICING = 40;
static constexpr double PRIMAL_TO_DUAL_GOAL_RATIO = 0.1;

class PossiblySimplified {
  public:
    PossiblySimplified(const ClauseDB& original, Var n_concrete)
        : m_original(original), m_original_concrete(n_concrete) {}

    ClauseDB& formula() noexcept {
        return m_simplified ? m_simplified->formula : m_original;
    }

    const ClauseDB& formula() const noexcept {
        return m_simplified ? m_simplified->formula : m_original;
    }

    Var n_concrete() const noexcept {
        return m_simplified ? m_simplified->num_concrete : m_original_concrete;
    }

    bool is_simplified() const noexcept {
        return static_cast<bool>(m_simplified);
    }

    void simplify(OutputObject& output) {
        m_simplification_stats.capture_before(m_original, m_original_concrete);
        m_simplifier.emplace(m_original, m_original_concrete);
        m_simplified.emplace(
            sammy::run_simplifier(*m_simplifier, m_simplification_stats));
        m_simplification_stats.capture_after(formula(), n_concrete());
        export_simplification(output, *m_simplifier, *m_simplified,
                              &statistics());
    }

    const SimplificationStats& statistics() const noexcept {
        return m_simplification_stats;
    }

    std::vector<std::vector<bool>>
    reconstruct_sample(const std::vector<std::vector<bool>>& s) const {
        if (!m_simplified)
            return s;
        return m_simplifier->reconstruct_sample(*m_simplified, s);
    }

    std::vector<Vertex> reconstruct_mes(const std::vector<Vertex>& mes) const {
        if (!m_simplified)
            return mes;
        return m_simplifier->reconstruct_lb(*m_simplified, mes);
    }

  private:
    ClauseDB m_original;
    Var m_original_concrete;
    std::optional<SimplifyDatastructure> m_simplifier;
    std::optional<SimplifiedInstance> m_simplified;
    SimplificationStats m_simplification_stats;
};

class PrimalDualDriver {
  public:
    PrimalDualDriver(EventRecorder& recorder, OutputObject& output,
                     const ExperimentFlagsConfig& config,
                     const ClauseDB& formula, Var n_concrete)
        : m_begin_time(Clock::now()), m_config(config), m_recorder(&recorder),
          m_output(&output), m_thread_pool(), m_instance(formula, n_concrete),
          m_inf_map((p_possibly_simplify(), m_instance.n_concrete())),
          m_fast_col_solver(&m_instance.formula(), m_instance.n_concrete(),
                            &m_thread_pool, &m_inf_map),
          m_fast_clq_solver(SharedDBPropagator(&m_instance.formula()),
                            &m_thread_pool),
          m_state(nullptr), m_last_primal_ended_at(m_begin_time) {
        m_fast_col_solver.set_quiet(true);
    }

    void run() {
        if (!m_state) {
            m_state = p_initial_state();
        }
        while (m_state) {
            m_state->main();
            auto next = m_state->leave();
            if (!next) {
                m_state.reset();
                break;
            } else {
                next->enter();
                m_state = std::move(next);
                if (m_state->is_final())
                    break;
            }
        }
    }

    bool mes_is_optimal() const noexcept {
        if (!m_best_mes || !m_best_sample)
            return false;
        if (m_best_mes->size() >= m_best_sample->size())
            return true;
        if (!m_state)
            return false;
        auto* last_state = m_state.get();
        if (dynamic_cast<LBOptimalState*>(last_state)) {
            return true;
        }
        return false;
    }

    std::size_t get_best_lb() const noexcept {
        std::size_t bound = 0;
        if (get_best_mes()) {
            bound = get_best_mes()->size();
        }
        if (m_best_sat_bound > bound) {
            bound = m_best_sat_bound;
        }
        return bound;
    }

    const std::optional<std::vector<Vertex>>& get_best_mes() const {
        return m_best_mes;
    }

    const std::optional<std::vector<std::vector<bool>>>&
    get_best_sample() const {
        return m_best_sample;
    }

    class Timeout : public std::exception {
        const char* what() const noexcept override { return "timeout"; }
    };

    bool update_best_sample(const std::vector<std::vector<bool>>& sample,
                            const char* source) {
        if (!m_best_sample || sample.size() < m_best_sample->size()) {
            m_recorder->store_event(
                "IMPROVED_SAMPLE",
                {{"size", sample.size()}, {"source", source}}, "size",
                "source");
            m_best_sample = sample;
            return true;
        }
        return false;
    }

    bool update_best_sample(const std::vector<SharedDBPropagator>& propagators,
                            const char* source) {
        if (m_best_sample && m_best_sample->size() <= propagators.size())
            return false;
        std::vector<std::vector<bool>> s;
        s.reserve(propagators.size());
        std::transform(propagators.begin(), propagators.end(),
                       std::back_inserter(s),
                       [](const SharedDBPropagator& prop) {
                           return prop.extract_assignment();
                       });
        return update_best_sample(s, source);
    }

    bool update_best_mes(const std::vector<Vertex>& mes, const char* source) {
        if (!m_best_mes || mes.size() > m_best_mes->size()) {
            m_recorder->store_event("IMPROVED_MUTUALLY_EXCLUSIVE_SET",
                                    {{"size", mes.size()},
                                     {"source", source},
                                     {"simplified_vertices", mes}},
                                    "size", "source");
            m_best_mes = mes;
            return true;
        }
        return false;
    }

    template <typename InputIterator>
    void update_spawners(InputIterator begin, InputIterator end) {
        auto old_size = m_all_spawners_set.size();
        m_all_spawners_set.insert(begin, end);
        if (m_all_spawners_set.size() != old_size) {
            m_all_spawners.assign(m_all_spawners_set.begin(),
                                  m_all_spawners_set.end());
        }
    }

    template <typename Container>
    void update_spawners(const Container& container) {
        using std::begin;
        using std::end;
        update_spawners(begin(container), end(container));
    }

    ClauseDB& formula() noexcept { return m_instance.formula(); }

    const ClauseDB& formula() const noexcept { return m_instance.formula(); }

    const std::vector<Vertex>& valid_interactions() {
        if (!m_valid_interactions) {
            m_recorder->store_event("BEGIN_COLLECTING_VALID_INTERACTIONS");
            m_valid_interactions.emplace(m_inf_map.collect_vertices(1'000'000));
            m_recorder->store_event("END_COLLECTING_VALID_INTERACTIONS",
                                    {{"count", m_valid_interactions->size()}},
                                    "count");
        }
        return *m_valid_interactions;
    }

    double time_remaining() const {
        return m_config.lb_mip_config.total_mip_timeout -
               seconds_between(m_begin_time, Clock::now());
    }

    bool rerun_primal_with_dual_info() {
        Clock::time_point before_prun = Clock::now();
        auto& lbsolver = *m_cnp_clq_solver;
        auto& ubsolver = m_fast_col_solver;
        const auto& best_mes = *m_best_mes;
        auto& rec = *m_recorder;
        double fraction =
            std::uniform_real_distribution<double>{0.0, 1.0}(sammy::rng());
        std::size_t max_num_constraints =
            std::size_t(fraction * (best_mes.size() - 1));
        rec.store_event("BEGIN_PRIMAL_WITH_DUAL_INFORMATION",
                        {{"max_configs_from_dual", max_num_constraints}},
                        "max_configs_from_dual");
        ubsolver.reset_coloring();
        lbsolver.export_highest_weights_to_primal(ubsolver,
                                                  max_num_constraints);
        std::vector<Vertex> vinitial = lbsolver.get_fractional_mes_support();
        vinitial.insert(vinitial.end(), best_mes.begin(), best_mes.end());
        std::sort(vinitial.begin(), vinitial.end());
        vinitial.erase(std::unique(vinitial.begin(), vinitial.end()),
                       vinitial.end());
        ubsolver.color_lazy(vinitial);
        update_spawners(ubsolver.class_spawners());
        update_best_sample(ubsolver.all_classes(),
                           "dual-initialized heuristic");
        std::size_t cuts_extracted = lbsolver.cuts_from_primal(
            ubsolver.all_classes(), MAX_CUTS_FROM_PRIMAL);
        rec.store_event("END_PRIMAL_WITH_DUAL_INFORMATION",
                        {{"solution_found", ubsolver.all_classes().size()},
                         {"best_known_solution", m_best_sample->size()},
                         {"cuts_extracted", cuts_extracted}},
                        "solution_found", "best_known_solution",
                        "cuts_extracted");
        Clock::time_point after_prun = Clock::now();
        m_last_primal_time = seconds_between(before_prun, after_prun);
        m_last_primal_ended_at = after_prun;
        return cuts_extracted > 0;
    }

  private:
    ClauseDB& p_possibly_simplify() {
        if (!m_config.dont_simplify && !m_instance.is_simplified()) {
            m_instance.simplify(*m_output);
        }
        return m_instance.formula();
    }

    class DriverState {
      public:
        explicit DriverState(PrimalDualDriver* driver)
            : driver(driver), successor(nullptr) {}

        void enter() { p_enter(); }

        virtual ~DriverState() = default;
        virtual void main() = 0;

        virtual bool is_final() const { return false; }

        virtual std::unique_ptr<DriverState> leave() {
            p_leave();
            std::unique_ptr<DriverState> result = std::move(successor);
            return result;
        }

      protected:
        PrimalDualDriver* driver;
        std::unique_ptr<DriverState> successor;

        // called with old state in driver->m_state
        virtual void on_enter() {}
        // called with current state in driver->m_state
        virtual void on_leave() {}

      private:
        void p_enter() {
            p_check_time();
            on_enter();
        }

        void p_leave() { on_leave(); }

        void p_check_time() {
            if (!std::isfinite(
                    driver->m_config.lb_mip_config.total_mip_timeout))
            {
                return;
            }
            auto tnow = seconds_between(driver->m_begin_time, Clock::now());
            if (tnow >= driver->m_config.lb_mip_config.total_mip_timeout) {
                throw Timeout();
            }
        }
    };

    /**
     * The typical initial state, where we run
     * the initial heuristic for the first time.
     * This requires some extra actions compared
     * to follow-up runs of the heuristic.
     */
    class InitialState : public DriverState {
      public:
        using DriverState::DriverState;

        void main() override {
            auto& driver = *this->driver;
            auto& rec = *driver.m_recorder;
            auto& fcol = driver.m_fast_col_solver;
            auto& fclq = driver.m_fast_clq_solver;
            rec.store_event("BEGIN_FIRST_SOLVE");
            fcol.initialize_feasibilities();
            fcol.color_lazy();
            rec.store_event("DONE_FIRST_SOLVE");
            driver.update_best_sample(fcol.all_classes(),
                                      "initial coloring heuristic");
            rec.store_event("BEGIN_FEASIBILITY_EXTRACTION");
            fcol.extract_feasibilities();
            rec.store_event("DONE_FEASIBILITY_EXTRACTION");
            rec.store_event("BEGIN_LEARN_INFEASIBILITIES");
            learn_infeasibilities(driver.formula(), &driver.m_inf_map);
            rec.store_event("DONE_LEARN_INFEASIBILITIES");
            driver.update_spawners(fcol.class_spawners());
            rec.store_event("BEGIN_FIRST_LB");
            auto r = driver.m_config.initial_heuristic_config
                         .random_clique_restarts_per_iteration;
            auto best_mes =
                fclq.random_multistart_best_clique(r, driver.m_all_spawners);
            driver.update_best_mes(best_mes, "initial clique heuristic");
            rec.store_event("DONE_FIRST_LB");
            ++driver.m_initial_iterations;
        }

        std::unique_ptr<DriverState> leave() override {
            return std::make_unique<InitialHeuristicState>(driver);
        }
    };

    /**
     * State where we repeat runs of the initial heuristic.
     */
    class InitialHeuristicState : public DriverState {
      public:
        using DriverState::DriverState;

        void main() override {
            auto& driver = *this->driver;
            auto& rec = *driver.m_recorder;
            auto& fcol = driver.m_fast_col_solver;
            auto& fclq = driver.m_fast_clq_solver;
            ++driver.m_initial_iterations;
            rec.store_event("BEGIN_NEXT_HEURISTIC_SAMPLE",
                            {{"iteration", driver.m_initial_iterations}},
                            "iteration");
            fcol.reset_coloring();
            fcol.color_lazy(*driver.m_best_mes);
            driver.update_best_sample(fcol.all_classes(),
                                      "repeat coloring heuristic");
            driver.update_spawners(fcol.class_spawners());
            rec.store_event("DONE_NEXT_HEURISTIC_SAMPLE",
                            {{"iteration", driver.m_initial_iterations}},
                            "iteration");
            rec.store_event("BEGIN_NEXT_HEURISTIC_MES",
                            {{"iteration", driver.m_initial_iterations}},
                            "iteration");
            auto r = driver.m_config.initial_heuristic_config
                         .random_clique_restarts_per_iteration;
            driver.update_best_mes(fclq.random_multistart_best_clique(
                                       r / 2, driver.m_all_spawners),
                                   "repeat clique heuristic (all spawners)");
            driver.update_best_mes(fclq.random_multistart_best_clique(
                                       r / 2, fcol.class_spawners()),
                                   "repeat clique heuristic (new spawners)");
            rec.store_event("DONE_NEXT_HEURISTIC_MES",
                            {{"iteration", driver.m_initial_iterations}},
                            "iteration");
        }

        std::unique_ptr<DriverState> leave() override {
            const double time_needed =
                seconds_between(driver->m_begin_time, Clock::now());
            auto goal_iter =
                driver->m_config.initial_heuristic_config.goal_iterations;
            if (time_needed >=
                    driver->m_config.initial_heuristic_config.max_time ||
                (time_needed >=
                     driver->m_config.initial_heuristic_config.min_time &&
                 driver->m_initial_iterations >= goal_iter))
            {
                if (driver->m_config.sat_coloring_config.initial_sat_coloring) {
                    return std::make_unique<InitialSATColoringState>(driver);
                }
                return std::make_unique<PrimalDualInitializationState>(driver);
            }
            return std::make_unique<InitialHeuristicState>(driver);
        }
    };

    /**
     * If we want to initialize using SAT coloring,
     * we enter this state.
     */
    class InitialSATColoringState : public DriverState {
      public:
        using DriverState::DriverState;

        void main() override {
            tbegin = Clock::now();
            successor = std::make_unique<PrimalDualInitializationState>(driver);
            auto& rec = *driver->m_recorder;
            auto& ubsolver = driver->m_fast_col_solver;
            auto& formula = driver->formula();
            std::vector<Vertex> initial_vertices(*driver->get_best_mes());
            add_and_make_unique(initial_vertices,
                                ubsolver.class_spawners().begin(),
                                ubsolver.class_spawners().end());
            UniverseSubgraph subgraph{&formula, &driver->m_thread_pool,
                                      &driver->m_inf_map, initial_vertices};
            subgraph.extend_matrix_by_propagation();
            p_update_clique_indices(*driver->get_best_mes(), subgraph);
            extra_constraints.update_vertex_count(subgraph);
            GurobiCliqueSolverG2 g2_clique_solver{
                &subgraph,
                formula.num_vars(),
                &rec,
                *driver->get_best_sample(),
                *driver->get_best_mes(),
                driver->m_config.lb_mip_config};
            p_run_satcoloring(subgraph, g2_clique_solver);
        }

      private:
        void p_run_satcoloring(UniverseSubgraph& subgraph,
                               GurobiCliqueSolverG2& clique_solver) {
            std::size_t best_bound = driver->get_best_mes()->size();
            std::size_t best_soln = driver->get_best_sample()->size();
            auto& recorder = *driver->m_recorder;
            auto& ubsolver = driver->m_fast_col_solver;
            while (best_bound < best_soln) {
                double trem = driver->m_config.sat_coloring_config
                                  .initial_sat_coloring_timeout;
                trem -= seconds_between(tbegin, Clock::now());
                if (trem <= 0) {
                    // timeout for initial sat coloring;
                    // return (successor already set)
                    driver->m_best_sat_bound = best_bound;
                    return;
                }
                SATKColoringSolver ksolve{&subgraph, driver->m_recorder,
                                          cindices, best_bound};
                ksolve.set_extra_constraints(&extra_constraints);
                auto ksolve_res = ksolve.solve(trem);
                if (!ksolve_res) {
                    // timeout for initial sat coloring;
                    // return (successor already set)
                    driver->m_best_sat_bound = best_bound;
                    return;
                }
                if (!*ksolve_res) {
                    best_bound += 1;
                    if (p_optimize_on_subgraph(clique_solver)) {
                        driver->update_best_mes(
                            clique_solver.get_best_mes(),
                            "LP rounding during initial SAT coloring");
                        const auto& mes = *driver->get_best_mes();
                        p_update_clique_indices(mes, subgraph);
                        if (mes.size() > best_bound)
                            best_bound = mes.size();
                    }
                } else {
                    std::vector<SharedDBPropagator> classes =
                        extra_constraints.coloring_to_classes(
                            subgraph.get_propagator(), subgraph,
                            ksolve.get_coloring(), best_bound);
                    if (classes.empty()) {
                        recorder.store_event("ADDED_CONSTRAINTS");
                        continue; // new extra constraints!
                    }
                    ubsolver.reset_coloring();
                    for (SharedDBPropagator& c : classes) {
                        ubsolver.add_color_class(std::move(c));
                    }
                    ubsolver.color_lazy(driver->m_all_spawners);
                    if (driver->update_best_sample(
                            ubsolver.all_classes(),
                            "extension of initial SAT coloring"))
                    {
                        best_soln = driver->get_best_sample()->size();
                    }
                    if (best_soln != best_bound) {
                        const auto& new_spawners = ubsolver.class_spawners();
                        clique_solver.add_new_vertices(new_spawners);
                        extra_constraints.update_vertex_count(subgraph);
                    }
                }
            }
        }

        bool p_optimize_on_subgraph(GurobiCliqueSolverG2& clique_solver) {
            std::size_t clique_size_before = driver->get_best_mes()->size();
            for (;;) {
                double trem = driver->m_config.sat_coloring_config
                                  .initial_sat_coloring_timeout;
                trem -= seconds_between(tbegin, Clock::now());
                switch (clique_solver.solve_full_relaxation(trem)) {
                case SolverState::OPTIMUM_ON_SUBGRAPH:
                case SolverState::TIMEOUT_IMPROVEMENT:
                case SolverState::TIMEOUT_NO_IMPROVEMENT:
                    return clique_size_before <
                           clique_solver.get_best_mes().size();

                default:
                    if (clique_solver.greedy_add_to_cutting_planes())
                        continue;
                    if (clique_solver.greedy_generate_cutting_planes())
                        continue;
                    return clique_size_before <
                           clique_solver.get_best_mes().size();
                }
            }
        }

        void p_update_clique_indices(const std::vector<Vertex>& clique,
                                     const UniverseSubgraph& subgraph) {
            cindices.clear();
            std::transform(clique.begin(), clique.end(),
                           std::back_inserter(cindices),
                           [&](Vertex v) { return subgraph.vertex_index(v); });
        }

        Clock::time_point tbegin;
        std::vector<std::size_t> cindices;
        ColoringExtraConstraints extra_constraints;
    };

    /**
     * State in which there is some solution and lower bound
     * from the initial heuristic and we initialize the primal/dual
     * information, build the subgraph, etc.
     */
    class PrimalDualInitializationState : public DriverState {
      public:
        using DriverState::DriverState;

        void main() override {
            auto& rec = *driver->m_recorder;
            rec.store_event("BEGIN_INITIAL_SUBGRAPH_CREATION");
            p_create_subgraph();
            auto& sg = *driver->m_subgraph;
            auto& formula = driver->formula();
            rec.store_event("END_INITIAL_SUBGRAPH_CREATION", {{"n", sg.n()}},
                            "n");
            rec.store_event("BEGIN_INITIAL_MODEL_CREATION");
            driver->m_cnp_clq_solver.emplace(
                &sg, formula.num_vars(), driver->m_recorder,
                *driver->m_best_sample, *driver->m_best_mes,
                driver->m_config.lb_mip_config);
            rec.store_event("END_INITIAL_MODEL_CREATION");
        }

        std::unique_ptr<DriverState> leave() override {
            return std::make_unique<PrimalDualBaseState>(driver);
        }

      private:
        std::vector<std::size_t>
        p_extract_coloring(const std::vector<std::vector<bool>>& sample,
                           const std::vector<Vertex>& vertices) {
            const auto n = vertices.size(), nc = sample.size();
            std::vector<std::size_t> result(n, 0);
            for (std::size_t i = 0; i < n; ++i) {
                Vertex v = vertices[i];
                for (std::size_t k = 0; k < nc; ++k) {
                    const std::vector<bool>& config = sample[k];
                    if (config[v.first] && config[v.second]) {
                        result[i] = k;
                        break;
                    }
                }
            }
            return result;
        }

        void p_run_satcoloring() {
            std::vector<Vertex> verts =
                driver->m_fast_col_solver.class_spawners();
            verts.insert(verts.end(), driver->m_best_mes->begin(),
                         driver->m_best_mes->end());
            std::sort(verts.begin(), verts.end());
            verts.erase(std::unique(verts.begin(), verts.end()), verts.end());
            UniverseSubgraph subgraph(&driver->formula(),
                                      &driver->m_thread_pool,
                                      &driver->m_inf_map, verts);
            subgraph.extend_matrix_by_propagation();
            SATColoringSolver solver{
                &subgraph, driver->m_recorder, *driver->m_best_mes,
                p_extract_coloring(*driver->m_best_sample, verts)};
            solver.optimize_coloring();
            auto [classes, spawners] = solver.coloring_to_classes();
            driver->m_recorder->store_event("CLASSES_FROM_COLORING",
                                            {{"num_classes", classes.size()}},
                                            "num_classes");
            auto& fcol = driver->m_fast_col_solver;
            fcol.reset_coloring();
            for (std::size_t ci = 0, cn = classes.size(); ci < cn; ++ci) {
                fcol.add_color_class(std::move(classes[ci]), spawners[ci]);
            }
            classes.clear();
            spawners.clear();
            fcol.color_lazy(driver->m_all_spawners);
            driver->update_best_sample(fcol.all_classes(),
                                       "extended SAT coloring");
            driver->update_spawners(fcol.class_spawners());
            subgraph.add_vertices(fcol.class_spawners());
            SATColoringSolver solver2{
                &subgraph, driver->m_recorder, *driver->m_best_mes,
                p_extract_coloring(*driver->m_best_sample,
                                   subgraph.vertex_set())};
            solver2.optimize_coloring();
            std::tie(classes, spawners) = solver.coloring_to_classes();
            fcol.reset_coloring();
            for (std::size_t ci = 0, cn = classes.size(); ci < cn; ++ci) {
                fcol.add_color_class(std::move(classes[ci]), spawners[ci]);
            }
            fcol.color_lazy(driver->m_all_spawners);
            driver->m_recorder->store_event(
                "LAZY_COLORING_RESULT",
                {{"num_classes", fcol.all_classes().size()}}, "num_classes");
        }

        void p_create_subgraph() {
            auto& opt_sg = driver->m_subgraph;
            auto& formula = driver->formula();
            std::vector<Vertex> initial = *driver->m_best_mes;
            const auto& last_spawners =
                driver->m_fast_col_solver.class_spawners();
            const auto& all_spawners = driver->m_all_spawners;
            initial.insert(initial.end(), last_spawners.begin(),
                           last_spawners.end());
            std::sort(initial.begin(), initial.end());
            initial.erase(std::unique(initial.begin(), initial.end()),
                          initial.end());
            if (initial.size() < 2000 && all_spawners.size() < 20000) {
                initial.insert(initial.end(), all_spawners.begin(),
                               all_spawners.end());
                std::sort(initial.begin(), initial.end());
                initial.erase(std::unique(initial.begin(), initial.end()),
                              initial.end());
            }
            opt_sg.emplace(&formula, &driver->m_thread_pool, &driver->m_inf_map,
                           std::move(initial));
            opt_sg->extend_matrix_by_propagation();
        }
    };

    /**
     * The base state of the primal/dual approach.
     * We enter this state if the relaxation needs
     * to be solved (or re-solved).
     */
    class PrimalDualBaseState : public DriverState {
      public:
        using DriverState::DriverState;

        void main() override {
            auto& driver = *this->driver;
            auto& lbsolver = *driver.m_cnp_clq_solver;
            double trem = driver.time_remaining();
            SolverState sol_state = lbsolver.solve_full_relaxation(trem);
            lbsolver.count_compressible_constraints(50);
            ++driver.m_full_relaxations_solved;
            if (sol_state == SolverState::TIMEOUT_IMPROVEMENT) {
                driver.update_best_mes(lbsolver.get_best_mes(), "LP rounding");
            }
            driver.m_relaxation_history.push_back(lbsolver.get_last_value());
            switch (sol_state) {
            default:
                throw Timeout{};

            case SolverState::IMPROVEMENT_FOUND:
                driver.update_best_mes(lbsolver.get_best_mes(), "LP rounding");
                successor =
                    std::make_unique<RerunPrimalHeuristicState>(&driver);
                return;

            case SolverState::NO_IMPROVEMENT_FOUND:
                if (driver.p_check_primal_goal_ratio()) {
                    successor =
                        std::make_unique<RerunPrimalHeuristicState>(&driver);
                } else {
                    successor =
                        std::make_unique<NotOptimalOnSubgraphState>(&driver);
                }
                return;

            case SolverState::OPTIMUM_ON_SUBGRAPH:
                driver.m_relaxation_history.clear();
                driver.m_have_optimum_on_subgraph = true;
                if (driver.update_best_mes(lbsolver.get_best_mes(),
                                           "LP rounding") ||
                    driver.p_check_primal_goal_ratio())
                {
                    successor =
                        std::make_unique<RerunPrimalHeuristicState>(&driver);
                } else {
                    successor =
                        std::make_unique<OptimalOnSubgraphState>(&driver);
                }
                return;
            }
        }
    };

    class RerunPrimalHeuristicState : public DriverState {
      public:
        using DriverState::DriverState;

        void main() override {
            auto& driver = *this->driver;
            bool primal_added_constraints =
                driver.rerun_primal_with_dual_info();
            if (primal_added_constraints) {
                successor = std::make_unique<PrimalDualBaseState>(&driver);
            } else if (!driver.m_have_optimum_on_subgraph) {
                successor =
                    std::make_unique<NotOptimalOnSubgraphState>(&driver);
            } else {
                successor = std::make_unique<OptimalOnSubgraphState>(&driver);
            }
        }
    };

    class NotOptimalOnSubgraphState : public DriverState {
      public:
        using DriverState::DriverState;

        void main() override {
            auto& driver = *this->driver;
            if (driver.m_relaxation_history.size() >=
                MAX_CUT_ROUNDS_PER_PRICING)
            {
                p_price_suboptimal();
                driver.m_relaxation_history.clear();
                successor = std::make_unique<PrimalDualBaseState>(&driver);
                return;
            }
            auto& lbsolver = *driver.m_cnp_clq_solver;
            if (lbsolver.greedy_add_to_cutting_planes() ||
                lbsolver.greedy_generate_cutting_planes() ||
                driver.rerun_primal_with_dual_info())
            {
                successor = std::make_unique<PrimalDualBaseState>(&driver);
                return;
            }
            successor = std::make_unique<CutFailedState>(&driver);
        }

      private:
        void p_price_suboptimal() {
            if (p_price_suboptimal(driver->m_fast_col_solver.class_spawners()))
                return;
            if (p_price_suboptimal(driver->m_all_spawners))
                return;
            p_price_suboptimal(
                driver->m_fast_col_solver.extract_coloring_order());
        }

        bool p_price_suboptimal(const std::vector<Vertex>& vs) {
            return driver->m_cnp_clq_solver->price_vertices(vs.begin(),
                                                            vs.end()) > 0;
        }
    };

    /**
     * A state we enter when we could not find any
     * violated cutting plane despite not having the
     * optimal clique on our current subgraph.
     * Currently, pricing is the only option in this case;
     * we could also try an exact separation procedure first.
     */
    class CutFailedState : public DriverState {
      public:
        using DriverState::DriverState;

        void main() override {
            const auto& valid = driver->valid_interactions();
            if (!driver->m_cnp_clq_solver->price_vertices(valid.begin(),
                                                          valid.end()))
            {
                successor = std::make_unique<CutAndPriceFailedState>(driver);
            } else {
                driver->m_relaxation_history.clear();
                successor = std::make_unique<PrimalDualBaseState>(driver);
            }
        }
    };

    class CutAndPriceFailedState : public DriverState {
      public:
        using DriverState::DriverState;

        void main() override {
            auto& rec = *driver->m_recorder;
            rec.store_event(
                "CUT_AND_PRICE_FAILED",
                {{"best_mes", driver->m_best_mes->size()},
                 {"lp_bound", driver->m_cnp_clq_solver->get_last_value()}},
                "best_mes", "lp_bound");
            successor = std::make_unique<CutAndPriceFailedState>(driver);
        }

        bool is_final() const override { return true; }
    };

    class OptimalOnSubgraphState : public DriverState {
      public:
        using DriverState::DriverState;

        void main() override {
            auto& driver = *this->driver;
            auto& lbsolver = *driver.m_cnp_clq_solver;
            auto& ubsolver = driver.m_fast_col_solver;
            const auto& subgraph = *driver.m_subgraph;
            const auto& spawners = ubsolver.class_spawners();
            double percent_of_graph =
                (std::max)(MIN_VERTICES_PER_PRICING,
                           MIN_RELATIVE_PER_PRICING * subgraph.n());
            std::size_t pos_vertices =
                lbsolver.price_vertices(spawners.begin(), spawners.end());
            if (pos_vertices >= percent_of_graph) {
                p_successor_found_new();
                return;
            }
            const auto& all_spawners = driver.m_all_spawners;
            pos_vertices += lbsolver.price_vertices(all_spawners.begin(),
                                                    all_spawners.end());
            if (pos_vertices >= percent_of_graph) {
                p_successor_found_new();
                return;
            }
            auto coloring_order = ubsolver.extract_coloring_order();
            pos_vertices += lbsolver.price_vertices(coloring_order.begin(),
                                                    coloring_order.end());
            if (pos_vertices >= percent_of_graph) {
                p_successor_found_new();
                return;
            }
            const auto& valid = driver.valid_interactions();
            pos_vertices += lbsolver.price_vertices(valid.begin(), valid.end());
            if (pos_vertices) {
                p_successor_found_new();
                return;
            }
            successor = std::make_unique<LBOptimalState>(&driver);
        }

      private:
        void p_successor_found_new() {
            driver->m_relaxation_history.clear();
            driver->m_have_optimum_on_subgraph = false;
            successor = std::make_unique<PrimalDualBaseState>(driver);
        }
    };

    class LBOptimalState : public DriverState {
      public:
        using DriverState::DriverState;

        bool is_final() const override { return true; }
        void main() override {
            auto& rec = *driver->m_recorder;
            rec.store_event(
                "OPTIMUM_CLIQUE_ON_G2",
                {{"lb", driver->m_best_mes->size()},
                 {"ub", driver->m_best_sample->size()},
                 {"lp_bound", driver->m_cnp_clq_solver->get_last_value()}},
                "lb", "ub", "lp_bound");
            successor = std::make_unique<LBOptimalState>(driver);
        }
    };

    std::unique_ptr<DriverState> p_initial_state() {
        std::unique_ptr<DriverState> result;
        if (!m_best_sample || !m_best_mes) {
            result = std::make_unique<InitialState>(this);
        }
        // TODO: state on re-run?
        result->enter();
        return result;
    }

    bool p_check_primal_goal_ratio() const {
        auto now = Clock::now();
        double tlast = m_last_primal_time;
        double tsince = seconds_between(m_last_primal_ended_at, now);
        double ftest = 1.0 / PRIMAL_TO_DUAL_GOAL_RATIO - 1.0;
        return tsince >= tlast * ftest;
    }

    Clock::time_point m_begin_time;
    ExperimentFlagsConfig m_config;
    EventRecorder* m_recorder;
    OutputObject* m_output;
    ThreadGroup<void> m_thread_pool;
    PossiblySimplified m_instance;
    PairInfeasibilityMap m_inf_map;
    ColoringHeuristicSolver m_fast_col_solver;
    ParallelFastCliqueBuilder m_fast_clq_solver;
    std::unique_ptr<DriverState> m_state;
    EdgeSet m_all_spawners_set;
    std::vector<Vertex> m_all_spawners;
    std::optional<std::vector<std::vector<bool>>> m_best_sample;
    std::optional<std::vector<Vertex>> m_best_mes;
    std::optional<UniverseSubgraph> m_subgraph;
    std::optional<GurobiCliqueSolverG2> m_cnp_clq_solver;
    std::optional<std::vector<Vertex>> m_valid_interactions;
    std::deque<double> m_relaxation_history;
    std::size_t m_initial_iterations = 0;
    std::size_t m_full_relaxations_solved = 0;
    std::size_t m_best_sat_bound = 0;
    double m_last_primal_time = 0.0;
    Clock::time_point m_last_primal_ended_at;
    bool m_have_optimum_on_subgraph = false;
};

} // namespace sammy

#endif
