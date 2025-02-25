#ifndef SAMMY_BARRAGE_H_INCLUDED_
#define SAMMY_BARRAGE_H_INCLUDED_

#include "experiment_flags.h"
#include "pair_infeasibility_map.h"
#include "initial_coloring_heuristic.h"
#include "simplification.h"
#include "output.h"
#include "thread_clauses.h"
#include "clique_storage.h"
#include "implied_vertices.h"
#include "initial_phase_result.h"
#include "barrage_lns_subproblem.h"
#include <optional>
#include <utility>
#include <any>
#include <queue>
#include <filesystem>


namespace sammy {

using EventMask = std::uint32_t;

enum class PortfolioEvent : EventMask {
    TIMEOUT = 1,
    BETTER_LOWER_BOUND = 2,
    BETTER_UPPER_BOUND = 4,
    BETTER_MES = 8,
    ALARM = 16,
    OPTIMALITY = 32
};

struct InterruptionCheckInfo {
    std::size_t best_mes;
    std::size_t best_lower_bound;
    std::size_t best_upper_bound;
};

class PortfolioSolver;
struct LNSTimeAndSuccessInfo {
    struct RemovedClassesInfo {
        std::size_t complete_try_successes = 0;
        std::size_t complete_tries_since_last_success = 0;
        std::size_t complete_tries_total = 0;
        double complete_tries_total_time = 0.0;
    };

    inline LNSTimeAndSuccessInfo(PortfolioSolver* solver) noexcept;

    std::map<std::size_t, RemovedClassesInfo> removed_classes_info;
    double current_goal_time = 30.0;
    double num_failure_threshold = 20.0;
    std::size_t universe_size;
    std::size_t min_num_tries_for_time = 10;

    void report_success(std::size_t removed_classes, double time) {
        auto& info = removed_classes_info[removed_classes];
        info.complete_tries_since_last_success = 0;
        info.complete_tries_total++;
        info.complete_try_successes++;
        info.complete_tries_total_time += time;
    }

    void report_failure(std::size_t removed_classes, double time) {
        auto& info = removed_classes_info[removed_classes];
        info.complete_tries_since_last_success++;
        info.complete_tries_total++;
        info.complete_tries_total_time += time;
    }

    std::size_t select_goal_num_removed() {
        if(removed_classes_info.empty()) {
            return p_empty_select_num_removed();
        }
        auto last_usable = removed_classes_info.end();
        for(auto it = removed_classes_info.begin(), e = removed_classes_info.end(); it != e; ++it) {
            const auto& info = it->second;
            if(info.complete_tries_since_last_success >= num_failure_threshold) {
                continue;
            }
            if(info.complete_tries_total < min_num_tries_for_time) {
                last_usable = it;
                continue;
            }
            if(info.complete_tries_total_time / info.complete_tries_total >= current_goal_time) {
                last_usable = it;
            }
        }
        if(last_usable == removed_classes_info.end()) {
            return p_enlarge_num_removed();
        }
        return last_usable->first;
    }

  private:
    std::size_t p_enlarge_num_removed() {
        return removed_classes_info.rbegin()->first + 1;
    }

    std::size_t p_empty_select_num_removed() {
        if(universe_size < 100'000) return 9;
        if(universe_size < 500'000) return 7;
        if(universe_size < 1'000'000) return 5;
        if(universe_size < 5'000'000) return 4;
        return 3;
    }
};

class PortfolioElement {
  public:
    explicit PortfolioElement(PortfolioSolver* solver) :
        solver(solver)
    {}

    /**
     * Destroy the portfolio element.
     * Joins if the worker has not yet been joined on.
     */
    virtual ~PortfolioElement() {
        join();
    }

    /**
     * Wait for the worker thread to end.
     * This does not notify the worker at all,
     * so the worker must be notified beforehand.
     */
    void join() {
        if(m_worker.joinable()) {
            m_worker.join();
        }
    }

    /**
     * Get a thread-local copy of the clauses.
     */
    inline ClauseDB& get_clauses() const;

    /**
     * Get the infeasibility map.
     */
    inline PairInfeasibilityMap &get_infeasibility_map() noexcept;

    /*
     * Routine that is called by the event
     * delivery mechanism (coordinator thread)
     * to check if this worker should be interrupted,
     * and interrupt it if necessary.
     * Called with the mutex held.
     */
    virtual void interrupt_if_necessary(const InterruptionCheckInfo& /*info*/) {}

    /**
     * Get an event recorder (if this element has any).
     */
    virtual const EventRecorder* get_recorder() const { return nullptr; }

    /**
     * If this element has any, synchronize the event recorder.
     */
    virtual void synchronize_recorder(const EventRecorder& /*other*/) {}

    /**
     * If this element has any, set the quiet flag of the event recorder.
     */
    virtual void set_recorder_quiet(bool quiet) {}

    /**
     * Get a description of this element (it it has any).
     */
    virtual std::string get_description() const { return ""; }

    /**
     * Called by the coordinator to
     * deliver new events to this worker.
     */
    inline void deliver_events(EventMask new_events, const InterruptionCheckInfo& info);

    /**
     * Called by the coordinator to deliver alarm events.
     */
    inline void deliver_alarm(std::size_t alarm_id, const InterruptionCheckInfo& info);

    /**
     * Set an alarm (this causes the coordinator to notify us,
     * generating an ALARM event at the requested time).
     * Automatically discards the currently set alarm (if there is any).
     */
    inline void set_alarm(double in_seconds);

    /**
     * Discard the currently set alarm.
     */
    inline void discard_alarm();

    /**
     * Can be called to obtain the raised events
     * and reset them to their unraised state.
     */
    EventMask consume_events() {
        std::unique_lock l{mutex};
        EventMask result = events;
        events = 0;
        return result;
    }

    /**
     * Called to start the worker thread.
     */
    void start() {
        if(m_worker.joinable()) throw std::logic_error("Trying to start running worker!");
        m_worker = std::thread{[&] () {this->main();}};
    }

  protected:
    /**
     * Method implementing the actual work.
     */
    virtual void main() = 0;

    /**
     * Mutex that locks the 'events' field,
     * and any data that can be accessed by
     * multiple threads, e.g., data accessed
     * in the 'interrupt_if_necessary' method.
     */
    std::mutex mutex;

    /**
     * Condition variable that can be 
     * waited on for events.
     */
    std::condition_variable condition;

    /**
     * The coordinating solver.
     */
    PortfolioSolver* solver;

    /**
     * The events currently raised/pending.
     */
    EventMask events = 0;

    /**
     * An atomic flag that we can check to see if
     * we should terminate ASAP. This flag is
     * already set on delivery (not consumption)
     * of the TIMEOUT event.
     */
    std::atomic<bool> should_terminate{false};

  private:
    std::thread m_worker; //< the worker thread
    std::optional<std::size_t> m_alarm{std::nullopt};
    std::size_t m_alarm_id_counter{0};
};

class PortfolioSolver {
  public:
    PortfolioSolver(ClausesTicket clauses, EventRecorder* global_recorder,
                    InitialPhaseResult&& initial_phase) :
        m_inf_map(std::move(initial_phase.inf_map)),
        m_implied_cache(&m_inf_map, initial_phase.universe_size),
        m_best_spawners(std::move(initial_phase.best_spawners)),
        m_all_spawners(std::move(initial_phase.all_spawners)),
        m_coloring_order(std::move(initial_phase.coloring_order)),
        m_universe_size(initial_phase.universe_size),
        m_ticket(clauses),
        m_global_recorder(global_recorder),
        m_best_mes(std::move(initial_phase.best_mutually_exclusive)),
        m_best_lb_vertex_set(m_best_mes),
        m_lower_bound(m_best_mes.size()),
        m_best_solution(
            local_clauses(clauses).num_vars(),
            &m_inf_map,
            initial_phase.best_solution.begin(),
            initial_phase.best_solution.end()),
        m_lns_info(this)
    {}

    void limited_reduce_universe(double time_limit) {
        if(!m_implied_cache.have_reduced_universe()) {
            m_global_recorder->store_event("BEGIN_LIMITED_IMPLIED_VERTEX_ELIMINATION");
            m_implied_cache.limited_reduce_universe(local_clauses(m_ticket), time_limit);
            p_post_reduce_universe();
        }
    }

    void reduce_universe() {
        if(!m_implied_cache.have_reduced_universe()) {
            m_global_recorder->store_event("BEGIN_IMPLIED_VERTEX_ELIMINATION");
            m_implied_cache.reduce_universe(local_clauses(m_ticket));
            p_post_reduce_universe();
        }
    }

    const ImpliedVertexCache &implied_cache() const noexcept {
        return m_implied_cache;
    }

    PartialSolution get_best_solution() const {
        std::unique_lock l{m_mutex};
        PartialSolution partial = m_best_solution;
        return partial;
    }

    const std::vector<Vertex>& get_best_spawners() const noexcept {
        return m_best_spawners;
    }

    const std::vector<Vertex>& get_all_spawners() const noexcept {
        return m_all_spawners;
    }

    const std::vector<Vertex>& get_coloring_order() const noexcept {
        return m_coloring_order;
    }

    std::size_t get_universe_size() const noexcept {
        return m_universe_size;
    }

    void set_alarm(PortfolioElement* element, double in_seconds, std::size_t alarm_id) {
        std::unique_lock l{m_mutex};
        AlarmTime time = Clock::now() + std::chrono::duration<double>(in_seconds);
        m_alarm_requests.emplace_back(AlarmRequest{time, alarm_id, element});
        if(m_alarm_requests.size() == 1) {
            m_condition.notify_one();
        }
    }

    ClauseDB& get_clauses() const {
        return local_clauses(m_ticket);
    }

    PairInfeasibilityMap& get_infeasibility_map() noexcept {
        return m_inf_map;
    }

    void run(double time_limit = 1000.0 * 365.0 * 24.0 * 60.0 * 60.0) {
        set_alarm(nullptr, time_limit, 0);
        p_main_loop();
    }

    void merge_events() {
        std::size_t worker_id = 0;
        for(const auto& e : m_elements) {
            const auto* events = e->get_recorder();
            if(events) {
                std::string description = e->get_description();
                if(description.empty()) description = std::to_string(worker_id);
                m_global_recorder->merge_events(*events, description);
            }
            ++worker_id;
        }
    }

    void add_element(std::unique_ptr<PortfolioElement> element) {
        std::unique_lock l{m_mutex};
        m_waiting_elements.emplace_back(std::move(element));
        if(m_waiting_elements.size() == 1) {
            m_condition.notify_one();
        }
    }

    template<typename ConcreteType, typename... Args>
    ConcreteType& emplace_element(Args&&... args)
    {
        auto element = std::make_unique<ConcreteType>(std::forward<Args>(args)...);
        ConcreteType& ref = *element;
        add_element(std::move(element));
        return ref;
    }

    /**
     * @brief Record a new event in the global event recorder.
     */
    template<typename Tag, typename... EventArgs>
    void report_global(Tag&& tag, OutputObject o, EventArgs&&... p)
    {
        std::unique_lock l{m_mutex};
        m_global_recorder->store_event(
            std::forward<Tag>(tag), std::move(o), std::forward<EventArgs>(p)...);
    }

    /**
     * @brief Obtain a clique cache view.
     */
    CliqueStorage::StorageView clique_cache_view() const noexcept {
        return m_clique_cache.obtain_view();
    }

    /**
     * @brief Report a new mutually exclusive set.
     */
    void report_mes(const std::vector<Vertex>& vertices, const char* source) {
        std::unique_lock l{m_mutex};
        if(m_best_mes.size() >= vertices.size()) return;
        m_best_mes = vertices;
        if(m_implied_cache.have_reduced_universe()) {
            m_implied_cache.replace_implied(m_best_mes);
            std::size_t old_size = m_best_mes.size();
            std::sort(m_best_mes.begin(), m_best_mes.end());
            m_best_mes.erase(std::unique(m_best_mes.begin(), m_best_mes.end()),
                             m_best_mes.end());
            if(m_best_mes.size() != old_size) {
                throw std::logic_error("MES size changed during implied vertex elimination!");
            }
        }
        EventMask event = static_cast<EventMask>(PortfolioEvent::BETTER_MES);
        m_global_recorder->store_event("IMPROVED_MES", 
            {{"size", m_best_mes.size()}, {"vertices", m_best_mes}, {"source", source}}, "size", "source");
        if(m_best_mes.size() > m_lower_bound) {
            m_lower_bound = m_best_mes.size();
            m_best_lb_vertex_set = m_best_mes;
            m_global_recorder->store_event("IMPROVED_LB", 
                {{"lb", m_best_mes.size()}, {"vertices", m_best_mes}, {"source", source}}, "lb", "source");
            event |= static_cast<EventMask>(PortfolioEvent::BETTER_LOWER_BOUND);
            if(m_lower_bound >= m_best_solution.size()) {
                m_global_recorder->store_event("OPTIMALITY_REACHED");
                event |= static_cast<EventMask>(PortfolioEvent::OPTIMALITY);
            }
        }
        p_raise_events(event);
    }

    /**
     * @brief Called by workers to report a new solution.
     * Returns true if the reported solution is better than the
     * current best and was stored in the portfolio solver.
     */
    bool report_solution(const PartialSolution& solution, const char* source) {
        std::unique_lock l{m_mutex};
        if(m_best_solution.size() <= solution.size()) return false;
        m_best_solution = solution;
        m_global_recorder->store_event("IMPROVED_SOLUTION", {{"size", m_best_solution.size()},
                                                             {"source", source},
                                                             {"lb", m_lower_bound}}, "size", "lb", "source");
        EventMask events = static_cast<EventMask>(PortfolioEvent::BETTER_UPPER_BOUND);
        if(m_best_solution.size() <= m_lower_bound) {
            m_global_recorder->store_event("OPTIMALITY_REACHED");
            events |= static_cast<EventMask>(PortfolioEvent::OPTIMALITY);
        }
        p_raise_events(events);
        return true;
    }

    /**
     * @brief Called by workers to report new lower bounds.
     */
    void report_lower_bound(std::size_t lower_bound, const std::vector<Vertex>& subgraph, const char* source) {
        std::unique_lock l{m_mutex};
        if(lower_bound <= m_lower_bound) return;
        m_best_lb_vertex_set = subgraph;
        m_lower_bound = lower_bound;
        m_global_recorder->store_event("IMPROVED_LB",
                                       {{"lb", m_lower_bound}, {"source", source}, 
                                        {"vertices", m_best_lb_vertex_set},
                                        {"ub", m_best_solution.size()}}, "lb", "source", "ub");
        EventMask events = static_cast<EventMask>(PortfolioEvent::BETTER_LOWER_BOUND);
        if(m_best_solution.size() <= m_lower_bound) {
            m_global_recorder->store_event("OPTIMALITY_REACHED");
            events |= static_cast<EventMask>(PortfolioEvent::OPTIMALITY);
        }
        p_raise_events(events);
    }

    /**
     * Add a clique to the cache.
     */
    void add_clique(const std::vector<Vertex>& vertices) {
        m_clique_cache.push_clique(vertices.begin(), vertices.end());
    }

    /**
     * Add a clique to the cache.
     */
    template<typename Iterator>
    void add_clique(Iterator begin, Iterator end) {
        m_clique_cache.push_clique(begin, end);
    }

    /**
     * Mark a clique as used.
     */
    void clique_was_used(CliqueStorage::StorageView& view, 
                         CliqueStorage::StorageView::Iterator iter) 
    {
        m_clique_cache.used_clique(view, iter);
    }

    /**
     * Mark a clique as used.
     */
    void clique_was_used(CliqueStorage::StorageView& view, Index index) 
    {
        m_clique_cache.used_clique(view, index);
    }

    /**
     * `Get the best mutually exclusive set currently available.
     */
    std::vector<Vertex> get_best_mes() const {
        std::unique_lock l{m_mutex};
        std::vector<Vertex> result = m_best_mes;
        return result;
    }

    /**
     * Get the best mutually exclusive set size currently available.
     */
    std::size_t get_best_mes_size() const noexcept {
        std::unique_lock l{m_mutex};
        return m_best_mes.size();
    }

    /**
     * Get the best lower bound currently available.
     */
    std::size_t get_best_lower_bound() const noexcept {
        std::unique_lock l{m_mutex};
        return m_lower_bound;
    }

    std::vector<Vertex> get_best_lower_bound_vertex_set() const {
        std::unique_lock l{m_mutex};
        return m_best_lb_vertex_set;
    }

    /**
     * Get the best upper bound currently available.
     */
    std::size_t get_best_solution_size() const noexcept {
        std::unique_lock l{m_mutex};
        return m_best_solution.size();
    }

    std::size_t lns_select_removed_class_count() {
        std::unique_lock l{m_mutex};
        return m_lns_info.select_goal_num_removed();
    }

    void lns_report_failure(std::size_t num_removed, double time) {
        std::unique_lock l{m_mutex};
        m_lns_info.report_failure(num_removed, time);
    }

    void lns_report_success(std::size_t num_removed, double time) {
        std::unique_lock l{m_mutex};
        m_lns_info.report_success(num_removed, time);
    }

    void enable_subproblem_reporting(const std::filesystem::path& out_dir) {
        if(m_subproblem_dir) {
            throw std::logic_error("Subproblem reporting already enabled!");
        }
        if(!exists(out_dir)) {
            create_directory(out_dir);
        }
        m_subproblem_dir = out_dir;
        auto universe_and_clauses_file = out_dir / "universe_and_clauses.json";
        nlohmann::json universe_and_clauses;
        auto universe = lit::externalize(m_implied_cache.get_universe());
        const auto& clauses = local_clauses(m_ticket);
        universe_and_clauses["infeasibility_map"] = m_inf_map.export_bits();
        universe_and_clauses["best_solution"] = m_best_solution.assignments_as<std::vector<bool>>();
        universe_and_clauses["best_spawners"] = lit::externalize(m_best_spawners);
        universe_and_clauses["best_mutually_exclusive"] = lit::externalize(m_best_mes);
        universe_and_clauses["all_spawners"] = lit::externalize(m_all_spawners);
        universe_and_clauses["coloring_order"] = lit::externalize(m_coloring_order);
        universe_and_clauses["universe_size"] = m_universe_size;
        universe_and_clauses["universe"] = std::move(universe);
        universe_and_clauses["clauses"] = clauses.export_all_clauses();
        universe_and_clauses["num_variables"] = clauses.num_vars();
        universe_and_clauses["num_concrete"] = m_inf_map.get_n_concrete();
        std::ofstream outfile;
        outfile.exceptions(std::ios::badbit | std::ios::failbit);
        outfile.open(universe_and_clauses_file, std::ios::out | std::ios::trunc);
        outfile << universe_and_clauses;
    }

    bool subproblem_reporting_enabled() const noexcept {
        return m_subproblem_dir.has_value();
    }

    std::uint64_t new_subproblem_id() noexcept {
        return m_subproblem_counter++;
    }

    void report_subproblem(const LNSSubproblem& subproblem,
                           const PartialSolution& remaining_configs,
                           std::size_t global_best_mes_size,
                           std::size_t global_best_lower_bound,
                           const char* subproblem_type_)
    {
        if(!subproblem_reporting_enabled()) return;
        auto id = new_subproblem_id();
        std::ostringstream name_fmt;
        std::string subproblem_type(subproblem_type_);
        auto name_filter = [] (const std::string& str) {
            std::string result;
            for(char c : str) {
                if(std::isalnum(c)) {
                    result.push_back(c);
                } else {
                    result.push_back('_');
                }
            }
            return result;
        };
        name_fmt << "subproblem-" << name_filter(subproblem_type) << "-"
                 << std::setw(5) << std::setfill('0') << id << ".json";
        auto subproblem_file = *m_subproblem_dir / name_fmt.str();
        nlohmann::json output;
        output["id"] = id;
        output["subproblem_type"] = subproblem_type;
        output["uncovered"] = lit::externalize(subproblem.uncovered_universe);
        output["initial_uncovered_mes"] = lit::externalize(subproblem.mutually_exclusive_set);
        output["num_nonremoved_configs"] = remaining_configs.size();
        output["num_removed_configs"] = subproblem.removed_configurations.size();
        output["lns_info"] = p_lns_info_to_json();
        output["removed_configs"] = bitsets_to_json(subproblem.removed_configurations);
        output["global_best_mes_size"] = global_best_mes_size;
        output["global_best_lower_bound"] = global_best_lower_bound;
        if(!remaining_configs.empty()) {
            output["remaining_config"] = std::vector<bool>(remaining_configs.get_assignment(0));
        } else {
            output["remaining_config"] = nlohmann::json{};
        }
        std::ofstream outfile;
        outfile.exceptions(std::ios::badbit | std::ios::failbit);
        outfile.open(subproblem_file, std::ios::out | std::ios::trunc);
        outfile << output;
    }

    void report_subproblem(const std::vector<Vertex>& uncovered,
                           const std::vector<Vertex>& initial_uncovered_mes,
                           const PartialSolution& remaining_configs,
                           const std::vector<DynamicBitset>& removed_configs,
                           std::size_t global_best_mes_size,
                           std::size_t global_best_lower_bound,
                           const char* subproblem_type)
    {
        if(!subproblem_reporting_enabled()) return;
        auto id = new_subproblem_id();
        std::ostringstream name_fmt;
        name_fmt << "subproblem-" << subproblem_type << "-" 
                 << std::setw(5) << std::setfill('0') << id << ".json";
        auto subproblem_file = *m_subproblem_dir / name_fmt.str();
        nlohmann::json subproblem;
        subproblem["id"] = id;
        subproblem["subproblem_type"] = subproblem_type;
        subproblem["uncovered"] = lit::externalize(uncovered);
        subproblem["initial_uncovered_mes"] = lit::externalize(initial_uncovered_mes);
        subproblem["num_nonremoved_configs"] = remaining_configs.size();
        subproblem["num_removed_configs"] = removed_configs.size();
        subproblem["lns_info"] = p_lns_info_to_json();
        subproblem["removed_configs"] = bitsets_to_json(removed_configs);
        subproblem["global_best_mes_size"] = global_best_mes_size;
        subproblem["global_best_lower_bound"] = global_best_lower_bound;
        if(!remaining_configs.empty()) {
            subproblem["remaining_config"] = std::vector<bool>(remaining_configs.get_assignment(0));
        } else {
            subproblem["remaining_config"] = nlohmann::json{};
        }
        std::ofstream outfile;
        outfile.exceptions(std::ios::badbit | std::ios::failbit);
        outfile.open(subproblem_file, std::ios::out | std::ios::trunc);
        outfile << subproblem;
    }

  private:
    void p_raise_events(EventMask events) {
        EventMask prev = m_raised;
        if((prev | events) == prev) return;
        m_raised |= events;
        m_condition.notify_one();
    }

    using AlarmTime = decltype(Clock::now() + std::chrono::duration<double>(1.0));
    struct AlarmRequest {
        AlarmTime time;
        std::size_t id;
        PortfolioElement* element;

        bool operator<(const AlarmRequest& other) const noexcept {
            return other.time < time;
        }
    };

    void p_post_reduce_universe() {
        m_global_recorder->store_event("DONE_IMPLIED_VERTEX_ELIMINATION",
                                       {{"original_size", m_implied_cache.original_universe_size()},
                                        {"reduced_size", m_implied_cache.reduced_universe_size()}},
                                       "original_size", "reduced_size");
        m_implied_cache.replace_implied(m_best_spawners);
        sort_unique(m_best_spawners);
        m_implied_cache.replace_implied(m_all_spawners);
        sort_unique(m_all_spawners);
        m_implied_cache.replace_implied(m_coloring_order);
        nosort_unique<PairHashSet<Vertex>>(m_coloring_order);
        m_implied_cache.replace_implied(m_best_mes);
        std::size_t old_mes_size = m_best_mes.size();
        sort_unique(m_best_mes);
        if(old_mes_size != m_best_mes.size()) {
            throw std::logic_error("MES size changed during implied vertex elimination!");
        }
    }

    void p_update_alarms() {
        for(const auto& r : m_alarm_requests) {
            m_alarm_queue.push(r);
        }
        m_alarm_requests.clear();
    }

    bool p_should_wake() const {
        return !m_alarm_requests.empty() || 
               !m_waiting_elements.empty() ||
               m_raised != 0;
    }

    void p_extract_due_alarms() {
        auto now = Clock::now();
        while(!m_alarm_queue.empty() && m_alarm_queue.top().time <= now) {
            m_due_alarms.push_back(m_alarm_queue.top());
            m_alarm_queue.pop();
        }
    }

    bool p_handle_due_alarms(const InterruptionCheckInfo& info) {
        for(const AlarmRequest& a : m_due_alarms) {
            if(!a.element) {
                p_raise_timeout();
                return true;
            } else {
                a.element->deliver_alarm(a.id, info);
            }
        }
        m_due_alarms.clear();
        return false;
    }

    void p_raise_timeout() {
        InterruptionCheckInfo info;
        {
            std::unique_lock l{m_mutex};
            m_shutting_down = true;
            info = p_get_info();
            m_global_recorder->store_event("PORTFOLIO_TIMEOUT");
        }
        p_deliver_events(static_cast<EventMask>(PortfolioEvent::TIMEOUT), info);
    }

    void p_deliver_events(EventMask events, const InterruptionCheckInfo& info) {
        for(auto &e : m_elements) {
            e->deliver_events(events, info);
        }
    }

    void p_update_elements() {
        if(!m_shutting_down) {
            for(auto& e : m_waiting_elements) {
                e->start();
                m_elements.emplace_back(std::move(e));
            }
            m_waiting_elements.clear();
        }
    }

    InterruptionCheckInfo p_main_locked(EventMask& new_events) {
        std::unique_lock l{m_mutex};
        p_update_elements();
        p_update_alarms();
        if(m_alarm_queue.empty()) {
            m_condition.wait(l, [&] () { return p_should_wake(); });
        } else {
            m_condition.wait_until(l, m_alarm_queue.top().time, [&] () { return p_should_wake(); });
        }
        new_events = m_raised;
        m_raised = 0;
        p_update_elements();
        p_update_alarms();
        return p_get_info();
    }

    OutputObject p_lns_info_to_json() const {
        std::unique_lock l{m_mutex};
        OutputObject result{
            {"current_goal_time", m_lns_info.current_goal_time},
            {"num_failure_threshold", m_lns_info.num_failure_threshold},
            {"min_num_tries_for_time", m_lns_info.min_num_tries_for_time},
            {"removed_classes_info", nlohmann::json{}}
        };
        auto& removed_classes_info = result["removed_classes_info"];
        for(const auto& [num_removed, info] : m_lns_info.removed_classes_info) {
            nlohmann::json info_json{
                {"num_removed", num_removed},
                {"complete_try_successes", info.complete_try_successes},
                {"complete_tries_since_last_success", info.complete_tries_since_last_success},
                {"complete_tries_total", info.complete_tries_total},
                {"complete_tries_total_time", info.complete_tries_total_time}
            };
            removed_classes_info.push_back(std::move(info_json));
        }
        return result;
    }

    bool p_main_unlocked(const InterruptionCheckInfo& info, EventMask new_events) {
        p_extract_due_alarms();
        if(p_handle_due_alarms(info)) {
            return true;
        }
        if(new_events) {
            p_deliver_events(new_events, info);
            if(new_events & static_cast<EventMask>(PortfolioEvent::OPTIMALITY)) {
                std::unique_lock l{m_mutex};
                m_shutting_down = true;
                return true;
            }
        }
        return false;
    }

    void p_await_completion() {
        for(auto& e : m_elements) {
            e->join();
        }
    }

    void p_main_loop() {
        for(;;) {
            EventMask new_events = 0;
            auto info = p_main_locked(new_events);
            if(p_main_unlocked(info, new_events)) {
                break;
            }
        }
        p_await_completion();
    }

    InterruptionCheckInfo p_get_info() const {
        return InterruptionCheckInfo{m_best_mes.size(), m_lower_bound, m_best_solution.size()};
    }

    // --- read-only elements ---
    /**
     * Pre-filled pair infeasibility map.
     */
    PairInfeasibilityMap m_inf_map;

    /**
     * ImpliedVertexCache for reducing the universe size.
     */
    ImpliedVertexCache m_implied_cache;

    /**
     * Class spawners for best sample from initial phase.
     */
    std::vector<Vertex> m_best_spawners;

    /**
     * Class spawners for all samples from the initial phase.
     */
    std::vector<Vertex> m_all_spawners;

    /**
     * Coloring order for the last heuristic run of the initial phase.
     */
    std::vector<Vertex> m_coloring_order;

    /**
     * The number of valid interactions.
     */
    std::size_t m_universe_size;

    /**
     * Ticket to access thread-local
     * mutable versions of the clause set.
     */
    ClausesTicket m_ticket;

    /**
     * Directory to which the subproblems are to be reported.
     */
    std::optional<std::filesystem::path> m_subproblem_dir;

    /**
     * Atomic counter for subproblem ids.
     */
    std::atomic<std::uint64_t> m_subproblem_counter{0};

    // --- events and coordination ---
    /**
     * Mutex for all global data and events.
     */
    mutable std::mutex m_mutex;

    /**
     * A flag that tracks if we are currently shutting down.
     */
    bool m_shutting_down{false};

    /**
     * Event recorder for output/logging purposes.
     */
    EventRecorder* m_global_recorder;

    /**
     * Condition variable to wait for events, alarm requests, ...
     */
    mutable std::condition_variable m_condition;

    /**
     * Pending alarm requests (new alarm requests go here).
     */
    std::vector<AlarmRequest> m_alarm_requests;

    /**
     * Unlocked container for alarms that have expired,
     * but have not yet been dispatched.
     */
    std::vector<AlarmRequest> m_due_alarms;

    /**
     * Priority queue of alarms. Contains a 'pseudo-alarm'
     * elements that signifies timeout expiry.
     */
    std::priority_queue<AlarmRequest> m_alarm_queue;

    /**
     * The set of raised events that have to be dispatched.
     */
    EventMask m_raised = 0;

    /**
     * The owning container of PortfolioElements that
     * are currently active.
     * It is not allowed to delete from this set while
     * the algorithm is running.
     */
    std::vector<std::unique_ptr<PortfolioElement>> m_elements;

    /**
     * List of elements that have been created and are waiting
     * to activate (and be started).
     */
    std::vector<std::unique_ptr<PortfolioElement>> m_waiting_elements;

    // --- best solution, LB vertex set, MES ---
    /**
     * Best mutually exclusive set.
     */
    std::vector<Vertex> m_best_mes;

    /**
     * Vertex set of best LB (these vertices require at least m_lower_bound configurations).
     * Need not be equal to m_best_mes.
     */
    std::vector<Vertex> m_best_lb_vertex_set;

    /**
     * Numerical LB value.
     */
    std::size_t m_lower_bound;

    /**
     * Best full solution.
     */
    PartialSolution m_best_solution;

    /**
     * Storage cache for cliques.
     * They are used for the destroy heuristic.
     */
    CliqueStorage m_clique_cache;

    /**
     * Information for deciding how many configurations
     * to remove in the destroy operation.
     */
    LNSTimeAndSuccessInfo m_lns_info;
};

void PortfolioElement::set_alarm(double in_seconds) {
    std::size_t alarm_id;
    {
        std::unique_lock l{mutex};
        if(m_alarm) {
            m_alarm.reset();
        }
        events &= ~static_cast<EventMask>(PortfolioEvent::ALARM);
        alarm_id = m_alarm_id_counter++;
        m_alarm = alarm_id;
    }
    solver->set_alarm(this, in_seconds, alarm_id);
}

void PortfolioElement::discard_alarm() {
    std::unique_lock l{mutex};
    if(m_alarm) {
        m_alarm.reset();
    }
    events &= ~static_cast<EventMask>(PortfolioEvent::ALARM);
}

void PortfolioElement::deliver_alarm(std::size_t alarm_id, const InterruptionCheckInfo& info) {
    std::unique_lock l{mutex};
    if(!m_alarm || alarm_id != *m_alarm) {
        return;
    }
    bool was_empty = (events == 0);
    events |= static_cast<EventMask>(PortfolioEvent::ALARM);
    interrupt_if_necessary(info);
    if(was_empty) {
        condition.notify_one();
    }
}

void PortfolioElement::deliver_events(EventMask new_events, const InterruptionCheckInfo& info) {
    std::unique_lock l{mutex};
    bool was_empty = (events == 0);
    events |= new_events;
    if(new_events & (static_cast<EventMask>(PortfolioEvent::TIMEOUT) | 
                     static_cast<EventMask>(PortfolioEvent::OPTIMALITY))) 
    {
        should_terminate.store(true);
    }
    interrupt_if_necessary(info);
    if(was_empty) {
        condition.notify_one();
    }
}

PairInfeasibilityMap &PortfolioElement::get_infeasibility_map() noexcept {
    return solver->get_infeasibility_map();
}

ClauseDB& PortfolioElement::get_clauses() const {
    return solver->get_clauses();
}

LNSTimeAndSuccessInfo::LNSTimeAndSuccessInfo(PortfolioSolver* solver) noexcept :
    universe_size(solver->get_universe_size())
{}

}

#endif
