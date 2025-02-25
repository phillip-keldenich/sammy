#ifndef SAMMY_LNS_DESTROY_H_INCLUDED_
#define SAMMY_LNS_DESTROY_H_INCLUDED_

#include "partial_solution.h"
#include "thread_interrupt.h"
#include "output.h"
#include "fast_clique.h"
#include "implied_vertices.h"
#include "barrage_lns_subproblem.h"

namespace sammy {

/**
 * When destroying a solution S by removal of
 * C configurations D, we call a destruction
 * 'doomed' if there is a MES of size C that
 * becomes uncovered by the removal.
 * We can argue that such MES M are uniquely covered,
 * i.e., consist only of interactions covered
 * by a single configuration in S.
 * This is due to the following:
 * - if a vertex of M is covered by two configurations in D,
 *   that leads to a contradiction with the fact that 
 *   the vertices in M are mutually exclusive and M has size C.
 * - if a vertex of M is covered by a configuration in D and
 *   another configuration not in D, that that vertex is not uncovered
 *   in the destroyed solution.
 * So, for avoidance of doomed destructions, looking for
 * uniquely covered MES is a good idea.
 */

/**
 * Info on a mutually exclusive set on
 * uniquely covered vertices
 * used to identify doomed destructions
 * during the LNS destroy procedure.
 */
struct LNSDestroyUniquelyCoveredMESInfo {
    LNSDestroyUniquelyCoveredMESInfo(std::vector<Vertex> vertices, const PartialSolution& current_solution) :
        m_uniquely_covered_mes(std::move(vertices))
    {
        for(const Vertex& v : m_uniquely_covered_mes) {
            std::optional<Index> cls;
            current_solution.find_covering_classes(v.first, v.second, [&] (Index i) {
                if(cls) {
                    throw std::logic_error("Uniquely covered MES vertex covered by two classes.");
                }
                cls.emplace(i);
            });
            if(!cls) {
                m_uncovered_vertices.push_back(v);
            } else {
                m_configurations_with_mes_vertex[*cls] = v;
            }
        }
    }

    /**
     * The list of all vertices in the MES; if the solution is a full solution,
     * these vertices are all uniquely covered.
     */
    std::vector<Vertex> m_uniquely_covered_mes;

    /**
     * The set of configuration indices which contain a vertex
     * from m_uniquely_covered_mes.
     */
    HashMap<Index, Vertex> m_configurations_with_mes_vertex;

    /**
     * If the solution is not a full solution, contains a list of
     * vertices from m_uniquely_covered_mes that are not currently covered.
     * Necessary to undo the removal of configurations.
     */
    std::vector<Vertex> m_uncovered_vertices;

    /**
     * On removal of configurations, update the indices.
     */
    void update_indices_on_removal(std::vector<Index> old_to_new) {
        HashMap<Index, Vertex> new_configs_with_vertex;
        for(const auto& entry : m_configurations_with_mes_vertex) {
            Index old_index = entry.first;
            if(old_to_new[old_index] == NIL) {
                m_uncovered_vertices.push_back(entry.second);
            } else {
                new_configs_with_vertex[old_to_new[old_index]] = entry.second;
            }
        }
        m_configurations_with_mes_vertex = std::move(new_configs_with_vertex);
    }

    /**
     * After the removal is undone on failure to find
     * an improved solution, restore the uncovered vertices.
     */
    void undo_removal(const PartialSolution& restored) {
        for(const Vertex& v : m_uncovered_vertices) {
            m_configurations_with_mes_vertex[restored.find_covering_class(v.first, v.second)] = v;
        }
        m_uncovered_vertices.clear();
    }

    /**
     * Check if this MES dooms the given destruction to fail.
     */
    bool dooms(const std::vector<Index>& destruction_indices) const {
        return std::all_of(destruction_indices.begin(), destruction_indices.end(), [this] (Index i) -> bool {
            return m_configurations_with_mes_vertex.count(i);
        });
    }

    void remove_included(std::vector<Index>& indices) const {
        indices.erase(std::remove_if(indices.begin(), indices.end(), [this] (Index i) {
            return m_configurations_with_mes_vertex.count(i);
        }), indices.end());
    }
};

/**
 * LNS destroy operation of different types.
 */
class LNSDestroy {
  public:
    LNSDestroy(EventRecorder* local_recorder, std::size_t worker_index, const ImpliedVertexCache* implied_cache, 
               PartialSolution full_solution, std::vector<Vertex> global_mes, SharedDBPropagator propagator) :
        m_local_recorder(local_recorder),
        m_worker_index(worker_index),
        m_implied_cache(implied_cache),
        m_clique_builder(std::move(propagator)),
        m_current_subproblem(std::move(full_solution)),
        m_best_global_mes(std::move(global_mes))
    {}

    std::size_t total_num_configs() const noexcept {
        return m_current_subproblem.size() + m_currently_removed.size();
    }

    void update_global_mes_if_better(const std::vector<Vertex>& global_mes) {
        if(global_mes.size() > m_best_global_mes.size()) {
            m_best_global_mes = global_mes;
        }
    }

    LNSSubproblem move_out_subproblem() const {
        Lit n_concrete = static_cast<Lit>(m_current_subproblem.get_n_concrete());
        return LNSSubproblem{
            std::move(m_filtered_uncovered), 
            std::move(m_best_initial_mes),
            std::move(m_currently_removed),
            n_concrete 
        };
    }

    /**
     * Get the current remaining assignments.
     */
    const PartialSolution& get_partial() const {
        return m_current_subproblem;
    }

    /**
     * Update the full solution that we want to destroy.
     * Among other side-effects, this invalidates our constraining MESs.
     */
    void update_full_solution_if_better(const PartialSolution& new_full_solution) {
        if(new_full_solution.size() < m_current_subproblem.size() + m_currently_removed.size()) {
            m_currently_removed.clear();
            m_current_subproblem = new_full_solution;
            m_constraining_unique_mes.clear();
            m_superset_constraints.clear();
        }
    }

    void update_full_solution_if_better(PartialSolution&& new_full_solution) {
        if(new_full_solution.size() < m_current_subproblem.size() + m_currently_removed.size()) {
            m_currently_removed.clear();
            m_current_subproblem = std::move(new_full_solution);
            m_constraining_unique_mes.clear();
            m_superset_constraints.clear();
        }
    }

    /**
     * Integrate a new MES into the constraints.
     * Must be a 'uniquely covered' MES.
     */
    void integrate_mes(const std::vector<Vertex>& mes) {
        m_constraining_unique_mes.emplace_back(mes, m_current_subproblem);
    }

    /**
     * Called on success to return the subproblem and the improved solution.
     */
    void return_subproblem_on_success(LNSSubproblem&& returned, PartialSolution&& improved) {
        m_currently_removed = std::move(returned.removed_configurations);
        m_filtered_uncovered = std::move(returned.uncovered_universe);
        m_best_initial_mes = std::move(returned.mutually_exclusive_set);
        m_currently_removed.clear();
        m_current_subproblem = std::move(improved);
        m_constraining_unique_mes.clear();
        m_superset_constraints.clear();
    }

    /**
     * Called on abort (or lost race against other solver) to return the subproblem 
     * and undo the destruction.
     */
    void return_subproblem_on_abort(LNSSubproblem&& returned) noexcept {
        m_currently_removed = std::move(returned.removed_configurations);
        m_filtered_uncovered = std::move(returned.uncovered_universe);
        m_best_initial_mes = std::move(returned.mutually_exclusive_set);
        p_undo_removal();
    }

    std::size_t worker_index() const noexcept {
        return m_worker_index;
    }

    /**
     * Called on impossible improvement to
     * return the subproblem and undo the destruction.
     */
    void improvement_impossible(LNSSubproblem&& returned, const std::vector<Vertex>& mes) {
        m_currently_removed = std::move(returned.removed_configurations);
        m_filtered_uncovered = std::move(returned.uncovered_universe);
        m_best_initial_mes = std::move(returned.mutually_exclusive_set);
        m_best_initial_mes.clear();
        if(mes.size() >= m_currently_removed.size()) {
            integrate_mes(mes);
        } else {
            std::vector<Index> indices(m_currently_removed.size(), 0);
            std::iota(indices.begin(), indices.end(), Index(m_current_subproblem.size()));
            m_superset_constraints.push_back(std::move(indices));
        }
        p_undo_removal();
    }

    /**
     * Trigger the destroy operation with a given number
     * of goal configurations to remove.
     * Throws an interrupt exception if the operation
     * detects the interruption flag being set.
     * After that exception, the class remains in a usable state.
     * Returns number of removed configurations;
     * return 0 iff the solution is provably optimal and cannot be improved.
     */
    std::size_t destroy(std::size_t goal_removed);

    /**
     * After finding a replacement for the destroyed
     * part of the solution, make a new full solution.
     */
    PartialSolution improve_destroyed(const std::vector<DynamicBitset>& replacement) {
        PartialSolution result = m_current_subproblem;
        for(const auto& bs : replacement) {
            result.add_assignment(bs);
        }
        return result;
    }

    const std::vector<Vertex>& best_global_mes() const noexcept {
        return m_best_global_mes;
    }

  private:
    /**
     * Event recorder for destruction events.
     */
    EventRecorder* m_local_recorder;

    /**
     * The worker index to add to events for identification.
     */
    std::size_t m_worker_index;

    /**
     * Probability for destruction by random + least-uniquely-covered removal.
     */
    double m_prob_random_least_unique = 0.5;

    /**
     * Probability for purely random destruction.
     */
    double m_prob_random_destruction = 0.2;

    /**
     * Probability for destruction avoiding doom.
     */
    //double m_prob_avoid_doom = 0.3;

    /**
     * Cache storing implied vertices.
     */
    const ImpliedVertexCache* m_implied_cache;

    /**
     * Local clique heuristic to build MESs.
     */
    FastCliqueBuilder m_clique_builder;

    /**
     * PartialSolution encoding the (remaining)
     * subproblem (before destruction, a solution
     * that we want to destroy).
     */
    PartialSolution m_current_subproblem;

    /**
     * List of currently removed configurations
     * (i.e., configurations we removed to obtain m_current_subproblem from a full solution).
     */
    std::vector<DynamicBitset> m_currently_removed;

    /**
     * Uniquely-covered MESs constrain our destructions
     * to avoid doomed destructions.
     */
    std::vector<LNSDestroyUniquelyCoveredMESInfo> m_constraining_unique_mes;

    /**
     * Filtered uncovered vertices after destruction.
     */
    std::vector<Vertex> m_filtered_uncovered;

    /**
     * Best initial MES found for the current destruction.
     */
    std::vector<Vertex> m_best_initial_mes;

    /**
     * The best global MES.
     */
    std::vector<Vertex> m_best_global_mes;

    /**
     * Sets for which we have established that no improvement is possible.
     */
    std::vector<std::vector<Index>> m_superset_constraints;

    /**
     * Attempt to destroy actual_to_remove classes purely at random
     * for at most num_tries tries.
     * Returns true if the destruction looks successful, i.e., does not
     * violate any constraints.
     */
    inline bool p_destroy_random(Index actual_to_remove, std::size_t num_tries);

    /**
     * Destroy mostly randomly but use the MESs to avoid doomed destructions
     * as much as possible.
     */
    inline bool p_destroy_avoid_doom(Index actual_to_remove, std::size_t num_tries);

    /**
     * Destroy actual_to_remove classes by a mix of random and least-unique-covered.
     */
    inline bool p_destroy_random_and_least_unique(Index actual_to_remove, std::size_t num_tries);

    /**
     * Check if the given list of indices is a doomed removal.
     */
    inline bool p_removal_doomed(const std::vector<Index>& removal_indices) const;

    /**
     * Check if the current removal is doomed by a MES;
     * if so, undo it and return true, else return false.
     */
    inline bool p_undo_if_doomed();

    /**
     * Actually perform the removal of the given indices.
     * This may be undone later using the m_currently_removed list,
     * but the indices may change.
     * Updates all indices of configurations accordingly.
     */
    inline void p_perform_removal(const std::vector<Index>& removal_indices, bool enum_uncovered = true);

    /**
     * Undo a removal if the removal did not lead to an improved solution.
     */
    inline void p_undo_removal();

    /**
     * Run the heuristic MES after removal of the given indices.
     */
    inline bool p_removal_heuristic_mes_dooms();

    /**
     * Consider the given MES (which may include covered interactions if all_known_uncovered = false).
     */
    inline bool p_heuristic_mes_from_dooms(std::vector<Vertex>& mes, 
                                           bool all_known_uncovered, bool integrate_unextended);

    template<typename... Args>
    inline void p_store_event(const char* event_name, OutputObject object, Args&&... args) {
        if(!m_local_recorder) return;
        object["worker_index"] = m_worker_index;
        m_local_recorder->store_event(event_name, std::move(object), std::forward<Args>(args)..., "worker_index");
    }
};

// -------------------- IMPLEMENTATION --------------------
std::size_t LNSDestroy::destroy(std::size_t goal_removed) {
    if(goal_removed < 3) {
        goal_removed = 3;
    }
    p_store_event("LNS_DESTROY_BEGIN", {{"goal_removed", goal_removed}}, "goal_removed");
    auto& rng = sammy::rng();
    if(!m_currently_removed.empty()) {
        throw std::logic_error("Entering destroy with partially-destroyed solution!");
    }
    if(goal_removed >= m_current_subproblem.size()) {
        goal_removed = m_current_subproblem.size();
        p_store_event("LNS_DESTROY_WHOLE_SOLUTION", {});
        std::vector<Index> indices = vector(range(Index(goal_removed)));
        if(p_removal_doomed(indices)) {
            p_store_event("LNS_DESTROY_SOLUTION_OPTIMAL", {});
            return 0;
        }
        p_perform_removal(indices);
        if(p_removal_heuristic_mes_dooms()) {
            p_store_event("LNS_DESTROY_SOLUTION_OPTIMAL", {});
            return 0;
        }
        p_store_event("LNS_DESTROY_WHOLE_SOLUTION_DONE", {});
        return goal_removed;
    }
    double x = std::uniform_real_distribution<double>(0.0, 1.0)(rng);
    if(x < m_prob_random_least_unique) {
        if(p_destroy_random_and_least_unique(goal_removed, 5)) {
            return goal_removed;
        }
    } else if(x < m_prob_random_destruction + m_prob_random_least_unique) {
        if(p_destroy_random(goal_removed, 10)) {
            return goal_removed;
        }
    }
    if(p_destroy_avoid_doom(goal_removed, 10)) {
        return goal_removed;
    }
    // TODO: current fallback is increasing the subproblem size
    return destroy(goal_removed + 1);
    //return p_destroy_fallback(goal_removed);
}

bool LNSDestroy::p_destroy_random(Index actual_to_remove, std::size_t num_tries) {
    for(std::size_t i = 0; i < num_tries; ++i) {
        throw_if_interrupted();
        p_store_event("LNS_DESTROY_RANDOM_BEGIN_TRIAL", {{"trial", i+1}, {"goal_removed", actual_to_remove}}, "goal_removed", "trial");
        std::vector<Index> indices = vector(range(actual_to_remove));
        std::shuffle(indices.begin(), indices.end(), sammy::rng());
        indices.resize(actual_to_remove);
        std::sort(indices.begin(), indices.end());
        if(p_removal_doomed(indices)) continue;
        p_perform_removal(indices);
        if(p_removal_heuristic_mes_dooms()) {
            p_undo_removal();
            continue;
        }
        p_store_event("LNS_DESTROY_RANDOM_DONE", {{"goal_removed", actual_to_remove}}, "goal_removed");
        return true;
    }
    p_store_event("LNS_DESTROY_RANDOM_FAILED", {{"goal_removed", actual_to_remove}, {"num_trials", num_tries}}, "goal_removed");
    return false;
}

bool LNSDestroy::p_destroy_avoid_doom(Index actual_to_remove, std::size_t num_tries) {
    auto& rng = sammy::rng();
    DynamicBitset available(m_current_subproblem.size(), true);
    DynamicBitset superset_buffer(m_current_subproblem.size(), false);
    std::vector<Index> buffer;
    std::vector<Index> current_removal;

    auto available_to_buffer = [&] (const auto& constraint) {
        buffer.clear();
        for(Index i : range(Index(m_current_subproblem.size()))) {
            if(available[i] && !constraint.m_configurations_with_mes_vertex.count(i)) {
                buffer.push_back(i);
            }
        }
    };

    auto superset_available_to_buffer = [&] (const auto& constraint) {
        buffer.clear();
        for(Index i : range(Index(m_current_subproblem.size()))) {
            if(available[i] && !superset_buffer[i]) {
                buffer.push_back(i);
            }
        }
    };
    
    auto superset_dooms = [&] (const auto& constraint) {
        superset_buffer.reset();
        for(Index i : constraint) {
            superset_buffer[i].set();
        }
        for(Index i : current_removal) {
            if(!superset_buffer[i]) {
                return false;
            }
        }
        return true;
    };

    auto random_extend = [&] () -> bool {
        if(buffer.empty()) {
            return false;
        }
        Index next_to_remove = buffer[std::uniform_int_distribution<Index>(0, buffer.size() - 1)(rng)];
        current_removal.push_back(next_to_remove);
        available[next_to_remove] = false;
        if(current_removal.size() > actual_to_remove) {
            return false;
        }
        return true;
    };

    auto run_try = [&] (std::size_t trial) -> bool {
        throw_if_interrupted();
        p_store_event("LNS_DESTROY_AVOID_DOOM_BEGIN_TRIAL", 
                      {{"goal_removed", actual_to_remove}, {"trial", trial}}, "goal_removed", "trial");
        current_removal.clear();
        available.set();
        for(const auto& constraint : m_constraining_unique_mes) {
            if(constraint.dooms(current_removal)) {
                available_to_buffer(constraint);
                if(!random_extend()) return false;
            }
        }
        for(const auto& constraint : m_superset_constraints) {
            if(superset_dooms(constraint)) {
                superset_available_to_buffer(constraint);
                if(!random_extend()) return false;
            }
        }
        std::size_t remaining_to_remove = actual_to_remove - current_removal.size();
        if(remaining_to_remove > 0) {
            buffer.clear();
            for(Index i : range(Index(m_current_subproblem.size()))) {
                if(available[i]) {
                    buffer.push_back(i);
                }
            }
            std::shuffle(buffer.begin(), buffer.end(), rng);
            buffer.resize(remaining_to_remove);
            current_removal.insert(current_removal.end(), buffer.begin(), buffer.end());
        }
        std::sort(current_removal.begin(), current_removal.end());
        p_store_event("LNS_DESTROY_AVOID_DOOM_DONE", 
                      {{"goal_removed", actual_to_remove}, {"trial", trial}}, "goal_removed", "trial");
        return true;
    };

    for(std::size_t i = 0; i < num_tries; ++i) {
        if(run_try(i + 1)) {
            p_perform_removal(current_removal);
            if(p_removal_heuristic_mes_dooms()) {
                p_undo_removal();
                continue;
            }
            return true;
        }
    }
    return false;
}

bool LNSDestroy::p_undo_if_doomed() {
    for(const auto& info : m_constraining_unique_mes) {
        if(info.m_uncovered_vertices.size() >= m_currently_removed.size()) {
            p_undo_removal();
            return true;
        }
    }
    return false;
}

bool LNSDestroy::p_destroy_random_and_least_unique(Index goal_removed, std::size_t num_tries) {
    std::size_t goal_removed_random = (std::min)(goal_removed, goal_removed / 3 + 1);
    std::size_t goal_removed_least_unique = goal_removed - goal_removed_random;
    if(goal_removed_least_unique == 0) {
        return p_destroy_random(goal_removed_random, num_tries);
    }
    p_store_event("BEGIN_DESTROY_RANDOM_AND_LEAST_UNIQUE", {
            {"goal_removed", goal_removed},
            {"random", goal_removed_random},
            {"least_unique", goal_removed_least_unique}
        }, "goal_removed");
    for(std::size_t t = 0; t < num_tries; ++t) {
        throw_if_interrupted();
        p_store_event("BEGIN_DESTROY_RANDOM_AND_LEAST_UNIQUE_TRIAL", {
                {"goal_removed", goal_removed},
                {"random", goal_removed_random},
                {"least_unique", goal_removed_least_unique},
                {"trial", t + 1}
            }, "goal_removed", "trial");
        std::vector<Index> indices = vector(range(Index(m_current_subproblem.size())));
        std::shuffle(indices.begin(), indices.end(), sammy::rng());
        indices.resize(goal_removed_random);
        std::sort(indices.begin(), indices.end());
        p_perform_removal(indices, /*enum_uncovered=*/false);
        std::vector<std::size_t> unique_per_class(m_current_subproblem.size(), 0);
        m_current_subproblem.iterate_all_uniquely_covered_with_class(
            [&] (Index i, Lit, Lit) {
                unique_per_class[i] += 1;
            });
        indices.clear();
        auto new_range = range(Index(m_current_subproblem.size()));
        indices.assign(new_range.begin(), new_range.end());
        std::nth_element(
            indices.begin(), indices.begin() + goal_removed_least_unique, indices.end(),
            [&] (Index a, Index b) { return unique_per_class[a] < unique_per_class[b]; }
        );
        indices.resize(goal_removed_least_unique);
        std::sort(indices.begin(), indices.end());
        p_perform_removal(indices, /*enum_uncovered=*/true);
        if(p_undo_if_doomed()) {
            continue;
        }
        if(p_removal_heuristic_mes_dooms()) {
            p_undo_removal();
            continue;
        }
        p_store_event("DONE_DESTROY_RANDOM_AND_LEAST_UNIQUE_DONE", {
                {"goal_removed", goal_removed},
                {"random", goal_removed_random},
                {"least_unique", goal_removed_least_unique},
                {"trial", t + 1},
                {"filtered_uncovered", m_filtered_uncovered.size()}
            }, "goal_removed", "trial", "filtered_uncovered");
        return true;
    }
    return false;
}

bool LNSDestroy::p_removal_doomed(const std::vector<Index>& removal_indices) const {
    if(std::any_of(m_constraining_unique_mes.begin(), m_constraining_unique_mes.end(), 
        [&] (const LNSDestroyUniquelyCoveredMESInfo& info) -> bool { 
            return info.dooms(removal_indices); 
        }
    )) {
        return true;
    }
    auto is_subseteq = [&] (const std::vector<Index>& superset) {
        return std::includes(superset.begin(), superset.end(), removal_indices.begin(), removal_indices.end());
    };
    return std::any_of(m_superset_constraints.begin(), m_superset_constraints.end(), is_subseteq);
}

void LNSDestroy::p_perform_removal(const std::vector<Index>& removal_indices, bool enum_uncovered) {
    std::vector<Index> old_to_new(m_current_subproblem.size(), NIL);
    auto current = removal_indices.begin(), end = removal_indices.end();
    Index new_index = 0;
    for(Index config_index : range(Index(m_current_subproblem.size()))) {
        if(config_index == *current) {
            if(++current == end) {
                for(Index c2 = config_index + 1, n = m_current_subproblem.size(); c2 < n; ++c2) {
                    old_to_new[c2] = new_index++;
                }
                break;
            }
            continue;
        } else {
            old_to_new[config_index] = new_index++;
        }
    }
    for(auto& info : m_constraining_unique_mes) {
        info.update_indices_on_removal(old_to_new);
    }
    // handle already-removed indices
    for(Index x = 0; x < m_currently_removed.size(); ++x) {
        old_to_new.push_back(new_index++);
    }
    for(Index i : removal_indices) {
        old_to_new[i] = new_index++;
    }
    for(auto& sup : m_superset_constraints) {
        std::transform(sup.begin(), sup.end(), sup.begin(), [&] (Index i) { return old_to_new[i]; });
        std::sort(sup.begin(), sup.end());
    }
    for(auto& r : m_current_subproblem.remove_assignments(removal_indices)) {
        m_currently_removed.push_back(std::move(r));
    }
    if(enum_uncovered) {
        try {
            m_filtered_uncovered.clear();
            m_current_subproblem.iterate_all_uncovered([this] (Lit lmin, Lit lmax) {
                m_filtered_uncovered.emplace_back(lmin, lmax);
            }, /*interruptible=*/true);
            m_implied_cache->remove_implied(m_filtered_uncovered, m_clique_builder.propagator());
        } catch(const InterruptError& e) {
            p_undo_removal();
            throw e;
        }
    }
}

void LNSDestroy::p_undo_removal() {
    m_filtered_uncovered.clear();
    for(auto& assignment : m_currently_removed) {
        m_current_subproblem.add_assignment(std::move(assignment));
    }
    m_currently_removed.clear();
    for(auto& info : m_constraining_unique_mes) {
        info.undo_removal(m_current_subproblem);
    }
}

bool LNSDestroy::p_heuristic_mes_from_dooms(std::vector<Vertex>& mes, bool all_known_uncovered, bool integrate_unextended) {
    if(!all_known_uncovered) {
        // if necessary, remove covered interactions
        auto is_covered = [this] (Vertex v) {
            return !std::binary_search(m_filtered_uncovered.begin(), m_filtered_uncovered.end(), v);
        };
        mes.erase(std::remove_if(mes.begin(), mes.end(), is_covered), mes.end());
    }
    if(mes.size() >= m_currently_removed.size()) {
        if(integrate_unextended) {
            // if this mes is not already in place, we can add it
            integrate_mes(mes);
        }
        return true;
    }
    mes = m_clique_builder.compute_clique_known_valid(mes.begin(), mes.end(), m_filtered_uncovered);
    if(mes.size() >= m_currently_removed.size()) {
        integrate_mes(mes);
        return true;
    }
    if(mes.size() > m_best_initial_mes.size()) {
        m_best_initial_mes = mes;
    }
    return false;
}

bool LNSDestroy::p_removal_heuristic_mes_dooms() {
    auto mes = m_clique_builder.random_multistart_best_clique_known_valid(10, m_filtered_uncovered);
    if(mes.size() >= m_currently_removed.size()) {
        integrate_mes(mes);
        p_store_event("INITIAL_MES_PRECLUDES_IMPROVEMENT",
                      {{"removed_configs", m_currently_removed.size()},
                       {"mes_size", mes.size()},
                       {"mes_source", "random_multistart_clique"}}, "removed_configs", "mes_size", "mes_source");
        return true;
    }
    m_best_initial_mes = mes;
    mes.clear();
    for(const auto& info : m_constraining_unique_mes) {
        if(info.m_uncovered_vertices.size() > mes.size()) {
            mes = info.m_uncovered_vertices;
        }
    }
    if(!mes.empty() && p_heuristic_mes_from_dooms(mes, true, false)) {
        p_store_event("UNIQUELY_COVERED_MES_PRECLUDES_IMPROVEMENT",
                      {{"removed_configs", m_currently_removed.size()},
                       {"mes_size", mes.size()},
                       {"mes_source", "cached"}}, "removed_configs", "mes_size", "mes_source");
        return true;
    }
    std::vector<Vertex> tmp = m_best_global_mes;
    if(!m_best_global_mes.empty() && p_heuristic_mes_from_dooms(tmp, false, true)) {
        p_store_event("UNIQUELY_COVERED_MES_PRECLUDES_IMPROVEMENT",
                      {{"removed_configs", m_currently_removed.size()},
                       {"mes_source", "global best locally extended"}}, "removed_configs", "mes_source");
        return true;
    }
    return false;
}

}

#endif
