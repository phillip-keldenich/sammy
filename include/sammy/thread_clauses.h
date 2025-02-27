#ifndef SAMMY_THREAD_CLAUSES_H_INCLUDED_
#define SAMMY_THREAD_CLAUSES_H_INCLUDED_

#include "clause_db.h"
#include "literals.h"
#include <mutex>
#include <optional>

namespace sammy {

template <typename T, typename Tag = void> class ThreadLocalManager {
  public:
    using Ticket = std::size_t;

    /**
     * Create a new 'slot' for a thread-local object of type T,
     * identified by the given Ticket value.
     * Threads can access a thread-local mutable copy of the object constructed
     * by the given arguments using the returned Ticket.
     */
    template <typename... Args>
    Ticket new_ticket(Args&&... unshared_object_args) {
        std::unique_lock l{m_mutex};
        Ticket result = m_original_objects.size();
        m_original_objects.emplace_back(
            std::forward<Args>(unshared_object_args)...);
        return result;
    }

    T& local(Ticket ticket) { return p_local(ticket); }

    const T& local(Ticket ticket) const { return p_local(ticket); }

  private:
    T& p_local(Ticket ticket) const {
        thread_local std::vector<std::optional<T>> objects;
        while (ticket > objects.size()) {
            objects.emplace_back(std::nullopt);
        }
        if (ticket == objects.size()) {
            std::unique_lock l{m_mutex};
            objects.emplace_back(m_original_objects.at(ticket));
        } else if (!objects[ticket]) {
            std::unique_lock l{m_mutex};
            objects[ticket] = m_original_objects.at(ticket);
        }
        return *objects[ticket];
    }

    mutable std::mutex m_mutex;
    std::vector<T> m_original_objects;
};

namespace detail {
inline ThreadLocalManager<ClauseDB>& get_tl_clause_manager() {
    static ThreadLocalManager<ClauseDB> manager;
    return manager;
}
} // namespace detail

using ClausesTicket = ThreadLocalManager<ClauseDB>::Ticket;

inline ClausesTicket publish_clauses(const ClauseDB& clause_db) {
    return detail::get_tl_clause_manager().new_ticket(clause_db);
}

inline ClauseDB& local_clauses(ClausesTicket ticket) {
    return detail::get_tl_clause_manager().local(ticket);
}

} // namespace sammy

#endif
