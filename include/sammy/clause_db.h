#ifndef HS_CLAUSE_DB_H_INCLUDED_
#define HS_CLAUSE_DB_H_INCLUDED_

#include "clause_reduction.h"
#include "literals.h"
#include "range.h"

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace sammy {

/**
 * Clause reference type.
 * An unsigned integer type indexing a clause in
 * a clause database.
 */
using CRef = std::uint32_t;

/**
 * A structure keeping a snapshot of 
 * the number of clauses in a clause database.
 */
struct ClauseCounts {
    CRef unary_clause_end;
    CRef binary_clause_end;
    CRef long_clause_end;
};

/**
 * @brief A database of clauses in internal format
 * (i.e., unsigned 0-based even-odd-encoding of literals)
 * that can be shared by multiple SharedDBPropagator instances.
 */
class ClauseDB {
  private:
    // the large array containing all literals of all clauses
    // that are not unary and not binary
    std::vector<Lit> m_literals;

    // array of all literals in unary clauses
    std::vector<Lit> m_unaries;

    // array of arrays of binary watches (other literal)
    std::vector<std::vector<Lit>> m_binaries;

    // array of binary clause literals:
    // duplicated to make keeping track of learnt binaries possible,
    // but separated from regular clauses because accesses
    // should be much rarer
    std::vector<std::pair<Lit, Lit>> m_binary_clauses;

    // number of variables
    std::uint32_t m_num_vars, m_num_clauses;

    // marked as frozen?
    // frozen: must not add new clauses.
    bool m_frozen = false;

    // turn an external clause into an internal one
    void p_internalize_external(const ExternalClause& c) {
        m_literals.push_back(c.size());
        for (auto l : c) {
            m_literals.push_back(lit::internalize(l));
        }
    }

    // turn an array of external clauses (that should be reduced)
    // into internal clauses
    void
    p_internalize_reduced_externals(const std::vector<ExternalClause>& cs) {
        m_num_clauses = cs.size();
        for (const auto& c : cs) {
            if (c.size() == 1) {
                m_unaries.push_back(lit::internalize(c[0]));
            } else if (c.size() == 2) {
                std::uint32_t l1 = lit::internalize(c[0]);
                std::uint32_t l2 = lit::internalize(c[1]);
                m_binaries[l1].push_back(l2);
                m_binaries[l2].push_back(l1);
                m_binary_clauses.emplace_back(l1, l2);
            } else {
                p_internalize_external(c);
            }
        }
    }

    // turn an external clause that we known nothing about
    // (might be tautological/contain duplicates) into an internal one
    void p_internalize_internal(const std::vector<std::vector<Lit>>& cs) {
        m_num_clauses = cs.size();
        std::vector<Lit> buffer;
        for(const auto& c : cs) {
            bool tautology = (std::find_if(c.begin(), c.end(), [&] (Lit l) {
                return std::find(c.begin(), c.end(), lit::negate(l)) != c.end();
            }) != c.end());
            if(tautology) {
                --m_num_clauses;
                continue;
            }
            buffer = c;
            std::sort(buffer.begin(), buffer.end());
            buffer.erase(std::unique(buffer.begin(), buffer.end()), buffer.end());
            if(c.size() == 1) {
                m_unaries.push_back(buffer[0]);
            } else if(c.size() == 2) {
                std::uint32_t l1 = buffer[0];
                std::uint32_t l2 = buffer[1];
                m_binaries[l1].push_back(l2);
                m_binaries[l2].push_back(l1);
                m_binary_clauses.emplace_back(l1, l2);
            } else {
                m_literals.push_back(c.size());
                for (auto l : c) {
                    m_literals.push_back(l);
                }
            }
        }
    }

  public:
    using LitIterator = const Lit*;
    using Lits = IteratorRange<LitIterator>;

    ClauseCounts get_clause_counts() const noexcept {
        return {num_unaries(), num_binaries(), literal_db_size()};
    }

    static ClauseDB import_from(std::istream& input) {
        std::uint32_t num_vars = 0;
        input >> num_vars;
        std::vector<ExternalClause> clauses;
        while (true) {
            ExternalClause cl;
            std::uint32_t s = 0;
            std::int32_t l;
            if (!(input >> s))
                break;
            cl.reserve(s);
            for (std::uint32_t i = 0; i < s; ++i) {
                input >> l;
                cl.push_back(l);
            }
            clauses.emplace_back(std::move(cl));
        }
        return ClauseDB(num_vars, clauses);
    }

    static ClauseDB import_from(const std::filesystem::path& path) {
        std::ifstream input;
        input.exceptions(std::ios::failbit | std::ios::badbit);
        input.open(path, std::ios::in);
        input.exceptions(std::ios::badbit);
        return import_from(input);
    }

    explicit ClauseDB(std::uint32_t num_vars) :
        m_binaries(2 * num_vars, std::vector<Lit>{}),
        m_num_vars(num_vars),
        m_num_clauses(0)
    {}

    explicit ClauseDB(std::uint32_t num_vars,
                      const std::vector<ExternalClause>& clauses)
        : m_binaries(2 * num_vars, std::vector<Lit>{}), m_num_vars(num_vars),
          m_num_clauses(0) {
        auto reduced = reduce_external_clauses(clauses);
        p_internalize_reduced_externals(reduced);
    }

    explicit ClauseDB(std::uint32_t num_vars,
                      const std::vector<std::vector<Lit>>& clauses)
        : m_binaries(2 * num_vars, std::vector<Lit>{}),
          m_num_vars(num_vars), m_num_clauses(0)
    {
        p_internalize_internal(clauses);
    }

    void export_to(std::ostream& output) const {
        output << m_num_vars << '\n';
        auto all_clauses = export_all_clauses();
        for (const auto& ec : all_clauses) {
            output << ec.size();
            for (auto l : ec) {
                output << ' ' << l;
            }
            output << '\n';
        }
    }

    /**
     * Mark the clause database as frozen.
     */
    void mark_frozen() noexcept {
        m_frozen = true;
    }

    /**
     * Check if the clause database is marked as frozen.
     */
    bool is_frozen() const noexcept {
        return m_frozen;
    }

    void export_to(const std::filesystem::path& path) const {
        std::ofstream output(path, std::ios::out | std::ios::trunc);
        output.exceptions(std::ios::badbit | std::ios::failbit);
        export_to(output);
    }

    void export_to_dimacs(std::ostream& output) const {
        output << "p cnf " << m_num_vars << ' ' << m_num_clauses << std::endl;
        auto all_clauses = export_all_clauses();
        for(const auto& ec : all_clauses) {
            for(auto l : ec) {
                output << l << ' ';
            }
            output << '0' << '\n';
        }
    }

    void export_to_dimacs(const std::filesystem::path& path) const {
        std::ofstream output(path, std::ios::out | std::ios::trunc);
        output.exceptions(std::ios::badbit | std::ios::failbit);
        export_to_dimacs(output);
    }

    std::vector<ExternalClause> export_all_clauses() const {
        std::vector<ExternalClause> result;
        result.reserve(m_num_clauses);
        for (Lit l : m_unaries) {
            result.emplace_back(
                std::initializer_list<ExternalLit>{lit::externalize(l)});
        }
        for (auto b : m_binary_clauses) {
            result.emplace_back(std::initializer_list<ExternalLit>{
                lit::externalize(b.first), lit::externalize(b.second)});
        }
        for (CRef current = 1, s = literal_db_size(); current < s;
             current = next_clause(current))
        {
            auto ls = lits_of(current);
            result.emplace_back();
            auto& out = result.back();
            out.reserve(ls.end() - ls.begin());
            for (Lit l : ls) {
                out.push_back(lit::externalize(l));
            }
        }
        return result;
    }

    template<typename Iterator>
    CRef add_clause(Iterator ibeg, Iterator iend)
    {
        if(m_frozen) {
            throw std::logic_error("Cannot add clauses to frozen ClauseDB!");
        }
        CRef res = NIL;
        auto length = std::distance(ibeg, iend);
        if(length == 1) {
            m_unaries.push_back(*ibeg);
        } else if(length == 2) {
            Lit l1 = ibeg[0];
            Lit l2 = ibeg[1];
            m_binaries[l1].push_back(l2);
            m_binaries[l2].push_back(l1);
            m_binary_clauses.emplace_back(l1, l2);
        } else {
            m_literals.push_back(length);
            res = m_literals.size();
            m_literals.insert(m_literals.end(), ibeg, iend);
        }
        ++m_num_clauses;
        return res;
    }

    CRef add_clause(Lits literals) {
        return add_clause(literals.begin(), literals.end());
    }

    CRef cref_of(Lits lits) const noexcept {
        return CRef(lits.begin() - m_literals.data());
    }

    Lits lits_of(CRef clause) const noexcept {
        LitIterator begin = m_literals.data() + clause;
        return {begin, begin + begin[-1]};
    }

    std::uint32_t clause_length(CRef clause) const noexcept {
        return m_literals[clause - 1];
    }

    CRef next_clause(CRef clause) const noexcept {
        return clause + m_literals[clause - 1] + 1;
    }

    std::uint32_t num_vars() const noexcept { return m_num_vars; }

    std::uint32_t num_clauses() const noexcept { return m_num_clauses; }

    std::uint32_t num_unaries() const noexcept { return m_unaries.size(); }

    std::uint32_t num_binaries() const noexcept {
        return m_binary_clauses.size();
    }

    std::size_t total_clause_size() const noexcept {
        return num_unaries() + 2 * num_binaries() + literal_db_size();
    }

    Lits unary_literals(std::uint32_t begin, std::uint32_t end) const noexcept {
        return {m_unaries.data() + begin, m_unaries.data() + end};
    }

    Lits unary_literals() const noexcept {
        return unary_literals(0, m_unaries.size());
    }

    using BinaryClauseIterator = const std::pair<Lit, Lit>*;
    using BinaryClauses = IteratorRange<BinaryClauseIterator>;
    BinaryClauses binary_clauses(std::uint32_t begin,
                                 std::uint32_t end) const noexcept {
        return {m_binary_clauses.data() + begin, m_binary_clauses.data() + end};
    }
    BinaryClauses binary_clauses() const noexcept {
        return binary_clauses(0, m_binary_clauses.size());
    }

    const std::vector<Lit>& binary_partners_of(Lit lit) const noexcept {
        return m_binaries[lit];
    }

    std::uint32_t literal_db_size() const noexcept { return m_literals.size(); }

    std::size_t total_memory_usage() const noexcept {
        std::size_t total_binary_usage = 0;
        for (const auto& v : m_binaries) {
            total_binary_usage += v.capacity();
        }
        total_binary_usage *= sizeof(Lit);
        return m_literals.capacity() * sizeof(Lit) +
               m_unaries.capacity() * sizeof(Lit) +
               m_binaries.capacity() * sizeof(std::vector<Lit>) +
               total_binary_usage +
               m_binary_clauses.capacity() * sizeof(std::pair<Lit, Lit>) +
               sizeof(ClauseDB);
    }
};

class ClauseDBView {
    ClauseDB* db;
    std::vector<std::uint32_t> seen_binary_partners_of;
    std::uint32_t seen_unary_count;
    std::uint32_t seen_binary_count;
    CRef seen_literal_end;

    template <typename NewClauseHandler>
    bool p_handle_new_unaries(NewClauseHandler& handler) {
        if (seen_unary_count != db->num_unaries()) {
            for (Lit l :
                 db->unary_literals(seen_unary_count, db->num_unaries()))
            {
                ++seen_unary_count;
                if (!handler.new_unary(l)) {
                    return false;
                }
            }
        }
        return true;
    }

    template <typename NewClauseHandler>
    bool p_handle_new_larger_clauses(NewClauseHandler& handler) {
        while (seen_literal_end != db->literal_db_size()) {
            CRef clause = seen_literal_end + 1;
            seen_literal_end = db->next_clause(clause) - 1;
            if (!handler.new_clause(db->lits_of(clause))) {
                return false;
            }
        }
        return true;
    }

    template <typename NewClauseHandler>
    bool p_handle_new_binaries(NewClauseHandler& handler) {
        if (seen_binary_count != db->num_binaries()) {
            for (const auto& bc :
                 db->binary_clauses(seen_binary_count, db->num_binaries()))
            {
                ++seen_binary_count;
                seen_binary_partners_of[bc.first]++;
                seen_binary_partners_of[bc.second]++;
                if (!handler.new_binary(bc.first, bc.second)) {
                    return false;
                }
            }
        }
        return true;
    }

  public:
    explicit ClauseDBView(ClauseDB* db) noexcept
        : db(db), seen_binary_partners_of(2 * db->num_vars(), 0),
          seen_unary_count(0), seen_binary_count(0), seen_literal_end(0) {}

    template <typename NewClauseHandler>
    void handle_new_unaries(NewClauseHandler& handler) {
        p_handle_new_unaries(handler);
    }

    template <typename NewClauseHandler>
    void handle_new_long_clauses(NewClauseHandler& handler) {
        p_handle_new_larger_clauses(handler);
    }

    template <typename NewClauseHandler>
    void handle_new_binary_clauses(NewClauseHandler& handler) {
        p_handle_new_binaries(handler);
    }

    template <typename NewClauseHandler>
    void handle_new_clauses(NewClauseHandler& handler) {
        if (!p_handle_new_unaries(handler))
            return;
        if (!p_handle_new_binaries(handler))
            return;
        if (!p_handle_new_larger_clauses(handler))
            return;
    }

    ClauseDB& database() noexcept { return *db; }

    const ClauseDB& database() const noexcept { return *db; }

    ClauseDB::Lits binary_partners_of(Lit l) const noexcept {
        std::size_t lim = seen_binary_partners_of[l];
        const auto& partners = db->binary_partners_of(l);
        return {partners.data(), partners.data() + lim};
    }

    std::size_t total_memory_usage() const noexcept {
        return 4 * seen_binary_partners_of.capacity() + sizeof(ClauseDBView);
    }
};

/**
 * Turn the given ClauseDB into a list of clauses (internal representation, 
 * i.e., 0-based indexing and even/odd encoding of true/false literals).
 */
template<typename ClauseType>
inline std::vector<ClauseType> to_clause_list(const ClauseDB& clauses)
{
    std::vector<ClauseType> result;
    for(Lit l : clauses.unary_literals()) {
        Lit a[1] = {l};
        result.emplace_back(+a, a+1);
    }
    for(auto [l1,l2] : clauses.binary_clauses()) {
        Lit a[2] = {l1,l2};
        result.emplace_back(+a, a+2);
    }
    for(CRef c = 1, ndb = clauses.literal_db_size(); c < ndb; c = clauses.next_clause(c)) {
        auto lits = clauses.lits_of(c);
        result.emplace_back(lits.begin(), lits.end());
    }
    return result;
}

} // namespace sammy

#endif
