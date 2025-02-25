#ifndef SAMMY_CLIQUE_STORAGE_H_INCLUDED_
#define SAMMY_CLIQUE_STORAGE_H_INCLUDED_

#include "literals.h"
#include "range.h"
#include "algorithm_ex.h"
#include <mutex>
#include <memory>
#include <boost/iterator/iterator_facade.hpp>

namespace sammy {
    
class CliqueStorage {
  public:
    /**
     * Create a clique storage that can take
     * cliques with a combined number of up to vertex_limit vertices.
     */
    explicit CliqueStorage(std::size_t vertex_limit = 10'000'000) :
        m_vertex_list(new Vertex[p_make_limit(vertex_limit)]),
        m_index_list(new std::size_t[vertex_limit / 2 + 1]),
        m_last_used(new std::size_t[vertex_limit / 2]),
        m_vertex_limit(vertex_limit)
    {
        m_index_list[0] = 0;
    }

    /**
     * A read-only view of the contents of the clique storage
     * at the point the view was acquired.
     * This allows us to iterate through the cliques without
     * holding a lock while other threads may be adding new cliques
     * (or even purging/moving around cliques when the storage is full).
     * Behaves like a container of CliqueView (IteratorRange<const Vertex*>).
     * The vertices are stored in a contiguous array.
     */
    class StorageView {
        std::shared_ptr<Vertex[]> vertex_list;
        std::shared_ptr<std::size_t[]> clique_index_list;
        std::shared_ptr<std::size_t[]> last_used;
        std::size_t clique_count;

        friend class CliqueStorage;

        StorageView(std::shared_ptr<Vertex[]> v,
                    std::shared_ptr<std::size_t[]> i,
                    std::shared_ptr<std::size_t[]> u,
                    std::size_t cc) noexcept :
            vertex_list(std::move(v)),
            clique_index_list(std::move(i)),
            last_used(std::move(u)),
            clique_count(cc)
        {}

      public:
        using CliqueView = IteratorRange<const Vertex*>;

        CliqueView operator[](std::size_t index) const noexcept {
            assert(index < clique_count);
            std::size_t ibeg = clique_index_list[index];
            std::size_t iend = clique_index_list[index+1];
            return CliqueView{vertex_list.get()+ibeg, vertex_list.get()+iend};
        }

        bool empty() const noexcept { return clique_count == 0; }

        std::size_t size() const noexcept { return clique_count; }

        class Iterator : 
            public boost::iterator_facade<Iterator, CliqueView,
                                          std::random_access_iterator_tag, CliqueView, std::ptrdiff_t> 
        {
          public:
            Iterator() = default;
            
            explicit Iterator(const Vertex* vbase, const std::size_t* index_ptr) noexcept :
                vbase(vbase), index_ptr(index_ptr)
            {}

          private:
            friend boost::iterator_core_access;
            friend class CliqueStorage;

            CliqueView dereference() const noexcept {
                std::size_t ibeg = index_ptr[0];
                std::size_t iend = index_ptr[1];
                return CliqueView{vbase + ibeg, vbase + iend};
            }

            bool equal(const Iterator& other) const noexcept {
                return index_ptr == other.index_ptr;
            }

            void increment() noexcept { ++index_ptr; }
            
            void decrement() noexcept { --index_ptr; }
            
            void advance(std::ptrdiff_t diff) noexcept {
                index_ptr += diff;
            }

            std::ptrdiff_t distance_to(const Iterator& other) const noexcept {
                return other.index_ptr - index_ptr;
            }

            const Vertex* vbase;
            const std::size_t* index_ptr;
        };

        Iterator begin() const noexcept {
            return Iterator{vertex_list.get(), clique_index_list.get()};
        }

        Iterator end() const noexcept {
            return Iterator{vertex_list.get(), clique_index_list.get() + clique_count};
        }
    };

    /**
     * Add a new clique to the CliqueStorage.
     * Two things can happen: 
     *  * The CliqueStorage has space for the new clique.
     *    In that case, the clique is simply added.
     *  * The CliqueStorage is full.
     *    In that case, the clique storage is reduced by
     *    throwing out a quarter of all cliques
     *    until there is space for the new clique,
     *    using a LRU pattern to decide which cliques are purged.
     */
    template<typename Iterator>
    void push_clique(Iterator begin, Iterator end)
    {
        std::size_t cs = std::distance(begin, end);
        if(cs <= 1) return;

        std::unique_lock l{m_lock};
        std::size_t* il = m_index_list.get();
        Vertex* vl = m_vertex_list.get();
        std::size_t current_size = il[m_clique_count];
        std::size_t new_end = current_size + cs;
        if(new_end >= m_vertex_limit) {
            if(2 * cs > m_vertex_limit) {
                m_vertex_limit = 2 * cs;
            }
            p_make_space(cs);
        }
        il[++m_clique_count] = new_end;
        std::copy(begin, end, &vl[current_size]);
        m_last_used[m_clique_count-1] = m_timestamp++;
    }

    /**
     * Notify the CliqueStorage that a certain clique was used.
     */
    void used_clique(StorageView& view, std::size_t index) {
        std::unique_lock l{m_lock};
        view.last_used[index] = m_timestamp++;
    }

    /**
     * Notify the CliqueStorage that a certain clique was used.
     */
    void used_clique(StorageView& view, StorageView::Iterator iterator) {
        auto index = std::size_t(iterator.index_ptr - view.clique_index_list.get());
        used_clique(view, index);
    }

    /**
     * Obtain a view of the clique list.
     * The part of the clique list that is returned is immutable
     * but will stay valid and does not require holding a lock
     * while walking the list of cliques (which may take long).
     */
    StorageView obtain_view() const {
        std::unique_lock l{m_lock};
        return StorageView{m_vertex_list, m_index_list, m_last_used, m_clique_count};
    }

  private:
    static std::size_t& p_make_limit(std::size_t& x) noexcept {
        if(x % 2) x += 1;
        if(x < 1000) x = 1000;
        return x;
    }

    std::size_t p_clique_size(std::size_t index) {
        return m_index_list[index + 1] - m_index_list[index];
    }

    void p_make_space(std::size_t required_space) {
        auto compare_lru = [&] (std::size_t i1, std::size_t i2) { 
            return m_last_used[i1] < m_last_used[i2];
        };
        std::vector<std::size_t> indices = vector(range(m_clique_count));
        std::size_t remove_cliques = (std::min)(1 + m_clique_count / 4, m_clique_count-1);
        auto removed_end = indices.begin() + remove_cliques;
        for(;;) {
            std::nth_element(indices.begin(), removed_end, 
                             indices.end(), compare_lru);
            std::size_t new_free_space = m_vertex_limit - m_index_list[m_clique_count];
            std::for_each(indices.begin(), removed_end, 
                [&] (std::size_t i) {new_free_space += p_clique_size(i);});
            if(new_free_space >= required_space) break;
            remove_cliques = (std::min)(2 * remove_cliques, m_clique_count-1);
            removed_end = indices.begin() + remove_cliques;
        }
        std::shared_ptr<Vertex[]> new_vertex_list(new Vertex[m_vertex_limit]);
        std::shared_ptr<std::size_t[]> new_index_list(new std::size_t[m_vertex_limit/2+1]);
        std::shared_ptr<std::size_t[]> new_last_used(new std::size_t[m_vertex_limit/2]);
        Vertex* new_out = new_vertex_list.get();
        std::size_t out_i = 0;
        new_index_list[0] = 0;
        auto copy_clique = [&] (std::size_t old_index) {
            std::size_t old_begin = m_index_list[old_index];
            std::size_t old_end = m_index_list[old_index + 1];
            new_out = std::copy(&m_vertex_list[old_begin], &m_vertex_list[old_end], new_out);
            new_index_list[out_i + 1] = std::size_t(new_out - new_vertex_list.get());
            new_last_used[out_i] = m_last_used[old_index];
            ++out_i;
        };
        std::for_each(removed_end, indices.end(), copy_clique);
        m_clique_count = out_i;
        m_vertex_list = std::move(new_vertex_list);
        m_index_list = std::move(new_index_list);
        m_last_used = std::move(new_last_used);
    }

    mutable std::mutex m_lock;
    std::shared_ptr<Vertex[]> m_vertex_list;
    std::shared_ptr<std::size_t[]> m_index_list;
    std::shared_ptr<std::size_t[]> m_last_used;
    std::size_t m_clique_count = 0;
    std::size_t m_vertex_limit;
    std::size_t m_timestamp = 0;
};

}

#endif
