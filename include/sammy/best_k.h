#ifndef SAMMY_BEST_K_INCLUDED_
#define SAMMY_BEST_K_INCLUDED_

#include <algorithm>
#include <utility>
#include <vector>

namespace sammy {

/**
 * A container for keeping track of the
 * best k elements of a set of elements.
 * Which elements are best is determined
 * using Compare::operator(), defaulting
 * to < on the elements (in which case,
 * the smallest elements are best).
 */
template <typename T, typename Compare = std::less<T>> class BestK {
  public:
    explicit BestK(std::size_t k, Compare&& compare = {})
        : m_k(k), m_compare(std::forward<Compare>(compare)) {
        m_elements.reserve(m_k);
    }

    void push(T&& element) { p_push(std::move(element)); }

    void push(const T& element) { p_push(element); }

    template <typename ComparableToElement>
    bool would_push(ComparableToElement&& score) const noexcept {
        if (m_elements.size() < m_k) {
            return true;
        }
        return !m_compare(m_elements[0], score);
    }

    template <typename... Args> void emplace(Args&&... args) {
        T new_element(std::forward<Args>(args)...);
        push(std::move(new_element));
    }

    void clear() noexcept { m_elements.clear(); }

    const std::vector<T>& elements() const noexcept { return m_elements; }

  private:
    template <typename TT> void p_push(TT&& element) {
        if (m_elements.size() < m_k) {
            m_elements.emplace_back(std::forward<TT>(element));
            if (m_elements.size() == m_k) {
                std::make_heap(m_elements.begin(), m_elements.end(), m_compare);
            }
            return;
        }
        if (!m_compare(element, m_elements[0])) {
            // we don't have to do anything for
            // elements that aren't better than
            // the worst of the best k elements so far
            return;
        }
        m_elements[0] = std::forward<TT>(element);
        p_sift_down();
    }

    void p_sift_down() {
        const std::size_t s = m_k;
        std::size_t ei = 0;
        while (2 * ei + 2 < s) {
            std::size_t c = 2 * ei + 1;
            if (m_compare(m_elements[c], m_elements[c + 1])) {
                // the left child is better,
                // so if any child, the right one moves up
                ++c;
            }
            // if the current element isn't better than its
            // worst child, we are done.
            if (!m_compare(m_elements[ei], m_elements[c]))
                return;
            std::swap(m_elements[ei], m_elements[c]);
            ei = c;
        }
        const std::size_t c = 2 * ei + 1;
        if (c < s) {
            // only one child; if we are better, swap
            if (m_compare(m_elements[ei], m_elements[c])) {
                std::swap(m_elements[ei], m_elements[c]);
            }
        }
    }

    std::size_t m_k;
    std::vector<T> m_elements;
    Compare m_compare;
};

} // namespace sammy

#endif
