#ifndef SAMMY_STAMP_SET_H_INCLUDED_
#define SAMMY_STAMP_SET_H_INCLUDED_

#include <cstdint>
#include <type_traits>
#include <limits>
#include <vector>

namespace sammy {

template<typename ValueType, typename StampType = std::uint32_t>
class StampSet
{
  private:
    static_assert(std::is_integral_v<ValueType>, "Value type must be integral");
    static_assert(std::is_integral_v<StampType> &&
                  std::is_unsigned_v<StampType>, "Stamp type must be unsigned");

    std::vector<StampType> m_stamps;
    StampType m_current_stamp;

  public:
    explicit StampSet(ValueType universe_size) :
        m_stamps(universe_size, StampType(0)),
        m_current_stamp(1)
    {}

    StampSet(const StampSet&) = default;
    StampSet &operator=(const StampSet&) = default;
    StampSet(StampSet&&) noexcept = default;
    StampSet &operator=(StampSet&&) noexcept = default;

    std::size_t universe_size() const noexcept {
        return m_stamps.size();
    }

    void clear() noexcept {
        if(++m_current_stamp == 0) {
            std::fill(m_stamps.begin(), m_stamps.end(), StampType(0));
            m_current_stamp = 1;
        }
    }

    template<typename ForwardIterator>
    void assign(ForwardIterator begin, ForwardIterator end) noexcept {
        clear();
        insert(begin, end);
    }

    template<typename ForwardIterator>
    void insert(ForwardIterator begin, ForwardIterator end) noexcept {
        std::for_each(begin, end, [&] (ValueType l) { insert(l); });
    }

    void insert(ValueType v) noexcept {
        m_stamps[v] = m_current_stamp;
    }

    void erase(ValueType v) noexcept {
        m_stamps[v] = 0;
    }

    bool check_insert(ValueType v) noexcept {
        StampType& s = m_stamps[v];
        bool result = (s != m_current_stamp);
        s = m_current_stamp;
        return result;
    }

    bool check_erase(ValueType v) noexcept {
        StampType& s = m_stamps[v];
        bool result = (s == m_current_stamp);
        s = 0;
        return result;
    }

    bool count(ValueType v) const noexcept {
        return m_stamps[v] == m_current_stamp;
    }

    bool contains(ValueType v) const noexcept {
        return count(v);
    }
};

}

#endif
