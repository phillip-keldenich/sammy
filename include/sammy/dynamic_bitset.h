#ifndef HS_DYNAMIC_BITSET_H_
#define HS_DYNAMIC_BITSET_H_

#if (defined(WIN32) || defined(_WIN32)) && defined(_MSC_VER) &&                \
    !defined(__clang__)
#include <intrin.h>
#pragma intrinsic(_BitScanForward64)
#pragma intrinsic(__popcnt64)
#endif

#include "range.h"
#include <cassert>
#include <climits>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <limits>
#include <vector>
#include <algorithm>

namespace sammy {

// the template stuff is necessary because the call might be ambiguous
// otherwise.
template <typename T, std::enable_if_t<std::is_integral_v<T> &&
                                           sizeof(T) == sizeof(std::uint32_t),
                                       int> = 0>
static inline auto count_set_bits(std::uint32_t x) noexcept {
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_popcount(x);
#elif defined(_MSC_VER)
    return __popcnt(x);
#else
    std::uint32_t count_2 =
        (x & UINT32_C(0x5555'5555)) + ((x & UINT32_C(0xAAAA'AAAA)) >> 1);
    std::uint32_t count_4 = (count_2 & UINT32_C(0x3333'3333)) +
                            ((count_2 & UINT32_C(0xCCCC'CCCC)) >> 2);
    std::uint32_t count_8 = (count_4 & UINT32_C(0x0F0F'0F0F)) +
                            ((count_4 & UINT32_C(0xF0F0'F0F0)) >> 4);
    std::uint32_t count_w = (count_8 & UINT32_C(0x00FF'00FF)) +
                            ((count_8 & UINT32_C(0xFF00'FF00)) >> 8);
    return (count_w & UINT32_C(0x0000'FFFF)) +
           ((count_w & UINT32_C(0xFFFF'0000)) >> 16);
#endif
}

// the template stuff is necessary because the call might be ambiguous
// otherwise.
template <typename T, std::enable_if_t<std::is_integral_v<T> &&
                                       sizeof(T) == sizeof(std::uint64_t), int> = 0>
static inline auto count_set_bits(T x) noexcept {
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_popcountll(x);
#elif defined(_MSC_VER)
    return __popcnt64(x);
#else
    std::uint64_t count_2 = (x & UINT64_C(0x5555'5555'5555'5555)) +
                            ((x & UINT64_C(0xAAAA'AAAA'AAAA'AAAA)) >> 1);
    std::uint64_t count_4 = (count_2 & UINT64_C(0x3333'3333'3333'3333)) +
                            ((count_2 & UINT64_C(0xCCCC'CCCC'CCCC'CCCC)) >> 2);
    std::uint64_t count_8 = (count_4 & UINT64_C(0x0F0F'0F0F'0F0F'0F0F)) +
                            ((count_4 & UINT64_C(0xF0F0'F0F0'F0F0'F0F0)) >> 4);
    std::uint64_t count_w = (count_8 & UINT64_C(0x00FF'00FF'00FF'00FF)) +
                            ((count_8 & UINT64_C(0xFF00'FF00'FF00'FF00)) >> 8);
    std::uint64_t count_d = (count_w & UINT64_C(0x0000'FFFF'0000'FFFF)) +
                            ((count_w & UINT64_C(0xFFFF'0000'FFFF'0000)) >> 16);
    return (count_d & UINT64_C(0x0000'0000'FFFF'FFFF)) +
           ((count_d & UINT64_C(0xFFFF'FFFF'0000'0000)) >> 32);
#endif
}

template <typename T, std::enable_if_t<std::is_integral_v<T> &&
                                       sizeof(T) == sizeof(std::uint32_t), int> = 0>
static inline auto count_trailing_zeros(T x) noexcept {
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_ctz(x);
#elif defined(_MSC_VER)
    unsigned long result;
    _BitScanForward(&result, x);
    return result;
#else
    std::size_t result = 0;
    for (std::uint32_t mask = 1; result < 32; ++result, mask <<= 1) {
        if (x & mask)
            break;
    }
    return result;
#endif
}

template <typename T, std::enable_if_t<std::is_integral_v<T> &&
                                       sizeof(T) == sizeof(std::uint64_t), int> = 0>
static inline auto count_trailing_zeros(T x) {
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_ctzll(x);
#elif defined(_MSC_VER)
    unsigned long result;
    _BitScanForward64(&result, x);
    return result;
#else
    std::size_t result = 0;
    for (std::uint64_t mask = 1; result < 64; ++result, mask <<= 1) {
        if (x & mask)
            break;
    }
    return result;
#endif
}

/**
 * @brief Since std::vector<bool> lacks efficient counting/iteration,
 *        and boost::dynamic_bitset does not do doubling growth
 *        unless std::vector::resize does, here goes another bitset.
 *        Also, we are hopefully more auto-vectorization-friendly
 *        than boost::dynamic_bitset, which struggles to get
 *        auto-vectorized (at least by GCC 11) for operations
 *        like -= and such.
 */
class DynamicBitset {
  public:
    using Block = std::size_t;
    static constexpr std::size_t bits_per_word = sizeof(Block) * CHAR_BIT;
    using const_reference = bool;
    using size_type = std::size_t;

    static constexpr std::size_t
    num_blocks_required(std::size_t n_bits) noexcept {
        return (n_bits / bits_per_word) + bool(n_bits % bits_per_word);
    }

    static constexpr Block broadcast(bool v) noexcept {
        // usually, the compiler is clever enough for this to
        // get compiled to sensible branch-free code
        return v ? ~Block(0) : Block(0);
    }

    // Constructors
    DynamicBitset() noexcept : m_num_bits(0) {}

    explicit DynamicBitset(std::size_t s, bool value = false)
        : m_blocks(num_blocks_required(s), broadcast(value)), m_num_bits(s) {
        p_zero_unused();
    }

    explicit DynamicBitset(const std::vector<bool>& vbool)
        : m_blocks(num_blocks_required(vbool.size()), broadcast(false)),
          m_num_bits(vbool.size())
    {
        for(std::size_t i = 0, s = vbool.size(); i != s; ++i) {
            if(vbool[i]) (*this)[i].set();
        }
    }

    bool operator==(const DynamicBitset& o) const noexcept {
        return size() == o.size() &&
            std::equal(m_blocks.begin(), m_blocks.end(), o.m_blocks.begin());
    }

    bool operator!=(const DynamicBitset& o) const noexcept {
        return !(*this == o);
    }

    DynamicBitset& operator|=(const DynamicBitset& o) noexcept {
        assert(m_num_bits == o.m_num_bits);
        std::transform(m_blocks.begin(), m_blocks.end(), o.m_blocks.begin(),
                       m_blocks.begin(),
                       [](Block b1, Block b2) { return b1 | b2; });
        return *this;
    }

    void binary_or(const DynamicBitset& o, std::size_t begin_offset) {
        assert(begin_offset < m_num_bits);
        assert(m_num_bits == o.m_num_bits);
        std::size_t bblock(begin_offset / bits_per_word);
        std::uint8_t boffs(begin_offset % bits_per_word);
        Block o_first = o.m_blocks[bblock];
        o_first &= ~Block((Block(1) << boffs) - 1);
        m_blocks[bblock] |= o_first;
        ++bblock;
        std::transform(m_blocks.begin() + bblock, m_blocks.end(),
                       o.m_blocks.begin() + bblock, m_blocks.begin() + bblock,
                       [](Block b1, Block b2) { return b1 | b2; });
    }

    void binary_subtract(const DynamicBitset& o, std::size_t begin_offset) {
        assert(begin_offset < m_num_bits);
        assert(m_num_bits <= o.m_num_bits);
        std::size_t bblock(begin_offset / bits_per_word);
        std::size_t boffs(begin_offset % bits_per_word);
        Block o_first = o.m_blocks[bblock];
        o_first &= ~Block((Block(1) << boffs) - 1);
        m_blocks[bblock] &= ~o_first;
        ++bblock;
        std::transform(m_blocks.begin() + bblock, m_blocks.end(),
                       o.m_blocks.begin() + bblock, m_blocks.begin() + bblock,
                       [] (Block b1, Block b2) { return b1 & ~b2; });
    }

    void binary_subtract(const std::uint32_t* bits_begin) {
        constexpr std::size_t blocks_per_u32 = sizeof(Block) / sizeof(std::uint32_t);
        const std::size_t nblocks = m_blocks.size();
        for(std::size_t i = 0; i < nblocks; ++i) {
            Block blk = p_read_block(bits_begin + i * blocks_per_u32);
            m_blocks[i] &= ~blk;
        }
    }

    void binary_or(const std::uint32_t* bits_begin) {
        constexpr std::size_t blocks_per_u32 = sizeof(Block) / sizeof(std::uint32_t);
        const std::size_t nblocks = m_blocks.size();
        for(std::size_t i = 0; i < nblocks; ++i) {
            Block blk = p_read_block(bits_begin + i * blocks_per_u32);
            m_blocks[i] |= blk;
        }
    }

    void binary_and(const std::uint32_t* bits_begin) {
        constexpr std::size_t blocks_per_u32 = sizeof(Block) / sizeof(std::uint32_t);
        const std::size_t nblocks = m_blocks.size();
        for(std::size_t i = 0; i < nblocks; ++i) {
            Block blk = p_read_block(bits_begin + i * blocks_per_u32);
            m_blocks[i] &= blk;
        }
    }

    void binary_and(const DynamicBitset& o, std::size_t begin_offset) {
        assert(begin_offset < m_num_bits);
        assert(m_num_bits <= o.m_num_bits);
        std::size_t bblock(begin_offset / bits_per_word);
        std::uint8_t boffs(begin_offset % bits_per_word);
        Block o_first = o.m_blocks[bblock];
        o_first |= Block((Block(1) << boffs) - 1);
        m_blocks[bblock] &= o_first;
        ++bblock;
        std::transform(m_blocks.begin() + bblock, m_blocks.end(),
                       o.m_blocks.begin() + bblock, m_blocks.begin() + bblock,
                       [](Block b1, Block b2) { return b1 & b2; });
    }

    DynamicBitset& operator&=(const DynamicBitset& o) noexcept {
        assert(m_num_bits <= o.m_num_bits);
        std::transform(m_blocks.begin(), m_blocks.end(), o.m_blocks.begin(),
                       m_blocks.begin(),
                       [](Block b1, Block b2) { return b1 & b2; });
        return *this;
    }

    DynamicBitset& operator^=(const DynamicBitset& o) noexcept {
        assert(m_num_bits == o.m_num_bits);
        std::transform(m_blocks.begin(), m_blocks.end(), o.m_blocks.begin(),
                       m_blocks.begin(),
                       [](Block b1, Block b2) { return b1 ^ b2; });
        return *this;
    }

    DynamicBitset& operator-=(const DynamicBitset& o) noexcept {
        assert(m_num_bits == o.m_num_bits);
        std::transform(m_blocks.begin(), m_blocks.end(), o.m_blocks.begin(),
                       m_blocks.begin(),
                       [](Block b1, Block b2) { return b1 & ~b2; });
        return *this;
    }

    DynamicBitset operator~() const {
        DynamicBitset result(*this);
        result.flip();
        return result;
    }

    // Pseudo-reference
    class reference {
      public:
        reference(const reference&) noexcept = default;

        /* implicit */ operator bool() const noexcept { return *m_ptr & m_msk; }

        bool operator~() const noexcept { return !(*m_ptr & m_msk); }

        bool operator!() const noexcept { return ~*this; }

        reference& flip() noexcept {
            *m_ptr ^= m_msk;
            return *this;
        }

        reference& set() noexcept {
            *m_ptr |= m_msk;
            return *this;
        }

        reference& set(bool v) noexcept { return *this = v; }

        reference& reset() noexcept {
            *m_ptr &= ~m_msk;
            return *this;
        }

        reference& operator=(bool v) noexcept {
            if (v)
                *m_ptr |= m_msk;
            else
                *m_ptr &= ~m_msk;
            return *this;
        }

        reference& operator=(const reference& o) noexcept {
            return *this = static_cast<bool>(o);
        }

        reference& operator|=(bool x) noexcept {
            if (x)
                *this = true;
            return *this;
        }

        reference& operator&=(bool x) noexcept {
            if (!x)
                *this = false;
            return *this;
        }

        reference& operator^=(bool x) noexcept {
            if (x)
                this->flip();
            return *this;
        }

        reference& operator-=(bool x) noexcept {
            if (x)
                this->reset();
            return *this;
        }

      private:
        void operator&() = delete;

        reference(Block* ptr, Block msk) noexcept : m_ptr(ptr), m_msk(msk) {}

        Block* m_ptr;
        Block m_msk;

        friend class DynamicBitset;
    };

    reference operator[](std::size_t idx) noexcept {
        std::size_t blk_i = idx / bits_per_word;
        std::uint8_t sub_i(idx % bits_per_word);
        return reference{m_blocks.data() + blk_i, Block(1) << sub_i};
    }

    bool operator[](std::size_t idx) const noexcept {
        std::size_t blk_i = idx / bits_per_word;
        std::uint8_t sub_i(idx % bits_per_word);
        return m_blocks[blk_i] & (Block(1) << sub_i);
    }

    std::size_t size() const noexcept { return m_num_bits; }

    void reserve(size_type bits) {
        m_blocks.reserve(num_blocks_required(bits));
    }

    void resize(size_type num_bits, bool value = false) {
        auto nreq = num_blocks_required(num_bits);
        if (num_bits < m_num_bits) {
            m_blocks.resize(nreq);
        } else if (num_bits > m_num_bits) {
            m_blocks.reserve(nreq);
            if (value)
                p_one_unused();
            m_blocks.resize(nreq, broadcast(value));
        }
        m_num_bits = num_bits;
        p_zero_unused();
    }

    void assign(size_type num_bits, bool value) {
        if(num_bits != size()) resize(num_bits, value);
        set(value);
    }

    void clear() noexcept {
        m_num_bits = 0;
        m_blocks.clear();
    }

    void push_back(bool bit) {
        std::size_t new_bit_subind = m_num_bits % bits_per_word;
        if (!new_bit_subind) {
            m_blocks.push_back(bit);
        } else {
            Block msk = bit;
            msk <<= new_bit_subind;
            m_blocks.back() |= msk;
        }
        ++m_num_bits;
    }

    void pop_back() noexcept {
        std::uint8_t bits_used_in_last(--m_num_bits % bits_per_word);
        if (bits_used_in_last == 0) {
            m_blocks.pop_back();
        } else {
            m_blocks.back() &= (Block(1) << bits_used_in_last) - 1;
        }
    }

    void set(bool value = true) noexcept {
        std::fill(m_blocks.begin(), m_blocks.end(), broadcast(value));
        p_zero_unused();
    }

    void reset() noexcept { set(false); }

    void flip() noexcept {
        std::transform(m_blocks.begin(), m_blocks.end(), m_blocks.begin(),
                       [](Block b) { return ~b; });
        p_zero_unused();
    }

    bool any() const noexcept {
        return std::any_of(m_blocks.begin(), m_blocks.end(),
                           [](Block b) { return b != 0; });
    }

    bool none() const noexcept { return !any(); }

    bool all() const noexcept {
        auto l = std::find_if(m_blocks.begin(), m_blocks.end(),
                              [](Block b) { return b != ~Block(0); });
        if (l == m_blocks.end())
            return true;
        if (l != m_blocks.end() - 1)
            return false;
        Block last = m_blocks.back();
        std::uint8_t bits_used_in_last = m_num_bits % bits_per_word;
        if (!bits_used_in_last)
            return false;
        return (last | ~Block((Block(1) << bits_used_in_last) - 1)) ==
               ~Block(0);
    }

    std::size_t count() const noexcept {
        return std::transform_reduce(m_blocks.begin(), m_blocks.end(),
                                     std::size_t(0), std::plus<>{},
                                     [](Block b) { return count_set_bits(b); });
    }

    std::size_t count_from(std::size_t begin_index) const noexcept {
        std::size_t word_idx = begin_index / bits_per_word;
        std::size_t sub_idx = begin_index % bits_per_word;
        Block partial = m_blocks[word_idx];
        partial &= ~Block((Block(1) << sub_idx) - 1);
        std::size_t initial = count_set_bits(partial);
        return std::transform_reduce(m_blocks.begin() + word_idx + 1, m_blocks.end(),
                                     initial, std::plus<>{},
                                     [](Block b) { return count_set_bits(b); });
    }

    class OnesIterator {
      public:
        using iterator_category = std::forward_iterator_tag;
        using difference_type = std::ptrdiff_t;
        using reference = std::size_t;
        using pointer = const std::size_t*;
        using value_type = std::size_t;

        bool operator==(const OnesIterator& o) const noexcept {
            return m_pos == o.m_pos && m_cnt == o.m_cnt;
        }

        bool operator!=(const OnesIterator& o) const noexcept {
            return m_pos != o.m_pos || m_cnt != o.m_cnt;
        }

        OnesIterator operator++(int) const noexcept {
            OnesIterator result(*this);
            ++result;
            return result;
        }

        OnesIterator& operator++() noexcept {
            if (--m_cnt == 0) {
                do {
                    m_offs += bits_per_word;
                    if (++m_pos == m_end) {
                        m_cnt = 0;
                        return *this;
                    }
                    m_buf = *m_pos;
                } while (!m_buf);
                m_cnt = count_set_bits(m_buf);
            } else {
                m_buf &= m_buf - 1;
            }
            return *this;
        }

        std::size_t operator*() const noexcept {
            return count_trailing_zeros(m_buf) + m_offs;
        }

        OnesIterator() noexcept = default;

      private:
        friend class DynamicBitset;

        explicit OnesIterator(const Block* end) noexcept
            : m_pos(end), m_end(end), m_buf(0), m_cnt(0), m_offs(0) {}

        explicit OnesIterator(const Block* begin, const Block* end) noexcept
            : m_pos(begin), m_end(end), m_buf(0), m_cnt(0), m_offs(0) {
            p_scroll();
        }

        explicit OnesIterator(const Block* begin, const Block* end,
                              std::uint8_t index_in_word,
                              std::size_t offs) noexcept
            : m_pos(begin), m_end(end), m_buf(0), m_cnt(0), m_offs(offs) {
            if (m_pos != m_end) {
                m_buf = *m_pos;
                m_buf &= ~Block((Block(1) << index_in_word) - 1);
                if (m_buf) {
                    m_cnt = count_set_bits(m_buf);
                    return;
                }
                ++m_pos;
                m_offs += bits_per_word;
                p_scroll();
            }
        }

        void p_scroll() {
            while (m_pos != m_end) {
                m_buf = *m_pos;
                if (m_buf) {
                    m_cnt = count_set_bits(m_buf);
                    return;
                }
                ++m_pos;
                m_offs += bits_per_word;
            }
        }

        const Block* m_pos;
        const Block* m_end;
        Block m_buf;
        Block m_cnt;
        std::size_t m_offs;
    };

    OnesIterator ones_begin() const noexcept {
        return OnesIterator(m_blocks.data(), m_blocks.data() + m_blocks.size());
    }

    OnesIterator ones_end() const noexcept {
        return OnesIterator(m_blocks.data() + m_blocks.size());
    }

    IteratorRange<OnesIterator> ones() const noexcept {
        return {ones_begin(), ones_end()};
    }

    OnesIterator ones_from_begin(std::size_t begin_index) const noexcept {
        std::size_t word_idx = begin_index / bits_per_word;
        std::size_t sub_idx = begin_index % bits_per_word;
        std::size_t offs = word_idx * bits_per_word;
        return OnesIterator(m_blocks.data() + word_idx,
                            m_blocks.data() + m_blocks.size(), sub_idx, offs);
    }

    IteratorRange<OnesIterator>
    ones_from(std::size_t begin_index) const noexcept {
        return {ones_from_begin(begin_index), ones_end()};
    }

    std::size_t bytes_used() const noexcept {
        return m_blocks.capacity() * sizeof(Block) * CHAR_BIT / 8;
    }

    DynamicBitset(const DynamicBitset&) = default;
    DynamicBitset& operator=(const DynamicBitset&) = default;
    ~DynamicBitset() = default;

    DynamicBitset(DynamicBitset&& o) noexcept
        : m_blocks(std::move(o.m_blocks)), m_num_bits(o.m_num_bits) {
        o.m_blocks.clear();
        o.m_num_bits = 0;
    }

    DynamicBitset& operator=(DynamicBitset&& o) noexcept {
        std::swap(m_num_bits, o.m_num_bits);
        m_blocks.swap(o.m_blocks);
        return *this;
    }

    explicit operator std::vector<bool>() const {
        const std::size_t n = size();
        std::vector<bool> result(n, false);
        for(std::size_t i = 0; i < n; ++i) {
            result[i] = bool((*this)[i]);
        }
        return result;
    }

    const std::vector<Block>& blocks() const noexcept { return m_blocks; }

  private:
    template<typename T, std::enable_if_t<sizeof(T) == sizeof(Block), int> = 0>
    Block p_read_block(const T* from) {
        return Block(*from);
    }

    template<typename T, std::enable_if_t<sizeof(T) * 2 == sizeof(Block) && 
                                          sizeof(Block) == sizeof(std::uint64_t), int> = 0>
    Block p_read_block(const T* from) {
        Block b1(*from);
        Block b2(*(from + 1));
        return (b2 << 32) | b1;
    }

    void p_zero_unused() {
        std::uint8_t bits_used_in_last = m_num_bits % bits_per_word;
        if (!bits_used_in_last)
            return;
        m_blocks.back() &= (Block(1) << bits_used_in_last) - 1;
    }

    void p_one_unused() {
        std::uint8_t bits_used_in_last = m_num_bits % bits_per_word;
        if (!bits_used_in_last)
            return;
        m_blocks.back() |= ~Block((Block(1) << bits_used_in_last) - 1);
    }

    std::vector<Block> m_blocks;
    std::size_t m_num_bits;
};

using Bitset = DynamicBitset;

}

#endif
