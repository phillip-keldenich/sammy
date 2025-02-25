#include <sammy/parallel_bit_filter.h>
#include <sammy/rng.h>
#include <doctest/doctest.h>
#include <random>
#include <iostream>

using namespace sammy;

static std::size_t generate_random_size(std::size_t min_size, std::size_t max_size) {
    auto& rng = sammy::rng();
    return std::uniform_int_distribution<std::size_t>(min_size, max_size)(rng);
}

static std::pair<DynamicBitset, std::vector<bool>> generate_random_set(std::size_t s, double p) {
    DynamicBitset result(s, false);
    std::vector<bool> result2(s, false);
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    auto& rng = sammy::rng();
    for(std::size_t i = 0; i < s; ++i) {
        if(dist(rng) < p) {
            result[i] = true;
            result2[i] = true;
        }
    }
    return {std::move(result), std::move(result2)};
}

bool same_set(const DynamicBitset& bs1, const std::vector<bool>& bs2) {
    if(bs1.size() != bs2.size()) return false;
    for(std::size_t i = 0, n = bs1.size(); i < n; ++i) {
        if(bool(bs1[i]) != bs2[i]) return false;
    }
    return true;
}

bool same_set(const DynamicBitset& bs1, const std::vector<std::size_t>& s2) {
    if(bs1.count() != s2.size()) return false;
    for(std::size_t x : s2) {
        if(!bs1[x]) return false;
    }
    return true;
}

TEST_CASE("[DynamicBitset] Iterating set bits") {
    DynamicBitset ones(12234, true);
    std::size_t expect_next = 0;
    for(std::size_t i : ones.ones()) {
        CHECK(expect_next == i);
        ++expect_next;
    }
    CHECK(expect_next == 12234);
    DynamicBitset zeros(1255, false);
    for(std::size_t i : zeros.ones()) {
        CHECK(i != i);
    }
    ones.resize(12667, false);
    expect_next = 0;
    for(std::size_t i : ones.ones()) {
        CHECK(expect_next == i);
        ++expect_next;
    }
    CHECK(expect_next == 12234);
}

TEST_CASE("[DynamicBitset] Offset bulk operations") {
    DynamicBitset s1(12234, false);
    DynamicBitset s2(12234, false);
    s1[0] = true;
    s1[1] = true;
    s1[65] = true;
    s1[76] = true;
    s1[79] = true;
    s1[1022] = true;
    s1[1025] = true;
    s1[12232] = true;
    s1[12233] = true;

    s2[0] = true;
    s2[1] = true;
    s2[65] = true;
    s2[76] = true;
    s2[78] = true;
    s2[1022] = true;
    s2[1026] = true;
    s2[12231] = true;
    s2[12232] = true;

    SUBCASE("binary_subtract0") {
        DynamicBitset tmp = s1;
        tmp.binary_subtract(s2, 0);
        CHECK(same_set(tmp, std::vector<std::size_t>{79, 1025, 12233}));
        tmp = s2;
        tmp.binary_subtract(s1, 0);
        CHECK(same_set(tmp, std::vector<std::size_t>{78, 1026, 12231}));
    }

    SUBCASE("binary_subtract76") {
        DynamicBitset tmp = s1;
        tmp.binary_subtract(s2, 76);
        CHECK(same_set(tmp, std::vector<std::size_t>{0, 1, 65, 79, 1025, 12233}));
        tmp = s2;
        tmp.binary_subtract(s1, 76);
        CHECK(same_set(tmp, std::vector<std::size_t>{0, 1, 65, 78, 1026, 12231}));
    }

	SUBCASE("binary_subtract_uint32_t") {
		std::uint32_t buffer[384];
		std::fill_n(+buffer, 384, 0);
		buffer[0] = 1;
		buffer[1] = 4495;
		buffer[2] = 2;
		buffer[119] = 995572;
		buffer[382] = (1 << 7);
		DynamicBitset tmp = s1;
		tmp.binary_subtract(buffer);
		CHECK(same_set(tmp, std::vector<std::size_t>{1, 76, 79, 1022, 1025, 12232, 12233}));
		tmp = s2;
		tmp.binary_subtract(buffer);
		CHECK(same_set(tmp, std::vector<std::size_t>{1, 76, 78, 1022, 1026, 12232}));
	}
}

TEST_CASE("[DynamicBitset] Unparallelized bitset bulk operations") {
    for(int i = 0; i < 100; ++i) {
        std::size_t rminsize = 512 * 1024;
        std::size_t rmaxsize = 2048 * 1024;
        std::size_t rsize = generate_random_size(rminsize, rmaxsize);
        auto [s1, s2] = generate_random_set(rsize, 0.5);
        auto [t1, t2] = generate_random_set(rsize, 0.5);
        CHECK(same_set(s1, s2));
        CHECK(same_set(t1, t2));
        auto c1 = s1;
        auto c2 = s2;
        CHECK(same_set(c1, c2));
        c1 &= t1;
        std::transform(c2.begin(), c2.end(), t2.begin(), c2.begin(), [] (bool b1, bool b2) { return b1 & b2; });
        CHECK(same_set(c1, c2));
        c1 = s1;
        c2 = s2;
        CHECK(same_set(c1, c2));
        c1 |= t1;
        std::transform(c2.begin(), c2.end(), t2.begin(), c2.begin(), [] (bool b1, bool b2) { return b1 | b2; });
        CHECK(same_set(c1, c2));
        c1 = s1;
        c2 = s2;
        CHECK(same_set(c1, c2));
        c1 ^= t1;
        std::transform(c2.begin(), c2.end(), t2.begin(), c2.begin(), [] (bool b1, bool b2) { return b1 ^ b2; });
        CHECK(same_set(c1, c2));
        c1 = s1;
        c2 = s2;
        CHECK(same_set(c1, c2));
        c1 -= t1;
        std::transform(c2.begin(), c2.end(), t2.begin(), c2.begin(), [] (bool b1, bool b2) { return b1 & !b2; });
        CHECK(same_set(c1, c2));
        c1.assign(1000, false);
        c2.assign(1000, false);
        CHECK(same_set(c1, c2));
        c1.assign(1000, true);
        c2.assign(1000, true);
        CHECK(same_set(c1, c2));
    }
}

TEST_CASE("[BitsetOperationsBuffer] Parallel bitset bulk operations") {
    ThreadGroup<void> tgroup;
    BitsetOperationsBuffer buffer{&tgroup};
    std::size_t num_sets = tgroup.num_threads() * 5 + 3;
    for(int i = 0; i < 100; ++i) {
        std::size_t rminsize = 512 * 1024;
        std::size_t rmaxsize = 2048 * 1024;
        std::size_t rsize = generate_random_size(rminsize, rmaxsize);
        std::vector<Bitset> sets;
        sets.reserve(num_sets);
        sets.emplace_back(rsize, true);
        for(std::size_t i = 0; i < num_sets; ++i) {
            sets.push_back(generate_random_set(rsize, 0.05).first);
        }
        Bitset initial_and_result(rsize, true);
        bitwise_filter(buffer, initial_and_result, sets.begin(), sets.end());
        for(std::size_t i = 0; i < rsize; ++i) {
            if(initial_and_result[i]) {
                CHECK(!std::any_of(sets.begin(), sets.end(), [i] (const Bitset& b) { return bool(b[i]); }));
            } else {
                CHECK(std::any_of(sets.begin(), sets.end(), [i] (const Bitset& b) { return bool(b[i]); }));
            }
        }
    }
}

