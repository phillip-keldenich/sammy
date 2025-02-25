#include <doctest/doctest.h>
#include <sammy/clique_storage.h>

using namespace sammy;

TEST_CASE("[CliqueStorage] Storage with few cliques") {
    CliqueStorage storage(1000);
    auto view = storage.obtain_view();
    CHECK(view.size() == 0);
    for(const auto& r : view) {
        CHECK(r.size() == 0);
        CHECK(false);
    }
    std::vector<Vertex> vertices{{1,2}, {2,3}, {3,4}, {4,5}, {5,6}};
    storage.push_clique(vertices.begin(), vertices.end());
    CHECK(view.size() == 0);
    view = storage.obtain_view();
    CHECK(view.size() == 1);
    CHECK(view[0].size() == 5);
    CHECK(view[0].begin()->first == 1);
    CHECK(view[0].begin()->second == 2);
    CHECK(std::distance(view.begin(), view.end()) == 1);
    for(const auto& r : view) {
        CHECK(r.size() == 5);
        for(const auto& [u, v] : r) {
            CHECK(u <= 5);
            CHECK(u == v - 1);
        }
    }
    storage.push_clique(vertices.begin(), vertices.end());
    CHECK(view.size() == 1);
    view = storage.obtain_view();
    CHECK(view.size() == 2);
    CHECK(std::distance(view.begin(), view.end()) == 2);
    CHECK(view[0].size() == 5);
    CHECK(view[1].size() == 5);
    for(const auto& r : view) {
        CHECK(r.size() == 5);
        for(const auto& [u, v] : r) {
            CHECK(u <= 5);
            CHECK(u == v - 1);
        }
    }
}
