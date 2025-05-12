#ifndef SAMMY_BARRAGE_LNS_SUBPROBLEM_H_INCLUDED_
#define SAMMY_BARRAGE_LNS_SUBPROBLEM_H_INCLUDED_

#include "dynamic_bitset.h"
#include "literals.h"
#include "output.h"
#include <memory>

namespace sammy {

/**
 * Structure for general representation of a subproblem
 * as encountered by our LNS.
 */
struct LNSSubproblem {
    std::vector<Vertex> uncovered_universe;
    std::vector<Vertex> mutually_exclusive_set;
    std::vector<DynamicBitset> removed_configurations;
    Lit num_concrete;

    /**
     * Dump to JSON output object,
     * to be used for internal serialization.
     * Used for debugging and testing purposes.
     */
    OutputObject dump() const {
        std::vector<std::vector<bool>> removed;
        for (const auto& config : removed_configurations) {
            removed.push_back(static_cast<std::vector<bool>>(config));
        }
        return OutputObject{{"uncovered", uncovered_universe},
                            {"mes", mutually_exclusive_set},
                            {"removed_configs", removed},
                            {"num_concrete", num_concrete}};
    }

    /**
     * Load from JSON output object,
     * to be used for internal deserialization.
     * Used for debugging and testing purposes.
     */
    static LNSSubproblem load(const OutputObject& obj) {
        std::vector<Vertex> uncovered =
            obj.at("uncovered").get<std::vector<Vertex>>();
        std::vector<Vertex> mes = obj.at("mes").get<std::vector<Vertex>>();
        std::vector<DynamicBitset> removed_configs;
        for (const auto& config : obj.at("removed_configs")) {
            removed_configs.emplace_back(config.get<std::vector<bool>>());
        }
        Lit num_concrete = obj.at("num_concrete").get<Lit>();
        return LNSSubproblem{std::move(uncovered), std::move(mes),
                             std::move(removed_configs), num_concrete};
    }
};

} // namespace sammy

#endif
