#ifndef SAMMY_GLOBAL_STAT_INFO_H_INCLUDED_
#define SAMMY_GLOBAL_STAT_INFO_H_INCLUDED_

#include "literals.h"

namespace sammy {

struct GlobalStatInfo {
    mutable std::mutex lock;
    HashMap<std::string, double> double_stats;
    HashMap<std::string, std::int64_t> int_stats;

    void double_stat_add(const std::string& name, double value) {
        std::unique_lock l{lock};
        auto it = double_stats.find(name);
        if (it == double_stats.end()) {
            double_stats[name] = value;
        } else {
            it->second += value;
        }
    }

    void int_stat_add(const std::string& name, std::int64_t value) {
        std::unique_lock l{lock};
        auto it = int_stats.find(name);
        if (it == int_stats.end()) {
            int_stats[name] = value;
        } else {
            it->second += value;
        }
    }
};

inline std::ostream& operator<<(std::ostream& o, const GlobalStatInfo& l) {
    std::unique_lock lck{l.lock};
    o << "GlobalStatInfo{\n";
    for (const auto& [k, v] : l.double_stats) {
        o << "\t" << k << ": " << v << ",\n";
    }
    for (const auto& [k, v] : l.int_stats) {
        o << "\t" << k << ": " << v << ",\n";
    }
    o << "}";
    return o;
}

inline GlobalStatInfo& get_global_stats() {
    static GlobalStatInfo stats;
    return stats;
}

} // namespace sammy

#endif
