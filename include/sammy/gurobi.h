#ifndef SAMMY_GUROBI_H_INCLUDED_
#define SAMMY_GUROBI_H_INCLUDED_

#include <gurobi_c++.h>
#include <mutex>

namespace sammy {

namespace detail {

inline void start_environment(GRBEnv& env, bool quiet) {
    static std::mutex start_environment_lock;
    std::unique_lock<std::mutex> lock{start_environment_lock};
    if (quiet) {
        env.set(GRB_IntParam_OutputFlag, 0);
        env.set(GRB_IntParam_LogToConsole, 0);
    }
    env.start();
}

} // namespace detail

inline GRBEnv& gurobi_environment(bool quiet) {
    try {
        thread_local bool initialized = false;
        thread_local GRBEnv environment{true};
        if (!initialized) {
            detail::start_environment(environment, quiet);
            initialized = true;
        }
        return environment;
    } catch (const GRBException& ex) {
        std::stringstream out;
        out << "Gurobi exception during initialization, likely "
               "license-related: "
            << ex.getMessage();
        throw std::runtime_error(out.str());
    }
}

} // namespace sammy

#endif
