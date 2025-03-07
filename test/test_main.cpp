#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest/doctest.h>
#include <sammy/external_sat_solver.h>

int main(int argc, char** argv) {
    auto sat_only = sammy::is_sat_only_call(argc, argv);
    if (sat_only) {
        return sammy::sat_only_entry_point(*sat_only);
    }

    doctest::Context context;
    context.applyCommandLine(argc, argv);
    int res = context.run(); // run
    if (context.shouldExit()) {
        // important - query flags (and --exit)
        // rely on the user doing this
        return res; // propagate the result of the tests
    }
    return res;
}
