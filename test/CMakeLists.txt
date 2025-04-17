find_package(doctest REQUIRED)

add_library(sammy_test_instances STATIC sammy_test_instance_soletta.cpp)
target_link_libraries(sammy_test_instances PUBLIC sammy)

# Option to disable tests if they can't be built due to system compatibility issues
option(SAMMY_SKIP_TESTS "Skip building tests if there are compatibility issues" OFF)

if(NOT SAMMY_SKIP_TESTS)
    add_executable(
        sammy_test
        test_main.cpp
        test_eliminate_subsumed.cpp
        test_detect_equalities.cpp
        test_bounded_variable_elimination.cpp
        test_initial_heuristic.cpp
        test_bitset_operations.cpp
        test_universe_subgraph.cpp
        test_partial_solution.cpp
        test_propagator.cpp
        test_satdsatur.cpp
        test_clique_storage.cpp
        test_implied_vertices.cpp
        test_clique_sat_dsatur.cpp
        test_cadical_binding.cpp
        test_cryptominisat_binding.cpp
        test_subproblem_solver_with_mes.cpp
        test_external_sat.cpp
    )

    if(UNIX)
        target_sources(sammy_test PRIVATE "../src/process_unix.cpp")
    endif()

    if(SAMMY_USE_CUDA)
        target_sources(sammy_test PRIVATE test_cuda.cpp)
    endif()

    if(EXISTS "${CMAKE_CURRENT_LIST_DIR}/sammy_test_data")
        target_compile_definitions(sammy_test PRIVATE "-DSAMMY_TEST_DATA=\"${CMAKE_CURRENT_LIST_DIR}/sammy_test_data\"")
    endif()

    target_link_libraries(sammy_test PRIVATE sammy doctest::doctest sammy_test_instances)

    if(LTO_SUPPORTED)
        set_property(TARGET sammy_test PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
        set_property(TARGET sammy_test PROPERTY INTERPROCEDURAL_OPTIMIZATION_DEBUG FALSE)
    endif()
    
    add_test(NAME run_sammy_test
             WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
             COMMAND sammy_test)
    
    message(STATUS "Tests enabled - set SAMMY_SKIP_TESTS=ON if you encounter linking errors")
else()
    message(STATUS "Tests disabled - SAMMY_SKIP_TESTS is ON")
endif()
