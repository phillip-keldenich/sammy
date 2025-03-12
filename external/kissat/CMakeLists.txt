set(KISSAT_LIB_SOURCES "allocate.c" "analyze.c"
                       "ands.c" "arena.c" "assign.c"
                       "averages.c" "backbone.c" "backtrack.c"
                       "build.c" "bump.c" "check.c" "classify.c" "clause.c"
                       "collect.c" "colors.c" "compact.c"
                       "config.c" "congruence.c" "decide.c" "deduce.c"
                       "definition.c" "dense.c" "dump.c"
                       "eliminate.c" "equivalences.c"
                       "error.c" "extend.c" "factor.c" "fastel.c" "file.c"
                       "flags.c" "format.c" "forward.c"
                       "gates.c" "heap.c" "ifthenelse.c"
                       "import.c" "internal.c" "kimits.c"
                       "kitten.c" "learn.c" "logging.c" "lucky.c"
                       "main.c" "minimize.c" "mode.c"
                       "options.c" "phases.c" "preprocess.c" "print.c"
                       "probe.c" "profile.c" "promote.c"
                       "proof.c" "propbeyond.c" "propdense.c" "propinitially.c"
                       "proprobe.c" "propsearch.c" "queue.c"
                       "reduce.c" "reluctant.c" "rephase.c"
                       "report.c" "reorder.c" "resize.c" "resolve.c"
                       "resources.c" "restart.c" "search.c"
                       "shrink.c" "smooth.c" "sort.c" "stack.c"
                       "statistics.c" "strengthen.c" "substitute.c"
                       "sweep.c" "terminate.c" "tiers.c" "trail.c" "transitive.c" "utilities.c"
                       "vector.c" "vivify.c" "walk.c" "warmup.c"
                       "watch.c" "weaken.c")
file(GENERATE OUTPUT "${CMAKE_BINARY_DIR}/kissat/build.h"
     CONTENT "#define VERSION \"0\"\n#define ID \"0\"\n#define COMPILER \"${CMAKE_C_COMPILER_ID}\"\n#define BUILD \"0\"\n")
list(TRANSFORM KISSAT_LIB_SOURCES PREPEND "${CMAKE_CURRENT_LIST_DIR}/kissat/src/")
add_library(kissat STATIC ${KISSAT_LIB_SOURCES})
target_compile_definitions(kissat PRIVATE QUIET=1)
add_library(kissat_binding STATIC "kissat_binding.cpp")

# kissat is currently broken with LTO;
# see https://github.com/arminbiere/kissat/issues/65
set_property(TARGET kissat PROPERTY INTERPROCEDURAL_OPTIMIZATION FALSE)
set_property(TARGET kissat_binding PROPERTY INTERPROCEDURAL_OPTIMIZATION FALSE)
#set_property(TARGET kissat PROPERTY INTERPROCEDURAL_OPTIMIZATION_DEBUG FALSE)
#set_property(TARGET kissat_binding PROPERTY INTERPROCEDURAL_OPTIMIZATION_DEBUG FALSE)

target_include_directories(kissat PRIVATE "${CMAKE_CURRENT_LIST_DIR}/kissat/src")
target_include_directories(kissat PRIVATE "${CMAKE_BINARY_DIR}/kissat")
target_include_directories(kissat_binding PRIVATE "${CMAKE_CURRENT_LIST_DIR}/kissat/src")
target_link_libraries(kissat_binding PUBLIC kissat)
target_include_directories(kissat_binding PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include> # include the public header directory in the source tree when building the lib
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)    # include the include directory when using the installed lib
target_compile_features(kissat_binding PUBLIC cxx_std_17)
