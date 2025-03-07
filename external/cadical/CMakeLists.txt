try_compile(CDC_HAVE_FLEXIBLE_ARRAYS SOURCES "${CMAKE_CURRENT_LIST_DIR}/cadical_feature_tests/flexible.cpp")
try_compile(CDC_HAVE_UNLOCKED_IO SOURCES "${CMAKE_CURRENT_LIST_DIR}/cadical_feature_tests/unlocked_io.cpp")
try_compile(CDC_HAVE_CLOSEFROM SOURCES "${CMAKE_CURRENT_LIST_DIR}/cadical_feature_tests/closefrom.cpp")

message("Cadical -- flexible arrays: ${CDC_HAVE_FLEXIBLE_ARRAYS}")
message("Cadical -- unlocked io: ${CDC_HAVE_UNLOCKED_IO}")
message("Cadical -- closefrom: ${CDC_HAVE_CLOSEFROM}")

add_library(cadical STATIC 
  "cadical/src/analyze.cpp"
  "cadical/src/arena.cpp"
  "cadical/src/assume.cpp"
  "cadical/src/averages.cpp"
  "cadical/src/backtrack.cpp"
  "cadical/src/backward.cpp"
  "cadical/src/bins.cpp"
  "cadical/src/block.cpp"
  #  "cadical/src/cadical.cpp"
  "cadical/src/ccadical.cpp"
  "cadical/src/checker.cpp"
  "cadical/src/clause.cpp"
  "cadical/src/collect.cpp"
  "cadical/src/compact.cpp"
  "cadical/src/condition.cpp"
  "cadical/src/config.cpp"
  "cadical/src/constrain.cpp"
  "cadical/src/contract.cpp"
  "cadical/src/cover.cpp"
  "cadical/src/decide.cpp"
  "cadical/src/decompose.cpp"
  "cadical/src/deduplicate.cpp"
  "cadical/src/drattracer.cpp"
  "cadical/src/elim.cpp"
  "cadical/src/ema.cpp"
  "cadical/src/extend.cpp"
  "cadical/src/external.cpp"
  "cadical/src/external_propagate.cpp"
  "cadical/src/file.cpp"
  "cadical/src/flags.cpp"
  "cadical/src/flip.cpp"
  "cadical/src/format.cpp"
  "cadical/src/frattracer.cpp"
  "cadical/src/gates.cpp"
  "cadical/src/idruptracer.cpp"
  "cadical/src/instantiate.cpp"
  "cadical/src/internal.cpp"
  "cadical/src/ipasir.cpp"
  "cadical/src/lidruptracer.cpp"
  "cadical/src/limit.cpp"
  "cadical/src/logging.cpp"
  "cadical/src/lookahead.cpp"
  "cadical/src/lratbuilder.cpp"
  "cadical/src/lratchecker.cpp"
  "cadical/src/lrattracer.cpp"
  "cadical/src/lucky.cpp"
  "cadical/src/message.cpp"
  "cadical/src/minimize.cpp"
  "cadical/src/mobical.cpp"
  "cadical/src/occs.cpp"
  "cadical/src/options.cpp"
  "cadical/src/parse.cpp"
  "cadical/src/phases.cpp"
  "cadical/src/probe.cpp"
  "cadical/src/profile.cpp"
  "cadical/src/proof.cpp"
  "cadical/src/propagate.cpp"
  "cadical/src/queue.cpp"
  "cadical/src/random.cpp"
  "cadical/src/reap.cpp"
  "cadical/src/reduce.cpp"
  "cadical/src/rephase.cpp"
  "cadical/src/report.cpp"
  "cadical/src/resources.cpp"
  "cadical/src/restart.cpp"
  "cadical/src/restore.cpp"
  "cadical/src/score.cpp"
  "cadical/src/shrink.cpp"
  "cadical/src/signal.cpp"
  "cadical/src/solution.cpp"
  "cadical/src/solver.cpp"
  "cadical/src/stats.cpp"
  "cadical/src/subsume.cpp"
  "cadical/src/terminal.cpp"
  "cadical/src/ternary.cpp"
  "cadical/src/transred.cpp"
  "cadical/src/util.cpp"
  "cadical/src/var.cpp"
  "cadical/src/veripbtracer.cpp"
  "cadical/src/version.cpp"
  "cadical/src/vivify.cpp"
  "cadical/src/walk.cpp"
  "cadical/src/watch.cpp"
)
set_property(TARGET cadical PROPERTY POSITION_INDEPENDENT_CODE On)

if(NOT CDC_HAVE_FLEXIBLE_ARRAYS)
  target_compile_definitions(cadical PRIVATE "-DNFLEXIBLE")
endif()
if(NOT CDC_HAVE_UNLOCKED_IO)
  target_compile_definitions(cadical PRIVATE "-DNUNLOCKED")
endif()
if(NOT CDC_HAVE_CLOSEFROM)
  target_compile_definitions(cadical PRIVATE "-DNCLOSEFROM")
endif()
target_compile_definitions(cadical PRIVATE "-DNCONTRACTS" "-DNTRACING" "-DNCONTRIB" "-DNIPASIR")
target_include_directories(cadical PRIVATE "${CMAKE_CURRENT_LIST_DIR}/cadical_build_header")
target_compile_definitions(cadical PRIVATE COMPILER="${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")

add_library(cadical_binding STATIC "cadical_binding.cpp")
target_compile_definitions(cadical_binding PRIVATE "-DNCONTRACTS" "-DNTRACING" "-DNCONTRIB" "-DNIPASIR")
target_link_libraries(cadical_binding PUBLIC cadical)
target_include_directories(cadical_binding PRIVATE "${CMAKE_CURRENT_LIST_DIR}/cadical_build_header"
                                                   "${CMAKE_CURRENT_LIST_DIR}/cadical/src")

# include the public header directory in the source tree when building the lib
# include the include directory when using the installed lib
target_include_directories(cadical_binding PUBLIC 
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include> 
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)    

if(LTO_SUPPORTED)
  set_property(TARGET cadical PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
  set_property(TARGET cadical PROPERTY INTERPROCEDURAL_OPTIMIZATION_DEBUG FALSE)
  set_property(TARGET cadical_binding PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
  set_property(TARGET cadical_binding PROPERTY INTERPROCEDURAL_OPTIMIZATION_DEBUG FALSE)
endif()
