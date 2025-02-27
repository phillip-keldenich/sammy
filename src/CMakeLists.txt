add_executable(sammy_solve solve.cpp)
target_link_libraries(sammy_solve PRIVATE sammy)

add_executable(generate_subproblem_mes_entry_points generate_subproblem_mes_entry_points.cpp)
target_link_libraries(generate_subproblem_mes_entry_points PRIVATE sammy)

if(LTO_SUPPORTED)
  set_property(TARGET PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
  set_property(TARGET PROPERTY INTERPROCEDURAL_OPTIMIZATION_DEBUG FALSE)
endif()

