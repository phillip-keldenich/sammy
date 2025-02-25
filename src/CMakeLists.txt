add_executable(sammy_solve solve.cpp)
target_link_libraries(sammy_solve PRIVATE sammy)

#add_executable(run_subproblem_solve run_subproblem_solve.cpp)
#target_link_libraries(run_subproblem_solve PRIVATE sammy nlohmann_json::nlohmann_json
#                      kissat_binding kissat gurobi::gurobi lingeling_binding lingeling
#                      cryptominisat5_binding cadical_binding cadical)

#add_executable(generate_subproblem_mes_entry_points generate_subproblem_mes_entry_points.cpp)
#target_link_libraries(generate_subproblem_mes_entry_points PRIVATE sammy)

if(LTO_SUPPORTED)
  set_property(TARGET PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
  set_property(TARGET PROPERTY INTERPROCEDURAL_OPTIMIZATION_DEBUG FALSE)
endif()

