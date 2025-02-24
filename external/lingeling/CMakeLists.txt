set(LINGELING_LIB_SOURCES "lglib.c" "lglopts.c" "lgldimacs.c")
list(TRANSFORM LINGELING_LIB_SOURCES PREPEND "${CMAKE_CURRENT_LIST_DIR}/lingeling/")
add_library(lingeling STATIC ${LINGELING_LIB_SOURCES})
target_compile_definitions(lingeling PRIVATE NLGLYALSAT=1 NLGLDRUPLIG=1 NLGLOG=1 NCHKSOL=1 NLGLFILES=1 NLGLDEMA=1)
add_library(lingeling_binding STATIC "lingeling_binding.cpp")
target_link_libraries(lingeling_binding PUBLIC lingeling)

if(LTO_SUPPORTED)
set_property(TARGET lingeling PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
set_property(TARGET lingeling_binding PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
set_property(TARGET lingeling PROPERTY INTERPROCEDURAL_OPTIMIZATION_DEBUG FALSE)
set_property(TARGET lingeling_binding PROPERTY INTERPROCEDURAL_OPTIMIZATION_DEBUG FALSE)
endif()

target_include_directories(lingeling PRIVATE "${CMAKE_CURRENT_LIST_DIR}/lingeling")
target_include_directories(lingeling_binding PRIVATE "${CMAKE_CURRENT_LIST_DIR}/lingeling")
target_include_directories(lingeling_binding PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include> # include the public header directory in the source tree when building the lib
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)    # include the include directory when using the installed lib
target_compile_features(lingeling_binding PUBLIC cxx_std_17)

