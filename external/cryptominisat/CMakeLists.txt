include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-Wno-deprecated-declarations" HAVE_DEPRECATION_WARNING_FLAG)
check_cxx_compiler_flag("-Wno-shadow" HAVE_SHADOW_WARNING_FLAG)
check_cxx_compiler_flag("-Wno-unused-variable" HAVE_UNUSED_VARIABLE_FLAG)
check_cxx_compiler_flag("-Wno-unused-but-set-variable" HAVE_UNUSED_BUT_SET_VARIABLE)
check_cxx_compiler_flag("-Wno-mismatched-tags" HAVE_NO_MISMATCHED_TAGS)

set(CMSAT5_SOURCE_NAMES 
    "cnf.cpp" "frat.cpp" "propengine.cpp" "varreplacer.cpp"
    "clausecleaner.cpp" "occsimplifier.cpp" "gatefinder.cpp" "subsumestrengthen.cpp"
    "clauseallocator.cpp" "sccfinder.cpp" "solverconf.cpp" "distillerlong.cpp"
    "distillerlitrem.cpp" "distillerbin.cpp" "distillerlongwithimpl.cpp"
    "str_impl_w_impl.cpp" "solutionextender.cpp" "completedetachreattacher.cpp"
    "searcher.cpp" "solver.cpp" "hyperengine.cpp" "subsumeimplicit.cpp"
    "datasync.cpp" "reducedb.cpp" "bva.cpp" "intree.cpp" "searchstats.cpp"
    "xorfinder.cpp" "cardfinder.cpp" "cryptominisat_c.cpp" "sls.cpp"
    "sqlstats.cpp" "vardistgen.cpp" "ccnr.cpp" "ccnr_cms.cpp"
    "lucky.cpp" "get_clause_query.cpp" "gaussian.cpp" "packedrow.cpp"
    "matrixfinder.cpp" "GitSHA1.cpp" "cryptominisat.cpp" "oracle/oracle.cpp"
    "picosat/picosat.c" "picosat/version.c"
)

# actual solver
list(TRANSFORM CMSAT5_SOURCE_NAMES PREPEND "${CMAKE_CURRENT_LIST_DIR}/cryptominisat5/src/")
add_library(cryptominisat5 STATIC ${CMSAT5_SOURCE_NAMES})
target_compile_features(cryptominisat5 PUBLIC cxx_std_17)

# binding library
add_library(cryptominisat5_binding STATIC "cmsat5_binding.cpp")
target_compile_features(cryptominisat5_binding PUBLIC cxx_std_17)
target_include_directories(cryptominisat5 PRIVATE "${CMAKE_CURRENT_LIST_DIR}/cryptominisat5/src")
target_include_directories(cryptominisat5_binding PRIVATE "${CMAKE_CURRENT_LIST_DIR}/cryptominisat5/src")
target_link_libraries(cryptominisat5_binding PUBLIC cryptominisat5)
target_include_directories(cryptominisat5_binding PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include> # include the public header directory in the source tree when building the lib
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)    # include the include directory when using the installed lib

if(HAVE_DEPRECATION_WARNING_FLAG)
    target_compile_options(cryptominisat5 PRIVATE "-Wno-deprecated-declarations")
    target_compile_options(cryptominisat5_binding PRIVATE "-Wno-deprecated-declarations")
endif()
if(HAVE_SHADOW_WARNING_FLAG)
    target_compile_options(cryptominisat5 PRIVATE "-Wno-shadow")
    target_compile_options(cryptominisat5_binding PRIVATE "-Wno-deprecated-declarations")
endif()
if(HAVE_UNUSED_VARIABLE_FLAG)
    target_compile_options(cryptominisat5 PRIVATE "-Wno-unused-variable")
    target_compile_options(cryptominisat5_binding PRIVATE "-Wno-unused-variable")
endif()
if(HAVE_UNUSED_BUT_SET_VARIABLE)
    target_compile_options(cryptominisat5 PRIVATE "-Wno-unused-but-set-variable")
    target_compile_options(cryptominisat5_binding PRIVATE "-Wno-unused-but-set-variable")
endif()
if(HAVE_NO_MISMATCHED_TAGS)
    target_compile_options(cryptominisat5 PRIVATE "-Wno-mismatched-tags")
    target_compile_options(cryptominisat5_binding PRIVATE "-Wno-mismatched-tags")
endif()
if(LTO_SUPPORTED)
    set_property(TARGET cryptominisat5 PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    set_property(TARGET cryptominisat5_binding PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    set_property(TARGET cryptominisat5 PROPERTY INTERPROCEDURAL_OPTIMIZATION_DEBUG FALSE)
    set_property(TARGET cryptominisat5_binding PROPERTY INTERPROCEDURAL_OPTIMIZATION_DEBUG FALSE)
endif()

