cmake_minimum_required(VERSION 3.21..3.31)

# detect conan toolchain files, if they are present
if(NOT DEFINED CMAKE_TOOLCHAIN_FILE) 
  if(EXISTS "${CMAKE_BINARY_DIR}/generators/conan_toolchain.cmake")
    set(CMAKE_TOOLCHAIN_FILE "${CMAKE_BINARY_DIR}/generators/conan_toolchain.cmake")
  elseif(EXISTS "${CMAKE_BINARY_DIR}/conan_toolchain.cmake")
    set(CMAKE_TOOLCHAIN_FILE "${CMAKE_BINARY_DIR}/conan_toolchain.cmake")
  endif()
endif()

project(sammy C CXX)

# options
option(SAMMY_NO_ORTOOLS "Do not use ortools" OFF)
option(SAMMY_NO_CUDA "Disable CUDA support; only relevant on CUDA platforms" OFF)
option(SAMMY_NO_LTO "Do not attempt link-time optimization" Off)

# check for CUDA support
if(NOT SAMMY_NO_CUDA)
    include(CheckLanguage)
    check_language(CUDA)
    if(CMAKE_CUDA_COMPILER)
      message(STATUS "CUDA compiler found, enabling CUDA support")
      set(SAMMY_USE_CUDA ON)
      enable_language(CUDA)
    else()
      message(STATUS "No CUDA compiler found, disabling CUDA support")
    endif()
endif()

# check for link-time optimization support
if(NOT SAMMY_NO_LTO)
	include(CheckIPOSupported)
  check_ipo_supported(RESULT LTO_SUPPORTED OUTPUT LTO_OUTPUT_DATA LANGUAGES C CXX)
else()
	set(LTO_SUPPORTED Off)
endif()

# find packages
find_package(Boost REQUIRED COMPONENTS iostreams program_options)
find_package(gurobi CONFIG REQUIRED)
find_package(nlohmann_json CONFIG REQUIRED)

if(NOT SAMMY_NO_ORTOOLS)
  find_package(ortools CONFIG REQUIRED)
endif()

if(SAMMY_USE_CUDA)
    add_library(sammycuda STATIC "src/cuda_iteration.cu")
    target_compile_features(sammycuda PUBLIC cxx_std_17)
    target_include_directories(sammycuda PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
    target_link_libraries(sammycuda PUBLIC nlohmann_json::nlohmann_json Boost::boost
	                        Boost::iostreams Boost::program_options)
    target_include_directories(sammy INTERFACE "${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}")
    target_compile_definitions(sammycuda PUBLIC "-DSAMMY_CUDA_SUPPORTED=1")
    target_link_libraries(sammy INTERFACE sammycuda)
endif()

# setup interface library target (sammy)
include(GNUInstallDirs)
add_library(sammy INTERFACE)
target_compile_features(sammy INTERFACE cxx_std_17)
target_include_directories(sammy INTERFACE 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include> # include the public header directory in the source tree when building the lib
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)    # include the include directory when using the installed lib
target_link_libraries(sammy INTERFACE nlohmann_json::nlohmann_json Boost::boost 
                                      Boost::iostreams Boost::program_options
                                      gurobi::gurobi)

add_subdirectory("external")

target_link_libraries(sammy INTERFACE cryptominisat5_binding cryptominisat5 
                                      kissat_binding kissat
                                      cadical_binding cadical
                                      lingeling_binding lingeling)

add_subdirectory("src")
enable_testing()
add_subdirectory("test")

