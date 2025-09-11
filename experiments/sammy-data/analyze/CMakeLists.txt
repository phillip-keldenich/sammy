cmake_minimum_required(VERSION 3.24) # 3.17 has been released March 2020

project(sample_analyzer CXX) # Feel free to add more project information
set(CMAKE_CXX_STANDARD 17) # Using C++-17 (most features with good support at
                           # 2022)
set(CMAKE_POSITION_INDEPENDENT_CODE ON) # The code needs to be compiled as PIC
                                        # to build the shared lib for python.
set(CMAKE_OSX_ARCHITECTURES x86_64)

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
endif()

# include(./cmake/CCache.cmake) # Magic to speed up compiling by caching.

# ~~~
# CPM Dependencies ~~~~~~~~~~~~~~~~
# CPM allows us to include some simple dependencies without much fuss.
# Great default, for everything more complicated, uses conan.
# Why not always use conan? CPM works without any interaction (conan needs
# to be called externally). CPM also has no problems with install PyBind11
# (which conan has at the time).
# If the dependency has complicated dependencies on its own or is slow to
# compile, better use conan.
# Check out https://github.com/cpm-cmake/CPM.cmake
include(./cmake/CPM.cmake) # Package manager for simple requirements.

add_subdirectory(src)

# ~~~
# Apps ~~~~~~~~~~~~~~~~~~~~~~~~~~~
if(NOT SKBUILD) # only make available if not using SKBUILD.
  # define your apps here
endif()