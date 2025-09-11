cpmaddpackage("gh:pybind/pybind11@2.13.6") # pybind11, essential


pybind11_add_module(
  _verify # target name with leading `_` as this will be a protected
  # module.
  sample_analyzer/verify/_verify.cpp) # the file with the API.
# enable compilation warnings
target_compile_options(
  _verify PRIVATE "$<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall>")
# This definition is important for skbuild to copy the binary to the right
# position in the package. !!! This part is critical and a wrong configuration
# leads to bad packages !!! without an error. Use a path relative to the project
# root (or !!! `cmake_install_dir`[setup.py], if used).
install(TARGETS _verify DESTINATION ./src/sample_analyzer/verify/)