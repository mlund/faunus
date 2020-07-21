cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

project(faunus)
include(cmake/StandardProjectSettings.cmake)
#include(cmake/PreventInSourceBuilds.cmake)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_EXTENSIONS YES)

# Link this 'library' to use the warnings specified in CompilerWarnings.cmake
add_library(project_options INTERFACE)
target_compile_features(project_options INTERFACE cxx_std_17)
add_library(project_warnings INTERFACE)
include(cmake/CompilerWarnings.cmake)
set_project_warnings(project_warnings)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(PythonInterp 3 REQUIRED)
find_package(Threads REQUIRED)
find_package(ZLIB REQUIRED)
enable_testing()

# Very basic PCH example
option(ENABLE_PCH "Enable Precompiled Headers" ON)
if(ENABLE_PCH)
  # This sets a global PCH parameter, each project will build its own PCH, which is a good idea if any #define's change
  #
  # consider breaking this out per project as necessary
  target_precompile_headers(
    project_options
    INTERFACE
    <vector>
    <string>
    <map>
    <set>
    <algorithm>
    <utility>
    <Eigen/Core>
    <Eigen/Geometry>
    <range/v3/view/filter.hpp>
    <range/v3/view/transform.hpp>
    <cereal/cereal.hpp>
    <cereal/types/base_class.hpp>
    <cereal/types/memory.hpp>
    <nlohmann/json.hpp>)
endif()

option(BUILD_STATIC "Build statically linked binaries" off)
option(ENABLE_FREESASA "Fetch 3rd-party SASA calculation software" on)
option(ENABLE_PYTHON "Try to compile python bindings (experimental!)" on)
option(ENABLE_LTO "Enable Link-Time Optimizations if available" off)
option(ENABLE_OPENMP "Try to use OpenMP parallisation" on)
option(ENABLE_MPI "Enable MPI code" off)

if (ENABLE_OPENMP)
  find_package(OpenMP)
  if (OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  endif()
endif()

if (ENABLE_MPI)
    find_package(MPI REQUIRED)
endif()

include(cmake/compilerflags.cmake)
include(cmake/externaltools.cmake)
include(cmake/gitversion.cmake)
 
add_subdirectory(src)
add_subdirectory(docs)

