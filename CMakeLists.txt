cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project("Faunus")

# ----------------------------
#   Locate external packages
# ----------------------------
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)
find_package(OpenMP)
find_package(OpenBabel2)
find_package(Boost COMPONENTS python)
find_package(PythonInterp)
find_package(PythonLibs)
find_package(Wget)
find_package(Subversion-patched)
find_package(CVS)
include(CTest)
include(CheckIncludeFileCXX)
enable_testing()

# ------------------------------
#   Default build is "Release"
# ------------------------------
include(OptimizationFlags)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release
    CACHE STRING "Build type: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()

# -----------------
#   User options
# -----------------
option(ENABLE_BABEL    "Try to use OpenBabel for file I/O" on)
option(ENABLE_OPENMP   "Try to use OpenMP parallization" off)
option(ENABLE_PYTHON   "Try to create a python interface" on)
option(ENABLE_FORTRAN  "Try to compile legacy fortran code" off)
option(ENABLE_TWISTER  "Enable Mersenne Twister random number generator" off)
option(ENABLE_SWIG     "Try to create python interface using SWIG" on)
mark_as_advanced( CLEAR CMAKE_VERBOSE_MAKEFILE CMAKE_CXX_COMPILER CMAKE_CXX_FLAGS )
mark_as_advanced( EXECUTABLE_OUTPUT_PATH LIBRARY_OUTPUT_PATH
                  CMAKE_OSX_ARCHITECTURES CMAKE_OSX_SYSROOT GCCXML DART_TESTING_TIMEOUT)

# ---------------------------
#   Add source directories
# ---------------------------
add_subdirectory(unittests)
add_subdirectory(legacy/fortran)
add_subdirectory(src)
add_subdirectory(doc EXCLUDE_FROM_ALL)
if(ENABLE_PYTHON)
  add_subdirectory(misc/python EXCLUDE_FROM_ALL)
endif()
if(ENABLE_SWIG)
  add_subdirectory(misc/swig EXCLUDE_FROM_ALL)
endif()

