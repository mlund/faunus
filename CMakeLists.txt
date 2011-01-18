cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project("Faunus")

# ----------------------------
#   Locate external packages
# ----------------------------
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)
find_package(OpenBabel2)
find_package(PythonInterp)
find_package(PythonLibs)
find_package(Subversion-patched)
include(CTest)
include(CheckIncludeFileCXX)
include(CheckCXXSourceCompiles)
enable_testing()

set(UNSIGNED_LONG_LONG_INT)
CHECK_CXX_SOURCE_COMPILES("int main() { unsigned long long int i; }" UNSIGNED_LONG_LONG_INT)

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
option(ENABLE_FORTRAN  "Try to compile legacy fortran code" off)
option(ENABLE_TWISTER  "Enable Mersenne Twister random number generator" off)
option(ENABLE_SWIG     "Try to create SWIG modules for python, tcl, ruby etc." on)
option(ENABLE_FAST_EXP "Approximate e^x function for certain pair potentials." off)
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

if(ENABLE_SWIG)
  add_subdirectory(misc/swig)
endif()

