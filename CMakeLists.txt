cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project("Faunus")

# ----------------------------
#   Locate external packages
# ----------------------------
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)
find_package(OpenBabel2)
find_package(Boost COMPONENTS python)
find_package(PythonInterp)
find_package(GCCXML)
find_package(GMX)
if (ENABLE_MACPORTS)
  find_package(PythonLibsNOFRAME)
else (ENABLE_MACPORTS)
  find_package(PythonLibs)
endif(ENABLE_MACPORTS)

#find_package(PythonLibsNOFRAME) # Fix linking of non-standard framework locations

# ------------------------------
#   Default build is "Release"
# ------------------------------
include(OptimizationFlags)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release
      CACHE STRING "Build type: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

# -----------------
#   User options
# -----------------
option(ENABLE_BABEL    "Try to use OpenBabel for file I/O" on)
option(ENABLE_OPENMP   "Try to use OpenMP parallization" off)
option(ENABLE_PYTHON   "Try to create a python interface" on)
option(ENABLE_FORTRAN  "Try to compile legacy fortran code" off)
option(ENABLE_GROMACS  "Try to enable Gromacs XTC trajectory output" off)
option(ENABLE_MACPORTS "Search Macports and Fink directories for Python frameworks" off)
option(ENABLE_LIBXTC   "Internal xtc support (experimental!)" off)
mark_as_advanced( CLEAR CMAKE_VERBOSE_MAKEFILE CMAKE_CXX_COMPILER CMAKE_CXX_FLAGS_RELEASE )
mark_as_advanced( EXECUTABLE_OUTPUT_PATH LIBRARY_OUTPUT_PATH
                  CMAKE_OSX_ARCHITECTURES CMAKE_OSX_SYSROOT GCCXML)

# ---------------------------
#   Add source directories
# ---------------------------
add_subdirectory(src)
add_subdirectory(legacy/fortran)
add_subdirectory(doc EXCLUDE_FROM_ALL)
if(ENABLE_PYTHON)
  add_subdirectory(misc/python EXCLUDE_FROM_ALL)
  #add_subdirectory(undone/python-test EXCLUDE_FROM_ALL)
endif(ENABLE_PYTHON)
