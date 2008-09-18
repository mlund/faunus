cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project("Faunus")

# ------------------------------
#   Default build is "Release"
# ------------------------------
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release
      CACHE STRING "Build type: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

# -----------------
#   User options
# -----------------
option(FAUNUS_OPENMP  "Enable OpenMP parallel execution on multicore CPU's" off)
option(FAUNUS_GROMACS "Enable Gromacs xtc trajectory output" off)
option(FAUNUS_PYTHON  "Create python interface" off)
mark_as_advanced(CLEAR CMAKE_CXX_COMPILER CMAKE_CXX_FLAGS_RELEASE
                 EXECUTABLE_OUTPUT_PATH LIBRARY_OUTPUT_PATH
                 CMAKE_OSX_ARCHITECTURES CMAKE_OSX_SYSROOT)

# ----------------------------
#   Locate external packages
# ----------------------------
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)
include(TestCXXAcceptsFlag)
include(OptimizationFlags)
find_package(OpenBabel2)
if(FAUNUS_PYTHON)
  find_package(Boost COMPONENTS python)
  find_package(PythonLibs)
  find_package(PythonInterp)
  find_package(GCCXML)
endif(FAUNUS_PYTHON)
if(FAUNUS_GROMACS)
  find_package(GMX)
  if(GMX_MOTIF)
    find_package(X11)
    find_package(Motif)
  endif(GMX_MOTIF)
endif(FAUNUS_GROMACS)

# ---------------------------
#   Add source directories
# ---------------------------
add_subdirectory(src)
add_subdirectory(src/examples)
add_subdirectory(legacy/fortran)
add_subdirectory(doc EXCLUDE_FROM_ALL)
add_subdirectory(undone/python-test EXCLUDE_FROM_ALL)
if(FAUNUS_PYTHON AND PYTHONLIBS_FOUND AND PYTHONINTERP_FOUND AND Boost_FOUND)
  #add_subdirectory(undone/python EXCLUDE_FROM_ALL)
  add_subdirectory(undone/python-proof-of-concept EXCLUDE_FROM_ALL)
endif(FAUNUS_PYTHON AND PYTHONLIBS_FOUND AND PYTHONINTERP_FOUND AND Boost_FOUND)
