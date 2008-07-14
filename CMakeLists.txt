cmake_minimum_required(VERSION 2.6)
enable_language(CXX)
project (FAUNUS CXX)

# --- Default build is "Release"
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif (NOT CMAKE_BUILD_TYPE)

# --- Build options
option (OPENMP "Enable OpenMP support")
option (GROMACS "Enable Gromacs support")

# --- Test platform
IF(UNIX)
  IF(APPLE)
  ELSE(APPLE)
  ENDIF(APPLE)
ENDIF(UNIX)

# --- Locate external packages
FIND_PACKAGE(DOXYGEN)
IF(DOXYGEN_FOUND)
ENDIF(DOXYGEN_FOUND)
#..add gromacs

# --- Compiler flags
if (CMAKE_COMPILER_IS_GNUCXX)
  set (CMAKE_CXX_FLAGS_RELEASE "-O3 -w -funroll-loops -DNDEBUG")
  if (OPENMP)
    #CMAKE_CXX_FLAGS_RELEASE
  endif(OPENMP)
endif (CMAKE_COMPILER_IS_GNUCXX)
if (CMAKE_CXX_COMPILER MATCHES "icpc")
  set (CMAKE_CXX_FLAGS_RELEASE "-O3 -w -DNDEBUG")
endif(CMAKE_CXX_COMPILER MATCHES "icpc")
if (CMAKE_CXX_COMPILER MATCHES "sunCC")
  set (CMAKE_CXX_FLAGS_RELEASE "-fast -w -DNDEBUG")
endif(CMAKE_CXX_COMPILER MATCHES "sunCC")
if (CMAKE_CXX_COMPILER MATCHES "pathCC")
  set (CMAKE_CXX_FLAGS_RELEASE "-Ofast -DNDEBUG")
endif(CMAKE_CXX_COMPILER MATCHES "pathCC")
if (CMAKE_CXX_COMPILER MATCHES "pgCC")
  set (CMAKE_CXX_FLAGS_RELEASE "-O3 -w -DNDEBUG")
endif(CMAKE_CXX_COMPILER MATCHES "pgCC")
MESSAGE("-- Setting compiler flags: "${CMAKE_CXX_FLAGS_RELEASE})

# --- Add source directories
add_subdirectory (classes build)

