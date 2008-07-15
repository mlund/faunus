cmake_minimum_required(VERSION 2.6)
enable_language(CXX)
project (FAUNUS CXX)

# ------------------------------
#   Default build is "Release"
# ------------------------------
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING
    "Build type: None Debug Release RelWithDebInfo MinSizeRel."
    FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

# -----------------
#   Build options
# -----------------
option (OPENMP "Enable OpenMP support")
option (GROMACS "Enable Gromacs file support" on)

# -----------------
#   Test platform
# -----------------
IF(UNIX)
  IF(APPLE)
  ELSE(APPLE)
  ENDIF(APPLE)
ENDIF(UNIX)

# ----------------------------
#   Locate external packages
# ----------------------------
set (CMAKE_MODULE_PATH classes)
find_package(DOXYGEN)

if (GROMACS)
  find_package(GMX)
  if(GMX_FOUND)
    message("-- Gromacs library: ${GMX_LIBRARIES}")
    add_definitions(-DGROMACS)
  endif(GMX_FOUND)
endif(GROMACS)

# ------------------
#   Compiler flags
# ------------------
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

#message("hej: ${CMAKE_INCLUDE_PATH}")

# ---------------------------
#   Add source directories
# ---------------------------
add_subdirectory (classes build)
add_subdirectory (doc build EXCLUDE_FROM_ALL)

