cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project("Faunus")

if (POLICY CMP0042)
    cmake_policy(SET CMP0042 OLD)
endif ()

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if (CMAKE_COMPILER_IS_GNUCC AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.9)
    message(FATAL "GCC 4.9 or higher required.")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.5)
        message(FATAL_ERROR "Clang 3.5 or higher required.")
    endif()
endif()

# ----------------------------
#   Locate external packages
# ----------------------------
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/scripts)
include(CheckIncludeFileCXX)
include(CheckCXXSourceCompiles)
enable_testing()

find_package(Eigen3 3.2.9)
if (NOT EIGEN3_FOUND)
    message(FATAL_ERROR "Eigen 3.1.2+ is required.")
endif ()
include_directories(${EIGEN3_INCLUDE_DIR})

# ------------------------------
#   Default build is "Release"
# ------------------------------
include(OptimizationFlags)
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RelWithDebInfo
            CACHE STRING "Build type: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif ()

# -----------------
#   User options
# -----------------
option(ENABLE_BABEL "Try to use OpenBabel for file I/O (experimental)" off)
option(ENABLE_OPENMP "Try to use OpenMP parallization" off)
option(ENABLE_MPI "Enable MPI code" off)
option(ENABLE_STATIC "Use static instead of dynamic linkage of faunus library" off)
option(ENABLE_PYTHON "Try to compile python bindings (experimental!)" on)
option(ENABLE_APPROXMATH "Use approximate math (Quake inverse sqrt, fast exponentials etc.)" off)
option(ENABLE_HASHTABLE "Use hash tables for bond bookkeeping - may be faster for big systems" off)
option(ENABLE_UNICODE "Use unicode characters in output" on)
option(ENABLE_POWERSASA "Fetch 3rd-party SASA calculation software" off)
mark_as_advanced(CLEAR CMAKE_VERBOSE_MAKEFILE CMAKE_CXX_COMPILER CMAKE_CXX_FLAGS)
mark_as_advanced(EXECUTABLE_OUTPUT_PATH LIBRARY_OUTPUT_PATH
        CMAKE_OSX_ARCHITECTURES CMAKE_OSX_SYSROOT GCCXML DART_TESTING_TIMEOUT)

# -------------------------------
#   Get git commit hash
# -------------------------------
find_package(Git)
if (GIT_FOUND)
    execute_process(
            COMMAND ${GIT_EXECUTABLE} log -1 --format="%h\ \(%ci\)"
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            OUTPUT_VARIABLE GIT_COMMIT_HASH
            OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if (GIT_COMMIT_HASH)
        add_definitions("-DGIT_COMMIT_HASH=${GIT_COMMIT_HASH}")
    endif ()
endif ()

# ---------------------------
#   Add source directories
# ---------------------------
add_subdirectory(src)
add_subdirectory(scripts)
add_subdirectory(doc EXCLUDE_FROM_ALL)

if (CMAKE_BUILD_TYPE MATCHES "RelWithDebInfo")
    string(REPLACE "NDEBUG" "_ndebug" CMAKE_CXX_FLAGS_RELWITHDEBINFO ${CMAKE_CXX_FLAGS_RELWITHDEBINFO})
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}"
            CACHE STRING "Flags used by the compiler during Release with Debug Info builds" FORCE)

    message(
            "...............................................................\n"
            "  Faunus is compiled in Debug mode to enable checks and\n"
            "  assertions highly recommended when developing and testing.\n"
            "  For production, switch to Release for optimal performance:\n\n"
            "  -DCMAKE_BUILD_TYPE=Release\n"
            "***************************************************************\n")
endif ()
if (CMAKE_BUILD_TYPE MATCHES "Release")
    message(
            "......................................................................\n"
            "  Faunus is running in Release mode recommended for production runs.\n"
            "  For setting up simulations it is HIGHLY recommended to use Debug\n"
            "  mode to enable run-time checks and assertions:\n\n"
            "  -DCMAKE_BUILD_TYPE=RelWithDebInfo\n"
            "**********************************************************************\n")
endif ()

