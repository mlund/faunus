cmake_minimum_required(VERSION 3.11.4 FATAL_ERROR)

project(faunus)

include(CheckIPOSupported)
include(ExternalProject)
include(FetchContent)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_EXTENSIONS NO)
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(PythonInterp 3 REQUIRED)
find_package(Threads REQUIRED)
find_package(ZLIB REQUIRED)
enable_testing()

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Default build is "RelWithDebInfo"
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RelWithDebInfo
        CACHE STRING "Build type: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif ()

if (CMAKE_BUILD_TYPE MATCHES "RelWithDebInfo")
    string(REPLACE "NDEBUG" "_ndebug" CMAKE_CXX_FLAGS_RELWITHDEBINFO ${CMAKE_CXX_FLAGS_RELWITHDEBINFO})
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}"
        CACHE STRING "Flags used by the compiler during Release with Debug Info builds" FORCE)
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

include("cmake/externaltools.cmake")
include("cmake/gitversion.cmake")
include("cmake/compilerflags.cmake")
 
# Add third-party headers to include path. Note this is done with SYSTEM
# to disable potential compiler warnings
include_directories(SYSTEM ${eigen_SOURCE_DIR} ${doctest_SOURCE_DIR} ${modernjson_SOURCE_DIR}/include ${rangev3_SOURCE_DIR}/include
    ${nanobench_SOURCE_DIR}/src/include
    ${Pybind11IncludeDir} ${DocoptIncludeDir} ${CppsidIncludeDir} ${XdrfileIncludeDir} ${SpdlogIncludeDir}
    ${ProgressTrackerIncludeDir} ${exprtk_SOURCE_DIR} ${cereal_SOURCE_DIR}/include ${zstr_SOURCE_DIR}/src)

add_subdirectory(src)
add_subdirectory(docs)

