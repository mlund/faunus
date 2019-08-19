cmake_minimum_required(VERSION 3.11.4 FATAL_ERROR)

project(faunus)

include(CheckIPOSupported)
include(ExternalProject)
include(FetchContent)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_EXTENSIONS NO)
find_package(PythonInterp 3 REQUIRED)
find_package(Threads REQUIRED)
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
option(ENABLE_POWERSASA "Fetch 3rd-party SASA calculation software" on)
option(ENABLE_PYTHON "Try to compile python bindings (experimental!)" on)
option(ENABLE_LTO "Enable Link-Time Optimizations if available" off)

option(ENABLE_APPROXMATH "Use approximate math" off)
if (ENABLE_APPROXMATH)
    add_definitions(-DFAU_APPROXMATH)
endif ()

option(ENABLE_OPENMP "Try to use OpenMP parallisation" on)
if (ENABLE_OPENMP)
  find_package(OpenMP)
  if (OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  endif()
endif()

option(ENABLE_MPI "Enable MPI code" off)
if (ENABLE_MPI)
    find_package(MPI REQUIRED)
endif()

option(ENABLE_SID "Enable SID emulation" off)

# EXPRTK
FetchContent_Declare(
    exprtk
    URL https://github.com/ArashPartow/exprtk/archive/e0e880c3797ea363d24782ba63fe362f7d94f89c.zip
    URL_HASH MD5=772293e80f8353961fcc8a2b337e8dec)
FetchContent_GetProperties(exprtk)
if(NOT exprtk_POPULATED)
    FetchContent_Populate(exprtk)
endif()
add_definitions("-Dexprtk_disable_string_capabilities")
add_definitions("-Dexprtk_disable_rtl_io_file")

# PROGRESS BAR
FetchContent_Declare(
    progressbar
    URL https://github.com/mlund/progress-cpp/archive/5ecc63ac263f7fa65b7c5559151328199f9fc5d5.zip
    URL_HASH MD5=788e8e35930ba66e4dc647890f9c0530)
FetchContent_GetProperties(progressbar)
if(NOT progressbar_POPULATED)
    FetchContent_Populate(progressbar)
endif()

# SPDLOG
ExternalProject_Add(
    project_spdlog
    PREFIX "${CMAKE_CURRENT_BINARY_DIR}/_deps"
    LOG_DOWNLOAD ON
    URL https://github.com/gabime/spdlog/archive/bb0f3839c1e52948f2ecb66e4e0aa79740a30205.zip
    URL_HASH SHA256=dece1767c33da47dbfdf60673af40b04a1557f22b6f4db4e1de8f99c2fb0a609
    CMAKE_ARGS -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -DSPDLOG_INSTALL=off -DCMAKE_POSITION_INDEPENDENT_CODE=on
    BUILD_COMMAND ${CMAKE_MAKE_PROGRAM} spdlog
    INSTALL_COMMAND ""
)
ExternalProject_Get_Property(project_spdlog source_dir)
ExternalProject_Get_Property(project_spdlog binary_dir)
set(SpdlogIncludeDir ${source_dir}/include)
add_library(spdlog STATIC IMPORTED)
set_property(TARGET spdlog PROPERTY IMPORTED_LOCATION ${binary_dir}/libspdlog.a)
add_dependencies(spdlog project_spdlog)

# MODERN JSON
FetchContent_Declare(
    modernjson
    URL "https://github.com/nlohmann/json/releases/download/v3.5.0/include.zip"
    URL_HASH SHA256=3564da9c5b0cf2e032f97c69baedf10ddbc98030c337d0327a215ea72259ea21)
FetchContent_GetProperties(modernjson)
if(NOT modernjson_POPULATED)
    FetchContent_Populate(modernjson)
endif()

# RANGE-V3
FetchContent_Declare(
    rangev3
    URL "https://github.com/ericniebler/range-v3/archive/0.5.0.tar.gz"
    URL_HASH MD5=42f932808ef51df9f1732a3d5b0ae06c)
FetchContent_GetProperties(rangev3)
if(NOT rangev3_POPULATED)
    FetchContent_Populate(rangev3)
endif()

# DOCOPT.CPP
ExternalProject_Add(project_docopt
    PREFIX "${CMAKE_CURRENT_BINARY_DIR}/_deps"
    CMAKE_ARGS -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -DCMAKE_POSITION_INDEPENDENT_CODE=on
    BUILD_COMMAND ${CMAKE_MAKE_PROGRAM} docopt_s
    INSTALL_COMMAND "" LOG_DOWNLOAD ON
    UPDATE_DISCONNECTED ON
    URL_MD5 79cfe2f0f057d8e2231ce6969e42d760
    URL "https://github.com/docopt/docopt.cpp/archive/v0.6.2.tar.gz")

ExternalProject_Get_Property(project_docopt binary_dir)
ExternalProject_Get_Property(project_docopt source_dir)
set(DocoptIncludeDir ${source_dir})
add_library(docopt STATIC IMPORTED GLOBAL)
add_dependencies(docopt project_docopt)
set_property(TARGET docopt PROPERTY IMPORTED_LOCATION ${binary_dir}/libdocopt.a)

# CPPSID
ExternalProject_Add(project_cppsid
    PREFIX "${CMAKE_CURRENT_BINARY_DIR}/_deps"
    CMAKE_ARGS -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -DCMAKE_POSITION_INDEPENDENT_CODE=on
    BUILD_COMMAND ${CMAKE_MAKE_PROGRAM} ""#cppsid player
    INSTALL_COMMAND "" LOG_DOWNLOAD ON
    UPDATE_DISCONNECTED ON
    URL_MD5 b420c4c114e00a147c2c9a974249f0d4
    URL "https://github.com/mlund/cppsid/archive/v0.2.1.tar.gz")
ExternalProject_Get_Property(project_cppsid binary_dir)
ExternalProject_Get_Property(project_cppsid source_dir)
set(CppsidIncludeDir ${source_dir}/include)
add_library(cppsid STATIC IMPORTED GLOBAL)
add_dependencies(cppsid project_cppsid)
set_property(TARGET cppsid PROPERTY IMPORTED_LOCATION ${binary_dir}/libcppsid.a)

# PYBIND11
FetchContent_Declare(
    pybind11
    URL https://github.com/pybind/pybind11/archive/v2.2.4.tar.gz
    URL_HASH MD5=c533a107bc95bd0c5525ff650b11bfa1)

# EIGEN
FetchContent_Declare(
    eigen
    URL "http://bitbucket.org/eigen/eigen/get/3.3.7.tar.gz"
    URL_HASH MD5=f2a417d083fe8ca4b8ed2bc613d20f07)
FetchContent_GetProperties(eigen)
if(NOT eigen_POPULATED)
    FetchContent_Populate(eigen)
endif()
 
# XRDFILE
ExternalProject_Add(
    project_xdrfile
    PREFIX "${CMAKE_CURRENT_BINARY_DIR}/_deps"
    URL "https://github.com/wesbarnett/libxdrfile/archive/2.1.2.tar.gz"
    URL_MD5 ee114404b4a01613b2f0167a2ad92536
    PATCH_COMMAND echo "add_library(xdrfile-static STATIC \${SRCS})" >> CMakeLists.txt
    BUILD_COMMAND ${CMAKE_MAKE_PROGRAM} xdrfile-static
    UPDATE_DISCONNECTED ON
    CMAKE_ARGS -Wno-dev -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_POSITION_INDEPENDENT_CODE=on
    LOG_DOWNLOAD ON INSTALL_COMMAND "")

ExternalProject_Get_Property(project_xdrfile source_dir)
ExternalProject_Get_Property(project_xdrfile binary_dir)
set(XdrfileIncludeDir ${source_dir}/include)
add_library(xdrfile STATIC IMPORTED)
set_property(TARGET xdrfile PROPERTY IMPORTED_LOCATION ${binary_dir}/libxdrfile-static.a)
add_dependencies(xdrfile project_xdrfile)
set_target_properties(xdrfile PROPERTIES POSITION_INDEPENDENT_CODE TRUE)

# DOCTEST
FetchContent_Declare(
    doctest
    URL "https://github.com/onqtam/doctest/archive/2.2.0.tar.gz"
    URL_MD5 3f91262877b792e545f971e54afb7602)
FetchContent_GetProperties(doctest)
if(NOT doctest_POPULATED)
    FetchContent_Populate(doctest)
    add_definitions(-DDOCTEST_CONFIG_DISABLE)
endif()

# POWERSASA
if (ENABLE_POWERSASA)
    foreach (powf array.hpp power_diagram.h power_sasa.h LICENSE)
        if (NOT EXISTS "${CMAKE_BINARY_DIR}/_deps/sasa/${powf}")
            file(DOWNLOAD
                "https://raw.githubusercontent.com/yesint/pteros/d0b65b7/thirdparty/powersasa/${powf}"
                "${CMAKE_BINARY_DIR}/_deps/sasa/${powf}" STATUS status)
        endif ()
    endforeach ()
    add_definitions("-DENABLE_POWERSASA")
    include_directories(SYSTEM "${CMAKE_BINARY_DIR}/_deps/sasa")
endif ()

# Download CoulombGalore header
foreach (file coulombgalore.h LICENSE)
    if (NOT EXISTS "${CMAKE_BINARY_DIR}/_deps/coulombgalore/${file}")
        file(DOWNLOAD
            "https://raw.githubusercontent.com/mlund/coulombgalore/52f2f99/${file}"
            "${CMAKE_BINARY_DIR}/_deps/coulombgalore/${file}" STATUS status)
    endif ()
endforeach ()
include_directories(SYSTEM "${CMAKE_BINARY_DIR}/_deps/coulombgalore")

# GIT
if (VERSION_STRING)
    set(GIT_LATEST_TAG ${VERSION_STRING})
else()
    find_package(Git)
    if (GIT_FOUND)
        execute_process(
            COMMAND ${GIT_EXECUTABLE} log -1 --format="%h\ \(%cd\)" --date short
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            OUTPUT_VARIABLE GIT_COMMIT_HASH
            OUTPUT_STRIP_TRAILING_WHITESPACE
            )
        if (GIT_COMMIT_HASH)
            add_definitions("-DGIT_COMMIT_HASH=${GIT_COMMIT_HASH}")
        endif ()
        execute_process(
            COMMAND ${GIT_EXECUTABLE} describe --abbrev=0 --tags
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            OUTPUT_VARIABLE GIT_LATEST_TAG
            OUTPUT_STRIP_TRAILING_WHITESPACE
            )
    endif()
endif()
if (GIT_LATEST_TAG)
    add_definitions("-DGIT_LATEST_TAG=${GIT_LATEST_TAG}")
endif ()

# Add third-party headers to include path. Note this is done with SYSTEM
# to disable potential compiler warnings
include_directories(SYSTEM ${eigen_SOURCE_DIR} ${doctest_SOURCE_DIR} ${modernjson_SOURCE_DIR} ${rangev3_SOURCE_DIR}/include
	${Pybind11IncludeDir} ${DocoptIncludeDir} ${CppsidIncludeDir} ${XdrfileIncludeDir} ${progressbar_SOURCE_DIR} ${SpdlogIncludeDir} ${exprtk_SOURCE_DIR})

#  Compiler specific flags
## GCC
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    add_compile_options(-Wall -Wextra -Wpedantic -Wunreachable-code -Wstrict-aliasing
        -Wno-sign-compare -Wno-unused-local-typedefs -Wno-unknown-pragmas)
## Clang
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    add_compile_options(-Wall -Wextra -Wpedantic -Wunreachable-code
        -Wstrict-aliasing -Wno-sign-compare -Wno-unused-local-typedef -Wno-unknown-pragmas)
endif()

# in Debug mode, all warnings are treated as errors and we want no optimisations
add_compile_options($<$<CONFIG:Debug>:-Werror>)
add_compile_options($<$<CONFIG:Debug>:-O0>)

add_subdirectory(src)
add_subdirectory(docs)

