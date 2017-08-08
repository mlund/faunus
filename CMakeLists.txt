cmake_minimum_required(VERSION 3.6.2 FATAL_ERROR)

project(neofaunus)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

enable_testing()

option(ENABLE_PYTHON "Try to compile python bindings (experimental!)" on)

# Download external header files (Eigen, JSON for Modern C++, Doctest, Range-V3)
include(ExternalProject)
ExternalProject_Add( modernjson
    PREFIX "${CMAKE_SOURCE_DIR}/3rdparty"
    URL "https://github.com/nlohmann/json/releases/download/v2.1.1/json.hpp"
    DOWNLOAD_NO_EXTRACT 0 DOWNLOAD_NO_PROGRESS ON
    CONFIGURE_COMMAND "" BUILD_COMMAND "" INSTALL_COMMAND "" )
ExternalProject_Get_Property(modernjson source_dir)
set(JSONIncludeDir ${source_dir}/../../src/)

ExternalProject_Add( range-v3
    PREFIX "${CMAKE_SOURCE_DIR}/3rdparty"
    URL "https://github.com/ericniebler/range-v3/archive/0.3.0.tar.gz"
    LOG_DOWNLOAD ON CONFIGURE_COMMAND "" BUILD_COMMAND "" INSTALL_COMMAND "" )
ExternalProject_Get_Property(range-v3 source_dir)
set(RangeV3IncludeDir ${source_dir}/include)

ExternalProject_Add( pybind11
    PREFIX "${CMAKE_SOURCE_DIR}/3rdparty"
    URL "https://github.com/pybind/pybind11/archive/v2.1.1.tar.gz"
    LOG_DOWNLOAD ON CONFIGURE_COMMAND "" BUILD_COMMAND "" INSTALL_COMMAND "" )
ExternalProject_Get_Property(pybind11 source_dir)
set(Pybind11IncludeDir ${source_dir}/include)

ExternalProject_Add( eigen
    PREFIX "${CMAKE_SOURCE_DIR}/3rdparty"
    URL "http://bitbucket.org/eigen/eigen/get/3.3.4.tar.gz"
    URL_MD5 1a47e78efe365a97de0c022d127607c3
    LOG_DOWNLOAD ON CONFIGURE_COMMAND "" BUILD_COMMAND "" INSTALL_COMMAND "" )
ExternalProject_Get_Property(eigen source_dir)
set(EigenIncludeDir ${source_dir})

ExternalProject_Add( libxdrfile
    PREFIX "${CMAKE_SOURCE_DIR}/3rdparty"
    URL "https://github.com/wesbarnett/libxdrfile/archive/2.1.2.tar.gz"
    CMAKE_ARGS "-Wno-dev -DCMAKE_INSTALL_PREFIX=${CMAKE_SOURCE_DIR}/3rdparty"
    LOG_DOWNLOAD ON CONFIGURE_COMMAND "" BUILD_COMMAND "" INSTALL_COMMAND "" )
ExternalProject_Get_Property(libxdrfile source_dir)
set(XDRIncludeDir "${CMAKE_SOURCE_DIR}/3rdparty/include")

ExternalProject_Add( doctest
    PREFIX "${CMAKE_SOURCE_DIR}/3rdparty"
    URL "https://raw.githubusercontent.com/onqtam/doctest/1.2.1/doctest/doctest.h"
    DOWNLOAD_NO_EXTRACT 0 DOWNLOAD_NO_PROGRESS ON
    CONFIGURE_COMMAND "" BUILD_COMMAND "" INSTALL_COMMAND "" )
ExternalProject_Get_Property(doctest source_dir)
set(DocTestIncludeDir ${source_dir}/../../src/)

ExternalProject_Add( catch
    PREFIX "${CMAKE_SOURCE_DIR}/3rdparty"
    URL "https://github.com/philsquared/Catch/releases/download/v1.9.6/catch.hpp"
    DOWNLOAD_NO_EXTRACT 0 DOWNLOAD_NO_PROGRESS ON
    CONFIGURE_COMMAND "" BUILD_COMMAND "" INSTALL_COMMAND "" )
ExternalProject_Get_Property(catch source_dir)
set(CatchIncludeDir ${source_dir}/../../src/)

set_target_properties(catch PROPERTIES EXCLUDE_FROM_ALL TRUE)

add_definitions(-DDOCTEST_CONFIG_DISABLE)
include_directories(${EigenIncludeDir} ${DocTestIncludeDir} ${JSONIncludeDir} ${RangeV3IncludeDir} ${XDRIncludeDir} ${Pybind11IncludeDir})

set(hdrs
        ${CMAKE_SOURCE_DIR}/src/average.h
        ${CMAKE_SOURCE_DIR}/src/auxiliary.h
        ${CMAKE_SOURCE_DIR}/src/core.h
        ${CMAKE_SOURCE_DIR}/src/energy.h
        ${CMAKE_SOURCE_DIR}/src/move.h
        ${CMAKE_SOURCE_DIR}/src/potentials.h
        )

add_executable(tests src/tests.cpp ${hdrs})
add_dependencies(tests modernjson doctest eigen range-v3 libxdrfile pybind11)
add_test(NAME tests COMMAND tests)

add_executable(example src/example.cpp ${hdrs})
add_dependencies(example modernjson doctest eigen range-v3)

if (ENABLE_PYTHON)
    INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/include)
    INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

    find_package(PythonLibs)
    find_package(PythonInterp)

    if (PYTHONLIBS_FOUND)
        add_library(pyfaunus SHARED src/pyfaunus.cpp)
        set_target_properties(pyfaunus PROPERTIES PREFIX "")
        set_target_properties(pyfaunus PROPERTIES SUFFIX ".so")
        target_link_libraries(pyfaunus ${PYTHON_LIBRARY})
        INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
    endif ()

endif ()
