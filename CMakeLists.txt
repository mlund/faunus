cmake_minimum_required(VERSION 3.9 FATAL_ERROR)

project(neofaunus)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)

enable_testing()

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

option(ENABLE_POWERSASA "Fetch 3rd-party SASA calculation software" on)
option(ENABLE_PYTHON "Try to compile python bindings (experimental!)" on)

option(ENABLE_APPROXMATH "Use approximate math" off)
if (ENABLE_APPROXMATH)
    add_definitions(-DFAU_APPROXMATH)
endif ()

option(ENABLE_OPENMP "Try to use OpenMP parallization" off)
if (ENABLE_OPENMP)
  find_package(OpenMP)
  if (OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  endif()
endif()

# MPI

option(ENABLE_MPI "Enable MPI code" off)
if (ENABLE_MPI)
    find_package(MPI)
    if (MPI_CXX_FOUND)
        set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_CXX_COMPILE_FLAGS})
        if (APPLE)
            set(CMAKE_CXX_LINK_FLAGS ${MPI_CXX_LINK_FLAGS})
        else ()
            set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_CXX_LINK_FLAGS})
        endif ()
        set(LINKLIBS ${LINKLIBS} ${MPI_CXX_LIBRARIES})
        include_directories(${MPI_CXX_INCLUDE_PATH})
        add_definitions(-DENABLE_MPI)
    endif ()
endif ()

include(ExternalProject)

# PROGRESS BAR

ExternalProject_Add( progressbar
    PREFIX "${CMAKE_SOURCE_DIR}/3rdparty"
    URL "https://raw.githubusercontent.com/prakhar1989/progress-cpp/master/ProgressBar.hpp"
    URL_MD5 0cdfc175b8a19e38a8b27c9e6ab81785
    UPDATE_DISCONNECTED ON
    DOWNLOAD_NO_EXTRACT 0 DOWNLOAD_NO_PROGRESS ON
    CONFIGURE_COMMAND "" BUILD_COMMAND "" INSTALL_COMMAND "" )
ExternalProject_Get_Property(progressbar source_dir)
set(ProgressIncludeDir ${source_dir}/../../src/)

# MODERN JSON

ExternalProject_Add( modernjson
    PREFIX "${CMAKE_SOURCE_DIR}/3rdparty"
    URL "https://github.com/nlohmann/json/releases/download/v3.0.1/json.hpp"
    URL_MD5 d00f62b5574bf20251b13d2ed0fd26db
    UPDATE_DISCONNECTED ON
    DOWNLOAD_NO_EXTRACT 0 DOWNLOAD_NO_PROGRESS ON
    CONFIGURE_COMMAND "" BUILD_COMMAND "" INSTALL_COMMAND "" )
ExternalProject_Get_Property(modernjson source_dir)
set(JSONIncludeDir ${source_dir}/../../src/)

# RANGE-V3

ExternalProject_Add( range-v3
    PREFIX "${CMAKE_SOURCE_DIR}/3rdparty"
    URL "https://github.com/ericniebler/range-v3/archive/0.3.5.tar.gz"
    URL_MD5 22f34c74f37a32cb9f2e76869cba27f5
    UPDATE_DISCONNECTED ON
    LOG_DOWNLOAD ON CONFIGURE_COMMAND "" BUILD_COMMAND "" INSTALL_COMMAND "" )
ExternalProject_Get_Property(range-v3 source_dir)
set(RangeV3IncludeDir ${source_dir}/include)

# DOCOPT.CPP

ExternalProject_Add(project_docopt
    PREFIX "${CMAKE_SOURCE_DIR}/3rdparty"
    CMAKE_ARGS -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    BUILD_COMMAND make docopt_s
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

# PYBIND11

ExternalProject_Add( pybind11
    PREFIX "${CMAKE_SOURCE_DIR}/3rdparty"
    URL "https://github.com/pybind/pybind11/archive/v2.2.1.tar.gz"
    URL_MD5 bab1d46bbc465af5af3a4129b12bfa3b
    UPDATE_DISCONNECTED ON
    INSTALL_DIR "${CMAKE_SOURCE_DIR}/3rdparty"
    CMAKE_ARGS -DPYBIND11_TEST=OFF -DPYBIND11_INSTALL=ON -DCMAKE_INSTALL_PREFIX=${CMAKE_SOURCE_DIR}/3rdparty/
    LOG_DOWNLOAD ON)

# EIGEN

ExternalProject_Add( eigen
    PREFIX "${CMAKE_SOURCE_DIR}/3rdparty"
    URL "http://bitbucket.org/eigen/eigen/get/3.3.4.tar.gz"
    URL_MD5 1a47e78efe365a97de0c022d127607c3
    UPDATE_DISCONNECTED ON
    LOG_DOWNLOAD ON CONFIGURE_COMMAND "" BUILD_COMMAND "" INSTALL_COMMAND "" )
ExternalProject_Get_Property(eigen source_dir)
set(EigenIncludeDir ${source_dir})

# XRDFILE

ExternalProject_Add( project_xdrfile
    PREFIX "${CMAKE_SOURCE_DIR}/3rdparty"
    URL "https://github.com/wesbarnett/libxdrfile/archive/2.1.2.tar.gz"
    URL_MD5 ee114404b4a01613b2f0167a2ad92536
    PATCH_COMMAND echo "add_library(xdrfile-static STATIC \${SRCS})" >> CMakeLists.txt
    BUILD_COMMAND make xdrfile-static
    UPDATE_DISCONNECTED ON
    CMAKE_ARGS -Wno-dev -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    LOG_DOWNLOAD ON INSTALL_COMMAND "")

ExternalProject_Get_Property(project_xdrfile source_dir)
ExternalProject_Get_Property(project_xdrfile binary_dir)
set(XdrfileIncludeDir ${source_dir}/include)
add_library(xdrfile STATIC IMPORTED)
set_property(TARGET xdrfile PROPERTY IMPORTED_LOCATION ${binary_dir}/libxdrfile-static.a)
add_dependencies(xdrfile project_xdrfile)

# DOCTEST

ExternalProject_Add( doctest
    PREFIX "${CMAKE_SOURCE_DIR}/3rdparty"
    URL "https://raw.githubusercontent.com/onqtam/doctest/1.2.7/doctest/doctest.h"
    URL_MD5 c521af6a674b3f507fd70b75510c0690
    UPDATE_DISCONNECTED ON
    DOWNLOAD_NO_EXTRACT 0 DOWNLOAD_NO_PROGRESS ON
    CONFIGURE_COMMAND "" BUILD_COMMAND "" INSTALL_COMMAND "" )

ExternalProject_Get_Property(doctest source_dir)
set(DocTestIncludeDir ${source_dir}/../../src/)

add_definitions(-DDOCTEST_CONFIG_DISABLE)
include_directories(${EigenIncludeDir} ${DocTestIncludeDir} ${JSONIncludeDir} ${RangeV3IncludeDir}
    ${Pybind11IncludeDir} ${DocoptIncludeDir} ${XdrfileIncludeDir} ${ProgressIncludeDir})

# POWERSASA

if (ENABLE_POWERSASA)
    foreach (powf array.hpp power_diagram.h power_sasa.h LICENSE)
        if (NOT EXISTS "${CMAKE_SOURCE_DIR}/include/faunus/sasa/${powf}")
            file(DOWNLOAD
                "https://raw.githubusercontent.com/yesint/pteros/49ab75b/thirdparty/sasa/${powf}"
                "${CMAKE_SOURCE_DIR}/include/faunus/sasa/${powf}" STATUS status)
        endif ()
    endforeach ()
    add_definitions("-DFAU_POWERSASA")
    include_directories("${CMAKE_SOURCE_DIR}/include/faunus/sasa")
endif ()

# GIT

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
endif ()

set(hdrs
    ${CMAKE_SOURCE_DIR}/src/analysis.h
    ${CMAKE_SOURCE_DIR}/src/average.h
    ${CMAKE_SOURCE_DIR}/src/auxiliary.h
    ${CMAKE_SOURCE_DIR}/src/core.h
    ${CMAKE_SOURCE_DIR}/src/energy.h
    ${CMAKE_SOURCE_DIR}/src/move.h
    ${CMAKE_SOURCE_DIR}/src/potentials.h
    ${CMAKE_SOURCE_DIR}/src/space.h
    )

add_executable(tests EXCLUDE_FROM_ALL src/tests.cpp ${hdrs})
add_dependencies(tests modernjson doctest eigen range-v3 xdrfile)
target_link_libraries(tests ${LINKLIBS})
add_test(NAME tests COMMAND tests)

add_executable(faunus src/faunus.cpp ${hdrs})
add_dependencies(faunus modernjson doctest eigen range-v3 docopt xdrfile progressbar pybind11)
target_link_libraries(faunus xdrfile docopt ${LINKLIBS})

add_executable(faunus_nopbc src/faunus.cpp ${hdrs})
add_dependencies(faunus_nopbc modernjson doctest eigen range-v3 docopt xdrfile progressbar pybind11)
target_link_libraries(faunus_nopbc xdrfile docopt ${LINKLIBS})
set_target_properties(faunus_nopbc PROPERTIES COMPILE_DEFINITIONS "CUBOID_NOPBC")

add_executable(faunus_sphere src/faunus.cpp ${hdrs})
add_dependencies(faunus_sphere modernjson doctest eigen range-v3 docopt xdrfile progressbar pybind11)
target_link_libraries(faunus_sphere xdrfile docopt ${LINKLIBS})
set_target_properties(faunus_sphere PROPERTIES COMPILE_DEFINITIONS "SPHERE")

if (ENABLE_PYTHON)
    INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/include)
    INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

    find_package(PythonLibs)
    find_package(PythonInterp)

    if (PYTHONINTERP_FOUND)
        execute_process(
            COMMAND ${PYTHON_EXECUTABLE} scripts/yason.py -h
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            ERROR_QUIET OUTPUT_QUIET
            RESULT_VARIABLE yasontest)
        if (${yasontest} EQUAL 1)
            message("missing required python modules 'yaml' or 'ruamel_yaml'")
        endif()
    endif()

    if (PYTHONLIBS_FOUND)
        find_package(pybind11 PATHS ${CMAKE_SOURCE_DIR}/3rdparty/share/cmake/pybind11/)
        if (pybind11_FOUND)
            message("Adding pyfaunus target.")
            pybind11_add_module(pyfaunus src/pyfaunus.cpp)
        endif()
    endif ()

endif ()

add_custom_target(
    manual
    COMMENT "Building PDF manual"
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/docs/_docs
    COMMAND cat
    header.md
    docs.md
    install.md
    running.md
    topology.md
    energy.md
    moves.md
    analysis.md |
    grep -v ".notice--" |
    grep -v "Github-Improve_this_page" |
    pandoc --self-contained --highlight-style=tango --template ../pandoc.tex
    -V date:"\\today"
    -V thanks:"Based on git revision "${GIT_COMMIT_HASH}
    -N --toc -o ${CMAKE_BINARY_DIR}/manual.pdf
    )
set_target_properties(manual PROPERTIES EXCLUDE_FROM_ALL TRUE)

# Install targets

install(TARGETS faunus DESTINATION bin)
install(TARGETS faunus_nopbc DESTINATION bin)
install(TARGETS faunus_sphere DESTINATION bin)

install(FILES ${CMAKE_BINARY_DIR}/manual.pdf DESTINATION share/faunus OPTIONAL)

install(FILES ${CMAKE_SOURCE_DIR}/scripts/yason.py DESTINATION bin
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

install(FILES
    ${CMAKE_SOURCE_DIR}/examples/README.md
    ${CMAKE_SOURCE_DIR}/examples/bulk.yml
    ${CMAKE_SOURCE_DIR}/examples/water.yml
    ${CMAKE_SOURCE_DIR}/examples/water.xyz
    ${CMAKE_SOURCE_DIR}/examples/membrane.yml
    ${CMAKE_SOURCE_DIR}/examples/membrane.state
    ${CMAKE_SOURCE_DIR}/examples/lipid.xyz
    ${CMAKE_SOURCE_DIR}/examples/penalty.yml
    ${CMAKE_SOURCE_DIR}/examples/isochoric.yml
    ${CMAKE_SOURCE_DIR}/examples/isochoric.xyz
    DESTINATION share/faunus/examples)

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

