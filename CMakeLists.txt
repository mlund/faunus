cmake_minimum_required(VERSION 3.11.4 FATAL_ERROR)

project(faunus)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_EXTENSIONS NO)
find_package(PythonInterp 3 REQUIRED)
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

option(BUILD_STATIC "Build statically linked binaries" off)
option(ENABLE_POWERSASA "Fetch 3rd-party SASA calculation software" on)
option(ENABLE_PYTHON "Try to compile python bindings (experimental!)" on)

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
if (ENABLE_SID)
   list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR})
   find_package(SDL2 REQUIRED)
endif()

include(ExternalProject)
include(FetchContent)

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

# CPPSID
ExternalProject_Add(project_cppsid
    PREFIX "${CMAKE_SOURCE_DIR}/3rdparty"
    CMAKE_ARGS -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    BUILD_COMMAND make cppsid player
    INSTALL_COMMAND "" LOG_DOWNLOAD ON
    UPDATE_DISCONNECTED ON
    URL_MD5 2d989c5dde42173ab4eba6374aabf042
    URL "https://github.com/mlund/cppsid//archive/8d18d8e4109518c15725326f90cef879a161e03f.zip")
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
    PREFIX "${CMAKE_SOURCE_DIR}/3rdparty"
    URL "https://github.com/wesbarnett/libxdrfile/archive/2.1.2.tar.gz"
    URL_MD5 ee114404b4a01613b2f0167a2ad92536
    PATCH_COMMAND echo "add_library(xdrfile-static STATIC \${SRCS})" >> CMakeLists.txt
    BUILD_COMMAND make xdrfile-static
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
        if (NOT EXISTS "${CMAKE_SOURCE_DIR}/include/faunus/sasa/${powf}")
            file(DOWNLOAD
                "https://raw.githubusercontent.com/yesint/pteros/d0b65b7/thirdparty/powersasa/${powf}"
                "${CMAKE_SOURCE_DIR}/include/faunus/sasa/${powf}" STATUS status)
        endif ()
    endforeach ()
    add_definitions("-DENABLE_POWERSASA")
    include_directories(SYSTEM "${CMAKE_SOURCE_DIR}/include/faunus/sasa")
endif ()

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
	${Pybind11IncludeDir} ${DocoptIncludeDir} ${CppsidIncludeDir} ${XdrfileIncludeDir} ${progressbar_SOURCE_DIR} ${exprtk_SOURCE_DIR})

# faunus cpp files
set(objs
    ${CMAKE_SOURCE_DIR}/src/functionparser.cpp
    ${CMAKE_SOURCE_DIR}/src/random.cpp
    ${CMAKE_SOURCE_DIR}/src/mpi.cpp
    ${CMAKE_SOURCE_DIR}/src/bonds.cpp
    ${CMAKE_SOURCE_DIR}/src/chainmove.cpp
    ${CMAKE_SOURCE_DIR}/src/clustermove.cpp
    ${CMAKE_SOURCE_DIR}/src/core.cpp
    ${CMAKE_SOURCE_DIR}/src/units.cpp
    ${CMAKE_SOURCE_DIR}/src/analysis.cpp
    ${CMAKE_SOURCE_DIR}/src/atomdata.cpp
    ${CMAKE_SOURCE_DIR}/src/energy.cpp
    ${CMAKE_SOURCE_DIR}/src/externalpotential.cpp
    ${CMAKE_SOURCE_DIR}/src/geometry.cpp
    ${CMAKE_SOURCE_DIR}/src/group.cpp
    ${CMAKE_SOURCE_DIR}/src/io.cpp
    ${CMAKE_SOURCE_DIR}/src/molecule.cpp
    ${CMAKE_SOURCE_DIR}/src/move.cpp
    ${CMAKE_SOURCE_DIR}/src/particle.cpp
    ${CMAKE_SOURCE_DIR}/src/penalty.cpp
    ${CMAKE_SOURCE_DIR}/src/potentials.cpp
    ${CMAKE_SOURCE_DIR}/src/reactioncoordinate.cpp
    ${CMAKE_SOURCE_DIR}/src/rotate.cpp
    ${CMAKE_SOURCE_DIR}/src/space.cpp
    ${CMAKE_SOURCE_DIR}/src/speciation.cpp
    ${CMAKE_SOURCE_DIR}/src/tensor.cpp)

set_source_files_properties(${objs} PROPERTIES LANGUAGE CXX)

# faunus header files
set(hdrs
    ${CMAKE_SOURCE_DIR}/src/analysis.h
    ${CMAKE_SOURCE_DIR}/src/average.h
    ${CMAKE_SOURCE_DIR}/src/atomdata.h
    ${CMAKE_SOURCE_DIR}/src/auxiliary.h
    ${CMAKE_SOURCE_DIR}/src/bonds.h
    ${CMAKE_SOURCE_DIR}/src/chainmove.h
    ${CMAKE_SOURCE_DIR}/src/clustermove.h
    ${CMAKE_SOURCE_DIR}/src/core.h
    ${CMAKE_SOURCE_DIR}/src/energy.h
    ${CMAKE_SOURCE_DIR}/src/externalpotential.h
    ${CMAKE_SOURCE_DIR}/src/functionparser.h
    ${CMAKE_SOURCE_DIR}/src/geometry.h
    ${CMAKE_SOURCE_DIR}/src/group.h
    ${CMAKE_SOURCE_DIR}/src/io.h
    ${CMAKE_SOURCE_DIR}/src/molecule.h
    ${CMAKE_SOURCE_DIR}/src/move.h
    ${CMAKE_SOURCE_DIR}/src/mpi.h
    ${CMAKE_SOURCE_DIR}/src/particle.h
    ${CMAKE_SOURCE_DIR}/src/penalty.h
    ${CMAKE_SOURCE_DIR}/src/potentials.h
    ${CMAKE_SOURCE_DIR}/src/reactioncoordinate.h
    ${CMAKE_SOURCE_DIR}/src/rotate.h
    ${CMAKE_SOURCE_DIR}/src/space.h
    ${CMAKE_SOURCE_DIR}/src/speciation.h
    ${CMAKE_SOURCE_DIR}/src/random.h
    ${CMAKE_SOURCE_DIR}/src/tensor.h
    ${CMAKE_SOURCE_DIR}/src/units.h
    )

# Compiler specific flags

## GCC
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    add_compile_options(-Wall -Wextra -Wpedantic -Wunreachable-code -Wstrict-aliasing
        -Wno-sign-compare -Wno-unused-local-typedefs -Wno-unknown-pragmas)
## Clang
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    add_compile_options(-Wall -Wextra -Wpedantic -Wunreachable-code
        -Wstrict-aliasing -Wno-sign-compare -Wno-unused-local-typedef -Wno-unknown-pragmas)
endif()

# target: unittests
add_executable(unittests src/unittests.cpp ${objs} ${hdrs})
target_link_libraries(unittests xdrfile)
add_test(NAME unittests COMMAND unittests)

# target: faunus

add_executable(faunus src/faunus.cpp ${objs} ${hdrs})
target_link_libraries(faunus PRIVATE xdrfile docopt)

if(ENABLE_SID)
    target_compile_definitions(faunus PRIVATE ENABLE_SID)
    target_link_libraries(faunus PRIVATE SDL2::Main cppsid)
endif()

if(BUILD_STATIC)
    target_link_libraries(faunus PRIVATE -static-libgcc -static-libstdc++ -static)
endif()

target_compile_definitions(faunus PRIVATE
    FAUNUS_TIPSFILE=\"${CMAKE_INSTALL_PREFIX}/share/faunus/tips.json\",\"${CMAKE_SOURCE_DIR}/tips.json\" )

if (MPI_CXX_FOUND)
    target_link_libraries(faunus PRIVATE MPI::MPI_CXX)
    target_compile_definitions(faunus PRIVATE ENABLE_MPI)
endif()

# python bindings
if (ENABLE_PYTHON)
    INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/include)
    INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
    find_package(PythonLibs 3 REQUIRED)

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
        FetchContent_GetProperties(pybind11)
        if(NOT pybind11_POPULATED)
            FetchContent_Populate(pybind11)
            add_subdirectory(${pybind11_SOURCE_DIR} ${pybind11_BINARY_DIR})
            pybind11_add_module(pyfaunus src/pyfaunus.cpp ${objs})
            target_link_libraries(pyfaunus PUBLIC xdrfile)
        endif()
    endif ()
endif ()

# MANUAL

find_program(PANDOC pandoc)
if (PANDOC)
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
        ${PANDOC} --self-contained --highlight-style=tango --template ../pandoc.tex
        -V date:${GIT_LATEST_TAG}""
        -N --toc -o ${CMAKE_BINARY_DIR}/manual.pdf
        )
    set_target_properties(manual PROPERTIES EXCLUDE_FROM_ALL TRUE)

    add_custom_target(
        manual_html
        COMMENT "Building html manual"
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
        ${PANDOC} --template ../GitHub.html5 -s --mathjax
        -M title:The\ Faunus\ User\ Guide
        -V date:${GIT_LATEST_TAG}""
        -N --toc -o ${CMAKE_BINARY_DIR}/manual.html
        )

    if (PYTHONINTERP_FOUND)
        add_custom_target(
            usagetips
            ALL
            COMMENT "extracting tips from manual"
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            COMMAND ${PYTHON_EXECUTABLE} scripts/extracttips.py)
        add_dependencies(usagetips manual_html)
    endif()
endif()

# TESTS

set(EXAMPLES_DIR ${CMAKE_SOURCE_DIR}/examples)
set(YASON ${CMAKE_SOURCE_DIR}/scripts/yason.py)
set(JSON_COMPARE ${CMAKE_SOURCE_DIR}/scripts/jsoncompare.py)
if (TARGET pyfaunus)
    add_test(
        NAME pythontest
        COMMAND sh -c "${PYTHON_EXECUTABLE} pythontest.py -v"
        WORKING_DIRECTORY ${EXAMPLES_DIR})
    set_tests_properties(pythontest PROPERTIES ENVIRONMENT "PYTHONPATH=$<TARGET_FILE_DIR:pyfaunus>")
endif()

add_test(
    NAME swapconf
    COMMAND sh -c "${PYTHON_EXECUTABLE} ${YASON} swapconf.yml\
    | $<TARGET_FILE:faunus> --nobar --state swapconf.state.json\
    ; ${PYTHON_EXECUTABLE} ${JSON_COMPARE} swapconf.out.json out.json --tol 0.05"
    WORKING_DIRECTORY ${EXAMPLES_DIR})

add_test(
    NAME doublelayer
    COMMAND sh -c "${PYTHON_EXECUTABLE} ${YASON} doublelayer.yml\
    | $<TARGET_FILE:faunus> --nobar --state doublelayer.state.json\
    ; ${PYTHON_EXECUTABLE} ${JSON_COMPARE} doublelayer.out.json out.json --tol 0.05"
    WORKING_DIRECTORY ${EXAMPLES_DIR})

add_test(
    NAME penalty
    COMMAND sh -c "jupyter nbconvert --execute penalty.ipynb --ExecutePreprocessor.timeout=None"
    WORKING_DIRECTORY ${EXAMPLES_DIR})
set_tests_properties(penalty PROPERTIES ENVIRONMENT FAUNUS_EXECUTABLE=$<TARGET_FILE:faunus>)

add_test(
    NAME minimal
    COMMAND sh -c "${PYTHON_EXECUTABLE} ${YASON} minimal.yml\
    | $<TARGET_FILE:faunus> --nobar --state minimal.state.json\
    ; ${PYTHON_EXECUTABLE} ${JSON_COMPARE} minimal.out.json out.json --tol 0.02"
    WORKING_DIRECTORY ${EXAMPLES_DIR})

add_test(
    NAME chain
    COMMAND sh -c "${PYTHON_EXECUTABLE} ${YASON} chain.yml\
    | $<TARGET_FILE:faunus> --nobar --state chain.state.json\
    ; ${PYTHON_EXECUTABLE} ${JSON_COMPARE} chain.out.json out.json --tol 0.05"
    WORKING_DIRECTORY ${EXAMPLES_DIR}/chain)

add_test(
    NAME polymers
    COMMAND sh -c "${PYTHON_EXECUTABLE} ${YASON} polymers.yml\
    | $<TARGET_FILE:faunus> --nobar --state polymers.state.json\
    ; ${PYTHON_EXECUTABLE} ${JSON_COMPARE} polymers.out.json out.json --tol 0.05"
    WORKING_DIRECTORY ${EXAMPLES_DIR})

add_test(
    NAME speciation
    COMMAND sh -c "${PYTHON_EXECUTABLE} ${YASON} speciation.yml\
    | $<TARGET_FILE:faunus> --quiet --state speciation.state.json\
    ; ${PYTHON_EXECUTABLE} ${JSON_COMPARE} speciation.out.json out.json --tol 0.02"
    WORKING_DIRECTORY ${EXAMPLES_DIR})

add_test(
    NAME phosphate
    COMMAND sh -c "${PYTHON_EXECUTABLE} ${YASON} phosphate.yml\
    | $<TARGET_FILE:faunus> --quiet --state phosphate.state.json\
    ; ${PYTHON_EXECUTABLE} ${JSON_COMPARE} phosphate.out.json out.json --tol 0.02"
    WORKING_DIRECTORY ${EXAMPLES_DIR})

add_test(
    NAME sodiumphosphate
    COMMAND sh -c "${PYTHON_EXECUTABLE} ${YASON} sodiumphosphate.yml\
    | $<TARGET_FILE:faunus> --quiet --state sodiumphosphate.state.json\
    ; ${PYTHON_EXECUTABLE} ${JSON_COMPARE} sodiumphosphate.out.json out.json --tol 0.05"
    WORKING_DIRECTORY ${EXAMPLES_DIR})

add_test(
    NAME cluster
    COMMAND sh -c "${PYTHON_EXECUTABLE} ${YASON} cluster.yml\
    | $<TARGET_FILE:faunus> --quiet --state cluster.state.json\
    ; ${PYTHON_EXECUTABLE} ${JSON_COMPARE} cluster.out.json out.json --tol 0.1"
    WORKING_DIRECTORY ${EXAMPLES_DIR})

add_test(
    NAME bulk
    COMMAND sh -c "${PYTHON_EXECUTABLE} ${YASON} bulk.yml\
    | $<TARGET_FILE:faunus> --quiet"
    WORKING_DIRECTORY ${EXAMPLES_DIR})

add_test(
    NAME seawater
    COMMAND sh -c "${PYTHON_EXECUTABLE} ${YASON} seawater.yml\
    | $<TARGET_FILE:faunus> --quiet"
    WORKING_DIRECTORY ${EXAMPLES_DIR})

add_test(
    NAME stockmayer
    COMMAND sh -c "${PYTHON_EXECUTABLE} ${YASON} stockmayer.yml\
    | $<TARGET_FILE:faunus> --quiet --state stockmayer.state.json\
    ; ${PYTHON_EXECUTABLE} ${JSON_COMPARE} stockmayer.out.json out.json --tol 0.065"
    WORKING_DIRECTORY ${EXAMPLES_DIR})

add_test(
    NAME membrane
    COMMAND sh -c "${PYTHON_EXECUTABLE} ${YASON} membrane.yml\
    | $<TARGET_FILE:faunus> --quiet --state membrane.state.json\
    ; ${PYTHON_EXECUTABLE} ${JSON_COMPARE} membrane.out.json out.json --tol 0.065"
    WORKING_DIRECTORY ${EXAMPLES_DIR})

# INSTALL TARGETS

install(TARGETS faunus DESTINATION bin)

if (TARGET pyfaunus)
    # find _relative_ path of python library and thus respect CMAKE_INSTALL_PREFIX.
    # could be replaced by FindPython() but req. cmake 3.12
    # https://stackoverflow.com/questions/1242904/finding-python-site-packages-directory-with-cmake
    execute_process(
        COMMAND "${PYTHON_EXECUTABLE}" -c "if True:
        from distutils import sysconfig as sc
        print(sc.get_python_lib(prefix='', plat_specific=True))"
        OUTPUT_VARIABLE PYTHON_SITE
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    install(TARGETS pyfaunus DESTINATION ${PYTHON_SITE} OPTIONAL)
endif()

install(FILES ${CMAKE_BINARY_DIR}/LICENSE DESTINATION share/faunus)
install(FILES ${CMAKE_BINARY_DIR}/tips.json DESTINATION share/faunus OPTIONAL)
install(FILES ${CMAKE_BINARY_DIR}/manual.pdf DESTINATION share/faunus OPTIONAL)
install(FILES ${CMAKE_BINARY_DIR}/manual.html DESTINATION share/faunus OPTIONAL)
install(FILES
    ${CMAKE_SOURCE_DIR}/scripts/yason.py
    ${CMAKE_SOURCE_DIR}/scripts/ob2aam.py
    ${CMAKE_SOURCE_DIR}/scripts/faunus-manual
    DESTINATION bin
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

install(FILES
    ${EXAMPLES_DIR}/README.md
    ${EXAMPLES_DIR}/bulk.yml
    ${EXAMPLES_DIR}/cluster.yml
    ${EXAMPLES_DIR}/cluster.agr
    ${EXAMPLES_DIR}/cluster.out.json
    ${EXAMPLES_DIR}/cluster.state.json
    ${EXAMPLES_DIR}/doublelayer.yml
    ${EXAMPLES_DIR}/doublelayer.agr
    ${EXAMPLES_DIR}/doublelayer.out.json
    ${EXAMPLES_DIR}/doublelayer.state.json
    ${EXAMPLES_DIR}/doublelayer.mfcorr.dat
    ${EXAMPLES_DIR}/isochoric.yml
    ${EXAMPLES_DIR}/isochoric.xyz
    ${EXAMPLES_DIR}/membrane.yml
    ${EXAMPLES_DIR}/membrane.out.json
    ${EXAMPLES_DIR}/membrane.state.json
    ${EXAMPLES_DIR}/lipid.xyz
    ${EXAMPLES_DIR}/minimal.yml
    ${EXAMPLES_DIR}/minimal.out.json
    ${EXAMPLES_DIR}/minimal.state.json
    ${EXAMPLES_DIR}/penalty.yml
    ${EXAMPLES_DIR}/penalty.ipynb
    ${EXAMPLES_DIR}/chain/chain.yml
    ${EXAMPLES_DIR}/chain/chain.xyz
    ${EXAMPLES_DIR}/chain/chain.state.json
    ${EXAMPLES_DIR}/chain/chain.out.json
    ${EXAMPLES_DIR}/polymers.yml
    ${EXAMPLES_DIR}/polymers.readme
    ${EXAMPLES_DIR}/polymers.pdf
    ${EXAMPLES_DIR}/polymers.state.json
    ${EXAMPLES_DIR}/polymers.out.json
    ${EXAMPLES_DIR}/pythontest.py
    ${EXAMPLES_DIR}/seawater.yml
    ${EXAMPLES_DIR}/speciation.yml
    ${EXAMPLES_DIR}/speciation.out.json
    ${EXAMPLES_DIR}/speciation.state.json
    ${EXAMPLES_DIR}/stockmayer.yml
    ${EXAMPLES_DIR}/stockmayer.out.json
    ${EXAMPLES_DIR}/stockmayer.state.json
    ${EXAMPLES_DIR}/swapconf.yml
    ${EXAMPLES_DIR}/swapconf.out.json
    ${EXAMPLES_DIR}/swapconf.state.json
    ${EXAMPLES_DIR}/swapconf.conformations.pqr
    ${EXAMPLES_DIR}/swapconf.weights.dat
    ${EXAMPLES_DIR}/water.yml
    ${EXAMPLES_DIR}/water.xyz
    DESTINATION share/faunus/examples)

# INFORMATION ABOUT DEBUG MODE

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
endif()
if (CMAKE_BUILD_TYPE MATCHES "Release")
    message(
        "......................................................................\n"
        "  Faunus is running in Release mode recommended for production runs.\n"
        "  For setting up simulations it is HIGHLY recommended to use Debug\n"
        "  mode to enable run-time checks and assertions:\n\n"
        "  -DCMAKE_BUILD_TYPE=RelWithDebInfo\n"
        "**********************************************************************\n")
endif()
