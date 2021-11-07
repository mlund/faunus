# ========== place libraries and executables in outermost dir ==========

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# ========== add to interface library (common to all faunus targets) ==========

target_link_libraries(project_options INTERFACE ZLIB::ZLIB)
target_link_libraries(project_options INTERFACE ${CMAKE_THREAD_LIBS_INIT})

# ========== pre-compiler system headers ==========

option(ENABLE_PCH "Enable Precompiled Headers" OFF)
if (ENABLE_PCH)
    target_precompile_headers(
            project_options INTERFACE <vector> <string> <map> <set> <algorithm> <utility> <random> <iostream>
            <limits> <memory>
            <Eigen/Core> <Eigen/Geometry>
            <range/v3/view/filter.hpp> <range/v3/view/transform.hpp> <range/v3/view/join.hpp>
            <cereal/cereal.hpp> <cereal/types/base_class.hpp> <cereal/types/memory.hpp>
            <nlohmann/json.hpp>
            <spdlog/fmt/fmt.h> <spdlog/spdlog.h>)
endif()

# ========== message parsing protocol ==========

option(ENABLE_MPI "Enable MPI code" off)
if (ENABLE_MPI)
    find_package(MPI REQUIRED)
    if (MPI_CXX_FOUND)
        target_link_libraries(project_options INTERFACE MPI::MPI_CXX)
        target_compile_definitions(project_options INTERFACE ENABLE_MPI)
        target_link_libraries(project_options INTERFACE mpl)
    endif()
endif()

# ========== support for free SASA ==========

if(ENABLE_FREESASA)
    target_link_libraries(project_options INTERFACE freesasa)
endif()

# ========== option to disable all test code ==========

option(ENABLE_TESTS "Build with unittests" on)
if (NOT ENABLE_TESTS)
    add_definitions(-DDOCTEST_CONFIG_DISABLE)
endif ()

# ========== target: functionparser - exprtk very slow, so split out to separate target ==========

add_library(functionparser STATIC ${CMAKE_SOURCE_DIR}/src/functionparser.cpp ${CMAKE_SOURCE_DIR}/src/functionparser.h)
set_target_properties(functionparser PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_compile_definitions(functionparser PRIVATE DOCTEST_CONFIG_DISABLE)
target_link_libraries(functionparser PRIVATE nlohmann_json exprtk)

# ========== faunus cpp and header files ==========

set(objs analysis.cpp average.cpp atomdata.cpp auxiliary.cpp bonds.cpp celllistimpl.cpp
	chainmove.cpp clustermove.cpp core.cpp forcemove.cpp units.cpp energy.cpp externalpotential.cpp
	geometry.cpp group.cpp io.cpp molecule.cpp montecarlo.cpp move.cpp mpicontroller.cpp
	particle.cpp penalty.cpp potentials.cpp random.cpp reactioncoordinate.cpp regions.cpp rotate.cpp sasa.cpp
        scatter.cpp space.cpp speciation.cpp tensor.cpp)

set(hdrs analysis.h average.h atomdata.h auxiliary.h bonds.h celllist.h celllistimpl.h
	chainmove.h clustermove.h core.h forcemove.h energy.h externalpotential.h geometry.h group.h io.h
	molecule.h montecarlo.h move.h mpicontroller.h particle.h penalty.h potentials.h
	reactioncoordinate.h rotate.h sasa.h space.h speciation.h random.h regions.h tensor.h units.h
	aux/eigen_cerealisation.h aux/eigensupport.h aux/iteratorsupport.h aux/multimatrix.h
	aux/eigen_cerealisation.h aux/eigensupport.h aux/equidistant_table.h aux/error_function.h
	aux/exp_function.h aux/invsqrt_function.h aux/iteratorsupport.h aux/legendre.h aux/multimatrix.h
	aux/pairmatrix.h aux/pow_function.h aux/table_1d.h aux/table_2d.h aux/timers.h aux/typeerasure.h aux/sparsehistogram.h)

set_source_files_properties(${objs} PROPERTIES LANGUAGE CXX)

# ========== target: libfaunus - main functionality of faunus ==========

add_library(libfaunus OBJECT ${objs} ${hdrs})
target_compile_definitions(libfaunus PRIVATE SPDLOG_COMPILED_LIB)
target_link_libraries(libfaunus PUBLIC
        project_options
        project_warnings
        cereal
        Eigen
        functionparser
        nlohmann_json
        pcg-cpp
        range-v3
        spdlog
        xdrfile
        zstr)

# ========== faunus (main executable) ==========

add_executable(faunus faunus.cpp)
target_link_libraries(faunus PRIVATE
    docopt_s
    Eigen
    libfaunus
    doctest
    progresstracker
    project_options
    range-v3
    spdlog)
target_compile_definitions(faunus PRIVATE SPDLOG_COMPILED_LIB)
set_target_properties(faunus PROPERTIES
        POSITION_INDEPENDENT_CODE ON
        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

if (ENABLE_SID)
    if (SDL2_FOUND)
        # fix buggy sdl2: https://stackoverflow.com/questions/45730098/cmake-leading-or-trailing-whitespace-policy-cmp0004
        string(STRIP ${SDL2_LIBRARIES} SDL2_LIBRARIES)
        include_directories(SYSTEM ${SDL2_INCLUDE_DIRS})
        target_link_libraries(faunus PRIVATE ${SDL2_LIBRARIES} cppsid)
        target_compile_definitions(faunus PRIVATE ENABLE_SID)
        if (PYTHONINTERP_FOUND)
            execute_process(
                    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/scripts/extractsids.py ${CMAKE_SOURCE_DIR}/scripts/music.yml
                    WORKING_DIRECTORY ${CMAKE_BINARY_DIR} ERROR_QUIET OUTPUT_QUIET)
        endif ()
    endif ()
endif ()

option(BUILD_STATIC "Build statically linked binaries" off)
if(BUILD_STATIC)
    target_link_libraries(faunus PRIVATE -static-libgcc -static-libstdc++ -static)
endif()

target_compile_definitions(faunus PRIVATE
    FAUNUS_TIPSFILE=\"${CMAKE_INSTALL_PREFIX}/share/faunus/tips.json\",\"${CMAKE_BINARY_DIR}/tips.json\")
target_compile_definitions(faunus PRIVATE FAUNUS_INSTALL_PREFIX=\"${CMAKE_INSTALL_PREFIX}\")
target_compile_definitions(faunus PRIVATE FAUNUS_BINARY_DIR=\"${CMAKE_BINARY_DIR}\")
install(TARGETS faunus DESTINATION bin)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/include)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# ========== python bindings ==========

option(ENABLE_PYTHON "Try to compile python bindings (experimental!)" off)
if (ENABLE_PYTHON)
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
            pybind11_add_module(pyfaunus pyfaunus.cpp ${objs})
            target_link_libraries(pyfaunus PRIVATE
                    project_options
                    cereal
                    Eigen
                    functionparser
                    nlohmann_json
                    range-v3
                    spdlog
                    xdrfile
                    zstr)
            target_compile_definitions(pyfaunus PUBLIC SPDLOG_COMPILED_LIB)
            target_compile_definitions(pyfaunus PRIVATE DOCTEST_CONFIG_DISABLE)
        endif()
    endif ()
endif ()

if (PYTHONINTERP_FOUND)
    # find _relative_ path of python library and thus respect CMAKE_INSTALL_PREFIX.
    # could be replaced by FindPython() but req. cmake 3.12
    # https://stackoverflow.com/questions/1242904/finding-python-site-packages-directory-with-cmake
    execute_process(
        COMMAND "${PYTHON_EXECUTABLE}" -c "if True:
        from distutils import sysconfig as sc
        print(sc.get_python_lib(prefix='', plat_specific=True))"
        OUTPUT_VARIABLE PYTHON_SITE
        OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()

if (TARGET pyfaunus)
    install(TARGETS pyfaunus DESTINATION ${PYTHON_SITE} OPTIONAL)
endif()

# ========== assorted installation targets ==========

install(FILES ${CMAKE_SOURCE_DIR}/LICENSE DESTINATION share/faunus)
install(FILES ${CMAKE_BINARY_DIR}/tips.json DESTINATION share/faunus OPTIONAL)
install(DIRECTORY ${CMAKE_BINARY_DIR}/sids DESTINATION share/faunus OPTIONAL)
install(FILES
    ${CMAKE_SOURCE_DIR}/scripts/yason.py
    ${CMAKE_SOURCE_DIR}/scripts/ob2aam.py
    ${CMAKE_SOURCE_DIR}/scripts/faunus-manual
    DESTINATION bin
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/scripts/faunustools DESTINATION ${PYTHON_SITE} OPTIONAL
    PATTERN ${CMAKE_SOURCE_DIR}/scripts/faunustool/__pycache__ EXCLUDE)
