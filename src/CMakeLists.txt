# faunus cpp files

# place libraries and executables in outermost dir
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

set(objs
    ${CMAKE_SOURCE_DIR}/src/random.cpp
    ${CMAKE_SOURCE_DIR}/src/mpicontroller.cpp
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
    ${CMAKE_SOURCE_DIR}/src/montecarlo.cpp
    ${CMAKE_SOURCE_DIR}/src/move.cpp
    ${CMAKE_SOURCE_DIR}/src/particle.cpp
    ${CMAKE_SOURCE_DIR}/src/penalty.cpp
    ${CMAKE_SOURCE_DIR}/src/potentials.cpp
    ${CMAKE_SOURCE_DIR}/src/reactioncoordinate.cpp
    ${CMAKE_SOURCE_DIR}/src/regions.cpp
    ${CMAKE_SOURCE_DIR}/src/rotate.cpp
    ${CMAKE_SOURCE_DIR}/src/space.cpp
    ${CMAKE_SOURCE_DIR}/src/speciation.cpp
    ${CMAKE_SOURCE_DIR}/src/tensor.cpp)

set_source_files_properties(${objs} PROPERTIES LANGUAGE CXX)

# faunus header files
set(hdrs
    ${CMAKE_SOURCE_DIR}/src/aux/eigensupport.h
    ${CMAKE_SOURCE_DIR}/src/aux/iteratorsupport.h
    ${CMAKE_SOURCE_DIR}/src/aux/multimatrix.h
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
    ${CMAKE_SOURCE_DIR}/src/geometry.h
    ${CMAKE_SOURCE_DIR}/src/group.h
    ${CMAKE_SOURCE_DIR}/src/io.h
    ${CMAKE_SOURCE_DIR}/src/molecule.h
    ${CMAKE_SOURCE_DIR}/src/montecarlo.h
    ${CMAKE_SOURCE_DIR}/src/move.h
    ${CMAKE_SOURCE_DIR}/src/mpicontroller.h
    ${CMAKE_SOURCE_DIR}/src/particle.h
    ${CMAKE_SOURCE_DIR}/src/penalty.h
    ${CMAKE_SOURCE_DIR}/src/potentials.h
    ${CMAKE_SOURCE_DIR}/src/reactioncoordinate.h
    ${CMAKE_SOURCE_DIR}/src/rotate.h
    ${CMAKE_SOURCE_DIR}/src/space.h
    ${CMAKE_SOURCE_DIR}/src/speciation.h
    ${CMAKE_SOURCE_DIR}/src/random.h
    ${CMAKE_SOURCE_DIR}/src/regions.h
    ${CMAKE_SOURCE_DIR}/src/tensor.h
    ${CMAKE_SOURCE_DIR}/src/units.h
    )

# faunus header files
set(tsts
    ${CMAKE_SOURCE_DIR}/src/atomdata_test.h
    ${CMAKE_SOURCE_DIR}/src/core_test.h
    ${CMAKE_SOURCE_DIR}/src/energy_test.h
    ${CMAKE_SOURCE_DIR}/src/group_test.h
    ${CMAKE_SOURCE_DIR}/src/potentials_test.h
    ${CMAKE_SOURCE_DIR}/src/tensor_test.h
    )

# target: functionparser - exprtk very slow, so split out to separate target
add_library(functionparser STATIC ${CMAKE_SOURCE_DIR}/src/functionparser.cpp ${CMAKE_SOURCE_DIR}/src/functionparser.h)
set_target_properties(functionparser PROPERTIES
    POSITION_INDEPENDENT_CODE ON)
    #LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# targer: libfaunus - main functionality of faunus
add_library(libfaunus STATIC ${objs} ${hdrs})
target_compile_definitions(libfaunus PRIVATE SPDLOG_COMPILED_LIB)
target_link_libraries(libfaunus PRIVATE ${CMAKE_THREAD_LIBS_INIT} xdrfile spdlog functionparser)
set_target_properties(libfaunus PROPERTIES
    POSITION_INDEPENDENT_CODE ON)
    #LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
if (MPI_CXX_FOUND)
    target_link_libraries(libfaunus PRIVATE MPI::MPI_CXX)
    target_compile_definitions(libfaunus PRIVATE ENABLE_MPI)
endif()
if(ENABLE_FREESASA)
    target_link_libraries(libfaunus PRIVATE freesasa)
endif()

# target: unittests
add_executable(unittests unittests.cpp ${tsts})
target_link_libraries(unittests libfaunus)
target_compile_definitions(unittests PRIVATE SPDLOG_COMPILED_LIB)
#set_target_properties(unittests PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
add_test(NAME unittests COMMAND unittests)

# target: faunus
add_executable(faunus faunus.cpp)
target_link_libraries(faunus PRIVATE libfaunus docopt progresstracker)
target_compile_definitions(faunus PRIVATE SPDLOG_COMPILED_LIB)
set_target_properties(faunus PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
if (MPI_CXX_FOUND)
    target_link_libraries(faunus PRIVATE MPI::MPI_CXX)
    target_compile_definitions(faunus PRIVATE ENABLE_MPI)
endif()

if (ENABLE_LTO)
    check_ipo_supported(RESULT result)
    if(result)
        set_property(TARGET faunus PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    endif()
endif()

if(ENABLE_SID)
    if (SDL2_FOUND)
        # fix buggy sdl2: https://stackoverflow.com/questions/45730098/cmake-leading-or-trailing-whitespace-policy-cmp0004
        string(STRIP ${SDL2_LIBRARIES} SDL2_LIBRARIES)
        include_directories(${SDL2_INCLUDE_DIRS})
        target_link_libraries(faunus PRIVATE ${SDL2_LIBRARIES} cppsid)
        target_compile_definitions(faunus PRIVATE ENABLE_SID)
        if (PYTHONINTERP_FOUND)
            execute_process(
                COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/scripts/extractsids.py ${CMAKE_SOURCE_DIR}/scripts/music.yml
                WORKING_DIRECTORY ${CMAKE_BINARY_DIR} ERROR_QUIET OUTPUT_QUIET)
        endif()
    endif()
endif()

if(BUILD_STATIC)
    target_link_libraries(faunus PRIVATE -static-libgcc -static-libstdc++ -static)
endif()

target_compile_definitions(faunus PRIVATE
    FAUNUS_TIPSFILE=\"${CMAKE_INSTALL_PREFIX}/share/faunus/tips.json\",\"${CMAKE_BINARY_DIR}/tips.json\")
target_compile_definitions(faunus PRIVATE FAUNUS_INSTALL_PREFIX=\"${CMAKE_INSTALL_PREFIX}\")
target_compile_definitions(faunus PRIVATE FAUNUS_BINARY_DIR=\"${CMAKE_BINARY_DIR}\")

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
            pybind11_add_module(pyfaunus pyfaunus.cpp)
            #target_link_libraries(pyfaunus PUBLIC ${CMAKE_THREAD_LIBS_INIT} xdrfile spdlog functionparser libfaunus)
            target_link_libraries(pyfaunus PRIVATE libfaunus)
            #set_target_properties(pyfaunus PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
        endif()
    endif ()
endif ()

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
    WORKING_DIRECTORY ${EXAMPLES_DIR}/swapconf)

add_test(
    NAME doublelayer
    COMMAND sh -c "${PYTHON_EXECUTABLE} ${YASON} doublelayer.yml\
    | $<TARGET_FILE:faunus> --nobar --state doublelayer.state.json\
    ; ${PYTHON_EXECUTABLE} ${JSON_COMPARE} doublelayer.out.json out.json --tol 0.05"
    WORKING_DIRECTORY ${EXAMPLES_DIR}/doublelayer)

add_test(
    NAME penalty
    COMMAND sh -c "jupyter nbconvert --execute penalty.ipynb --ExecutePreprocessor.timeout=None"
    WORKING_DIRECTORY ${EXAMPLES_DIR}/penalty)
set_property(TEST penalty APPEND PROPERTY ENVIRONMENT FAUNUS_EXECUTABLE=$<TARGET_FILE:faunus>)
set_property(TEST penalty APPEND PROPERTY ENVIRONMENT YASON_EXECUTABLE=${YASON})

add_test(
    NAME minimal
    COMMAND sh -c "${PYTHON_EXECUTABLE} ${YASON} minimal.yml\
    | $<TARGET_FILE:faunus> --nobar --state minimal.state.json\
    ; ${PYTHON_EXECUTABLE} ${JSON_COMPARE} minimal.out.json out.json --tol 0.02"
    WORKING_DIRECTORY ${EXAMPLES_DIR}/minimal)

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
    WORKING_DIRECTORY ${EXAMPLES_DIR}/polymers)

add_test(
    NAME speciation
    COMMAND sh -c "${PYTHON_EXECUTABLE} ${YASON} speciation.yml\
    | $<TARGET_FILE:faunus> --quiet --state speciation.state.json\
    ; ${PYTHON_EXECUTABLE} ${JSON_COMPARE} speciation.out.json out.json --tol 0.02"
    WORKING_DIRECTORY ${EXAMPLES_DIR}/speciation)

add_test(
    NAME phosphate
    COMMAND sh -c "${PYTHON_EXECUTABLE} ${YASON} phosphate.yml\
    | $<TARGET_FILE:faunus> --quiet --state phosphate.state.json\
    ; ${PYTHON_EXECUTABLE} ${JSON_COMPARE} phosphate.out.json out.json --tol 0.02"
    WORKING_DIRECTORY ${EXAMPLES_DIR}/phosphate)

add_test(
    NAME sodiumphosphate
    COMMAND sh -c "${PYTHON_EXECUTABLE} ${YASON} sodiumphosphate.yml\
    | $<TARGET_FILE:faunus> --quiet --state sodiumphosphate.state.json\
    ; ${PYTHON_EXECUTABLE} ${JSON_COMPARE} sodiumphosphate.out.json out.json --tol 0.05"
    WORKING_DIRECTORY ${EXAMPLES_DIR}/sodiumphosphate)

add_test(
    NAME cluster
    COMMAND sh -c "${PYTHON_EXECUTABLE} ${YASON} cluster.yml\
    | $<TARGET_FILE:faunus> --quiet --state cluster.state.json\
    ; ${PYTHON_EXECUTABLE} ${JSON_COMPARE} cluster.out.json out.json --tol 0.1"
    WORKING_DIRECTORY ${EXAMPLES_DIR}/cluster)

add_test(
    NAME bulk
    COMMAND sh -c "${PYTHON_EXECUTABLE} ${YASON} bulk.yml\
    | $<TARGET_FILE:faunus> --quiet"
    WORKING_DIRECTORY ${EXAMPLES_DIR}/bulk)

add_test(
    NAME seawater
    COMMAND sh -c "${PYTHON_EXECUTABLE} ${YASON} seawater.yml\
    | $<TARGET_FILE:faunus> --quiet"
    WORKING_DIRECTORY ${EXAMPLES_DIR}/seawater)

add_test(
    NAME stockmayer
    COMMAND sh -c "${PYTHON_EXECUTABLE} ${YASON} stockmayer.yml\
    | $<TARGET_FILE:faunus> --quiet --state stockmayer.state.json\
    ; ${PYTHON_EXECUTABLE} ${JSON_COMPARE} stockmayer.out.json out.json --tol 0.065"
    WORKING_DIRECTORY ${EXAMPLES_DIR}/stockmayer)

add_test(
    NAME membrane
    COMMAND sh -c "${PYTHON_EXECUTABLE} ${YASON} membrane.yml\
    | $<TARGET_FILE:faunus> --quiet --state membrane.state.json\
    ; ${PYTHON_EXECUTABLE} ${JSON_COMPARE} membrane.out.json out.json --tol 0.065"
    WORKING_DIRECTORY ${EXAMPLES_DIR}/membrane)

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

install(FILES ${CMAKE_SOURCE_DIR}/LICENSE DESTINATION share/faunus)
install(FILES ${CMAKE_BINARY_DIR}/tips.json DESTINATION share/faunus OPTIONAL)
install(FILES ${CMAKE_BINARY_DIR}/manual.pdf DESTINATION share/faunus OPTIONAL)
install(FILES ${CMAKE_BINARY_DIR}/manual.html DESTINATION share/faunus OPTIONAL)
install(DIRECTORY ${CMAKE_BINARY_DIR}/sids DESTINATION share/faunus OPTIONAL)
install(FILES
    ${CMAKE_SOURCE_DIR}/scripts/yason.py
    ${CMAKE_SOURCE_DIR}/scripts/ob2aam.py
    ${CMAKE_SOURCE_DIR}/scripts/faunus-manual
    DESTINATION bin
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

install(FILES
    ${EXAMPLES_DIR}/README.md
    ${EXAMPLES_DIR}/bulk/bulk.yml
    ${EXAMPLES_DIR}/cluster/cluster.yml
    ${EXAMPLES_DIR}/cluster/cluster.agr
    ${EXAMPLES_DIR}/cluster/cluster.out.json
    ${EXAMPLES_DIR}/cluster/cluster.state.json
    ${EXAMPLES_DIR}/cluster/water.xyz
    ${EXAMPLES_DIR}/doublelayer/doublelayer.yml
    ${EXAMPLES_DIR}/doublelayer/doublelayer.agr
    ${EXAMPLES_DIR}/doublelayer/doublelayer.out.json
    ${EXAMPLES_DIR}/doublelayer/doublelayer.state.json
    ${EXAMPLES_DIR}/doublelayer/doublelayer.mfcorr.dat
    ${EXAMPLES_DIR}/isochoric.yml
    ${EXAMPLES_DIR}/isochoric.xyz
    ${EXAMPLES_DIR}/membrane/membrane.yml
    ${EXAMPLES_DIR}/membrane/membrane.out.json
    ${EXAMPLES_DIR}/membrane/membrane.state.json
    ${EXAMPLES_DIR}/lipid.xyz
    ${EXAMPLES_DIR}/minimal/minimal.yml
    ${EXAMPLES_DIR}/minimal/minimal.out.json
    ${EXAMPLES_DIR}/minimal/minimal.state.json
    ${EXAMPLES_DIR}/penalty/penalty.yml
    ${EXAMPLES_DIR}/penalty/penalty.ipynb
    ${EXAMPLES_DIR}/chain/chain.yml
    ${EXAMPLES_DIR}/chain/chain.xyz
    ${EXAMPLES_DIR}/chain/chain.state.json
    ${EXAMPLES_DIR}/chain/chain.out.json
    ${EXAMPLES_DIR}/polymers/polymers.yml
    ${EXAMPLES_DIR}/polymers/polymers.readme
    ${EXAMPLES_DIR}/polymers/polymers.pdf
    ${EXAMPLES_DIR}/polymers/polymers.state.json
    ${EXAMPLES_DIR}/polymers/polymers.out.json
    ${EXAMPLES_DIR}/pythontest.py
    ${EXAMPLES_DIR}/seawater/seawater.yml
    ${EXAMPLES_DIR}/speciation/speciation.yml
    ${EXAMPLES_DIR}/speciation/speciation.out.json
    ${EXAMPLES_DIR}/speciation/speciation.state.json
    ${EXAMPLES_DIR}/stockmayer/stockmayer.yml
    ${EXAMPLES_DIR}/stockmayer/stockmayer.out.json
    ${EXAMPLES_DIR}/stockmayer/stockmayer.state.json
    ${EXAMPLES_DIR}/swapconf/swapconf.yml
    ${EXAMPLES_DIR}/swapconf/swapconf.out.json
    ${EXAMPLES_DIR}/swapconf/swapconf.state.json
    ${EXAMPLES_DIR}/swapconf/swapconf.conformations.pqr
    ${EXAMPLES_DIR}/swapconf/swapconf.weights.dat
    ${EXAMPLES_DIR}/water.yml
    ${EXAMPLES_DIR}/water.xyz
    DESTINATION share/faunus/examples)


