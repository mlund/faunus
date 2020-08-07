option(BUILD_STATIC "Build statically linked binaries" off)
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

# place libraries and executables in outermost dir
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# target: functionparser - exprtk very slow, so split out to separate target
add_library(functionparser STATIC ${CMAKE_SOURCE_DIR}/src/functionparser.cpp ${CMAKE_SOURCE_DIR}/src/functionparser.h)
target_compile_definitions(functionparser PRIVATE project_options project_warnings)
set_target_properties(functionparser PROPERTIES POSITION_INDEPENDENT_CODE ON)

# faunus cpp files
set(objs analysis.cpp atomdata.cpp bonds.cpp chainmove.cpp clustermove.cpp core.cpp
        forcemove.cpp units.cpp energy.cpp externalpotential.cpp geometry.cpp group.cpp
        io.cpp molecule.cpp montecarlo.cpp move.cpp mpicontroller.cpp particle.cpp
        penalty.cpp potentials.cpp random.cpp reactioncoordinate.cpp regions.cpp rotate.cpp
        space.cpp speciation.cpp tensor.cpp)

# faunus header files
set(hdrs analysis.h average.h atomdata.h auxiliary.h bonds.h chainmove.h clustermove.h core.h
        forcemove.h energy.h externalpotential.h geometry.h group.h io.h molecule.h montecarlo.h
        move.h mpicontroller.h particle.h penalty.h potentials.h reactioncoordinate.h rotate.h
        space.h speciation.h random.h regions.h tensor.h units.h
        aux/eigen_cerealisation.h aux/eigensupport.h aux/iteratorsupport.h aux/multimatrix.h
        aux/eigen_cerealisation.h aux/eigensupport.h aux/equidistant_table.h aux/error_function.h
        aux/exp_function.h aux/invsqrt_function.h aux/iteratorsupport.h aux/legendre.h aux/multimatrix.h
        aux/pairmatrix.h aux/pow_function.h aux/table_1d.h aux/table_2d.h aux/timers.h)

# faunus test header files
set(tsts atomdata_test.h auxiliary_test.h bonds_test.h core_test.h energy_test.h geometry_test.h
        externalpotential_test.h group_test.h io_test.h molecule_test.h particle_test.h potentials_test.h
        scatter_test.h space_test.h tensor_test.h forcemove_test.h)

#
# target: libfaunus - main functionality of faunus
#
set_source_files_properties(${objs} PROPERTIES LANGUAGE CXX)
add_library(libfaunus STATIC ${objs} ${hdrs})
target_compile_definitions(libfaunus PRIVATE SPDLOG_COMPILED_LIB)
target_link_libraries(libfaunus
    PRIVATE ${CMAKE_THREAD_LIBS_INIT} project_options project_warnings xdrfile spdlog functionparser ZLIB::ZLIB)
set_target_properties(libfaunus PROPERTIES POSITION_INDEPENDENT_CODE ON)
if (MPI_CXX_FOUND)
    target_link_libraries(libfaunus PRIVATE MPI::MPI_CXX)
    target_compile_definitions(libfaunus PRIVATE ENABLE_MPI)
endif()
if(ENABLE_FREESASA)
    target_link_libraries(libfaunus PRIVATE freesasa)
endif()

#
# target: unittests
#
add_executable(unittests EXCLUDE_FROM_ALL unittests.cpp ${tsts})
target_link_libraries(unittests libfaunus)
target_compile_definitions(unittests PRIVATE SPDLOG_COMPILED_LIB)
add_test(
    NAME unittests
    COMMAND sh -c "$<TARGET_FILE:unittests> --test-case-exclude=*Benchmark*")

#
# target: faunus (main executable)
#
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
	include_directories(SYSTEM ${SDL2_INCLUDE_DIRS})
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
install(TARGETS faunus DESTINATION bin)

#
# python bindings
#
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
install(DIRECTORY ${CMAKE_BINARY_DIR}/sids DESTINATION share/faunus OPTIONAL)
install(FILES
    ${CMAKE_SOURCE_DIR}/scripts/yason.py
    ${CMAKE_SOURCE_DIR}/scripts/ob2aam.py
    ${CMAKE_SOURCE_DIR}/scripts/faunus-manual
    DESTINATION bin
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
