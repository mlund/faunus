include_directories(${CMAKE_SOURCE_DIR}/include)

# -------------------------------
#   Faunus library object files
# -------------------------------
set(hdrs
        ${CMAKE_SOURCE_DIR}/include/faunus/average.h
        ${CMAKE_SOURCE_DIR}/include/faunus/analysis.h
        ${CMAKE_SOURCE_DIR}/include/faunus/common.h
        ${CMAKE_SOURCE_DIR}/include/faunus/auxiliary.h
        ${CMAKE_SOURCE_DIR}/include/faunus/ewald.h
        ${CMAKE_SOURCE_DIR}/include/faunus/externalpotential.h
        ${CMAKE_SOURCE_DIR}/include/faunus/faunus.h
        ${CMAKE_SOURCE_DIR}/include/faunus/json.h
        ${CMAKE_SOURCE_DIR}/include/faunus/point.h
        ${CMAKE_SOURCE_DIR}/include/faunus/move.h
        ${CMAKE_SOURCE_DIR}/include/faunus/inputfile.h
        ${CMAKE_SOURCE_DIR}/include/faunus/io.h
        ${CMAKE_SOURCE_DIR}/include/faunus/geometry.h
        ${CMAKE_SOURCE_DIR}/include/faunus/group.h
        ${CMAKE_SOURCE_DIR}/include/faunus/energy.h
        ${CMAKE_SOURCE_DIR}/include/faunus/species.h
        ${CMAKE_SOURCE_DIR}/include/faunus/mpi.h
        ${CMAKE_SOURCE_DIR}/include/faunus/mcloop.h
        ${CMAKE_SOURCE_DIR}/include/faunus/molecule.h
        ${CMAKE_SOURCE_DIR}/include/faunus/physconst.h
        ${CMAKE_SOURCE_DIR}/include/faunus/potentials.h
        ${CMAKE_SOURCE_DIR}/include/faunus/range.h
        ${CMAKE_SOURCE_DIR}/include/faunus/slump.h
        ${CMAKE_SOURCE_DIR}/include/faunus/scatter.h
        ${CMAKE_SOURCE_DIR}/include/faunus/space.h
        ${CMAKE_SOURCE_DIR}/include/faunus/spherocylinder.h
        ${CMAKE_SOURCE_DIR}/include/faunus/tabulate.h
        ${CMAKE_SOURCE_DIR}/include/faunus/textio.h
        ${CMAKE_SOURCE_DIR}/include/faunus/titrate.h
        ${CMAKE_SOURCE_DIR}/include/faunus/topology.h
        )
set(objs
        geometry potentials slump species analysis)
set_source_files_properties(${objs} PROPERTIES LANGUAGE CXX)

add_subdirectory(xdrfile)

if (ENABLE_NOPBC)
    add_definitions(-DFAU_NOPBC)
endif ()

if (ENABLE_HASHTABLE)
    add_definitions(-DFAU_HASHTABLE)
endif ()

if (NOT ENABLE_UNICODE)
    add_definitions(-DAVOID_UNICODE)
endif ()

# -------------------------------------
#   Use approximate math funtions?
# -------------------------------------
if (ENABLE_APPROXMATH)
    add_definitions(-DFAU_APPROXMATH)
endif ()

# ----------------------------
#  Fetch 3rd-party sasa class
# ----------------------------
if (ENABLE_POWERSASA)
    foreach (powf array.hpp power_diagram.h power_sasa.h LICENSE)
        if (NOT EXISTS "${CMAKE_SOURCE_DIR}/include/faunus/sasa/${powf}")
            file(DOWNLOAD
                    "http://sourceforge.net/p/pteros/code/ci/master/tree/src/core/sasa/${powf}?format=raw"
                    "${CMAKE_SOURCE_DIR}/include/faunus/sasa/${powf}" STATUS status)
        endif ()
    endforeach ()
    add_definitions("-DFAU_POWERSASA")
    include_directories("${CMAKE_SOURCE_DIR}/include/faunus/sasa")
endif ()

# -----------------------
#   Link with openbabel
# -----------------------
if (ENABLE_BABEL)
    find_package(OpenBabel2)
    if (OPENBABEL_FOUND)
        set(objs ${objs} iobabel)
        include_directories(${OPENBABEL_INCLUDE_DIR})
        set(LINKLIBS ${LINKLIBS} ${OPENBABEL_LIBRARIES})
        add_definitions(-DBABEL)
    endif ()
endif ()

# -----------------------
#   Link with MPI
# -----------------------
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

# --------------------
#   Faunus libraries
# --------------------
if (ENABLE_STATIC)
    add_library(libfaunus STATIC ${objs} ${hdrs})
else ()
    add_library(libfaunus SHARED ${objs} ${hdrs})
endif ()
set_target_properties(libfaunus PROPERTIES OUTPUT_NAME faunus)
target_link_libraries(libfaunus xdrfile ${LINKLIBS})
export(TARGETS xdrfile libfaunus FILE libfaunus.cmake)
#install(TARGETS libfaunus LIBRARY DESTINATION lib)

#----- header install target -----
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/faunus
        DESTINATION include
        FILES_MATCHING PATTERN "*.h"
        PATTERN ".git" EXCLUDE
        PATTERN "CMakeFiles" EXCLUDE)

# -----------------------------------------
#   Add "cppcheck" static analyzer target
# -----------------------------------------
find_program(CPPCHECK_EXE "cppcheck")
if (CPPCHECK_EXE)
    add_custom_target(
            cppcheck
            COMMAND ${CPPCHECK_EXE}
            -I ${CMAKE_SOURCE_DIR}/include -j2 --quiet --enable=warning
            #--suppress=noExplicitConstructor --suppress=uninitMemberVar
            ${CMAKE_SOURCE_DIR}/src/*.cpp ${CMAKE_SOURCE_DIR}/src/moves/*.cpp
            COMMENT "Running C++ static analyzer..."
    )
    set_target_properties(cppcheck PROPERTIES EXCLUDE_FROM_ALL TRUE)
endif ()

# -----------------------------------------
#   Function to add a generic C++ example
# -----------------------------------------
function(fau_example tname tdir tsrc)
    add_executable(${tname} "${tdir}/${tsrc}")
    set_source_files_properties("${tdir}/${tsrc}" PROPERTIES LANGUAGE CXX)
    set_target_properties(${tname}
            PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${tdir}/")
    target_link_libraries(${tname} libfaunus)
    install(DIRECTORY "${tdir}"
            DESTINATION "share/faunus/examples"
            PATTERN ".svn" EXCLUDE)
endfunction(fau_example)

#add_subdirectory(playground EXCLUDE_FROM_ALL)
add_subdirectory(examples)

if (EXISTS ${MYPLAYGROUND})
    add_subdirectory(${MYPLAYGROUND} ${MYPLAYGROUND})
endif ()

