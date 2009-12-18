include_directories(${CMAKE_SOURCE_DIR}/include)

#----- object files for libfaunus -----
set(objs slump physconst inputfile point io titrate container
    hardsphere group particles analysis species widom histogram
    mcloop ensemble faunus
    moves/base moves/translational moves/rotational moves/charge
    moves/miscmove moves/rosenbluth moves/saltbath moves/clustermove
    moves/crankshaft moves/membrane) 
set_source_files_properties(${objs} PROPERTIES LANGUAGE CXX)

#----- xdrfile static library for gromacs output -----
add_subdirectory(xdrfile)

#----- link with openbabel ------
if(ENABLE_BABEL AND OPENBABEL_FOUND)
  set(objs ${objs} iobabel)
  include_directories(${OPENBABEL_INCLUDE_DIR})
  set(LINKLIBS ${LINKLIBS} ${OPENBABEL_LIBRARIES})
  add_definitions(-DBABEL)
endif(ENABLE_BABEL AND OPENBABEL_FOUND)

#----- link with fortran library -----
if(ENABLE_FORTRAN AND CMAKE_Fortran_COMPILER_WORKS)
  set(LINKLIBS ${LINKLIBS} fortranfunc)
endif(ENABLE_FORTRAN AND CMAKE_Fortran_COMPILER_WORKS)

#----- static faunus library -----
add_library(libfaunus.a STATIC EXCLUDE_FROM_ALL ${objs})
set_target_properties(libfaunus.a PROPERTIES OUTPUT_NAME faunus)
target_link_libraries(libfaunus.a xdrfile ${LINKLIBS})

#----- dynamic faunus library -----
add_library(libfaunus SHARED ${objs})
set_target_properties(libfaunus PROPERTIES
                      OUTPUT_NAME faunus
                      COMPILE_DEFINITIONS "__SVN_REV__=${SVNREVISION}"
                      )

target_link_libraries(libfaunus xdrfile ${LINKLIBS})
install(TARGETS libfaunus
        LIBRARY DESTINATION lib)

#----- dynamic faunus library - HYPERSPHERE version -----
set (compdef "__SVN_REV__=${SVNREVISION}" "HYPERSPHERE")
add_library(libhyperfaunus SHARED ${objs})
set_target_properties(libhyperfaunus PROPERTIES
                      OUTPUT_NAME hyperfaunus
                      COMPILE_DEFINITIONS "HYPERSPHERE"
                      EXCLUDE_FROM_ALL TRUE
                      )

target_link_libraries(libhyperfaunus xdrfile ${LINKLIBS})

#----- header install target -----
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/faunus
        DESTINATION include
        FILES_MATCHING PATTERN "*.h"
        PATTERN ".svn" EXCLUDE
        PATTERN "CMakeFiles" EXCLUDE)

#----- Add example programs -----
add_subdirectory(examples)

