include_directories(${CMAKE_SOURCE_DIR}/include)

#----- object files for libfaunus -----
set(objs inputfile io titrate point physconst slump container
    hardsphere group particles analysis species widom histogram
    moves/base moves/translational moves/rotational) 
set_source_files_properties(${objs} PROPERTIES LANGUAGE CXX)

#----- link with openbabel ------
if(ENABLE_BABEL AND OPENBABEL_FOUND)
  set(objs ${objs} iobabel)
  include_directories(${OPENBABEL_INCLUDE_DIR})
  set(LINKLIBS ${LINKLIBS} ${OPENBABEL_LIBRARIES})
endif(ENABLE_BABEL AND OPENBABEL_FOUND)

#---- link with libxtc -----
if (ENABLE_LIBXTC)
  add_subdirectory(libxtc)
  add_definitions(-DLIBXTC)
  include_directories("${CMAKE_SOURCE_DIR}/src/libxtc/")
  set(LINKLIBS ${LINKLIBS} libxtc)
endif (ENABLE_LIBXTC)

#----- link with gromacs -----
if (ENABLE_GROMACS AND GMX_FOUND)
  add_definitions(-DGROMACS)
  include_directories(${GMX_INCLUDE_DIR})
  set(LINKLIBS ${LINKLIBS} ${GMX_LIBRARIES})
endif (ENABLE_GROMACS AND GMX_FOUND)

#----- static faunus library -----
add_library(libfaunus.a STATIC EXCLUDE_FROM_ALL ${objs})
set_target_properties(libfaunus.a PROPERTIES OUTPUT_NAME faunus)
target_link_libraries(libfaunus.a ${LINKLIBS})

#----- dynamic faunus library -----
add_library(libfaunus SHARED ${objs})
set_target_properties(libfaunus PROPERTIES OUTPUT_NAME faunus)
target_link_libraries(libfaunus ${LINKLIBS})
install(TARGETS libfaunus
        LIBRARY DESTINATION lib)

#----- header install target -----
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/faunus
        DESTINATION include
        FILES_MATCHING PATTERN "*.h"
        PATTERN ".svn" EXCLUDE
        PATTERN "CMakeFiles" EXCLUDE)

#----- Add example programs -----
add_subdirectory(examples)

