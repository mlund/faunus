include_directories(${CMAKE_SOURCE_DIR}/include)

#----- object files for libfaunus -----
set(objs inputfile io titrate point physconst slump container
    potentials hardsphere group particles analysis species) 

#----- locate openbabel ------
if (OPENBABEL_FOUND)
  set(objs ${objs} iobabel)
  include_directories(${OPENBABEL_INCLUDE_DIR})
  set(LINKLIBS ${LINKLIBS} ${OPENBABEL_LIBRARIES})
endif (OPENBABEL_FOUND)

#----- locate gromacs -----
if (GROMACS AND GMX_FOUND)
  include_directories(${GMX_INCLUDE_DIR})
  if (NOT GMX_MOTIF)
    add_definitions(-DGROMACS)
    set(LINKLIBS ${LINKLIBS} ${GMX_LIBRARIES})
  endif (NOT GMX_MOTIF)
  if (GMX_MOTIF AND X11_FOUND AND MOTIF_FOUND)
    add_definitions(-DGROMACS)
    set(LINKLIBS ${LINKLIBS} ${GMX_LIBRARIES} ${X11_Xt_LIB} ${MOTIF_LIBRARIES})
  endif (GMX_MOTIF AND X11_FOUND AND MOTIF_FOUND)
endif(GROMACS AND GMX_FOUND)

#----- add libfaunus targets -----
set_source_files_properties(${objs} PROPERTIES LANGUAGE CXX)
add_library (fns STATIC ${objs})
target_link_libraries(fns ${LINKLIBS})
add_library (fnsdyn SHARED ${objs})
target_link_libraries(fnsdyn ${LINKLIBS})

#----- set install targets -----
install(TARGETS fns ARCHIVE DESTINATION lib)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/faunus
  DESTINATION include/faunus FILES_MATCHING PATTERN "*.h"
  PATTERN "CMakeFiles" EXCLUDE)

