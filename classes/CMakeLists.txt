include_directories(${CMAKE_SOURCE_DIR}/classes)

set(objs inputfile io titrate point physconst slump container
    potentials hardsphere group particles analysis species) 
set_source_files_properties(${objs} PROPERTIES LANGUAGE CXX)
add_library (fns STATIC ${objs})

if (GROMACS)
  if(GMX_FOUND) 
    include_directories(${GMX_INCLUDE_DIR})
    if (GMX_MOTIF AND X11_FOUND AND MOTIF_FOUND)
      add_definitions(-DGROMACS)
      target_link_libraries(fns ${GMX_LIBRARIES} ${X11_Xt_LIB} ${MOTIF_LIBRARIES})
      message("-- Gromacs MOTIF/X dependency - yes")
    endif (GMX_MOTIF AND X11_FOUND AND MOTIF_FOUND)
    if (NOT GMX_MOTIF)
      add_definitions(-DGROMACS)
      target_link_libraries(fns ${GMX_LIBRARIES})
    endif (NOT GMX_MOTIF)
    message("-- Gromacs library: ${GMX_LIBRARIES}")
  endif(GMX_FOUND)
endif(GROMACS)

install(TARGETS fns ARCHIVE DESTINATION lib)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/classes/
  DESTINATION include/faunus FILES_MATCHING PATTERN "*.h"
  PATTERN "CMakeFiles" EXCLUDE
  PATTERN "examples" EXCLUDE)
