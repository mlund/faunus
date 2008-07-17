include_directories($faunus_source_dir/classes)
set(objs inputfile io titrate point physconst slump container
    potentials hardsphere group particles analysis species) 
set_source_files_properties(${objs} PROPERTIES LANGUAGE CXX)

add_library ("fns" SHARED ${objs})

if (GROMACS)
  if(GMX_FOUND AND X11_FOUND AND MOTIF_FOUND)
    message("-- Gromacs library: ${GMX_LIBRARIES}")
    include_directories(${GMX_INCLUDE_DIR})
    target_link_libraries("fns" ${GMX_LIBRARIES} ${X11_Xt_LIB} ${MOTIF_LIBRARIES})
    add_definitions(-DGROMACS)
  endif(GMX_FOUND AND X11_FOUND AND MOTIF_FOUND)
endif(GROMACS)


