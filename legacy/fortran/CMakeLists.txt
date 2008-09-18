enable_language(Fortran OPTIONAL)

# Function to add a generic fortran target
function( fortran_target tname tdir tsrc )
  set_source_files_properties( "${tsrc}" PROPERTIES LANGUAGE Fortran)
  add_executable( ${tname} "${tsrc}" )
  set_target_properties(${tname}
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${tdir}/")
  install (DIRECTORY "${tdir}" DESTINATION "share/faunus/examples"
    PATTERN ".svn" EXCLUDE
    PATTERN "*.f" EXCLUDE
    PATTERN "*.inc" EXCLUDE)
endfunction( fortran_target )

# Add targets
if (CMAKE_Fortran_COMPILER_WORKS)
  fortran_target(bulk bulk "bulk/bulk.f")
endif (CMAKE_Fortran_COMPILER_WORKS)
