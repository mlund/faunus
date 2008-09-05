enable_language(Fortran OPTIONAL)
include_directories(${CMAKE_SOURCE_DIR}/legacy/fortran)

# Function to add a generic fortran target
function( fortran_target tname tdir tsrc )
  set_source_files_properties( "${tsrc}" PROPERTIES LANGUAGE Fortran)
  add_executable( ${tname} "${tsrc}" )
  set_target_properties(${tname}
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${tdir}/")
  install (TARGETS ${tname} RUNTIME DESTINATION "bin")
endfunction( fortran_target )

# Add targets:
fortran_target(bulk bulk "bulk/bulk.f")


