if (ENABLE_FORTRAN)
  enable_language(Fortran)

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

  if (CMAKE_Fortran_COMPILER_WORKS)
    # Add exe targets
    fortran_target(bulk bulk "bulk/bulk.f")

    # Add fortran library
    set(fobj functions/vscoul2pc.f functions/gauss.f)
    set_source_files_properties(${fobj}, PROPERTIES LANGUAGE Fortran)
    add_library( fortranfunc SHARED EXCLUDE_FROM_ALL ${fobj}) 
  endif (CMAKE_Fortran_COMPILER_WORKS)

endif (ENABLE_FORTRAN)
