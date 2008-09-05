include_directories(${CMAKE_SOURCE_DIR}/classes)

#----- function to add a generic C++ example -----
function( fau_example tname tdir tsrc )
  add_executable( ${tname} "${tdir}/${tsrc}" )
  set_source_files_properties( "${tdir}/${tsrc}" PROPERTIES LANGUAGE CXX)
  set_target_properties(${tname}
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${tdir}/")
  target_link_libraries(${tname} fns)
  install (TARGETS ${tname} RUNTIME DESTINATION "bin")
endfunction( fau_example )

#----- add examples to build process -----
fau_example(pka titration pka.C)
fau_example(GCpka titration GCpka.C)
fau_example(twobody twobody twobody.C)
fau_example(twobody-hof twobody-hofmeister twobody-hof.C)
fau_example(widom widom widom.C)
fau_example(manybody manybody manybody.C)
fau_example(isobaric isobaric isobaric.C)
fau_example(binding binding binding.C)

#----- examples where openbabel is needed -----
if (OPENBABEL_FOUND)
  include_directories(${OPENBABEL_INCLUDE_DIR})
  fau_example(babel2aam tools babel2aam.C)
endif (OPENBABEL_FOUND)
