include_directories(${CMAKE_SOURCE_DIR}/classes)
include_directories(${OPENBABEL_INCLUDE_DIR})

# Function to add a generic C++ example
function( fau_example tname tdir tsrc )
  add_executable( ${tname} "${tdir}/${tsrc}" )
  set_source_files_properties( "${tdir}/${tsrc}" PROPERTIES LANGUAGE CXX)
  set_target_properties(${tname}
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${tdir}/")
  target_link_libraries(${tname} fns ${GMX_LIBRARIES})
  install (TARGETS ${tname} RUNTIME DESTINATION "bin")
endfunction( fau_example )

fau_example(pka titration pka.C)
fau_example(GCpka titration GCpka.C)
fau_example(twobody twobody twobody.C)
fau_example(twobody-hof twobody-hofmeister twobody-hof.C)
fau_example(widom widom widom.C)
fau_example(manybody manybody manybody.C)
fau_example(isobaric isobaric isobaric.C)

if (OPENBABEL_FOUND)
  fau_example(babel2aam tools babel2aam.C)
endif (OPENBABEL_FOUND)
