include_directories(${CMAKE_SOURCE_DIR}/classes)

# Function to add a generic C++ example
function( fau_example tname tsrc )
  add_executable( ${tname} ${tsrc} )
  set_source_files_properties( ${tsource} PROPERTIES LANGUAGE CXX)
  target_link_libraries(${tname} fns ${GMX_LIBRARIES})
  install (TARGETS ${tname} RUNTIME DESTINATION "bin")
endfunction( fau_example )

fau_example(pka titration/pka.C)
fau_example(GCpka titration/GCpka.C)
fau_example(twobody twobody/twobody.C)
fau_example(twobody-hof twobody-hofmeister/twobody-hof.C)
fau_example(widom widom/widom.C)
fau_example(manybody manybody/manybody.C)
fau_example(isobaric isobaric/isobaric.C)

