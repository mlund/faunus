include_directories(${CMAKE_SOURCE_DIR}/include)

#----- function to add a generic C++ example -----
function( fau_example tname tdir tsrc )
  add_executable( ${tname} "${tdir}/${tsrc}" )
  set_source_files_properties( "${tdir}/${tsrc}" PROPERTIES LANGUAGE CXX)
  set_target_properties(${tname}
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${tdir}/")
  target_link_libraries(${tname} libfaunus)
  install (DIRECTORY "${tdir}"
           DESTINATION "share/faunus/examples"
           PATTERN ".svn" EXCLUDE)
endfunction( fau_example )

#----- add examples to build process -----
fau_example(pka titration pka.cpp)
fau_example(pkaGC titration pka.cpp)
fau_example(twobody twobody twobody.cpp)
fau_example(twobody-hof twobody-hofmeister twobody-hof.cpp)
fau_example(widom_cube widom widom.cpp)
fau_example(widom_sphere widom widom.cpp)
fau_example(manybody manybody manybody.cpp)
fau_example(isobaric isobaric isobaric.cpp)
fau_example(binding binding binding.cpp)

#----- target specific build options -----
set_target_properties( pkaGC PROPERTIES COMPILE_DEFINITIONS "GCPKA")
set_target_properties( widom_sphere PROPERTIES COMPILE_DEFINITIONS "WIDOM_SPHERE")

#----- examples where openbabel is needed -----
if (ENABLE_BABEL AND OPENBABEL_FOUND)
  include_directories(${OPENBABEL_INCLUDE_DIR})
  fau_example(babel2aam tools babel2aam.cpp)
endif (ENABLE_BABEL AND OPENBABEL_FOUND)

