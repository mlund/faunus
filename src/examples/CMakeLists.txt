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
# Arguments: (target directory source)
fau_example(pka titration pka.cpp)
fau_example(pkaGC titration pka.cpp)
#fau_example(pkaDH titration pkaDH.cpp)
fau_example(twobody twobody twobody.cpp)
fau_example(twobodyDH twobody twobodyDH.cpp)
fau_example(twobody-hof twobody-hofmeister twobody-hof.cpp)
fau_example(widom_cube widom widom.cpp)
fau_example(widom_sphere widom widom.cpp)
fau_example(manybody manybody manybody.cpp)
fau_example(manybodyMP manybody manybody.cpp)
fau_example(isobaric isobaric isobaric.cpp)
fau_example(isobaricMP isobaric isobaric.cpp)
fau_example(isobaricPLANE isobaric isobaric.cpp)
#fau_example(isobaricFDH isobaric isobaric.cpp)
fau_example(binding binding binding.cpp)
fau_example(printpotential tools printpotential.cpp)
fau_example(random tools random.cpp)
fau_example(cgp1 tools cg1.cpp)
fau_example(cgp2 tools cg2.cpp)
fau_example(pmf2muex tools pmf2muex.cpp)
fau_example(g_widom tools g_widom.cpp)
fau_example(spce spce spce.cpp)
fau_example(spcecavity spce spcecavity.cpp)
fau_example(acid spce spceacid.cpp)
fau_example(spceMI spce spce.cpp)
fau_example(rosenbluth rosenbluth rb.cpp)
fau_example(membrane membrane mem.cpp)
fau_example(memp wallprotein memp.cpp)

#----- target specific build options -----
set_target_properties( pkaGC PROPERTIES COMPILE_DEFINITIONS "GCPKA")
set_target_properties( manybodyMP PROPERTIES COMPILE_DEFINITIONS "MONOPOLE")
set_target_properties( isobaricMP PROPERTIES COMPILE_DEFINITIONS "MONOPOLE")
set_target_properties( isobaricPLANE PROPERTIES COMPILE_DEFINITIONS "XYPLANE")
#set_target_properties( isobaricFDH PROPERTIES COMPILE_DEFINITIONS "FASTDH")
set_target_properties( widom_sphere PROPERTIES COMPILE_DEFINITIONS "WIDOM_SPHERE")
set_target_properties( spceMI PROPERTIES COMPILE_DEFINITIONS "SPCE_MINIMAGE")

#----- examples where openbabel is needed -----
if (ENABLE_BABEL AND OPENBABEL_FOUND)
  include_directories(${OPENBABEL_INCLUDE_DIR})
  fau_example(babel2aam tools babel2aam.cpp)
  fau_example(wallprotein wallprotein wp.cpp)
endif (ENABLE_BABEL AND OPENBABEL_FOUND)

