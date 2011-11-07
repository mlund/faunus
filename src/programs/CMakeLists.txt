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

#----- function to add a generic C++ example (CIGARPARTICLES) -----
function( fau_cigarexample tname tdir tsrc )
  add_executable( ${tname} "${tdir}/${tsrc}" )
  set_source_files_properties( "${tdir}/${tsrc}" PROPERTIES LANGUAGE CXX)
  set_target_properties(${tname}
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${tdir}/"
    COMPILE_DEFINITIONS "CIGARPARTICLE"
    EXCLUDE_FROM_ALL TRUE )
  target_link_libraries(${tname} libcigarfaunus)
  install (DIRECTORY "${tdir}"
           DESTINATION "share/faunus/examples"
           PATTERN ".svn" EXCLUDE)
endfunction( fau_cigarexample )

#----- add examples to build process -----
# Arguments: (target directory source)
#fau_example(gc "./" gc.cpp)
fau_example(bulk "./" bulk.cpp)
fau_example(polymer-npt "./" polymer-npt.cpp)
fau_cigarexample(cigar "./" cigar.cpp)

#----- target specific build options -----
set_target_properties( cigar PROPERTIES COMPILE_DEFINITIONS "CIGARPARTICLE" EXCLUDE_FROM_ALL TRUE)

#----- examples where openbabel is needed -----
if (ENABLE_BABEL AND OPENBABEL_FOUND)
  include_directories(${OPENBABEL_INCLUDE_DIR})
  #fau_example(babel2aam tools babel2aam.cpp)
endif (ENABLE_BABEL AND OPENBABEL_FOUND)

