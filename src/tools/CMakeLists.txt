include_directories(${CMAKE_SOURCE_DIR}/include)

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
#fau_example(polymer-npt-sphere "./" polymer-npt.cpp)
#fau_example(polymer-npt-cuboid "./" polymer-npt.cpp)
#fau_example(titrate_implicit "./" pka_implicit.cpp)
#fau_example(manybody "./" manybody.cpp)

if(ENABLE_MPI AND MPI_CXX_FOUND)
  fau_example(manybodyMPI "./" manybodyMPI.cpp)
  fau_example(manybodyMPI_slit "./" manybodyMPI.cpp)
  set_target_properties( manybodyMPI_slit PROPERTIES COMPILE_DEFINITIONS "SLIT" EXCLUDE_FROM_ALL TRUE)
endif()

#----- target specific build options -----
#set_target_properties( polymer-npt-cuboid PROPERTIES COMPILE_DEFINITIONS "CUBOID" EXCLUDE_FROM_ALL TRUE)

#----- examples where openbabel is needed -----
if (ENABLE_BABEL AND OPENBABEL_FOUND)
  include_directories(${OPENBABEL_INCLUDE_DIR})
  #fau_example(babel2aam tools babel2aam.cpp)
endif (ENABLE_BABEL AND OPENBABEL_FOUND)

