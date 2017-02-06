find_package(Doxygen)
if(DOXYGEN_FOUND)
  add_custom_target(manual
    COMMAND ${DOXYGEN_EXECUTABLE} doc/Doxyfile
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR} ) 
  install( DIRECTORY "${CMAKE_SOURCE_DIR}/doc/doxyhtml" DESTINATION share/doc/faunus )
endif(DOXYGEN_FOUND)

