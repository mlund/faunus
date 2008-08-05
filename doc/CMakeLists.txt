find_package(Doxygen)
if(DOXYGEN_FOUND AND DOXYGEN_DOT_EXECUTABLE)
  add_custom_target(manual
    COMMAND ${DOXYGEN_EXECUTABLE} doc/Doxyfile
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR} ) 
  add_custom_target(manualul
    COMMAND scp -rC doc/html/* mikaellund@shell.sourceforge.net:/home/groups/f/fa/faunus/htdocs/
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR} ) 
  add_custom_target(docclean
    COMMAND rm -fR doc/html doc/latex 
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR} ) 
  install( DIRECTORY "${CMAKE_SOURCE_DIR}/doc/html" DESTINATION share/doc/faunus )
endif(DOXYGEN_FOUND AND DOXYGEN_DOT_EXECUTABLE)

