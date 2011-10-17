find_package(Doxygen)
if(DOXYGEN_FOUND AND DOXYGEN_DOT_EXECUTABLE)
  add_custom_target(manual
    COMMAND ${DOXYGEN_EXECUTABLE} doc/Doxyfile
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR} ) 
  add_custom_target(upload_homepage_and_manual
    COMMAND scp -rC doc/homepage/* mikaellund,faunus@web.sourceforge.net:htdocs
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR} ) 
  add_custom_target(upload_homepage
    COMMAND scp -rC doc/homepage/*.html doc/homepage/images mikaellund,faunus@web.sourceforge.net:htdocs
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR} ) 
  add_custom_target(docclean
    COMMAND rm -fR doc/homepage/doxyhtml doc/latex 
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR} ) 
  install( DIRECTORY "${CMAKE_SOURCE_DIR}/doc/homepage" DESTINATION share/doc/faunus )
endif(DOXYGEN_FOUND AND DOXYGEN_DOT_EXECUTABLE)

