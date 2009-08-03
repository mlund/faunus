if(PYTHONLIBS_FOUND AND PYTHONINTERP_FOUND AND Boost_FOUND)

  if (CVS_FOUND AND WGET_FOUND)

    # GCCXML dl&install target
    add_custom_target( gccxml
      COMMAND bash install_gccxml.sh ${CVS_EXECUTABLE}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )

    # Py++ and pygccxml dl&install target
    add_custom_target( pyplusplus
      COMMAND bash install_pyplusplus.sh ${PYTHON_EXECUTABLE} ${WGET_EXECUTABLE}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
    add_dependencies( pyplusplus gccxml )

    # Target for building boost definitions using pyplusplus and
    # pygccxml
    add_custom_target( python-generate
      COMMAND ${PYTHON_EXECUTABLE} generate.py
      COMMAND touch CMakeLists.txt
      COMMAND cmake ../../
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
    add_dependencies( python-generate pyplusplus )

  endif (CVS_FOUND AND WGET_FOUND)

  if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/generated/" )
    # Target for building the python module
    include_directories( ${CMAKE_SOURCE_DIR}/include
      ${Boost_INCLUDE_DIR}
      ${PYTHON_INCLUDE_PATH}
      ${CMAKE_CURRENT_SOURCE_DIR}
      )
    file(GLOB_RECURSE objs generated/ *.cpp)
    #aux_source_directory("generated/" objs)
    set_source_files_properties(${objs} PROPERTIES LANGUAGE CXX)
    add_library(python SHARED ${objs})
    target_link_libraries(python libfaunus ${Boost_PYTHON_LIBRARY} ${PYTHON_LIBRARY})
    set_target_properties(python PROPERTIES PREFIX ""
      LINK_FLAGS "-F/opt/local/Library/Frameworks/"
      OUTPUT_NAME "faunus" SUFFIX ".so")

    # Cleanup python bindings
    add_custom_target( python-clean
      COMMAND rm -fR generated
      COMMAND cmake ../../
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}  )

  endif (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/generated/" )

endif(PYTHONLIBS_FOUND AND PYTHONINTERP_FOUND AND Boost_FOUND)
