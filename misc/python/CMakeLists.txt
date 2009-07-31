if(PYTHONLIBS_FOUND AND PYTHONINTERP_FOUND AND Boost_FOUND)

  if (CVS_FOUND AND WGET_FOUND)

    # GCCXML dl&install target
    add_custom_target( gccxml_install
      COMMAND bash install_gccxml.sh ${CVS_EXECUTABLE}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )

    # Py++ and pygccxml dl&install target
    add_custom_target( pyplusplus_install
      COMMAND bash install_pyplusplus.sh ${PYTHON_EXECUTABLE} ${WGET_EXECUTABLE}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
    add_dependencies( pyplusplus_install gccxml_install )

    # Target for building boost definitions using pyplusplus and
    # pygccxml
    add_custom_target( python-generate
      COMMAND ${PYTHON_EXECUTABLE} generate.py
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
    add_dependencies( python-generate pyplusplys_install gccxml_install )

  endif (CVS_FOUND AND WGET_FOUND)

  # Target for building the python module
  # (only if we have generated python bindings)
  if (EXISTS "generated")
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
      OUTPUT_NAME "faunus" SUFFIX ".so")
  endif (EXISTS "generated")

endif(PYTHONLIBS_FOUND AND PYTHONINTERP_FOUND AND Boost_FOUND)
