if(PYTHONLIBS_FOUND AND PYTHONINTERP_FOUND AND Boost_FOUND)

  if (GCCXML AND WGET_FOUND)
    # Download and install Py++ and pygccxml in misc/python/
    add_custom_target( pyplusplus_install
      COMMAND bash install_pyplusplus.sh ${PYTHON_EXECUTABLE} ${WGET_EXECUTABLE}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )

    # Target for building boost definitions using pyplusplys and
    # pygccxml (required)
    add_custom_target( pyplusplus_rerun
      COMMAND bash generate.sh ${PYTHON_EXECUTABLE}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
  endif (GCCXML AND WGET_FOUND)

  # Target for building the python module
  include_directories( ${CMAKE_SOURCE_DIR}/include
                       ${Boost_INCLUDE_DIR}
                       ${PYTHON_INCLUDE_PATH}
                       ${CMAKE_CURRENT_SOURCE_DIR}
                     )
  file(GLOB_RECURSE objs generated/ *.cpp)
  set_source_files_properties(${objs} PROPERTIES LANGUAGE CXX)
  add_library(python SHARED ${objs})
  target_link_libraries(python libfaunus ${Boost_PYTHON_LIBRARY} ${PYTHON_LIBRARY})
  set_target_properties(python PROPERTIES PREFIX ""
                        OUTPUT_NAME "faunus" SUFFIX ".so")
endif(PYTHONLIBS_FOUND AND PYTHONINTERP_FOUND AND Boost_FOUND)
