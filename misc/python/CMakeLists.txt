if(PYTHONLIBS_FOUND AND PYTHONINTERP_FOUND AND Boost_FOUND)
  # Target for building boost definitions using pyplusplys and
  # pygccxml (required)
  # Perhaps add "generate.sh" that sets environment and
  # calls generate.py (python_exe passes as argument)?
  add_custom_target( updatepython
                     COMMAND ${PYTHON_EXECUTABLE} generate.py
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                   )

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
