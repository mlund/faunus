if(PYTHONLIBS_FOUND AND PYTHONINTERP_FOUND AND Boost_FOUND)
  include_directories( ${CMAKE_SOURCE_DIR}/include
                       ${Boost_INCLUDE_DIR}
                       ${PYTHON_INCLUDE_PATH} ) 
  file(GLOB_RECURSE objs generated/ *.cpp)
  set_source_files_properties(${objs} PROPERTIES LANGUAGE CXX)
  add_library(pyfaunus SHARED ${objs})
  target_link_libraries(pyfaunus libfaunus ${Boost_PYTHON_LIBRARY} ${PYTHON_LIBRARY})
  set_target_properties(pyfaunus PROPERTIES PREFIX "")
endif(PYTHONLIBS_FOUND AND PYTHONINTERP_FOUND AND Boost_FOUND)
